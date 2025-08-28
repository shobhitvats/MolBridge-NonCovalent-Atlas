"""
REST API endpoints for Protein Interaction Explorer.
Provides programmatic access to analysis capabilities.
"""

from fastapi import FastAPI, HTTPException, BackgroundTasks, UploadFile, File, Depends, Query
from fastapi.responses import StreamingResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from pydantic import BaseModel, Field
from typing import List, Dict, Any, Optional, Union
import asyncio
import io
import uuid
import json
from datetime import datetime
from pathlib import Path

from analysis.batch_processor import BatchProcessor
from utils.config import AppConfig, InteractionConfig
from utils.pdb_handler import PDBHandler
from utils.cache import CacheManager
from reporting.report_generator import ReportGenerator

# Initialize app
app = FastAPI(
    title="Protein Interaction Explorer API",
    description="REST API for comprehensive noncovalent interaction analysis",
    version="1.0.0"
)

# Configure CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Global configuration
config = AppConfig()
pdb_handler = PDBHandler(config)
cache_manager = CacheManager(config)
batch_processor = BatchProcessor(config)
report_generator = ReportGenerator(config)

# Job tracking
active_jobs: Dict[str, Dict[str, Any]] = {}

# Pydantic models
class AnalysisRequest(BaseModel):
    pdb_ids: List[str] = Field(..., description="List of PDB identifiers")
    interaction_types: List[str] = Field(
        default=["hydrogen_bond", "halogen_bond", "pi_pi", "ionic", "hydrophobic"],
        description="Types of interactions to analyze"
    )
    use_cache: bool = Field(default=True, description="Use cached results if available")
    force_refresh: bool = Field(default=False, description="Force refresh of cached data")
    include_metadata: bool = Field(default=True, description="Include structure metadata")
    config_preset: str = Field(default="literature_default", description="Configuration preset")

class ConfigurationUpdate(BaseModel):
    interaction_config: Optional[Dict[str, Any]] = None
    processing_config: Optional[Dict[str, Any]] = None

class ReportRequest(BaseModel):
    pdb_ids: List[str]
    format: str = Field(default="pdf", description="Report format: pdf, excel, powerpoint, csv, latex")
    include_metadata: bool = Field(default=True)
    include_methodology: bool = Field(default=True)
    include_visualizations: bool = Field(default=False)

class AnalysisResponse(BaseModel):
    job_id: str
    status: str
    message: str
    started_at: datetime
    pdb_ids: List[str]

class JobStatus(BaseModel):
    job_id: str
    status: str
    progress: float
    message: str
    started_at: datetime
    completed_at: Optional[datetime] = None
    results: Optional[Dict[str, Any]] = None
    error: Optional[str] = None

# API Endpoints

@app.get("/", tags=["General"])
async def root():
    """Root endpoint with API information."""
    return {
        "name": "Protein Interaction Explorer API",
        "version": "1.0.0",
        "description": "REST API for comprehensive noncovalent interaction analysis",
        "endpoints": {
            "analysis": "/analyze",
            "status": "/jobs/{job_id}",
            "results": "/results/{job_id}",
            "reports": "/reports",
            "config": "/config"
        }
    }

@app.get("/health", tags=["General"])
async def health_check():
    """Health check endpoint."""
    return {
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "services": {
            "cache": cache_manager.get_cache_stats(),
            "config": "active"
        }
    }

@app.post("/analyze", response_model=AnalysisResponse, tags=["Analysis"])
async def start_analysis(
    request: AnalysisRequest,
    background_tasks: BackgroundTasks
):
    """
    Start protein interaction analysis for multiple structures.
    
    This endpoint initiates a background analysis job and returns immediately
    with a job ID. Use the /jobs/{job_id} endpoint to check progress and
    retrieve results.
    """
    # Generate unique job ID
    job_id = str(uuid.uuid4())
    
    # Validate input
    if not request.pdb_ids:
        raise HTTPException(status_code=400, detail="At least one PDB ID is required")
    
    if len(request.pdb_ids) > 100:
        raise HTTPException(status_code=400, detail="Maximum 100 structures per request")
    
    # Update configuration if needed
    if request.config_preset in ["conservative", "literature_default", "exploratory"]:
        config.set_preset(request.config_preset)
    
    # Initialize job tracking
    active_jobs[job_id] = {
        "job_id": job_id,
        "status": "queued",
        "progress": 0.0,
        "message": "Analysis queued",
        "started_at": datetime.now(),
        "completed_at": None,
        "pdb_ids": request.pdb_ids,
        "request": request.dict(),
        "results": None,
        "error": None
    }
    
    # Start background analysis
    background_tasks.add_task(
        run_analysis_job,
        job_id,
        request
    )
    
    return AnalysisResponse(
        job_id=job_id,
        status="queued",
        message="Analysis started",
        started_at=active_jobs[job_id]["started_at"],
        pdb_ids=request.pdb_ids
    )

@app.get("/jobs/{job_id}", response_model=JobStatus, tags=["Analysis"])
async def get_job_status(job_id: str):
    """Get status and progress of an analysis job."""
    if job_id not in active_jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job_info = active_jobs[job_id]
    
    return JobStatus(**job_info)

@app.get("/results/{job_id}", tags=["Analysis"])
async def get_analysis_results(job_id: str):
    """Get complete results from a finished analysis job."""
    if job_id not in active_jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job_info = active_jobs[job_id]
    
    if job_info["status"] == "completed":
        return job_info["results"]
    elif job_info["status"] == "failed":
        raise HTTPException(status_code=500, detail=f"Analysis failed: {job_info['error']}")
    else:
        raise HTTPException(status_code=202, detail="Analysis still in progress")

@app.post("/reports", tags=["Reports"])
async def generate_report(request: ReportRequest):
    """Generate and download analysis reports in various formats."""
    
    # Validate format
    supported_formats = ["pdf", "excel", "powerpoint", "csv", "latex"]
    if request.format not in supported_formats:
        raise HTTPException(
            status_code=400, 
            detail=f"Unsupported format. Supported: {supported_formats}"
        )
    
    # Get analysis results for requested structures
    analysis_results = {}
    for pdb_id in request.pdb_ids:
        # Try to get from cache or recent jobs
        result = None
        for job_info in active_jobs.values():
            if (job_info["status"] == "completed" and 
                pdb_id in job_info["pdb_ids"] and 
                job_info["results"]):
                result = job_info["results"].get(pdb_id)
                break
        
        if result is None:
            # Run quick analysis for this structure
            try:
                results = batch_processor.process_batch([pdb_id])
                result = results.get(pdb_id)
            except Exception as e:
                raise HTTPException(
                    status_code=500,
                    detail=f"Failed to analyze {pdb_id}: {str(e)}"
                )
        
        analysis_results[pdb_id] = result
    
    # Generate report
    try:
        if request.format == "pdf":
            content = report_generator.generate_pdf_report(
                request.pdb_ids,
                analysis_results,
                request.dict()
            )
            media_type = "application/pdf"
            filename = f"protein_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pdf"
            
        elif request.format == "excel":
            content = report_generator.generate_excel_report(
                request.pdb_ids,
                analysis_results
            )
            media_type = "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            filename = f"protein_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx"
            
        elif request.format == "powerpoint":
            content = report_generator.generate_powerpoint_report(
                request.pdb_ids,
                analysis_results
            )
            media_type = "application/vnd.openxmlformats-officedocument.presentationml.presentation"
            filename = f"protein_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.pptx"
            
        elif request.format == "csv":
            if len(request.pdb_ids) == 1:
                content = report_generator.generate_csv_report(
                    request.pdb_ids[0],
                    analysis_results[request.pdb_ids[0]]
                ).encode('utf-8')
            else:
                raise HTTPException(
                    status_code=400,
                    detail="CSV format only supported for single structures"
                )
            media_type = "text/csv"
            filename = f"protein_analysis_{request.pdb_ids[0]}_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv"
            
        elif request.format == "latex":
            content = report_generator.generate_latex_export(
                request.pdb_ids,
                analysis_results
            ).encode('utf-8')
            media_type = "text/plain"
            filename = f"protein_analysis_{datetime.now().strftime('%Y%m%d_%H%M%S')}.tex"
        
        # Return as streaming response
        return StreamingResponse(
            io.BytesIO(content),
            media_type=media_type,
            headers={"Content-Disposition": f"attachment; filename={filename}"}
        )
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Report generation failed: {str(e)}")

@app.get("/config", tags=["Configuration"])
async def get_configuration():
    """Get current analysis configuration."""
    return {
        "interaction_config": config.interaction_config.__dict__,
        "processing_config": config.processing_config.__dict__,
        "visualization_config": config.visualization_config.__dict__,
        "version": config.version,
        "available_presets": ["conservative", "literature_default", "exploratory"]
    }

@app.post("/config", tags=["Configuration"])
async def update_configuration(update: ConfigurationUpdate):
    """Update analysis configuration."""
    try:
        if update.interaction_config:
            for key, value in update.interaction_config.items():
                if hasattr(config.interaction_config, key):
                    setattr(config.interaction_config, key, value)
        
        if update.processing_config:
            for key, value in update.processing_config.items():
                if hasattr(config.processing_config, key):
                    setattr(config.processing_config, key, value)
        
        return {"message": "Configuration updated successfully"}
    
    except Exception as e:
        raise HTTPException(status_code=400, detail=f"Configuration update failed: {str(e)}")

@app.post("/upload", tags=["Data"])
async def upload_structure(file: UploadFile = File(...)):
    """Upload and analyze custom PDB structure file."""
    if not file.filename.endswith(('.pdb', '.cif')):
        raise HTTPException(status_code=400, detail="Only PDB and CIF files are supported")
    
    try:
        # Save uploaded file temporarily
        content = await file.read()
        temp_path = Path(f"/tmp/{uuid.uuid4()}_{file.filename}")
        temp_path.write_bytes(content)
        
        # Analyze the structure
        structure = pdb_handler.load_structure_from_file(str(temp_path))
        if structure is None:
            raise HTTPException(status_code=400, detail="Invalid structure file")
        
        # Run analysis
        results = batch_processor.analyze_single_structure(structure, str(temp_path))
        
        # Clean up
        temp_path.unlink()
        
        return {
            "filename": file.filename,
            "analysis_results": results,
            "message": "File analyzed successfully"
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")

@app.get("/structures/{pdb_id}/info", tags=["Data"])
async def get_structure_info(pdb_id: str):
    """Get basic information about a PDB structure."""
    try:
        structure = pdb_handler.load_structure(pdb_id)
        if structure is None:
            raise HTTPException(status_code=404, detail="Structure not found")
        
        # Extract basic information
        info = pdb_handler._process_structure(structure)
        
        return {
            "pdb_id": pdb_id,
            "structure_info": info,
            "message": "Structure information retrieved successfully"
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get structure info: {str(e)}")

@app.delete("/jobs/{job_id}", tags=["Analysis"])
async def cancel_job(job_id: str):
    """Cancel a running analysis job."""
    if job_id not in active_jobs:
        raise HTTPException(status_code=404, detail="Job not found")
    
    job_info = active_jobs[job_id]
    
    if job_info["status"] in ["completed", "failed", "cancelled"]:
        return {"message": "Job already finished"}
    
    # Mark as cancelled
    job_info["status"] = "cancelled"
    job_info["message"] = "Job cancelled by user"
    job_info["completed_at"] = datetime.now()
    
    return {"message": "Job cancelled successfully"}

@app.get("/cache/stats", tags=["System"])
async def get_cache_stats():
    """Get cache statistics and status."""
    return cache_manager.get_cache_stats()

@app.post("/cache/clear", tags=["System"])
async def clear_cache():
    """Clear all cached data."""
    try:
        cache_manager.clear_cache()
        return {"message": "Cache cleared successfully"}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Cache clear failed: {str(e)}")

# Background task functions

async def run_analysis_job(job_id: str, request: AnalysisRequest):
    """Background task to run protein interaction analysis."""
    try:
        # Update job status
        active_jobs[job_id]["status"] = "running"
        active_jobs[job_id]["message"] = "Analysis in progress"
        
        # Run analysis
        results = batch_processor.process_batch(
            request.pdb_ids,
            interaction_types=request.interaction_types,
            use_cache=request.use_cache,
            force_refresh=request.force_refresh,
            include_metadata=request.include_metadata,
            progress_callback=lambda progress, msg: update_job_progress(job_id, progress, msg)
        )
        
        # Update job with results
        active_jobs[job_id]["status"] = "completed"
        active_jobs[job_id]["progress"] = 100.0
        active_jobs[job_id]["message"] = "Analysis completed successfully"
        active_jobs[job_id]["completed_at"] = datetime.now()
        active_jobs[job_id]["results"] = results
        
    except Exception as e:
        # Update job with error
        active_jobs[job_id]["status"] = "failed"
        active_jobs[job_id]["message"] = "Analysis failed"
        active_jobs[job_id]["completed_at"] = datetime.now()
        active_jobs[job_id]["error"] = str(e)

def update_job_progress(job_id: str, progress: float, message: str):
    """Update job progress."""
    if job_id in active_jobs:
        active_jobs[job_id]["progress"] = progress
        active_jobs[job_id]["message"] = message

# Exception handlers
@app.exception_handler(404)
async def not_found_handler(request, exc):
    return JSONResponse(
        status_code=404,
        content={"detail": "Resource not found"}
    )

@app.exception_handler(500)
async def internal_error_handler(request, exc):
    return JSONResponse(
        status_code=500,
        content={"detail": "Internal server error"}
    )

if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=8000)
