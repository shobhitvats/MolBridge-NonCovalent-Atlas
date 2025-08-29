"""
Report generation for Protein Interaction Explorer.
Generates PDF, PowerPoint, and other export formats.
"""

import io
import json
from datetime import datetime
from typing import Dict, Any, List, Optional
from pathlib import Path
import base64

import pandas as pd
from reportlab.lib.pagesizes import letter, A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, PageBreak
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_LEFT
from openpyxl import Workbook
from openpyxl.styles import Font, PatternFill, Alignment
from pptx import Presentation
from pptx.util import Inches
from pptx.enum.text import PP_ALIGN

from utils.config import AppConfig

class ReportGenerator:
    """Generates comprehensive reports from interaction analysis results."""
    
    def __init__(self, config: AppConfig):
        self.config = config
    
    def _get_interaction_property(self, interaction: Any, property_name: str, default_value: Any = None) -> Any:
        """
        Safely get a property from an interaction object, whether it's a dictionary or an object with attributes.

        Args:
            interaction: The interaction object (dict or custom object)
            property_name: Name of the property to retrieve
            default_value: Default value if property is not found

        Returns:
            The property value or default_value
        """
        if hasattr(interaction, property_name):
            # Object-style interaction (e.g., HydrogenBond dataclass)
            return getattr(interaction, property_name, default_value)
        elif isinstance(interaction, dict):
            # Dictionary-style interaction
            return interaction.get(property_name, default_value)
        else:
            # Unknown format, return default
            return default_value
    
    def generate_pdf_report(self,
                           pdb_ids: List[str],
                           analysis_results: Dict[str, Dict[str, Any]],
                           options: Dict[str, Any]) -> bytes:
        """
        Generate comprehensive PDF report.
        
        Args:
            pdb_ids: List of PDB identifiers
            analysis_results: Analysis results for each PDB
            options: Report generation options
            
        Returns:
            PDF report as bytes
        """
        buffer = io.BytesIO()
        doc = SimpleDocTemplate(buffer, pagesize=A4)
        styles = getSampleStyleSheet()
        story = []
        
        # Custom styles
        title_style = ParagraphStyle(
            'CustomTitle',
            parent=styles['Heading1'],
            fontSize=24,
            spaceAfter=30,
            alignment=TA_CENTER
        )
        
        heading_style = ParagraphStyle(
            'CustomHeading',
            parent=styles['Heading2'],
            fontSize=16,
            spaceAfter=12,
            spaceBefore=20
        )
        
        # Cover page
        story.extend(self._create_cover_page(pdb_ids, title_style, styles))
        story.append(PageBreak())
        
        # Table of contents
        if options.get('include_metadata', True):
            story.extend(self._create_table_of_contents(pdb_ids, analysis_results, heading_style))
            story.append(PageBreak())
        
        # Executive summary
        story.extend(self._create_executive_summary(pdb_ids, analysis_results, heading_style, styles))
        story.append(PageBreak())
        
        # Methodology section
        if options.get('include_methodology', True):
            story.extend(self._create_methodology_section(heading_style, styles))
            story.append(PageBreak())
        
        # Individual structure reports
        for pdb_id in pdb_ids:
            if pdb_id in analysis_results and analysis_results[pdb_id]:
                story.extend(self._create_structure_report(
                    pdb_id, 
                    analysis_results[pdb_id], 
                    heading_style, 
                    styles,
                    options
                ))
                story.append(PageBreak())
        
        # Comparative analysis (if multiple structures)
        if len(pdb_ids) > 1:
            story.extend(self._create_comparative_analysis(
                pdb_ids,
                analysis_results,
                heading_style,
                styles
            ))
            story.append(PageBreak())
        
        # Appendix
        story.extend(self._create_appendix(pdb_ids, analysis_results, heading_style, styles))
        
        # Build PDF
        doc.build(story)
        buffer.seek(0)
        return buffer.getvalue()
    
    def _create_cover_page(self, pdb_ids: List[str], title_style, styles) -> List:
        """Create cover page elements."""
        story = []
        
        # Title
        story.append(Paragraph("Protein Interaction Analysis Report", title_style))
        story.append(Spacer(1, 0.5*inch))
        
        # Structure list
        if len(pdb_ids) == 1:
            story.append(Paragraph(f"Structure: {pdb_ids[0]}", styles['Heading2']))
        else:
            story.append(Paragraph(f"Structures Analyzed: {len(pdb_ids)}", styles['Heading2']))
            structure_list = ", ".join(pdb_ids[:10])
            if len(pdb_ids) > 10:
                structure_list += f" and {len(pdb_ids) - 10} more"
            story.append(Paragraph(structure_list, styles['Normal']))
        
        story.append(Spacer(1, 0.5*inch))
        
        # Metadata
        story.append(Paragraph(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}", styles['Normal']))
        story.append(Paragraph(f"Software: Protein Interaction Explorer v{self.config.version}", styles['Normal']))
        story.append(Spacer(1, 1*inch))
        
        # Summary logo/image placeholder
        story.append(Paragraph("🧬 Comprehensive Noncovalent Interaction Analysis", styles['Heading3']))
        
        return story
    
    def _create_table_of_contents(self, pdb_ids: List[str], analysis_results: Dict, heading_style) -> List:
        """Create table of contents."""
        story = []
        story.append(Paragraph("Table of Contents", heading_style))
        story.append(Spacer(1, 0.3*inch))
        
        toc_items = [
            "Executive Summary",
            "Methodology",
        ]
        
        # Add individual structures
        for pdb_id in pdb_ids:
            toc_items.append(f"Structure Analysis: {pdb_id}")
        
        if len(pdb_ids) > 1:
            toc_items.append("Comparative Analysis")
        
        toc_items.append("Appendix")
        
        for i, item in enumerate(toc_items, 1):
            story.append(Paragraph(f"{i}. {item}", getSampleStyleSheet()['Normal']))
        
        return story
    
    def _create_executive_summary(self, pdb_ids: List[str], analysis_results: Dict, heading_style, styles) -> List:
        """Create executive summary."""
        story = []
        story.append(Paragraph("Executive Summary", heading_style))
        
        # Calculate overall statistics
        total_interactions = 0
        interaction_type_counts = {}
        successful_analyses = 0
        
        for pdb_id in pdb_ids:
            if pdb_id in analysis_results and analysis_results[pdb_id]:
                successful_analyses += 1
                interactions = analysis_results[pdb_id].get('interactions', {})
                
                for int_type, int_list in interactions.items():
                    count = len(int_list)
                    total_interactions += count
                    interaction_type_counts[int_type] = interaction_type_counts.get(int_type, 0) + count
        
        # Summary text
        summary_text = f"""
        This report presents a comprehensive analysis of noncovalent interactions in {len(pdb_ids)} protein structure(s).
        The analysis identified a total of {total_interactions} interactions across {len(interaction_type_counts)} different 
        interaction types. {successful_analyses} of {len(pdb_ids)} structures were successfully analyzed.
        """
        
        story.append(Paragraph(summary_text, styles['Normal']))
        story.append(Spacer(1, 0.3*inch))
        
        # Summary table
        if interaction_type_counts:
            summary_data = [['Interaction Type', 'Total Count', 'Average per Structure']]
            
            for int_type, count in sorted(interaction_type_counts.items(), key=lambda x: x[1], reverse=True):
                display_name = self._get_interaction_display_name(int_type)
                avg_count = count / successful_analyses if successful_analyses > 0 else 0
                summary_data.append([display_name, str(count), f"{avg_count:.1f}"])
            
            table = Table(summary_data)
            table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.grey),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 12),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.beige),
                ('GRID', (0, 0), (-1, -1), 1, colors.black)
            ]))
            
            story.append(table)
        
        return story
    
    def _create_methodology_section(self, heading_style, styles) -> List:
        """Create methodology section."""
        story = []
        story.append(Paragraph("Methodology", heading_style))
        
        methodology_text = f"""
        <b>Analysis Parameters:</b><br/>
        This analysis was performed using the Protein Interaction Explorer (version {self.config.version}),
        which implements state-of-the-art algorithms for detecting noncovalent interactions in protein structures.
        
        <br/><br/>
        <b>Interaction Types Analyzed:</b><br/>
        • Hydrogen bonds (conventional, low-barrier, C5-type)<br/>
        • Halogen bonds (Cl, Br, I, F interactions)<br/>
        • π-π stacking interactions<br/>
        • Ionic interactions<br/>
        • Hydrophobic contacts<br/>
        • C-H···π interactions<br/>
        • Chalcogen bonds<br/>
        • Pnictogen bonds<br/>
        • Tetrel bonds<br/>
        • Anion-π interactions<br/>
        • n→π* interactions<br/>
        • London dispersion forces (heuristic)<br/>
        
        <br/>
        <b>Detection Criteria:</b><br/>
        All interactions were detected using literature-based geometric criteria including distance and angle cutoffs.
        The analysis includes both intra- and inter-molecular interactions, with specific handling for different
        assembly types (biological assembly vs. asymmetric unit).
        
        <br/><br/>
        <b>Quality Control:</b><br/>
        Structures underwent validation checks including completeness assessment and clash detection.
        Missing hydrogen atoms were added using standard geometric algorithms where applicable.
        """
        
        story.append(Paragraph(methodology_text, styles['Normal']))
        
        return story
    
    def _create_structure_report(self, 
                               pdb_id: str, 
                               result: Dict[str, Any], 
                               heading_style, 
                               styles,
                               options: Dict[str, Any]) -> List:
        """Create report section for individual structure."""
        story = []
        story.append(Paragraph(f"Structure Analysis: {pdb_id}", heading_style))
        
        # Structure information
        structure_info = result.get('structure_info', {})
        metadata = result.get('metadata', {})
        
        info_text = f"""
        <b>PDB ID:</b> {pdb_id}<br/>
        <b>Analysis Time:</b> {metadata.get('analysis_time', 0):.2f} seconds<br/>
        <b>Total Chains:</b> {len(structure_info.get('chains', []))}<br/>
        <b>Total Residues:</b> {structure_info.get('total_residues', 0)}<br/>
        <b>Total Atoms:</b> {structure_info.get('total_atoms', 0)}<br/>
        """
        
        story.append(Paragraph(info_text, styles['Normal']))
        story.append(Spacer(1, 0.2*inch))
        
        # Interaction summary
        interactions = result.get('interactions', {})
        total_interactions = sum(len(int_list) for int_list in interactions.values())
        
        story.append(Paragraph(f"Total Interactions Detected: {total_interactions}", styles['Heading3']))
        
        # Detailed interaction table
        if total_interactions > 0:
            interaction_data = [['Interaction Type', 'Count', 'Average Distance (Å)']]
            
            for int_type, int_list in interactions.items():
                if int_list:
                    display_name = self._get_interaction_display_name(int_type)
                    count = len(int_list)
                    
                    distances = [self._get_interaction_property(i, 'distance', 0) for i in int_list if self._get_interaction_property(i, 'distance', None) is not None]
                    avg_distance = sum(distances) / len(distances) if distances else 'N/A'
                    avg_distance_str = f"{avg_distance:.2f}" if isinstance(avg_distance, float) else avg_distance
                    
                    interaction_data.append([display_name, str(count), avg_distance_str])
            
            if len(interaction_data) > 1:
                table = Table(interaction_data)
                table.setStyle(TableStyle([
                    ('BACKGROUND', (0, 0), (-1, 0), colors.lightblue),
                    ('TEXTCOLOR', (0, 0), (-1, 0), colors.black),
                    ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                    ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                    ('FONTSIZE', (0, 0), (-1, 0), 10),
                    ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                    ('BACKGROUND', (0, 1), (-1, -1), colors.lightgrey),
                    ('GRID', (0, 0), (-1, -1), 1, colors.black)
                ]))
                
                story.append(table)
        
        # Hotspots
        hotspots = result.get('hotspots', [])
        if hotspots:
            story.append(Spacer(1, 0.2*inch))
            story.append(Paragraph("Interaction Hotspots", styles['Heading3']))
            
            hotspot_text = "Top interaction hotspots (residues with highest interaction density):<br/>"
            for i, hotspot in enumerate(hotspots[:5], 1):
                hotspot_text += f"{i}. {hotspot['residue']} ({hotspot['interaction_count']} interactions)<br/>"
            
            story.append(Paragraph(hotspot_text, styles['Normal']))
        
        return story
    
    def _create_comparative_analysis(self, 
                                   pdb_ids: List[str],
                                   analysis_results: Dict[str, Dict[str, Any]],
                                   heading_style,
                                   styles) -> List:
        """Create comparative analysis section."""
        story = []
        story.append(Paragraph("Comparative Analysis", heading_style))
        
        # Create comparison table
        comparison_data = [['PDB ID', 'Total Interactions'] + [
            self._get_interaction_display_name(int_type) 
            for int_type in ['hydrogen_bond', 'halogen_bond', 'pi_pi', 'ionic', 'hydrophobic']
        ]]
        
        for pdb_id in pdb_ids:
            if pdb_id in analysis_results and analysis_results[pdb_id]:
                interactions = analysis_results[pdb_id].get('interactions', {})
                total = sum(len(int_list) for int_list in interactions.values())
                
                row = [pdb_id, str(total)]
                for int_type in ['hydrogen_bond', 'halogen_bond', 'pi_pi', 'ionic', 'hydrophobic']:
                    count = len(interactions.get(int_type, []))
                    row.append(str(count))
                
                comparison_data.append(row)
        
        if len(comparison_data) > 1:
            table = Table(comparison_data)
            table.setStyle(TableStyle([
                ('BACKGROUND', (0, 0), (-1, 0), colors.darkblue),
                ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),
                ('ALIGN', (0, 0), (-1, -1), 'CENTER'),
                ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
                ('FONTSIZE', (0, 0), (-1, 0), 9),
                ('BOTTOMPADDING', (0, 0), (-1, 0), 12),
                ('BACKGROUND', (0, 1), (-1, -1), colors.lightsteelblue),
                ('GRID', (0, 0), (-1, -1), 1, colors.black)
            ]))
            
            story.append(table)
        
        return story
    
    def _create_appendix(self, pdb_ids: List[str], analysis_results: Dict, heading_style, styles) -> List:
        """Create appendix with detailed data."""
        story = []
        story.append(Paragraph("Appendix", heading_style))
        
        appendix_text = f"""
        <b>Software Information:</b><br/>
        Protein Interaction Explorer v{self.config.version}<br/>
        Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br/>
        
        <br/>
        <b>Data Provenance:</b><br/>
        Structures downloaded from RCSB Protein Data Bank<br/>
        Analysis performed using literature-based geometric criteria<br/>
        
        <br/>
        <b>Citation:</b><br/>
        If you use this analysis in your research, please cite:<br/>
        Protein Interaction Explorer: Comprehensive Analysis of Noncovalent Interactions in Protein Structures<br/>
        """
        
        story.append(Paragraph(appendix_text, styles['Normal']))
        
        return story
    
    def generate_csv_report(self, pdb_id: str, result: Dict[str, Any]) -> str:
        """Generate CSV report for a single structure."""
        interactions = result.get('interactions', {})
        
        # Combine all interactions into a single list
        all_interactions = []
        
        for interaction_type, interaction_list in interactions.items():
            for interaction in interaction_list:
                interaction_data = {
                    'PDB_ID': pdb_id,
                    'Interaction_Type': self._get_interaction_display_name(interaction_type),
                    'Residue_1': self._get_interaction_property(interaction, 'residue1', ''),
                    'Chain_1': self._get_interaction_property(interaction, 'chain1', ''),
                    'Residue_2': self._get_interaction_property(interaction, 'residue2', ''),
                    'Chain_2': self._get_interaction_property(interaction, 'chain2', ''),
                    'Distance_A': self._get_interaction_property(interaction, 'distance', ''),
                    'Angle_deg': self._get_interaction_property(interaction, 'angle', ''),
                    'Strength': self._get_interaction_property(interaction, 'strength', ''),
                    'Atom_1': self._get_interaction_property(interaction, 'atom1', ''),
                    'Atom_2': self._get_interaction_property(interaction, 'atom2', '')
                }
                all_interactions.append(interaction_data)
        
        # Convert to DataFrame and then CSV
        df = pd.DataFrame(all_interactions)
        return df.to_csv(index=False)
    
    def generate_excel_report(self, 
                            pdb_ids: List[str],
                            analysis_results: Dict[str, Dict[str, Any]]) -> bytes:
        """Generate Excel report with multiple sheets."""
        buffer = io.BytesIO()
        
        with pd.ExcelWriter(buffer, engine='openpyxl') as writer:
            # Summary sheet
            summary_data = []
            for pdb_id in pdb_ids:
                if pdb_id in analysis_results and analysis_results[pdb_id]:
                    interactions = analysis_results[pdb_id].get('interactions', {})
                    total = sum(len(int_list) for int_list in interactions.values())
                    
                    row = {'PDB_ID': pdb_id, 'Total_Interactions': total}
                    for int_type, int_list in interactions.items():
                        display_name = self._get_interaction_display_name(int_type).replace(' ', '_')
                        row[display_name] = len(int_list)
                    
                    summary_data.append(row)
            
            if summary_data:
                summary_df = pd.DataFrame(summary_data)
                summary_df.to_excel(writer, sheet_name='Summary', index=False)
            
            # Individual sheets for each structure
            for pdb_id in pdb_ids:
                if pdb_id in analysis_results and analysis_results[pdb_id]:
                    csv_data = self.generate_csv_report(pdb_id, analysis_results[pdb_id])
                    df = pd.read_csv(io.StringIO(csv_data))
                    
                    if not df.empty:
                        sheet_name = f"Interactions_{pdb_id}"[:31]  # Excel sheet name limit
                        df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        buffer.seek(0)
        return buffer.getvalue()
    
    def generate_powerpoint_report(self,
                                 pdb_ids: List[str],
                                 analysis_results: Dict[str, Dict[str, Any]]) -> bytes:
        """Generate PowerPoint presentation."""
        prs = Presentation()
        
        # Title slide
        title_slide_layout = prs.slide_layouts[0]
        slide = prs.slides.add_slide(title_slide_layout)
        title = slide.shapes.title
        subtitle = slide.placeholders[1]
        
        title.text = "Protein Interaction Analysis"
        subtitle.text = f"Analysis of {len(pdb_ids)} Structure(s)\n{datetime.now().strftime('%Y-%m-%d')}"
        
        # Summary slide
        bullet_slide_layout = prs.slide_layouts[1]
        slide = prs.slides.add_slide(bullet_slide_layout)
        shapes = slide.shapes
        
        title_shape = shapes.title
        body_shape = shapes.placeholders[1]
        
        title_shape.text = 'Analysis Summary'
        
        tf = body_shape.text_frame
        tf.text = f'Structures Analyzed: {len(pdb_ids)}'
        
        # Calculate totals
        total_interactions = 0
        for pdb_id in pdb_ids:
            if pdb_id in analysis_results and analysis_results[pdb_id]:
                interactions = analysis_results[pdb_id].get('interactions', {})
                total_interactions += sum(len(int_list) for int_list in interactions.values())
        
        p = tf.add_paragraph()
        p.text = f'Total Interactions: {total_interactions}'
        
        # Individual structure slides
        for pdb_id in pdb_ids:
            if pdb_id in analysis_results and analysis_results[pdb_id]:
                self._add_structure_slide(prs, pdb_id, analysis_results[pdb_id])
        
        # Save to buffer
        buffer = io.BytesIO()
        prs.save(buffer)
        buffer.seek(0)
        return buffer.getvalue()
    
    def _add_structure_slide(self, prs, pdb_id: str, result: Dict[str, Any]):
        """Add slide for individual structure."""
        slide_layout = prs.slide_layouts[1]
        slide = prs.slides.add_slide(slide_layout)
        
        title_shape = slide.shapes.title
        body_shape = slide.placeholders[1]
        
        title_shape.text = f'Structure: {pdb_id}'
        
        tf = body_shape.text_frame
        interactions = result.get('interactions', {})
        total = sum(len(int_list) for int_list in interactions.values())
        
        tf.text = f'Total Interactions: {total}'
        
        for int_type, int_list in interactions.items():
            if int_list:
                p = tf.add_paragraph()
                display_name = self._get_interaction_display_name(int_type)
                p.text = f'{display_name}: {len(int_list)}'
    
    def generate_latex_export(self,
                            pdb_ids: List[str],
                            analysis_results: Dict[str, Dict[str, Any]]) -> str:
        """Generate LaTeX snippets for publication."""
        latex_content = []
        
        # Document header
        latex_content.append(r"\section{Protein Interaction Analysis}")
        latex_content.append("")
        
        # Summary table
        latex_content.append(r"\subsection{Interaction Summary}")
        latex_content.append(r"\begin{table}[h]")
        latex_content.append(r"\centering")
        latex_content.append(r"\caption{Summary of detected noncovalent interactions}")
        latex_content.append(r"\begin{tabular}{|l|c|}")
        latex_content.append(r"\hline")
        latex_content.append(r"\textbf{Structure} & \textbf{Total Interactions} \\")
        latex_content.append(r"\hline")
        
        for pdb_id in pdb_ids:
            if pdb_id in analysis_results and analysis_results[pdb_id]:
                interactions = analysis_results[pdb_id].get('interactions', {})
                total = sum(len(int_list) for int_list in interactions.values())
                latex_content.append(f"{pdb_id} & {total} \\\\")
        
        latex_content.append(r"\hline")
        latex_content.append(r"\end{tabular}")
        latex_content.append(r"\label{tab:interaction_summary}")
        latex_content.append(r"\end{table}")
        latex_content.append("")
        
        # Methods section
        latex_content.append(r"\subsection{Methods}")
        latex_content.append("Noncovalent interactions were analyzed using the Protein Interaction Explorer ")
        latex_content.append(f"software (version {self.config.version}). The analysis included detection of ")
        latex_content.append("hydrogen bonds, halogen bonds, $\\pi$-$\\pi$ stacking interactions, ionic ")
        latex_content.append("interactions, and hydrophobic contacts using literature-based geometric criteria.")
        latex_content.append("")
        
        return "\n".join(latex_content)
    
    def _get_interaction_display_name(self, interaction_type: str) -> str:
        """Get human-readable name for interaction type."""
        display_names = {
            'hydrogen_bond': 'Hydrogen Bonds',
            'halogen_bond': 'Halogen Bonds',
            'pi_pi': 'π-π Stacking',
            'ionic': 'Ionic Interactions',
            'hydrophobic': 'Hydrophobic Contacts',
            'ch_pi': 'C-H···π Interactions',
            'chalcogen_bond': 'Chalcogen Bonds',
            'pnictogen_bond': 'Pnictogen Bonds',
            'tetrel_bond': 'Tetrel Bonds',
            'anion_pi': 'Anion-π Interactions',
            'n_pi_star': 'n→π* Interactions',
            'dispersion': 'London Dispersion'
        }
        return display_names.get(interaction_type, interaction_type.replace('_', ' ').title())
