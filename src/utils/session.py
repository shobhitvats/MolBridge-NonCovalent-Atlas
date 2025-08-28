"""
Session management for Protein Interaction Explorer.
Handles saving/loading analysis sessions and state management.
"""

import json
import uuid
from datetime import datetime, timedelta
from pathlib import Path
from typing import Dict, Any, Optional, List
import streamlit as st
import pickle
from dataclasses import asdict

class SessionManager:
    """Manages user sessions and analysis state."""
    
    def __init__(self, sessions_dir: Optional[Path] = None):
        self.sessions_dir = sessions_dir or Path("cache/sessions")
        self.sessions_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize session state
        if 'session_id' not in st.session_state:
            st.session_state.session_id = str(uuid.uuid4())
        
        if 'session_data' not in st.session_state:
            st.session_state.session_data = self._create_empty_session()
    
    def _create_empty_session(self) -> Dict[str, Any]:
        """Create an empty session structure."""
        return {
            "session_id": st.session_state.get('session_id', str(uuid.uuid4())),
            "created_at": datetime.now().isoformat(),
            "updated_at": datetime.now().isoformat(),
            "version": "1.0.0",
            "settings": {
                "interaction_cutoffs": {},
                "visualization_settings": {},
                "selected_interactions": [],
                "preset": "literature_default"
            },
            "analysis_results": {},
            "uploaded_pdbs": [],
            "bookmarks": [],
            "user_notes": {},
            "metadata": {
                "total_structures": 0,
                "total_interactions": 0,
                "analysis_time": 0.0
            }
        }
    
    def save_session(self, session_name: Optional[str] = None) -> str:
        """Save current session to disk."""
        session_data = st.session_state.session_data.copy()
        session_data["updated_at"] = datetime.now().isoformat()
        
        if session_name:
            session_data["name"] = session_name
        
        session_file = self.sessions_dir / f"{session_data['session_id']}.json"
        
        try:
            with open(session_file, 'w') as f:
                json.dump(session_data, f, indent=2, default=str)
            
            return session_data['session_id']
        except Exception as e:
            st.error(f"Failed to save session: {e}")
            return ""
    
    def load_session(self, session_id: str) -> bool:
        """Load session from disk."""
        session_file = self.sessions_dir / f"{session_id}.json"
        
        if not session_file.exists():
            st.error(f"Session {session_id} not found")
            return False
        
        try:
            with open(session_file, 'r') as f:
                session_data = json.load(f)
            
            st.session_state.session_data = session_data
            st.session_state.session_id = session_id
            
            return True
        except Exception as e:
            st.error(f"Failed to load session: {e}")
            return False
    
    def list_sessions(self) -> List[Dict[str, Any]]:
        """List all available sessions."""
        sessions = []
        
        for session_file in self.sessions_dir.glob("*.json"):
            try:
                with open(session_file, 'r') as f:
                    session_data = json.load(f)
                
                sessions.append({
                    "session_id": session_data["session_id"],
                    "name": session_data.get("name", "Unnamed Session"),
                    "created_at": session_data["created_at"],
                    "updated_at": session_data["updated_at"],
                    "total_structures": session_data["metadata"]["total_structures"],
                    "total_interactions": session_data["metadata"]["total_interactions"]
                })
            except Exception:
                continue
        
        # Sort by updated time (most recent first)
        sessions.sort(key=lambda x: x["updated_at"], reverse=True)
        return sessions
    
    def delete_session(self, session_id: str) -> bool:
        """Delete a session."""
        session_file = self.sessions_dir / f"{session_id}.json"
        
        try:
            if session_file.exists():
                session_file.unlink()
            return True
        except Exception as e:
            st.error(f"Failed to delete session: {e}")
            return False
    
    def export_session(self, session_id: str) -> Optional[bytes]:
        """Export session as downloadable file."""
        session_file = self.sessions_dir / f"{session_id}.json"
        
        if not session_file.exists():
            return None
        
        try:
            return session_file.read_bytes()
        except Exception:
            return None
    
    def import_session(self, session_data: bytes) -> bool:
        """Import session from uploaded file."""
        try:
            session_dict = json.loads(session_data.decode())
            
            # Validate session structure
            required_keys = ["session_id", "settings", "analysis_results"]
            if not all(key in session_dict for key in required_keys):
                st.error("Invalid session file format")
                return False
            
            # Generate new session ID to avoid conflicts
            session_dict["session_id"] = str(uuid.uuid4())
            session_dict["imported_at"] = datetime.now().isoformat()
            
            st.session_state.session_data = session_dict
            st.session_state.session_id = session_dict["session_id"]
            
            return True
        except Exception as e:
            st.error(f"Failed to import session: {e}")
            return False
    
    def add_bookmark(self, pdb_id: str, residue_info: Dict[str, Any], note: str = ""):
        """Add a bookmark for a specific residue."""
        bookmark = {
            "id": str(uuid.uuid4()),
            "pdb_id": pdb_id,
            "residue_info": residue_info,
            "note": note,
            "created_at": datetime.now().isoformat()
        }
        
        if 'bookmarks' not in st.session_state.session_data:
            st.session_state.session_data['bookmarks'] = []
        
        st.session_state.session_data['bookmarks'].append(bookmark)
    
    def remove_bookmark(self, bookmark_id: str):
        """Remove a bookmark."""
        if 'bookmarks' in st.session_state.session_data:
            st.session_state.session_data['bookmarks'] = [
                b for b in st.session_state.session_data['bookmarks'] 
                if b['id'] != bookmark_id
            ]
    
    def add_note(self, pdb_id: str, note: str):
        """Add a note for a PDB structure."""
        if 'user_notes' not in st.session_state.session_data:
            st.session_state.session_data['user_notes'] = {}
        
        st.session_state.session_data['user_notes'][pdb_id] = {
            "note": note,
            "updated_at": datetime.now().isoformat()
        }
    
    def get_note(self, pdb_id: str) -> str:
        """Get note for a PDB structure."""
        return st.session_state.session_data.get('user_notes', {}).get(pdb_id, {}).get('note', '')
    
    def update_metadata(self, **kwargs):
        """Update session metadata."""
        if 'metadata' not in st.session_state.session_data:
            st.session_state.session_data['metadata'] = {}
        
        st.session_state.session_data['metadata'].update(kwargs)
        st.session_state.session_data['updated_at'] = datetime.now().isoformat()
    
    def compare_sessions(self, session_id_1: str, session_id_2: str) -> Dict[str, Any]:
        """Compare two sessions and return differences."""
        session_1 = self._load_session_data(session_id_1)
        session_2 = self._load_session_data(session_id_2)
        
        if not session_1 or not session_2:
            return {}
        
        comparison = {
            "session_1_id": session_id_1,
            "session_2_id": session_id_2,
            "structures_comparison": {
                "session_1_count": session_1["metadata"]["total_structures"],
                "session_2_count": session_2["metadata"]["total_structures"],
                "difference": session_2["metadata"]["total_structures"] - session_1["metadata"]["total_structures"]
            },
            "interactions_comparison": {
                "session_1_count": session_1["metadata"]["total_interactions"], 
                "session_2_count": session_2["metadata"]["total_interactions"],
                "difference": session_2["metadata"]["total_interactions"] - session_1["metadata"]["total_interactions"]
            },
            "new_structures": [],
            "removed_structures": [],
            "interaction_differences": {}
        }
        
        # Compare structure lists
        structures_1 = set(session_1.get("uploaded_pdbs", []))
        structures_2 = set(session_2.get("uploaded_pdbs", []))
        
        comparison["new_structures"] = list(structures_2 - structures_1)
        comparison["removed_structures"] = list(structures_1 - structures_2)
        
        return comparison
    
    def _load_session_data(self, session_id: str) -> Optional[Dict[str, Any]]:
        """Load session data without affecting current session."""
        session_file = self.sessions_dir / f"{session_id}.json"
        
        if not session_file.exists():
            return None
        
        try:
            with open(session_file, 'r') as f:
                return json.load(f)
        except Exception:
            return None
    
    def cleanup_old_sessions(self, days: int = 30):
        """Clean up sessions older than specified days."""
        cutoff_date = datetime.now() - timedelta(days=days)
        
        for session_file in self.sessions_dir.glob("*.json"):
            try:
                with open(session_file, 'r') as f:
                    session_data = json.load(f)
                
                updated_at = datetime.fromisoformat(session_data["updated_at"])
                
                if updated_at < cutoff_date:
                    session_file.unlink()
            except Exception:
                continue
    
    def get_session_share_url(self, session_id: str) -> str:
        """Generate a shareable URL for a session."""
        # In a production environment, this would generate a secure token
        # For now, we'll use a simple approach
        base_url = st.get_option("server.baseUrlPath") or "http://localhost:8501"
        return f"{base_url}?session={session_id}"
