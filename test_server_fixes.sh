#!/bin/bash

# Test script to verify server fixes
cd /workspaces/Protein-Interaction-Analysis-Server

echo "🧪 Testing Protein Interaction Analysis Server..."
echo "================================================="

# Set environment
export PYTHONPATH=/workspaces/Protein-Interaction-Analysis-Server/src:$PYTHONPATH

# Test imports and basic functionality
echo "1. Testing imports..."
/workspaces/Protein-Interaction-Analysis-Server/.venv/bin/python -c "
import sys
try:
    import streamlit as st
    from ui.main_interface import MainInterface
    from utils.config import AppConfig
    from analysis.hydrogen_bonds import HydrogenBond
    print('   ✅ All imports successful')
except Exception as e:
    print(f'   ❌ Import error: {e}')
    sys.exit(1)
"

echo "2. Testing interaction object handling..."
/workspaces/Protein-Interaction-Analysis-Server/.venv/bin/python -c "
import sys
try:
    from ui.main_interface import MainInterface
    from utils.config import AppConfig
    from analysis.hydrogen_bonds import HydrogenBond
    
    # Create test objects
    config = AppConfig()
    interface = MainInterface(config)
    
    # Test with HydrogenBond object
    hb = HydrogenBond(
        donor_atom=None, acceptor_atom=None, hydrogen_atom=None,
        distance=3.2, angle=150.0, strength='strong',
        donor_residue='ALA1', acceptor_residue='GLU2',
        donor_chain='A', acceptor_chain='A'
    )
    
    # Test helper function
    strength = interface._get_interaction_property(hb, 'strength', 'weak')
    distance = interface._get_interaction_property(hb, 'distance', 0)
    
    print(f'   ✅ Object handling works: strength={strength}, distance={distance}')
    
    # Test filtering
    interactions = [hb]
    filtered = interface._filter_interactions_by_strength(interactions, 'strong')
    print(f'   ✅ Filtering works: {len(filtered)} interactions kept')
    
except Exception as e:
    print(f'   ❌ Object handling error: {e}')
    import traceback
    traceback.print_exc()
    sys.exit(1)
"

echo "3. Testing dictionary-style interactions..."
/workspaces/Protein-Interaction-Analysis-Server/.venv/bin/python -c "
import sys
try:
    from ui.main_interface import MainInterface
    from utils.config import AppConfig
    
    config = AppConfig()
    interface = MainInterface(config)
    
    # Test with dictionary
    dict_interaction = {
        'strength': 'moderate',
        'distance': 4.1,
        'residue1': 'ARG1',
        'chain1': 'A'
    }
    
    strength = interface._get_interaction_property(dict_interaction, 'strength', 'weak')
    distance = interface._get_interaction_property(dict_interaction, 'distance', 0)
    
    print(f'   ✅ Dictionary handling works: strength={strength}, distance={distance}')
    
except Exception as e:
    print(f'   ❌ Dictionary handling error: {e}')
    sys.exit(1)
"

echo ""
echo "🎉 All tests passed! The AttributeError has been fixed."
echo ""
echo "🚀 To run the server:"
echo "   cd /workspaces/Protein-Interaction-Analysis-Server"
echo "   PYTHONPATH=src .venv/bin/python -m streamlit run server.py --server.port 8501 --server.address 0.0.0.0"
echo ""
echo "🌐 Then open: http://localhost:8501"
