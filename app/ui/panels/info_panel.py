"""Info panel extracted from monolithic MainInterface for modularity."""
from __future__ import annotations
import streamlit as st


def render():
    st.header("ℹ️ Interaction Detection Criteria")
    st.write("This section details the criteria used to detect each type of noncovalent interaction, along with their sources.")

    with st.expander("Hydrogen Bonds", expanded=True):
        st.markdown(
            """
            **Criteria:**  
            - Distance: 3.5 Å (donor–acceptor)  
            - Angle: 120° (D–H···A minimum)  
            **Source:** IUPAC Recommendations 2011  
            [doi.org/10.1351/PAC-REP-10-01-01](https://doi.org/10.1351/PAC-REP-10-01-01)
            """
        )
    with st.expander("Halogen Bonds", expanded=True):
        st.markdown(
            """
            **Criteria:**  
            - Distance: ≤ Σ(vdW radii)  
            - Angle: ~160° preferred (C–X···A)  
            **Source:** IUPAC Recommendations 2013  
            [doi.org/10.1351/PAC-REC-12-05-10](https://doi.org/10.1351/PAC-REC-12-05-10)
            """
        )
    with st.expander("Pnictogen Bonds", expanded=True):
        st.markdown(
            """
            **Criteria:**  
            - Distance: ≤ Σ(vdW radii)  
            - Angle: ≥150° (σ-hole alignment)  
            **Source:** IUPAC 2023  
            [doi.org/10.1515/pac-2020-1002](https://doi.org/10.1515/pac-2020-1002)
            """
        )
    with st.expander("Chalcogen Bonds", expanded=True):
        st.markdown(
            """
            **Criteria:**  
            - Distance (Chalcogen atom → acceptor): ≤ **4.0 Å** (typical exploration 3.2–4.5 Å)  
            - θ alignment (C–S···A or substituent–S···A axis): **115° – 155°**  
            - |φ| out-of-plane deviation: **≤ 50°** (planarity / approach coplanarity)  
            - Directionality arises from σ-hole on S/Se/Te; tighter θ window → higher specificity  
            **Source:** Review of unconventional noncovalent (chalcogen) interactions in protein structure  
            [doi.org/10.1017/qrd.2023.3](https://doi.org/10.1017/qrd.2023.3)
            """
        )
    with st.expander("Other Interactions (π-π, C-H···π, Anion-π, etc.)", expanded=True):
        st.markdown(
            """
            **Criteria (Representative):**  
            - **π-π Stacking:** Centroid–centroid ≤ **5.5 Å** (face-to-face often 3.8–4.1 Å); interplanar angle ≤ **30°** for parallel; T-shaped recognized at ~60°–120°.  
            - **C–H···π:** C (projected donor carbon) to ring centroid 2.0–4.8 Å; C–H···centroid angle ≥ **90°**; perpendicular offset ≤ **2.5 Å**.  
            - **Anion–π:** Anion center to ring centroid ≤ **5.0 Å** (enrichment <4.3 Å); angle between anion→centroid vector and ring normal ≤ **30°** (favors approach along positive quadrupole axis).  
            - **n→π*** (carbonyl n→π*): O(i)→C(i+1) distance 2.9–3.6 Å (≤3.5 Å default); approach angle (O lone pair vector proxy) **95°–140°**.  
            - **Hydrophobic Contacts:** C···C (nonpolar sidechain atoms) ≤ **5.0 Å** (exploratory upper range 5.5–6.0 Å).  
            - **London Dispersion (Heuristic):** vdW window lower bound ~3.3 Å; upper sampling cut ≤ **5.0 Å** – emphasis on cumulative clustering rather than single strength.  
            - **Sulfur–π:** S (CYS/MET) to ring centroid ≤ **6.0 Å** (strong <5.0 Å); perpendicular offset ≤ **3.0 Å**.  
            **Source:** Comprehensive review of unconventional noncovalent interactions in proteins.  
            [doi.org/10.1021/acsomega.3c00205](https://doi.org/10.1021/acsomega.3c00205)
            """
        )
