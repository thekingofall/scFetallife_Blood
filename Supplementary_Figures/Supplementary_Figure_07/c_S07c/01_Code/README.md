# Supplementary Figure 7c: NK-cell temporal transitions

`01_compute_nk_temporal_transitions.py` reads `00_Global_Data/Shared_Inputs/Supplementary_Figures/Supplementary_Figure_07/01_NK_Trajectories/Fetal_NK_CellRank_Input.h5ad`. It solves the MOSCOT temporal problem using the object's post-conception age field (`day`), with epsilon 0.001, tau_a 0.95 and mean cost scaling. CellRank combines the transport and connectivity matrices using a connectivity weight of 0.2.

The workflow generates panels b–d together: CXCR6+ NK, CD56highCD16low NK and CX3CR1+ NK, respectively. Each panel receives PDF and PNG files in its own `02_Figures` folder. The shared transition matrix and cell metadata are saved in `b_S07b/03_PlotData`.

The pcw axis is labelled at 11.6, 21.7, 26.9, 32.4 and 37.9 weeks. Tick positions follow the continuous age scale.
