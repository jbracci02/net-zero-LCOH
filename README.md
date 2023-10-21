# net-zero-LCOH
This repository contains the primary code and all data inputs used to generate results in the scientific paper titled "Leveling the Playing Field: A Cost Comparison of Various Hourly-Reliable, Net-Zero Hydrogen Production Pathways"

This repository contains the following files:
1. ND_final_CA.jl: Techno-economic and emission assessment model used to generate results for this paper. This model is constructed in Julia.
2. NDwEB_final_CA.jl: Version of the H2 production cost model used to create the error bounds for Figure 4 in the main text.
3. lowEB_CA.jld2: lower LCOH error bound data for Figure 4
4. highEB_CA.jld2: higher LCOH error bound data for Figure 4
5. LCOHallEBmatrix_CA.jld2: Output LCOH data for each run of the Monte Carlo simulation performed in the TEA model found in NDwEB_final_CA.jl.
6. CF_Data_SAM.mat: MATLAB data file containing the hourly capacity factor data for solar PV in each location analyzed in this paper.
7. CA_grid_data_cambium.mat: MATLAB data file containing the hourly grid electricity price and emissions data pulled from NREL's Cambium dataset.
8. H2TEA_SI_final (nature_comm): Supporting Information to the main text. Contains model input data and additional figures not shown in the main text.
9. Source Data.xlsx: Contains all input data for the modeling performed in this study. Also contains raw results data for each plot shown in the main text.
