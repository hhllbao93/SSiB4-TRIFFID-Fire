# SSiB4/TRIFFID-Fire
The directory contains the code of SSiB4/TRIFFID-Fire (The Simplified Simple Biosphere Model coupled with the Top-down Representation of Interactive Foliage and Flora Including Dynamics Model-coupled with fire model from Li et al. (2012;2013))

The code will be used in the submission of the manuscript, "Modeling long-term fire impact on ecosystem characteristics and surface energy using the dynamic global vegetation model SSiB4/TRIFFID-Fire".

1. In v1.0 directory, runglob.sh and configure_ssib4.ys.IN are scripts for compilation. The ssib4_v1 directory includes the complete code for running the SSiB4/TRIFFID-Fire
2. In ssib4_src/ssib4_v1, the Fortran code related to the fire model includes: 
(1) ssib4_module_fireseason.F90: The code to calculate burned area
(2) ssib4_module_fireimp.F90: The code to calculate carbon loss due to biomass combustion and post-fire mortality
(3) ssib4_module_fireparms.F90: The parameter table for fire model
(4) ssib4_module_sitriffid.F90: To update fire impact in vegetation dynamics in TRIFFID
(5) ssib4_module_sfctrif.F90: Interface for fire model, TRIFFID, and SSiB4
