# SSiB4/TRIFFID-Fire
The directory contains the code of SSiB4/TRIFFID-Fire (The Simplified Simple Biosphere Model coupled with the Top-down Representation of Interactive Foliage and Flora Including Dynamics Model-coupled with fire model from Li et al. (2012;2013))


The Fortran code related to the fire model includes: 
(1) ssib4_module_fireseason.F90: The code to calculate burned area
(2) ssib4_module_fireimp.F90: The code to calculate carbon loss due to biomass combustion and post-fire mortality
(3) ssib4_module_fireparms.F90: The parameter table for fire model
(4) ssib4_module_sitriffid.F90: To update fire impact in vegetation dynamics in TRIFFID
(5) ssib4_module_sfctrif.F90: Interface for fire model, TRIFFID, and SSiB4


Version 1.0: the fire modeling on global scale. It has ben used in the submission of the manuscript, "Modeling long-term fire impact on ecosystem characteristics and surface energy using the dynamic global vegetation model SSiB4/TRIFFID-Fire" (Huang et al. 2020)

Version 1.1: the fire model application in Southern Africa. Updates in the model includes: 
1) A constant agricultural fraction has been replaced by the annual-updated agricuture data. As they have a different crop fraction and spatial distribution in tropical regions, we have recalibrated the parameters for fire spread, fuel combustibility, and carbon combustion to reproduce the observed magnitude and temporal variations of burned area and carbon emission in satellite data.
2) The seasonality of vegetation productivity (GPP) is underestimated in Huang et al. (2020) as the model overestimated GPP in the dry season. We have adjusted paramters in the calculation of root-zone soil moisture potential factor f(θ) to reflect the effects of soil water deficit on transpiration. 
More details can be referred to Huang et al. (2021) 

Ref: 
1. Huang, H., Xue, Y., Li, F., and Liu, Y.: Modeling long-term fire impact on ecosystem characteristics and surface energy using a process-based vegetation–fire model SSiB4/TRIFFID-Fire v1.0, Geosci. Model Dev., 13, 6029-6050, 10.5194/gmd-13-6029-2020, 2020.
1. Huang, H., Xue, Y., Liu. Y., Li, F., and Okin, G.: Modeling the short-term fire effects on vegetation dynamics and surface energy in Southern Africa using the improved SSiB4/TRIFFID-Fire model, Submitted.
