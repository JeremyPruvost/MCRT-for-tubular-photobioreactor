# MCRT-for-tubular-photobioreactor
Matlab code for light transfer modeling in tubular photo(bio)reactors using Monte Carlo Ray Tracing

IMPORTANT NOTE (17Sept2025): The Full code will be made available once the related publication will be accepted (under review)

%%%% J.Hoengies, J.Pruvost (Nantes University, UCLA)
%%%% Version : 17 sept. 2025

%%%% The Matlab version must support the use of parfor for parallel
%%%% calculation. If not the case, parfor will have to be replace by common
%%%% for loops.


This code enables calculation of the light radiation field inside tubular photo(bio)reactors geometries. It has set for solar conditions, considering possibilities of non-normal incidence, various inclination and orientation toward the sun, and direct and diffuse light.

Main codes are: 
•	Main_MCRT_TubePBR_Collimated.m : main code for MCRT for light transfer prediction in tubes - Case of collimated (direct) part of the sunlight
•	Main_MCRT_TubePBR_Diffuse.m: main code for MCRT for light transfer prediction in tubes - Case of diffuse part of the sunlight
Other codes are subroutines. 


As outputs, the code gives : 
-	Local (LRPA) and mean (MRPA) values of the Rate of Photons Absorption 
-	Fluence rate field (G)
-	Light fraction (gamma) as defined by the value of compensation point of photosynthesis (Ac)
Those values allows further calculation such as:
-	the LRPA or MRPA determination corresponding to a sunlight composed of both direct and diffuse light (LRPA and MRPA obtained by each code have just to be added), 
-	Prediction of resulting photosynthetic growth rate (an appropriate kinetic growth model is however requested), 
-	extension to locations/periods of the year (by using appropriate values of sunlight received using for example solar database). 
Note that calculations are made for 3 grids refinements (the mean value being defined in Grid), to check dependency of the prediction to the number of rays (Nrays) used in MCRT method.

For more details, see ????
