# MCRT-for-tubular-photobioreactor
Matlab code for light transfer modeling in tubular photo(bio)reactors using Monte Carlo Ray Tracing

%%%% J.Hoengies, J.Pruvost (Nantes University, UCLA)
%%%% Version : 23 March 2026

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

For more details, see 

J. Hoeniges, L. Pilon, J. Dauchet, J. Pruvost, Light transfer investigation in solar tubular photobioreactors using Monte Carlo ray tracing, Chemical Engineering Journal Advances, Volume 26,
2026. https://doi.org/10.1016/j.ceja.2026.101081.
Abstract: Tubular photobioreactors (PBRs) are widely used for microalgae solar culture. In such geometry, light transfer is rendered complex by the effects of refractive index mismatches across curved air/glass/culture interfaces. This study aims to develop an open source code implementing the Monte Carlo ray-tracing (MCRT) method and to identify the design rules to optimize light transfer in tubular PBRs. The light transfer model accounted for the tube wall thickness, the angle of the collimated and diffuse incident solar radiation, the multiple reflections and refractions at interfaces, and anisotropic scattering and absorption by the microalgae cells. The importance of the dimensionless optical thickness of the suspension and of the tube thickness normalized by the outer tube radius in determining the radiation field was demonstrated. The concave tube wall was responsible for light concentration and a normalized tube thickness in the range 0.2–0.4 was found to be optimum to increase photon flux in the culture by up to 42 % and the biomass productivity by 60 % compared to very thin tube thickness. Furthermore, due to the effect of curved optical surfaces, the conditions for maximum growth rate were different from those for full light absorption without dark volume, by contrast with flat panel PBRs. Predictions for several PBR configurations confirmed the usefulness of our MCRT method as a generic tool to predict performances of solar tubular PBRs. For example, an increase of 60 % in maximum biomass concentration was obtained with vertical inclination of PBR tubes compared to horizontal inclination, but with an increase in light stressing conditions.
Keywords: Modeling; Photobioreactors; Tubular; Microalgae; Light transfer; Ray-tracing

