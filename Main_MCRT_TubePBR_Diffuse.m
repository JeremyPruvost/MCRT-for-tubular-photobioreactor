%%%% Main code for MCRT for light transfer prediction in tubes - Case of
%%%% diffuse part of the sunlight
%%%% J.Hoengies, J.Pruvost (Nantes University, UCLA)
%%%% Version : 17 sept. 2025

%%%% The Matlab version must support the use of parfor for parallel
%%%% calculation (available in MatlabR2024b at least). If not the case, parfor will have to be replace by common
%%%% for loops.

%%%%%%% Microalgae related parameters
Ea = 198.6;     % Absorption cross section [kg/m2]
Es = 2873;      % Scattering cross section [kg/m2]

g = 0.974;      % Microalgae cell asymmetry factor []

Cx=0.5;         % Dry-weight biomass concentration, DW [g/l]
Ac_ref = 2800;  % Compensation point of photosynthesis (µmole/g/s)

%%%%%% Tube geometry related parameters
r_i = 0.03;             %   tube inner radius [m]
r_o=0.04;               %   tube outer radius [m]

%%%%%% Light illuminatation related parameters
q_sol = 200;            % PFD [µmol/m2/s]
theta_in = 0;           % Incident angle of light with the tube centerline [°]

%%%%%% Manually defined solar position
% theta_z - solar azimuth angle defined with respect to z-axis
% gamma_s - angle between projection of ray onto the x-y plane and the
% (due south direction, east is positive, west is negative)

gamma_s = 0;            % degrees
theta_z=90-theta_in;    % degrees 

%%%%%% Manually define tube orientation
% theta_t - angle between tube centerline and positive z-axis
% gamma_t - angle between projection of tube centerline onto the x-y plane
% and the positive x-axis, counter clockwise is positive

theta_t = 0;        % degrees
gamma_t = 0;        % degrees


%%%%%% Ray tracing related parameters
Nrays = 1e6;        %Number of rays
Grid = 100;         % Averaged size of the grid used in tube spatial discretization


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Calculation section%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CxR = Cx*r_i;       % Areal concentration [kg/m²]
rr = r_i/r_o;       % radius ratio (r_inner/r_outer) []
t_tube=r_o-r_i;     %tube wall thickness [m]
t_norm=t_tube/r_o;  %Normalized tube wall thickness []
albedo=Es/(Ea+Es);  %Albedo of microalgal culture []

%%%%%%%%%%%%% Ray tracing section
INPUTS=[Cx q_sol rr r_i Ac_ref Ea Es albedo Nrays theta_z Grid g gamma_s theta_t gamma_t];
[gamma,LRPA,LRPA_gamma,q_PFD] = MCRT_Diffuse_fun(INPUTS);

%%%%%%%% Display of general information before calculation
disp('PFDin (q''in) [µmol photons/m²/s]');disp(q_PFD)
disp('Areal concentration [kg/m²]');disp(CxR)
disp('Albedo []');disp(albedo)

disp('radius ratio (r_inner/r_outer) []');disp(rr)
disp('Tube wall thickness [m]');disp(t_tube)
disp('Normalized tube wall thickness []');disp(t_norm)
alight=1/r_o; %% Tube with one side illumination
disp('Specific illuminated surface [m²/m3]');disp(alight)




%%%%%%%%%%%%% Grid effect checking
disp('Grid refinement:');disp([0.5*Grid Grid 1.5*Grid])
disp('Light fraction (gamma) prediction as a function of grid refinement');disp(gamma)

%%%%%%%%%%%%% Analysis 
gamma_i=gamma(2);
LRPA_i=LRPA{2};
LRPA_gamma_i=LRPA_gamma{2};

%%%%%%%%%%%% Photons mass balance checking
MPRA__S_integral=q_PFD*alight/Cx;
ind_fin=find(isfinite(LRPA_i)&LRPA_i>0);

MRPA_V_integral=mean(LRPA_i(ind_fin));
disp('MRPA estimatation from surface balance (µmole photons/kg/s)');disp(MPRA__S_integral)
disp('MRPA from volume integral on LRPA values (µmole photons/kg/s');disp(MRPA_V_integral)

r_trace=linspace(-r_i,r_i,Grid);
ind_nul=LRPA_i==0;
ind_int=isfinite(LRPA_i);
sum(ind_nul,'all')/sum(ind_int,'all');
disp('% of grid cells without ray absorption');disp(sum(ind_nul,'all')/sum(ind_int,'all'))
disp('Light fraction (gamma) obtained from MRCT preduction');disp(gamma_i)




%%%%%%%%%%%%% Graphic section
figure(1)
clf
[X,Y]=meshgrid(r_trace,r_trace);
plot(X(ind_fin),Y(ind_fin),'.')
title('Grid cells where absorption occured')


DefineTurboColorMap;

figure(2)
clf
title('LRPA field')
surf(r_trace,r_trace,LRPA_i/1000,'FaceColor','interp','EdgeColor','none')
view(0,90)
colormap(turboColorMap)
colorbar
circle(0,0,r_i);
circle(0,0,r_o);
axis([-r_o r_o -r_o r_o])
title('LRPA field (in µmole/g/s)')
grid off



figure(3)
clf
surf(r_trace,r_trace,LRPA_gamma_i,'FaceColor','interp','EdgeColor','none')
view(0,90)
colormap(turboColorMap)
%colorbar
hold on
circle(0,0,r_i);
circle(0,0,r_o);
axis([-r_o r_o -r_o r_o])
title('Map of culture locations where LRPA>Ac')
grid off



figure(4)
clf
title('Fluence rate field')
surf(r_trace,r_trace,LRPA_i/Ea,'FaceColor','interp','EdgeColor','none')
view(0,90)
colormap(turboColorMap)
colorbar
circle(0,0,r_i);
circle(0,0,r_o);
axis([-r_o r_o -r_o r_o])
title('Fluence rate field (in µmole/m²/s)')
grid off













