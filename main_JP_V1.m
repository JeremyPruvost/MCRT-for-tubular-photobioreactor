cd C:\Recherche\Modelisation\MCRT_Tube_Jack\CodesPerso\single_simulation
addpath C:\Recherche\Modelisation\MCRT_Tube_Jack\CodesPerso\CommonSubRoutines

% INPUTS = [Cx, w_pig, q_sol, rr, r_i, Ac]
clear all


%%%%%%% General parameter
calculate_rx=1;     %Boolean. Specify if growth kinetics (1) or if only light transfer (0) are calculated
NameFileKineticsParameters='ParametreCinetiqueChlorelle_POMoy_V4_Ac2800.mat';  %Give here the name of the mat file containing kinetics parameters

%%%%%%% Microalgae related parameters
Ea = 198.6;   % absorption cross section [kg/m2]
Es = 2873;   % scattering cross section [kg/m2]
%Es=0;

g = 0.974;  % microalgae cell asymmetry factor []

Cx=0.5; % Biomass concentration, DW [g/l]
Ac_ref = 2800; % compensation rate of photon absorption


%%%%%% Tube geometry related parameters
t_norm=0.5;
%t_norm_i=[0.05 0.2 0.3 0.35 0.4 0.45 0.5 0.6 0.8 0.9 0.92]; 
%rr = 0.9;           % ratio of inner to outer radius
r_i = 0.02;         % tube inner radius [m]

r_o=r_i/(1-t_norm);    %outer radius [m]
%r_o=0.02;    %outer radius [m]
%r_o = r_i/rr;       % tube outer radius [m]
t = r_o - r_i;      % tube thickness [m]
Cx=0.01/r_i;
%Cx=5;

%%%%%% Light illuminatation related parameters
q_sol = 200;  % PFD [µmol/m2/s]
theta_in = 45; % Incident angle of light with the tube centerline [°]

% theta_z - solar azimuth angle defined with respect to z-axis
% gamma_s - angle between projection of ray onto the x-y plane and the
% due south direction, east is positive, west is negative

% Manually define solar position
gamma_s = 0;    % degrees
theta_z=theta_in;    % degrees 

% Define tube orientation
% theta_t - angle between tube centerline and positive z-axis
% gamma_t - angle between projection of tube centerline onto the x-y plane
% and the positive x-axis, counter clockwise is positive
theta_t = 0;        % degrees
gamma_t = 0;        % degrees


%%%%%% Ray tracing related parameters
Nrays = 1e6;    %Number of rays
Grid = 100; % Size of the grid used in tube spatial discretization


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%Calculation section%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CxR = Cx*r_i; % Areal concentration [kg/m²]
rr = r_i/r_o; % radius ratio (r_inner/r_outer) []
t_tube=r_o-r_i;  %tube wall thickness [m]
t_norm=t_tube/r_o; %Normalized tube wall thickness []
albedo=Es/(Ea+Es); %Albedo []

%%%%%%%%%%%%% Ray tracing section

INPUTS=[Cx q_sol rr r_i Ac_ref Ea Es albedo Nrays theta_z Grid g calculate_rx gamma_s theta_t gamma_t];
[gamma,rx,rx_mean,LRPA,LRPA_gamma,q_PFD] = tube_ray_tracing_JP_fun(INPUTS,NameFileKineticsParameters);

%%%%%%%% General information
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
disp('Volumetric growth rate (rx) prediction as a function of grid refinement');disp(rx_mean)



%%%%%%%%%%%%% Post-treatment section (for the most refined grid -->
%%%%%%%%%%%%% ind_trace=3)
ind_trace=2;

LRPA_gamma=flipud(LRPA_gamma{ind_trace});
rx=flipud(rx{ind_trace});
LRPA=flipud(LRPA{ind_trace});
% LRPA_gamma=(LRPA_gamma{ind_trace});
% rx=(rx{ind_trace});
% LRPA=(LRPA{ind_trace});

r_trace=linspace(-r_i,r_i,Grid);

% [X,Y]=meshgrid(r_trace,r_trace);
% R=(X.^2+Y.^2).^0.5;
% ind_ext=find(R>r_i);
% rx(ind_ext)=1/0;
% LRPA_gamma(ind_ext)=1/0;
% LRPA(ind_ext)=1/0;

disp('Gamma value');disp(gamma(ind_trace));
disp('Volumetric growth rate (g/l/d)');disp(rx_mean(ind_trace));
disp('Areal growth rate (g/m²/d)');disp(1000*rx_mean(ind_trace)/alight);

verif_rx=mean(rx(isfinite(rx)))

%%%%%%%%%%%%% Graph section
DefineTurboColorMap;

figure(1)
clf
title('LRPA field')
surf(r_trace,r_trace,LRPA,'FaceColor','interp','EdgeColor','none')
view(0,90)
colormap(turboColorMap)
colorbar
circle(0,0,r_i);
circle(0,0,r_o);
axis([-r_o r_o -r_o r_o])
title('LRPA field')
grid off

figure(2)
clf
surf(r_trace,r_trace,LRPA_gamma,'FaceColor','interp','EdgeColor','none')
view(0,90)
colormap(turboColorMap)
%colorbar
hold on
circle(0,0,r_i);
circle(0,0,r_o);
axis([-r_o r_o -r_o r_o])
title('Map of culture locations where LRPA>Ac')
grid off

figure(3)
clf
surf(r_trace,r_trace,rx,'FaceColor','interp','EdgeColor','none')
view(0,90)
colormap(turboColorMap)
colorbar
hold on
circle(0,0,r_i);
[h,xunit,yunit]=circle(0,0,r_o);
axis([-r_o r_o -r_o r_o])
title('Local growth rate distribution')
grid off


%%%%%%%%%%%% Photons mass balance checking
MPRA__S_integral=q_PFD*alight/Cx;
ind_fin=isfinite(LRPA)|LRPA>0;
MRPA_V_integral=mean(LRPA(ind_fin));
disp('MRPA from surface integrale (µmole photons/kg/s)');disp(MPRA__S_integral)
disp('MRPA from volume integrale (µmole photons/kg/s');disp(MRPA_V_integral)

toc

%mean(rx(find(R<r_i)))



%%%%%%%%%%%% Photons mass balance checking

MPRA__S_integral=q_PFD*alight/Cx;
ind_fin=find(isfinite(LRPA)&LRPA>0);

MRPA_V_integral=mean(LRPA(ind_fin));
disp('MRPA from surface integrale (µmole photons/kg/s)');disp(MPRA__S_integral)
disp('MRPA from volume integrale (µmole photons/kg/s');disp(MRPA_V_integral)

r_trace=linspace(-r_i,r_i,1.5*Grid);
ind_nul=LRPA==0;
ind_int=isfinite(LRPA);
sum(ind_nul,'all')/sum(ind_int,'all');
disp("% Cells without ray absorption");disp(sum(ind_nul,'all')/sum(ind_int,'all'))




save SauveSimulRT_1_valid_theta



