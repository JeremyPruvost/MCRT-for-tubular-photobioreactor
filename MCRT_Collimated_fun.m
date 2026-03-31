function [gamma, LRPA, LRPA_gamma, q_PFD] = tube_ray_tracing_fun(INPUTS)

%% INPUTS
% INPUTS = [Cx, w_pig, q_sol, rr, r_i, Ac]

%%%%%| OPTICAL PROPERTIES | %%%%%
n1 = 1;             % air - refractive index 
k1 = 0;             % air - absorption index 
n2 = 1.5;           % tube - refractive index 
k2 = 0;             % tube - absorption index 
n3 = 1.333;         % medium - refractive index 
k3 = 0.000000;      % medium - absorption index 

%%%%%| MICROALGAE CELLS | %%%%%
% Calculate based on pigment content (for Chlorella vulgaris)
Cx = INPUTS(1);           % biomass concentration [g/L] 
g=INPUTS(12);

%%%%%| SOLAR | %%%%%
% Solar intensity
q_sol = INPUTS(2);        % PAR photon flux density [umol/m^{-2}s^{-1}]

% theta_z - solar azimuth angle defined with respect to z-axis
% gamma_s - angle between projection of ray onto the x-y plane and the
% due south direction, east is positive, west is negative

% Manually define solar position
theta_s = INPUTS(9);   % degrees
gamma_s = INPUTS(13);    % degrees

% Define solar direction
sun = [-sind(theta_s).*cosd(gamma_s),-sind(theta_s).*sind(gamma_s),-cosd(theta_s)];

%%%%%| TUBE | %%%%%
% Define inner radius and ratio
rr = INPUTS(3);           % ratio of inner to outer radius
r_i = INPUTS(4);         % tube inner radius [m]
r_o = r_i/rr;       % tube outer radius [m]
t = r_o - r_i;      % tube thickness [m]

% Define inner radius and tube thickenss
a_light = 1/r_i;    % specific illuminated area S_light/Volume [1/m]

% Define tube orientation
% theta_t - angle between tube centerline and positive z-axis
% gamma_t - angle between projection of tube centerline onto the x-y plane
% and the positive x-axis, counter clockwise is positive
theta_t = INPUTS(14);        % degrees
gamma_t = INPUTS(15);        % degrees

% Define tube centerline
tube = [-sind(theta_t)*cosd(gamma_t),-sind(theta_t)*sind(gamma_t),-cosd(theta_t)];
theta = acosd(dot(tube,sun,2));  % Angle between the sun and tube centerline


q_PFD = q_sol*cosd(90-theta); % use angle between sunlight and tube orthagonal vector


%%%%%| RAY TRACING | %%%%%
Grid=INPUTS(11);
Grid=[0.5*Grid Grid 1.5*Grid]; %create coarsen and more refined grid for grid influence checking

el_len = (2*r_i)./Grid;
el_vol = ((2*r_i)./Grid).^2; % element height of 1 m
A_illum = pi*r_o;   %%% Half tube perimeter (illuminated surface, considering external diameter)
V_tube = pi*r_i^2; % element height of 1 m

Ac = INPUTS(5);
Nrays = INPUTS(9);

%% START
 Ea_moy = INPUTS(6);
 Es_moy = INPUTS(7);

 
kappa = Cx*Ea_moy;
sigma = Cx*Es_moy;
beta = kappa + sigma;
albedo = sigma/beta;

E_photon = (q_PFD*A_illum)/Nrays;

%  0 - UNKNOWN
%  1 - air/tube
%  2 - tube/culture
%  3 - culture/tube
%  4 - tube/air
%  5 - no interface - within culture medium
%  matrix containing optical

n_infc = [n1+sqrt(-1)*k1 n2+sqrt(-1)*k2; 
          n2+sqrt(-1)*k2 n3+sqrt(-1)*k3; 
          n3+sqrt(-1)*k3 n2+sqrt(-1)*k2; 
          n2+sqrt(-1)*k2 n1+sqrt(-1)*k1;
          n3+sqrt(-1)*k3 n3+sqrt(-1)*k3];

%% Plot cylinders
% for loop added to make plot section easily collapsable
for i = 1
asdf = 1;

if asdf == 1
    [X_o,Y_o,Z_o] = cylinder(r_o,100);
    [X_i,Y_i,Z_i] = cylinder(r_i,100);
    Z_o = (Z_o-0.5)*5e-3;
    Z_i = (Z_i-0.5)*5e-3;

    rot_z = [cosd(gamma_t), -sind(gamma_t), 0 ;
        sind(gamma_t), cosd(gamma_t), 0 ;
        0, 0, 1];
    rot_y = [cosd(-theta_t), 0, sind(-theta_t) ;
        0, 1, 0 ;
        -sind(-theta_t), 0, cosd(-theta_t)];

    cyl_bot = [transpose(X_o(1,:)), transpose(Y_o(1,:)), transpose(Z_o(1,:))];
    cyl_top = [transpose(X_o(2,:)), transpose(Y_o(2,:)), transpose(Z_o(2,:))];

    cyl_bot = cyl_bot*rot_y;
    cyl_bot = cyl_bot*rot_z;
    cyl_top = cyl_top*rot_y;
    cyl_top = cyl_top*rot_z;

    X_o = [transpose(cyl_bot(:,1)) ; transpose(cyl_top(:,1)) ];
    Y_o = [transpose(cyl_bot(:,2)) ; transpose(cyl_top(:,2)) ];
    Z_o = [transpose(cyl_bot(:,3)) ; transpose(cyl_top(:,3)) ];

    cyl_bot = [transpose(X_i(1,:)), transpose(Y_i(1,:)), transpose(Z_i(1,:))];
    cyl_top = [transpose(X_i(2,:)), transpose(Y_i(2,:)), transpose(Z_i(2,:))];

    cyl_bot = cyl_bot*rot_y;
    cyl_bot = cyl_bot*rot_z;
    cyl_top = cyl_top*rot_y;
    cyl_top = cyl_top*rot_z;

    X_i = [transpose(cyl_bot(:,1)) ; transpose(cyl_top(:,1)) ];
    Y_i = [transpose(cyl_bot(:,2)) ; transpose(cyl_top(:,2)) ];
    Z_i = [transpose(cyl_bot(:,3)) ; transpose(cyl_top(:,3)) ];
end
end

%% Start ray tracing
% Gives vector orthagonal to plane defined by ray and tube
plane_int = cross(sun, tube)/norm(cross(sun,tube));

% Sets desired distance of emission "screen" away from tube
screen = 1.1*(r_o)/sind(theta);

% Defines emission points on plane_int with offset of 'screen'
sc1 = plane_int*r_o + (-sun)*screen;
sc2 = plane_int*(-r_o) + (-sun)*screen;

screen = sc1 - sc2;

t = linspace(0,1,Nrays);

% initialize
Nref = 0;
Nabs = 0; 
Ntube = 0;
N_refl_it = 0;
N_touch_tube = 0; 
count = 1;
counter = 1;
rho_count = 0;
abs_loc = [];

parfor i = 1:Nrays
    disp('Number of rays remaining to calculate');disp(i);
    fin = 0;
    interface = 1;
    dist = 0;
    ray1 = sun;
    Ninterface=0;
    while fin == 0
        % calculate intersection location and incidence angle with tube
        if interface == 1
            R = r_o;
            emit = sc2 + rand*screen; % random emission
            rand_beta = rand(1);
            l_beta = -log(rand_beta)/beta;
            [int, ang_i] = circ_intrsct_fun(R, tube, emit, ray1, theta_t, gamma_t,interface);
            
        elseif interface == 0
        

            [int_i, ang_i_i] = circ_intrsct_fun(r_i, tube, emit, ray1, theta_t, gamma_t,2);
            [int_o, ang_i_o] = circ_intrsct_fun(r_o, tube, emit, ray1, theta_t, gamma_t,4);
            dist_i = norm(int_i - int);
            dist_o = norm(int_o - int);

            if ang_i_o>pi/2
                ang_i_o=pi/2;
            end
            if ang_i_i>pi/2
                ang_i_i=pi/2;
            end
            
            if dist_i < dist_o
                int = int_i;
                ang_i = ang_i_i;
                interface = 2; % tube/culture
                N_touch_tube = N_touch_tube + 1; 
                
            else
                int = int_o;
                ang_i = ang_i_o;
                interface = 4; % tube/air           

            end             
        else
            [int, ang_i] = circ_intrsct_fun(R, tube, emit, ray1, theta_t, gamma_t,interface);            
        end
        ni = real(n_infc(interface,1));
        ki = imag(n_infc(interface,1));
        nt = real(n_infc(interface,2));
        kt = imag(n_infc(interface,2));

        
        if interface ~= 5
            if ang_i>pi/2
                ang_i=pi/2;
            end
            [rho,theta_2, Ni, Nt] = refl(ni,ki,nt,kt,ang_i);
        else
            rho = 0;
        end

        if rho > 1
            rho_count = rho_count + 1; 
        end
        
        rando = rand;
        if isnan(ang_i)
            interface = 1;
            rho = 1;
        end
        
        if rando <= rho % reflection occurs 
            switch interface                
                case 1 % air/tube
                    % NO POSSIBLE INTERSECTIONS WITH TUBES
                    Nref = Nref + 1;
                    fin = 1;
                case 2 % tube/culture
                    % MUST INTERSECT OUTER CYL
                    norml  = normal_vector(int,tube);
                    ray2 = refl_dir(norml,ray1);
                    R = r_o;
                    emit = int;
                    ray1 = ray2;
                    interface = 4;                    
                case 3 % culture/tube
                    % MUST INTERSECT INNER CYL
                    Ninterface=0;
                    N_refl_it = N_refl_it + 1;
                    interface = 3;
                    norml  = normal_vector(int,tube);
                    norml = -norml;
                    ray2 = refl_dir(norml,ray1);
                    emit = int;
                    ray1 = ray2;
                    R = r_i;
                    [int, ang_i] = circ_intrsct_fun(R, tube, emit, ray1, theta_t, gamma_t,interface);
                    inc = norm(emit-int);
                    if dist + inc > l_beta
                        inc_attn = l_beta - dist; % distance traveled current inc before attenuation
                        abs_or_sca = rand();
                        if abs_or_sca < albedo % ray is scattered
                            emit = emit + ray1*inc_attn;    % set attenuation location as new "emission" location
                            if (emit(1)^2 + emit(2)^2)^0.5 > r_i
                                emit;
                            end
                            [ray1(1),ray1(2),ray1(3)]=scatter_hg(g,ray1(1),ray1(2),ray1(3)); % use HG approx to determine new direction
                            l_beta = -log(rand())/beta; % generate new l_beta since attenuation took place
                            dist = 0;   % update distance traveled since last attenuation event
%                             ray_path = [ray_path; emit];
                            interface = 5;
                        else % ray is absorbed
                            Nabs = Nabs + 1;
                            location = emit + ray1*inc_attn;
                            if (location(1)^2 + location(2)^2)^0.5 < r_i
                                abs_loc(i,:) = location;
                            else
                                location;
                            end
%                             ray_path = [ray_path; abs_loc(i,:)];
                            fin = 1;
                        end

                    else % no attenuation, ray continues to culture/tube interface
                        interface = 3; % culture/tube
                        dist = dist + inc;
%                         ray_path = [ray_path; int];
                    end

                case 4 % tube/air
                    % CAN INTERSECT OUTER OR INNER CYL                     
                    norml  = normal_vector(int,tube);                    

                    ray2 = refl_dir(norml,ray1);
                    emit = int;
                    ray1 = ray2;
                    interface = 0;
                                    Ninterface=Ninterface+1;
                                    if Ninterface>100
                                        disp('Photon trapped in tube wall due to multiple reflections : destroyed')
                                    fin=1;
                                    end
            end
        elseif rando > rho % refraction occurs            
            switch interface
                case 1 % refraction at air/tube interface
                    % CAN INTERSECT OUTER OR INNER TUBE
                    norml  = normal_vector(int,tube);
                    norml = -norml;
                    ray2  = refr_dir(norml,ray1,Ni,Nt);
                    emit = int;
                    ray1 = ray2;
                    interface = 0;
                case 2 % refraction at tube/culture interface
                    % MUST INTERSECT INNER TUBE
                    interface = 3;
                    norml  = normal_vector(int,tube);
                    norml = -norml;
                    ray2  = refr_dir(norml,ray1,Ni,Nt);
                    emit = int;
                    ray1 = ray2;
                    R = r_i;
                    [int, ang_i] = circ_intrsct_fun(R, tube, emit, ray1, theta_t, gamma_t,interface);
                    inc = norm(emit-int);
                    if dist + inc > l_beta
                        inc_attn = l_beta - dist; % distance traveled current inc before attenuation
                        abs_or_sca = rand();
                        if abs_or_sca < albedo % ray is scattered
                            emit = emit + ray1*inc_attn;    % set attenuation location as new "emission" location
                            if (emit(1)^2 + emit(2)^2)^0.5 > r_i
                                emit;
                            end
                            [ray1(1),ray1(2),ray1(3)]=scatter_hg(g,ray1(1),ray1(2),ray1(3)); % use HG approx to determine new direction
                            l_beta = -log(rand())/beta; % generate new l_beta since attenuation took place
                            dist = 0;   % update distance traveled since last attenuation event
                            interface = 5; 
                        else % ray is absorbed
                            Nabs = Nabs + 1;
                            location = emit + ray1*inc_attn;
                            if (location(1)^2 + location(2)^2)^0.5 - r_i < 1e-15
                                abs_loc(i,:) = location;
                            else 
                                location;
                            end
                            fin = 1;
                        end

                    else % no attenuation, ray continues to culture/tube interface
                        interface = 3; % culture/tube
                        dist = dist + inc;
                    end

                case 3 % refraction at culture/tube interface
                    % MUST INTERSECT OUTER TUBE
                    norml  = normal_vector(int,tube);
                    ray2  = refr_dir(norml,ray1,Ni,Nt);
                    emit = int;
                    ray1 = ray2;
                    R = r_o;
                    interface = 4;
                case 4 % refraction at tube/air interface
                    % NO POSSIBLE INTERSECTIONS WITH TUBES
                    norml  = normal_vector(int,tube);
                    ray2  = refr_dir(norml,ray1,Ni,Nt);
                    emit = int;
                    ray1 = ray2;
                    fin = 1;
                    interface = 1;
                case 5 % propagating in culture
                    % MUST INTERSECT INNER TUBE
                    R = r_i;
                    inc = norm(emit-int);
                    if dist + inc > l_beta
                        inc_attn = l_beta - dist; % distance traveled current inc before attenuation
                        abs_or_sca = rand();
                        if abs_or_sca < albedo % ray is scattered
                            emit = emit + ray1*inc_attn;
                            if (emit(1)^2 + emit(2)^2)^0.5 > r_i
                                emit;
                            end
                            [ray1(1),ray1(2),ray1(3)]=scatter_hg(g,ray1(1),ray1(2),ray1(3)); % use HG approx to determine new direction
                            l_beta = -log(rand())/beta; % generate new l_beta since attenuation took place
                            dist = 0;   % update distance traveled since last attenuation event
                            interface = 5;
                        else % ray is absorbed
                            Nabs = Nabs + 1;
                            location = emit + ray1*inc_attn;
                            if (location(1)^2 + location(2)^2)^0.5 < r_i
                                abs_loc(i,:) = location;
                            else
                                location;
                            end
                            fin = 1;
                        end

                    else % no attenuation, ray continues to culture/tube interface
                        interface = 3; % culture/tube
                        dist = dist + inc;
                    end

            end % refraction case structure
        end % reflection refraction loop
    end % current ray
    
    interface = 1;
    dist = 0;
    fin = 0;
    count = 1;

%    Save Temp_Rays
end % current "bin"

rho = Nref/Nrays;
alpha = Nabs/Nrays; 

%min(abs_loc,[],'all')
disp('% of rays abosrbed');disp(alpha)
disp('% of rays reflected');disp(rho)

%%%% Validation that number of complex number due to problem in angle
%%%% calculation remains negligible

for i = 1:length(Grid)
ind_complex=find(imag(abs_loc(:,i))~=0); %%% identification of complex location (if any)
if length(ind_complex)/Nrays>0.01
    disp('WARNING : Too much incoherent calculations leading to excessive complex locations - due certainly to given angle calculations')
    pause
end
NumberofcomplexLocations(i)=length(ind_complex);
abs_loc(:,i)=real(abs_loc(:,i)); %%% Correction
end



clear Z G_loc LRPA
%%% set minimum acceptable Grid size, then calculate number of photons
%%% necessary to hit desired resolution (Ac/"res")
for i = 1:length(Grid)
    Z{i} = zeros(Grid(i),Grid(i));
    Z{i} = hist3([nonzeros(abs_loc(:,1)),nonzeros(abs_loc(:,2))],'Nbins',[Grid(i),Grid(i)]);
    G_loc{i} = (Z{i}.*E_photon)./(el_vol(i)*Cx*Ea_moy);
    LRPA{i} = (Z{i}.*E_photon)./(el_vol(i)*Cx);
end

gamma = 0;


for k = 1:length(Grid)
    r_mesh = size(LRPA{k},1)/2;
    gamma_num = 0;
    gamma_denom = 0;
    for i = 1:size(LRPA{k},1)
        for j = 1:size(Z{k},1)
            x = i-r_mesh;
            y = j-r_mesh;
            r = sqrt(x^2 + y^2);
            if r > r_mesh
                LRPA_gamma{k}(i,j) = -1;
            elseif LRPA{k}(i,j) > Ac 
                gamma_denom = gamma_denom + 1;
                gamma_num = gamma_num + 1;

                LRPA_gamma{k}(i,j) = 1;
            else
                gamma_denom = gamma_denom + 1;
            end
        end
    end

    gamma(k) = gamma_num/gamma_denom;


% Calculate values on the tube volume (finite values, orherwise set equal to -inf), excluding positions where no absorption occured (LRPA>0) 
r_trace=linspace(-r_i,r_i,Grid(k));
[X,Y]=meshgrid(r_trace,r_trace);
R=(X.^2+Y.^2).^0.5;
ind_ext=find(R>r_i);
LRPA_gamma{k}(ind_ext)=-1/0;
LRPA{k}(ind_ext)=-1/0;


end