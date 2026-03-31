function [int, ang_i] = circ_intrsct_fun(R, tube, p1, ray, theta_t, gamma_t, interface, same_circ)
% INPUTS
% R - tube radius
% tube - vector defining tube orientation
% p1 - emission location
% ray - ray direction vector (does not need to be unit vector)
% OUTPUTS
% int - point of intersection with tube
% ang_i - [rad] angle of incidence between ray and tube surface

% Calculations for determine ray/tube intersection are based on the
% equations found on this webpage:
% https://mathworld.wolfram.com/Circle-LineIntersection.html

% theta_t = 90;        % degrees
% gamma_t = 0;         % degrees

% Rotate everything such that tube vector is positive x-axis
% Normalize ray/circle to consider unit circle
ray = ray./sqrt(ray(1)^2 + ray(2)^2 + ray(3)^2);

if abs((p1(1)^2 + p1(2)^2)^0.5 - R) < 1e-5
    same_circ = 1;
else
    same_circ = 0;
end 

if (p1(1)^2 + p1(2)^2)^0.5 > R &&  interface == 5
    p1;
end 

p1 = p1./R;
r = 1;

rot_z = [cosd(gamma_t), -sind(gamma_t), 0 ;
    sind(gamma_t), cosd(gamma_t), 0 ;
    0, 0, 1];
rot_y = [cosd(theta_t), 0, sind(theta_t) ;
    0, 1, 0 ;
    -sind(theta_t), 0, cosd(theta_t)];

p1 = p1*rot_y;
p1 = p1*rot_z;

ray = ray*rot_y;
ray = ray*rot_z;

% % Projection of "ray" onto "tube"
% proj = (dot(ray,tube)*tube)/norm(tube)^2;
% % Projection of "ray" onto plane defined by normal vector "tube"
% ray_proj = ray - proj;

ray_proj = [ray(1),ray(2),0];

% Calculate a second point along ray direction
p2 = p1 + ray_proj.*1;

% Calculation to determine circle/ray intersection
x = zeros(2,1);
y = x;
check = x;
dx = p2(1) - p1(1);
dy = p2(2) - p1(2);
dr = sqrt(dx^2 + dy^2);
D = p1(1)*p2(2) - (p2(1)*p1(2));

delta = r^2*dr^2-D^2;

if delta <= 0
    int = [NaN, NaN, NaN];
    t = [NaN, NaN];
else
    if sign(dy) == 0
        x(1,1) = (D*dy+dx*sqrt(delta))/dr^2;
        x(2,1) = (D*dy-dx*sqrt(delta))/dr^2;
    else
        x(1,1) = (D*dy+sign(dy)*dx*sqrt(delta))/dr^2;
        x(2,1) = (D*dy-sign(dy)*dx*sqrt(delta))/dr^2;
    end
    y(1,1) = (-D*dx+abs(dy)*sqrt(delta))/dr^2;
    y(2,1) = (-D*dx-abs(dy)*(delta^0.5))/dr^2;

    ty = (y-p1(2))./ray(2);
    tx = (x-p1(1))./ray(1);
    
    diff_y = zeros(2,2);
    diff_x = diff_y;

    diff_y(:,1) = x - (p1(:,1) + ray(1)*ty);
    diff_y(:,2) = y - (p1(:,2) + ray(2)*ty);
    diff_x(:,1) = x - (p1(:,1) + ray(1)*tx);
    diff_x(:,2) = y - (p1(:,2) + ray(2)*tx);


    for i = 1:size(diff_x,1)*size(diff_x,2)
        if isnan(diff_x(i))
            diff_x(i) = Inf;
        end
    end

    for i = 1:size(diff_y,1)*size(diff_y,2)
        if isnan(diff_y(i))
            diff_y(i) = Inf;
        end
    end


    max_diff_x = max(abs(diff_x),[],'all');
    max_diff_y = max(abs(diff_y),[],'all');

    
    if max_diff_x < max_diff_y
        t = tx;
    else
        t = ty;
    end



    z = p1(3) + ray(3)*t;
end


if t(1) <= 0
    t(1) = 0;
end
if t(2) <= 0
    t(2) = 0;
end


p1 = [p1; p1];

int_x = p1(:,1) + ray(1)*t;
int_y = p1(:,2) + ray(2)*t;
int_z = p1(:,3) + ray(3)*t;
int = [int_x, int_y, int_z];

check = sqrt((p1(:,1)-int_x).^2 + (p1(:,2)-int_y).^2 + (p1(:,3)-int_z).^2);

if same_circ == 1
    if check(1) >= check(2)
        int = int(1,:);
    elseif check(2) > check(1)
        int = int(2,:);
    end 
elseif same_circ == 0
    if check(1) < 1e-7
        int = int(2,:);
    elseif check(2) <= 1e-7
        int = int(1,:);
    elseif check(1) > check(2)
        int = int(2,:);
    elseif check(1) <= check(2)
        int = int(1,:);
    end
end

if interface == 1 || interface == 2
    n_cyl = [-int(1), -int(2), 0]; % INWARD pointing normal vector
elseif interface == 3 || interface == 4 || interface == 5
    n_cyl = [int(1), int(2), 0]; % OUTWARD pointing normal vector
end

ang_i = acos(dot(n_cyl,ray)/(norm(ray)*norm(n_cyl)));

if isnan(ang_i)
    ang_i;
    int = [NaN, NaN, NaN];
elseif ang_i > pi/2
    ang_i;
end

if size(int,1) == 2
    check; 
end 


int = int.*R;

rot_z = [cosd(-gamma_t), -sind(-gamma_t), 0 ;
    sind(-gamma_t), cosd(-gamma_t), 0 ;
    0, 0, 1];
rot_y = [cosd(-theta_t), 0, sind(-theta_t) ;
    0, 1, 0 ;
    -sind(-theta_t), 0, cosd(-theta_t)];

% second check
if abs((int(1)^2 + int(2)^2)^0.5)- R > 1e-10
    int;
end



int = int*rot_z;
int = int*rot_y;