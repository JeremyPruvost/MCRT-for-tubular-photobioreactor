function [norml] = normal_vector(ray, tube)
% gives OUTWARD pointing normal vector

proj = (dot(ray,tube)/norm(tube))*tube;
norml = ray-proj; 
norml = norml/norm(norml);
end