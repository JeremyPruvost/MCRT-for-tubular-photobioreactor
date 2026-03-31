function [new_ray] = refl_dir(norml,ray)
%refl - calculate direction of reflected ray 
proj=dot(ray,norml)*norml;
new_ray = ray - 2*proj; 
new_ray=new_ray/norm(new_ray);
end