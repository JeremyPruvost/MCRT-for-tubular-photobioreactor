function [refrray] = refr_dir(norml, ray, N1, N2)

theta_i = acos(dot(norml,ray)/(norm(ray)*norm(norml)));

full_reflec_angle=asin(N2/N1); %% Calculate limite angle for full reflection (part added by JP in 2025 to consider extreme angle cases)

if theta_i>full_reflec_angle 
refrray = refl_dir(norml,ray);
disp('Full reflexion occurs during refraction calculation due to too large incident angle'); 
else
theta_t = asin(N1.*sin(theta_i)./N2);

ray = ray/sqrt(dot(ray,ray));
norml = norml/sqrt(dot(norml,norml));
refrray=N1*ray/N2+(cos(theta_t)-N1*cos(theta_i)/N2)*norml;
refrray=refrray/sqrt(dot(refrray,refrray));

end

end