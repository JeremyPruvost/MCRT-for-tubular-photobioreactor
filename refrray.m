% the fuction of refractive ray
function yy=refrray(ni,ki,nt,kt,incray,normal)
% incray is the vector of the incident ray
% normal is the normal of the interface pointing to transmitted medium
% N1 and N2 are effective refractive index and absorpton index

% effective refractive index
% Ni=eri1(ni,ki,incray,normal);
% 
% Nt=eri2(nt,kt,Ni,incray,normal);
Ni = ni;
Nt = nt;
% incident angle
thetai=acos(dot(incray,normal));
% transmitted angle
thetat=asin(Ni*sin(thetai)/Nt);
yy=Ni*incray/Nt+(cos(thetat)-Ni*cos(thetai)/Nt)*normal;
yy=yy/sqrt(dot(yy,yy));