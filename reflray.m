function yy=reflray(incray,normal)
% incray is the vector of the incident ray
% normal is the normal of the interface pointing to transmitted medium
yy=incray-2*dot(incray,normal)*normal;
yy=yy/sqrt(dot(yy,yy));
