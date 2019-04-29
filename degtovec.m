function zs=degtovec(long, lat)


degtorad=3.14159265/180;
theta=lat*degtorad;
phi=long*degtorad;
x=cos(theta)*cos(phi);
y=cos(theta)*sin(phi);
z=sin(theta);
zs=[x y z];

