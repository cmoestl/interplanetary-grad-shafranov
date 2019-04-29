function [long,lat]=vectodeg(zs)

lat=asin(zs(3))*180/3.14159265;
long=acos(zs(1)/sqrt(zs(1)^2+zs(2)^2))*180/3.14159265;
if zs(2) < 0
	long=360-long;
end

