function [bxh, byh, bzh]=rotatey(bxh1, byh1, bzh1, a)

%a=input('Euler angle alpha (0-90) =');
b=0;
g=0;

a=a*pi/180;
b=b*pi/180;
g=g*pi/180;    

C=[cos(g)*cos(b)*cos(a)-sin(g)*sin(a)  cos(g)*cos(b)*sin(a)+sin(g)*cos(a) -cos(g)*sin(b); ...
  -sin(g)*cos(b)*cos(a)-cos(g)*sin(a) -sin(g)*cos(b)*sin(a)+cos(g)*cos(a)  sin(g)*sin(b); ...
   sin(b)*cos(a) sin(b)*sin(a) cos(b)];

%3D rotation

%rrz=[bxh1; byh1; bzh1];
rrz=[bzh1; bxh1; byh1];

nrrz=C*rrz;

bxh=nrrz(2,:);
byh=nrrz(3,:);
bzh=nrrz(1,:);

%cross(bxh,byh)-bzh
%xangle=acos(dot(bxh,bxh1))*180/pi
%yangle=acos(dot(byh,byh1))*180/pi
%zangle=acos(dot(bzh,bzh1))*180/pi
