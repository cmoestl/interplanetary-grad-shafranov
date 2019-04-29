%hu12n.m
%fitting the Pt(A) function and numerical integration of the GS equation
%calls hu34n.m


format short;
format compact;

 %adjp=1	%  0 for force-free and fit Bz; 1.5 for additional electron
 %pressure
 %order=1;	% order of the polynomial -- ptol = exp ( poly(A) )
 %adjp taken from CASE1.par

 nden=1e6;  
 nv=1e3; 
 nb=1e-9; 
 k=1.38*1e-23; 
 miu=4.0*pi*1e-7;
 ec=1.6*1e-19;

f0=fopen('ALLCASE.par', 'rt');
%dirstr=fscanf(f0, '%s\n')
dirstr=fscanf(f0, '%c', 14);
ctrl=fscanf(f0, '%i', 6);
porder=ctrl(2); nx=ctrl(3); ny=ctrl(4); dmid=ctrl(5); get_Ab=ctrl(6);
ctrl2=fscanf(f0, '%f', 2);
dAl0=ctrl2(1); dAr0=ctrl2(2);

fclose(f0);


%if strcmp(deblank(pwd), 'C:\Qiang\matlab\2001-2002')
cd(dirstr);

disp('--------------- START hu12n.m ------------------')
diary('output_gs_solver.txt') % in Output stehen dann alle Paramete

disp('ALLCASE.par')
dirstr
ctrl
ctrl2


parfile=['CASE1.par'];

fpar=fopen(parfile, 'rt');
str=fscanf(fpar, '%c', 21);
while feof(fpar) ~= 1
i12=fscanf(fpar,'%i %i %i\n', 3); 
end
fclose(fpar);

disp('CASE1.par')
str
i12

cd ..
load(str);
cd(dirstr)

CASE1=clouddata;
CASE1=double(CASE1);
i1=i12(1);
i2=i12(2);
adjp=abs(i12(3));
%i1=input('i1 ?')
%i2=input('i2 ?')

rn=i2-i1+1;  %Anzahl der Datenpunkte
%nx=21; %21;

%load(str);
%CASE1=eval(str(9:(length(str)-4)));

%[ltime,lT]=issireform(str);
Dtime=CASE1(:,7);
ltime=cumsum(Dtime);

time=ltime(i1:i2)-ltime(i1);
bx=CASE1(i1:i2,8);
by=CASE1(i1:i2,9);
bz=CASE1(i1:i2,10);

vx=CASE1(i1:i2,13);
vy=CASE1(i1:i2,14);
vz=CASE1(i1:i2,15);
vv=sqrt(vx.^2+vy.^2+vz.^2);

if size(CASE1,2)>15
 Ne=CASE1(i1:i2, 16)*nden;
 Te=CASE1(i1:i2, 17);
 pe=adjp*Ne.*Te*k*1e6; 
else
    pe=0
end

den=CASE1(i1:i2,11);
%T=CASE1(i1:i2,6);
%t=lT(i1:i2);
t=CASE1(i1:i2,12); %

dt(1)=0;
for i=2:rn
	dt(i)=time(i)-time(i-1);
end

n=den.*nden;
pp=adjp*n.*t.*k*1e6; 	%temperature is in 10^6 K
%%%%%%%%%%%%%%%%
p=pe+pp;
%%%%%%%%%%%%%%%%

%%%%
Y=CASE1(i1:i2,3); M=CASE1(i1:i2,1); D=CASE1(i1:i2,2); hh=CASE1(i1:i2,4);
mmm=CASE1(i1:i2,5); ssm=CASE1(i1:i2,6);
Datem=datenum(Y,M,D,hh,mmm,ssm);


%figure;
%subplot(211);
%plot(Datem,bx,Datem,by,Datem,bz);ylabel('B');legend('B_x','B_y','B_z');
%datetick('x',13);grid;
%subplot(212);
%plot(Datem,den,Datem,t);legend('density','Temp.');
%datetick('x',13);grid;
%%%%

% load bxh.dat
% load byh.dat
% load bzh.dat
load dht.dat    %deHoffmann-Teller velocity
% xhg=bxh;
% yhg=byh;% minimum invariance direction
% zhg=bzh;
vhtg=dht;% deHoffmann-Teller velocity

%figure;
%subplot(511)
%plot(bx);ylabel('B_x');axis([0 250 -inf inf]);grid;title('Case 1');
%subplot(512)
%plot(by);ylabel('B_y');axis([0 250 -inf inf]);grid;
%subplot(513)
%plot(bz);ylabel('B_z');axis([0 250 -inf inf]);grid;
%subplot(514)
%plot(n);ylabel('N');axis([0 250 -inf inf]);grid;
%subplot(515)
%plot(t);ylabel('T');axis([0 250 -inf inf]);grid;
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%hu2.m

%xhn=sqrt(xhg(1).*xhg(1)+xhg(2).*xhg(2)+xhg(3).*xhg(3));
%yhn=sqrt(yhg(1).*yhg(1)+yhg(2).*yhg(2)+yhg(3).*yhg(3));
%zhn=sqrt(zhg(1).*zhg(1)+zhg(2).*zhg(2)+zhg(3).*zhg(3));
%xh=xhg/xhn; yh=yhg/yhn; zh=zhg/zhn;

vht=vhtg.*nv;
vhtn=sqrt(vht(1).*vht(1)+vht(2).*vht(2)+vht(3).*vht(3));


% find the invariant axis on (xh,zh) plane
% counter-clockwise rotation by theta

% theta=input('theta=');
% 
%     delta=theta*pi/180;
%     C=[cos(delta) sin(delta); -sin(delta) cos(delta)];
% 
%     ya=yh;
%     XZ=[xh;zh];
%     XZ=C*XZ;
 %    xa=XZ(1,:);
%     za=XZ(2,:); %invariant axis

%ya=yh;
%xa=xh;
%za=zh;

%%%xa=[.8700 .4266 .2471];	%10/18/95
%%%ya=[.00971 -.5159 .8566];
%%%za=[.4929 -.7428 -.4530];
%
%xa=[0.88479237957879   0.42384433015091   0.19364511054565];%10/18/95
%ya=[0.01072130535050  -0.43394666284666   0.90082162606980];
%za=[0.46588576481140  -0.79494195744691  -0.38848680249907];

%%%%%%%%%%%%%%%%% satellite coord.
 
%zs=za;

%alpha1=input('Rotate around X_2 by (degrees):');
%alpha2=input('Rotate around X_3 by (degrees):');

%[xs, ys, zs]=rotatez(xa, ya, zs, alpha1); %rotate around z by 0-360 deg
%[xs, ys, zs]=rotatex(xs, ys, zs, alpha2); %rotate around x by 0-90 deg

load zs.dat; 
%zs=zsg;
%zs=bzh;

ZlatGSE=asin(zs(3))*180/3.14159265
ZlongGSE=acos(zs(1)/sqrt(zs(1)^2+zs(2)^2))*180/3.14159265
if zs(2) < 0
	ZlongGSE=360-ZlongGSE
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vza=vht*zs'; 
vzav=vza*zs;
vhtsv=vht-vzav;
vhts=sqrt(vhtsv*vhtsv');
vhtsn=vhtsv/vhts;
%xs=-vhtsn;		% -x = V_HT - dot(V_HT, z)z
%xs=-vhtsn;  
xs=-vhtsn;  %
%xs=+vhtsn;
%ys(1)=zs(2)*xs(3)-zs(3)*xs(2);
%ys(2)=zs(3)*xs(1)-zs(1)*xs(3);
%ys(3)=zs(1)*xs(2)-zs(2)*xs(1);

ys=cross(zs,xs);

save xs.dat xs -ascii;
save ys.dat ys -ascii;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% alternative approach to determine xs ys zs
%%%%% by setting VHT_1 \dot xs = VHT_2 \dot xs
%%%%% load xs.dat
%%%%%if dot(xs,byh) > 0
%%%%% zs=cross(xs,byh);
%%%%%else
%%%%% zs=-cross(xs,byh);
%%%%%end
%%%%% zs=zs/norm(zs);
%%%%% ys=cross(zs,xs);

%[xs, ys, zs]=rotate(xs, ys, zs);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ACEpos=[243 8.0 24];      % ACE position in GSE (Re)

%ACExyz=ACEpos*[xs' ys' zs'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 vxs=-vht*xs';
%  vxs=+vht*xs';

dxa=vxs*dt;%spatial interval

bgse=[bx by bz]; V=[vx*nv-vht(1) vy*nv-vht(2) vz*nv-vht(3)];

b=bgse.*nb;
bb=sqrt(b(:,1).^2+b(:,2).^2+b(:,3).^2);
bmax=max(bb);

% project b v into the satellite coord
bxs=b*xs'; Vxs=V*xs';
bys=b*ys'; Vys=V*ys';
bzs=b*zs'; Vzs=V*zs'; 

Vt=sqrt(Vxs.^2+Vys.^2);       % m/s
VA=bb./sqrt(miu*n*1.67e-27);  % m/s

% integrate A along the satellite trajectory


dxa(1)=0;
 dA(1)=0;
 for i=2:(rn)

    dA(i)=-(bys(i)+bys(i-1)).*dxa(i).*0.5;
 end
 A1=[cumsum(dA)];

 A1=A1';


    bzn=bzs;
    bz2=(bzn.*bzn)./(2*miu);
    ptol=bz2+p;
	PTol=bb.^2/(2*miu)+p;
	PB=bb.^2/(2*miu);
	beta=p./PB;

Xa=cumsum(dxa');
Xa=Xa'; Ya=zeros(rn,1);
xL=Xa(rn)-Xa(1);

%[A1,In]=sort(A1);
%ptol=ptol(In);

%%%%%%%% display A ~ x, ptol ~ x, ptol versus A
%	try the exponential fit of ptol ~ A

%porder=3;       %6 
lnptol=log(ptol);
[P,S]=polyfit(A1,ptol,porder);
%[P,S]=polyfit(A1,lnptol,order);		%
residue0=S.normr

pf=polyval(P, A1);
%pf=exp(polyval(P,A1));

figure;
subplot(411);plot(Xa, A1,'--');xlabel('x');ylabel('A');grid on;
title('Grössen im Invarianz-System')
%title(['theta=',num2str(theta)]);
subplot(412);plot(Xa, bxs,'-',Xa, bys,'--', Xa, bzs, '-.');
xlabel('x');grid on;legend('B_x', 'B_y', 'B_z');
subplot(414);plot(Xa, beta,'-');xlabel('x');ylabel('\beta');grid on;
subplot(413);plot(Xa, PTol,'-', Xa, PB, '--', Xa, p, '-.', Xa, ptol, ':');
xlabel('x');grid on;
legend('P_T', 'P_B', 'p', 'P_t');


%figure;
%subplot(211);plot(Xa, A1);xlabel('x');ylabel('A');
%title(['theta=',num2str(theta)]);
%subplot(212);plot(Xa, ptol);xlabel('x');ylabel('P_T');
figure; subplot(221)
plot(A1,ptol,'--',A1,pf);
line(A1(1), ptol(1), 'marker','d','markerfacecolor','r');
line(A1(length(A1)), ptol(length(ptol)), 'marker','x','markersize',12);
xlabel('A');ylabel('P_t');
%title(['\theta=',num2str(theta)]);

%%%%%%%% display B_z versus A and p vs A
%PBz=polyfit(A1,bzn,3); Pp=polyfit(A1,p,3);
%pf2=polyval(PBz, A1); pf3=polyval(Pp, A1);
%figure;
subplot(222);
plot(A1,bzn,'-');xlabel('A');ylabel('B_z');
subplot(223);
plot(A1,p,'-');xlabel('A');ylabel('p');
subplot(224);
plot(A1,bz2,'-');xlabel('A');ylabel('B_z^2 /2\mu_0');
%axis([-inf inf 0 1.2e-9]);


%figure;
%plot(A1,p,'-');xlabel('A');ylabel('p');


Vht_t=norm(vht-dot(vht,zs)*zs);
%figure;
%quiver([Xa Xa(3)], [zeros(size(Xa)) Xa(2)], [Vxs' Vht_t], [Vys' 0], .3, 'k-');
%hold on; quiver(Xa, zeros(size(Xa)), bxs, bys, .3, 'r-');axis equal;
ave_beta=mean(beta)
ave_Vt=mean(Vt)
ave_VA=mean(VA)
Mta=ave_Vt/ave_VA
Mta2=mean(Vt./VA)

cd ..;

pwd

disp('--------------- ENDE hu12n.m ------------------')

  

hu34n;

%hu34n_multi; %multi spacecraft analysis

diary off;



