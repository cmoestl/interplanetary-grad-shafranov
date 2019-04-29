function II=hu12n_opt(adjp, zs)
 % called by optcloud.m
 %Modified by C. Möstl Aug 2006
 
 %adjp=1.3; %adjust p (only for oct 19, 1984)
 %adjp=1.6;	% based on the assumption of the total pressure equalibrium
 %across a TD	-- IRM
 %adjp=1	%  0 for force-free and fit Bz; 1.5 for additional electron pressure
 %%Wind location in Earth Radii, GSE; 1 Re=6378 km
  
 %xs,ys,zs Basisvektoren des Invarianz-System in GSE BASIS
 %bxs,bys,bzs B-Vektoren in Invarianz-Basis
 %Vxs,Vys,Vzs plasma v-Vektoren in Invarianz-Basis
 %weiters selbsterklärend:
 %ave_beta=mean(beta)
 %ave_Vt=mean(Vt)
 %ave_VA=mean(VA)
 %mean_Vt_to_Va=ave_Vt/ave_VA
 %ZlatGSE
 %ZlongGSE
 %Zdelta_Hu_GSE wie in Hu 8/2005 und 3/2004
 %Zlamba_Hu_GSE
 % 2 figures als Output
 

 disp('--------------- ab hier hu12n_opt.m: -------------------')
  
 %*************** PARAMETERS **********************
 porder=3   % order of the polynomial -- ptol = exp ( poly(A) )
 wi=[0 0 0]' %%Wind location in Earth Radii, GSE; 1 Re=6378 km
 %***************************************************
 %Anzahl der Datenpunkte ist 'rn', es wird dA=-By dx usw. über alle
 %verfügbaren echten Datenpunkte integriert
 
 
 format short;
 format compact;
 
 nden=1e6;  
 nv=1e3; 
 nb=1e-9; 
 k=1.38*1e-23; 
 miu=4.0*pi*1e-7;
 ec=1.6*1e-19;

cd ..
f0=fopen('ALLCASE.par', 'rt');
%dirstr=fscanf(f0, '%s\n')
dirstr=deblank(fgets(f0));
fclose(f0);
cd(dirstr)

%if strcmp(deblank(pwd), 'C:\chris\Matlab\GS-code')
%cd('C:\chris\Matlab\GS-code')
%cd(dirstr);
%end

parfile=['CASE1.par'];

fpar=fopen(parfile, 'rt');
str=fscanf(fpar, '%c', 21);
while feof(fpar) ~= 1
i12=fscanf(fpar,'%i %i %i\n', 3); 
end
fclose(fpar);
cd ..
load(str);
cd(dirstr)
CASE1=clouddata;
i1=i12(1);
i2=i12(2);
% adjp=abs(i12(3));
%i1=input('i1 ?')
%i2=input('i2 ?')

rn=i2-i1+1;
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
pe=abs(adjp)*Ne.*Te*k*1e6; 
%pe=Ne.*Te*k*1e6; 
else
    pe=0;
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
pp=abs(adjp)*n.*t.*k*1e6; 	%temperature is in 10^6 K
%%%%%%%%%%%%%%%%
p=pe+pp;
%%%%%%%%%%%%%%%%

%%%%
Y=CASE1(i1:i2,3); M=CASE1(i1:i2,1); D=CASE1(i1:i2,2); hh=CASE1(i1:i2,4);
mmm=CASE1(i1:i2,5);
ssm=CASE1(i1:i2,6);
Datem=datenum(double(Y),double(M),double(D),double(hh),double(mmm),double(ssm));


%figure;
%subplot(211);
%plot(Datem,bx,Datem,by,Datem,bz);ylabel('B');legend('B_x','B_y','B_z');
%datetick('x',13);grid;
%subplot(212);
%plot(Datem,den,Datem,t);legend('density','Temp.');
%datetick('x',13);grid;
%%%%

load bxh.dat
load byh.dat
load bzh.dat
load dht.dat
xhg=bxh;
yhg=byh;% minimum invariance direction
zhg=bzh;
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



%hu2.m

%xhn=sqrt(xhg(1).*xhg(1)+xhg(2).*xhg(2)+xhg(3).*xhg(3));
%yhn=sqrt(yhg(1).*yhg(1)+yhg(2).*yhg(2)+yhg(3).*yhg(3));
%zhn=sqrt(zhg(1).*zhg(1)+zhg(2).*zhg(2)+zhg(3).*zhg(3));
%xh=xhg/xhn; yh=yhg/yhn; zh=zhg/zhn;
xh=xhg; yh=yhg; zh=zhg;

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

%ya=yh
%xa=xh
%za=zh

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

%zs=zsg;
%zs=bzh;















%%%%%%%%%%%%%% INVARIANCE AXIS IN DIFFERENT COORDSYS %%%%%%%%%%%%%%%%%%


ZlatGSE=asin(zs(3))*180/3.14159265
ZlongGSE=acos(zs(1)/sqrt(zs(1)^2+zs(2)^2))*180/3.14159265;
if zs(2) < 0
	ZlongGSE=360-ZlongGSE
else ZlongGSE
end


%****** convert to Hu 8/2005 **********
%delta measured from X-axis

Zdelta_Hu_GSE=asin(zs(1))*180/3.14159265;

if Zdelta_Hu_GSE < 0
        Zdelta_Hu_GSE=90+abs(Zdelta_Hu_GSE)
else Zdelta_Hu_GSE=90-Zdelta_Hu_GSE
end

%Lamda wie GSElatlong, nur Koordsys drehen 
%x->y, y->z, z-> x;
Zlambda_Hu_GSE=acos(zs(2)/sqrt(zs(2)^2+zs(3)^2))*180/3.14159265;
if zs(3) < 0
	Zlambda_Hu_GSE=360-Zlambda_Hu_GSE
else Zlambda_Hu_GSE
end

%****** AXIS FROM MVAB (MVUB) *****************************
load mvaby.dat;

LepZlatGSE=asin(mvaby(3))*180/3.14159265
LepZlongGSE=acos(mvaby(1)/sqrt(mvaby(1)^2+mvaby(2)^2))*180/3.14159265;
if mvaby(2) < 0
	LepZlongGSE=360-LepZlongGSE
else
 LepZlongGSE   
end

%****** again convert to Hu 8/2005 **********
%delta measured from X-axis
%lambda from Y towards Z-axis

LepZdelta_Hu_GSE=asin(mvaby(1))*180/3.14159265;

if LepZdelta_Hu_GSE < 0
        LepZdelta_Hu_GSE=90+abs(LepZdelta_Hu_GSE)
else LepZdelta_Hu_GSE=90-LepZdelta_Hu_GSE
end

%Lambda wie GSElatlong, nur Koordsys drehen 
%x->y, y->z, z-> x;
LepZlambda_Hu_GSE=acos(mvaby(2)/sqrt(mvaby(2)^2+mvaby(3)^2))*180/3.14159265;
if mvaby(3) < 0
	LepZlambda_Hu_GSE=360-LepZlambda_Hu_GSE
else LepZlambda_Hu_GSE
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%










vza=vht*zs'; 
vzav=vza*zs;
vhtsv=vht-vzav;
vhts=sqrt(vhtsv*vhtsv');
vhtsn=vhtsv/vhts;
xs=-vhtsn;		% -x = V_HT - dot(V_HT, z)z
%xs=+vhtsn;
%ys(1)=zs(2)*xs(3)-zs(3)*xs(2);
%ys(2)=zs(3)*xs(1)-zs(1)*xs(3);
%ys(3)=zs(1)*xs(2)-zs(2)*xs(1);

ys=cross(zs,xs);

%save xs.dat xs -ascii;
%save ys.dat ys -ascii;
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


%%%%%%%%%%%%%%%%%%%%%%%%%%% FITTING AND DISPLAYING %%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% (1) FIT ptol (transverse pressure) 
%%	try the exponential fit of ptol ~ A
%lnptol=log(ptol);
[P,S]=polyfit(A1,ptol,porder);
%[P,S]=polyfit(A1,lnptol,order);		%
residue0=S.normr;
pf=polyval(P, A1);
%pf=exp(polyval(P,A1));


%%%%%%%%%%%%%%%%% (2) FIT bz2 (magnetic transverse pressure) 
[P1,S1]=polyfit(A1,bz2,porder);
pf1=polyval(P1, A1);

%%%%%%%%%%%%%%%%% (3) FIT bz  

[P2,S2]=polyfit(A1,bzn,porder);
pf2=polyval(P2, A1);



%%%%%%%%%%%%%%%%%% calculate associated fitting residues %%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%% equation (4) in Hu 3/2004

N=length(A1);
Rfp=sqrt((1/N)*sum((ptol-pf).^2))/(max(ptol)-min(ptol));
Rfb=sqrt((1/N)*sum((bz2-pf1).^2))/(max(bz2)-min(bz2));
Rfbz=sqrt((1/N)*sum((bzn-pf2).^2))/(max(bzn)-min(bzn));
save Rf Rfp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

















%%%%%%%%%%%%%%%%%%%%%% 4 PLOTS untereinander
figure;
subplot(411);plot(Xa, A1,'--');xlabel('x');ylabel('A');grid on;
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


%%%%%%%%%%% PLOT  Pt, Bz, Bz2, p vs A1 %%%%%%%%%%%%%%%%%%%%

figure; subplot(221)
plot(A1,ptol,'--',A1,pf);
line(A1(1), ptol(1), 'marker','d','markerfacecolor','r');
line(A1(length(A1)), ptol(length(ptol)), 'marker','o','markersize',6, 'markerfacecolor','g');
xlabel('A (red=start)');ylabel('P_t');
%title(['\theta=',num2str(theta)]);

%%%%%%%% display B_z versus A and p vs A
%PBz=polyfit(A1,bzn,3); Pp=polyfit(A1,p,3);
%pf2=polyval(PBz, A1); pf3=polyval(Pp, A1);
%figure;
subplot(222);
plot(A1,bzn,'--',A1,pf2);xlabel('A');ylabel('B_z');
line(A1(1), bzn(1), 'marker','d','markerfacecolor','r');
line(A1(length(A1)), bzn(length(bzn)), 'marker','o','markersize',6, 'markerfacecolor','g');
subplot(223);
plot(A1,p,'-');xlabel('A');ylabel('p');
subplot(224);
plot(A1,bz2,'--',A1,pf1);xlabel('A');ylabel('B_z^2 /2\mu_0');
line(A1(1), bz2(1), 'marker','d','markerfacecolor','r');
line(A1(length(A1)), bz2(length(bz2)), 'marker','o','markersize',6, 'markerfacecolor','g');

%axis([-inf inf 0 1.2e-9]);

Vht_t=norm(vht-dot(vht,zs)*zs)

%%%%% Pfeile entlang Trajektorie
%figure;
%quiver([Xa Xa(3)], [zeros(size(Xa)) Xa(2)], [Vxs' Vht_t], [Vys' 0], .3, 'k-');
%hold on; quiver(Xa, zeros(size(Xa)), bxs, bys, .3, 'r-');axis equal;

ave_beta=mean(beta)
ave_Vt=mean(Vt)
ave_VA=mean(VA)
mean_Vt_to_Va=ave_Vt/ave_VA
%Mta2=mean(Vt./VA)


%PARAMETER ausgeben

disp('Basisvektoren des Invarianz-Systems in GSE:')
xs
ys
zs
ZlatGSE
ZlongGSE
Rfp
Rfb
Rfbz
disp('ratio of V_transversal to local Alfven speed:')
disp(mean_Vt_to_Va)



bb=sqrt(bxs.^2+bys.^2+bzs.^2);
bxs=bxs./bb;
bys=bys./bb;
bzs=bzs./bb;

vv=sqrt(Vxs.^2+Vys.^2+Vzs.^2);
Vxs=Vxs./vv;
Vys=Vys./vv;
Vzs=Vzs./vv;

disp('ratio of remaining V to B: X Y Z')
rx=mean(Vxs./bxs);
ry=mean(Vys./bys);
rz=mean(Vzs./bzs);
rvb=[rx ry rz]
rvbdev=[std(Vxs./bxs) std(Vys./bys) std(Vzs./bzs)]

% %%%%%%%%%%%%%%%%%%%%%%%%%%% Basis visualisieren %%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %%Wind location in Earth Radii, GSE; 1 Re=6378 km
% dht=vht/vhtn
% 
% %x=zeros(2,1);
% %y=zeros(2,1);
% %z=zeros(2,1);
% 
% x=[wi(1) wi(1)]';
% y=[wi(2) wi(2)]';
% z=[wi(3) wi(3)]';
% 
% u=[xs(1) ys(1)]';
% v=[xs(2) ys(2)]';
% w=[xs(3) ys(3)]';
% 
% figure;
% h1=quiver3(x,y,z,u,v,w);
% set(h1,'LineWidth',1.8);    
% hold on;
% h2=quiver3(wi(1),wi(2),wi(3),dht(1),dht(2),dht(3));
% set(h2,'LineWidth',1.8);   
% h3=quiver3(wi(1)-zs(1),wi(2)-zs(2),wi(3)-zs(3),2*zs(1),2*zs(2),2*zs(3));
% set(h3,'LineWidth',1.8);       
% 
% h4=quiver3(wi(1),wi(2),wi(3),mvaby(1),mvaby(2),mvaby(3));
% set(h4,'LineWidth',1.8);       
% 
% xlabel('X (GSE)');ylabel('Y (GSE)');zlabel('Z (GSE)');
% legend([h4 h3 h2 h1],'MVUB y-Achse','Invarianz-Achse z','deHoffmann-Teller velocity','x,y')
% axis([wi(1)-1 wi(1)+1 wi(2)-1 wi(2)+1 wi(3)-1 wi(3)+1])
% 
% %axis equal;
% axis vis3d;
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd ..

%%%%%%%%%%%%%%%%%%%%%%%%
%II=input('Is the optimal invariance direction found? y/n [n]:', 's');
%if isempty(II)
%	II = 'n';
%end
%if II=='y'
           hgsave(2, str(1:(length(str)-4)));
  cd(dirstr)         
    parid=fopen('CASE1_good.par', 'at');
    fprintf(parid, '%s\n', str);
    fprintf(parid, '%s\n', [int2str(i12(1:2)') '  ' int2str(adjp) '  %optcloud good']);
    fprintf(parid, '%f %f %f\n', zs');
    fclose(parid);
    
%       parid=fopen('CASE1.par', 'at');
%     fprintf(parid, '%i %i %3d\n', [i12(1:2); 0]);
%     fclose(parid);
%     %cd('..');
II='y';
 disp('--------------- Ende hu12n_opt.m: -------------------')
 
 save 'hu12n_opt_var.mat' ;
 
end