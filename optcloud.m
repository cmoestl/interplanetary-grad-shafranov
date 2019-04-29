%function optcloud(adjp)

%optcloud.m is a function that does the optimal search for flux rope axis. 
%(here, by default optcloud is equivalent to optcloud(0)  which optimizes
%the single-variable function
%$B_z^2/2 \mu_0$. One can also specify optcloud(1) and optcloud(-1) which utilize $P_t(A)$ and $B_z(A)$, 
%respectively)
% flux rope axis determination - cloud events
% y axis obtained from MVAB with constraint - all points
% no fitting - residue from point-wise subtraction
% fitting in hu12n_opt.m
% auf residue map links klicken-> zeigt für alternative direction alle Parameter
% rechts klicken -> direction passt und speichern 

% Remarks Christian Möstl Aug. 2006 IWF GRAZ
% INPUT: ALLCASE.par (oberste Zeilen, Directory Name des Events), 
%        CASE1.par (gibt an wo die Daten sind (Filename) und dann die Indizes der Wolke Anfang-Ende, dann 'adjp'   
%        dht.dat  (deHoffmann-Teller velocity)  
%
% OUTPUT:
% Zzp ist die Invarianz-Richtung in GSE
% Zzpg ist die alternative Invarianz-Richtung in GSE (durch Mausklick mit
% links)
% alle Invarianz-Basisvektoren stehen in 'hu12n_opt.m' als xs,ys,zs
% bgse: Magnetfeld-Messdaten in GSE
% [xa,ya,za]: aus MVAB mit constraint x1*vht=0
% min_direction: Winkel [alpha1,alpha2] im [xa,ya,za] System
% Zp ist Invarianz-axis im [xa,ya,za] System
% Zmin ist alternative Invarianz-axis im [xa,ya,za] System
% Bz=bgse*Zzp' soll grösser 0 sein, dann stimmt die Richtung
% optres (matrix) enthält die Residuen R für alle Richtungen

% FILES:
% 'zs.dat' enthält Invarianzrichtung

%callse
%mvab2.m, mvab2D.m, rotatexo.m, rotatey.m, rotatezo.m, mvab_perp.m, BengtHu3.m, hu12n_opt.m

clear
format short;




%*****************PARAMETERS**************************************

NL=0.00 %noise level
latsteps=5;
longsteps=10;

adjp=1 % p=abs(adjp)*NkT; if adjp<=-1 Pt=Bz, also Bz wird optimised 
        % adjp=0 -> B_z^2/2 /mu0 wird optimised
points=101; %101 % as in Hu original; Hu's SR was .05; %.1;         %.3


%rn=i2-i1+1; Anzahl der gemessenen Datenpunkte
%SR=points/rn; Verhältnis grid zu messung, ungefähr ein Drittel??
%.05; %.1;         %.3
%nx=round(SR*rn); 	%for resampling if nx~=rn
%nx=min([rn nx])    %d.h. resampling maximal mit rn


% order=1;	% order of the polynomial -- ptol = exp ( poly(A) )

%*******************************************************************



 nden=1e6;  
 nv=1e3; 
 nb=1e-9; 
 k=1.38*1e-23; 
 mu0=4.0*pi*1e-7;
 ec=1.6*1e-19;

f0=fopen('ALLCASE.par', 'rt');
%dirstr=fscanf(f0, '%s\n')
dirstr=deblank(fgets(f0));
fclose(f0);

cd(dirstr);
%end
diary('output_optcloud.txt') % in Output stehen dann alle Parameter!

disp('******************************************************************')
disp('------------------------- START optcloud.m ------------------')


parfile=['CASE1.par'];

fpar=fopen(parfile, 'rt');
str=fscanf(fpar, '%c', 21);
while feof(fpar) ~= 1
i12=fscanf(fpar,'%i %i %i\n', 3); 
end
fclose(fpar);



%****************************** Daten laden*****************************
cd ..
load(str);
cd(dirstr)
CASE1=clouddata;
i1=i12(1);
i2=i12(2);
%if nargin==0
%adjp=i12(3);
%end
%adjp
%i1=input('i1 ?')
%i2=input('i2 ?')

vhtfile=['dht.dat'];


rn=i2-i1+1;
SR=points/rn; %.05; %.1;         %.3
nx=round(SR*rn); 	%for resampling nx~=rn
nx=min([rn nx])


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

if size(CASE1,2)>15
 Ne=CASE1(i1:i2, 16)*nden;
 Te=CASE1(i1:i2, 17);
 pe=abs(adjp)*Ne.*Te*k*1e6; 
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

%************************** Ergebnisse aus prefr.m laden ***************

load(vhtfile);

vhtg=dht;	% deHoffmann-Teller velocity

vht=vhtg.*nv; 
vhtn=sqrt(vht(1).*vht(1)+vht(2).*vht(2)+vht(3).*vht(3)); %Betrag vht -> 'vhtn'

%%%%
bgse=[bx by bz]; 
V=[vx*nv-vht(1) vy*nv-vht(2) vz*nv-vht(3)];

b=bgse.*nb;
bb=sqrt(b(:,1).^2+b(:,2).^2+b(:,3).^2);
bmax=max(bb);
%%%%


%%%%%%%%%%%%%%%%  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%xa=[1 0 0]; % ;
realy=[0 1 0];			%% old lepping axis
load bzh.dat;
realz=bzh;				%% old axis X(2) from MVAB %%% mit <Bn>=0 Anm. C.M.

 	nvht=-vht;
	vht_unit=nvht/norm(nvht);


%********************* (x', y', z') *************************************
Xperp=mvab_perp(b, 0, vht_unit);       %%%%%%%%%% MVAB mit x1*vht=0


ya=Xperp(:,1)';
xa=vht_unit;
za=cross(xa,ya);


angleZmin=100;

nw=0;
















%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SCHLEIFE ********************************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% geht bis ~Zeile 850 %%%%%%%%%%%


%***********************************************************************

while abs(angleZmin-90)>=3 & nw<1
    
    
    
 b=b*[xa' ya' za'];
 nvht=nvht*[xa' ya' za'];
 
 
%%************************** Lepping's axis ***************
 [XX, Lamerr]=mvab2(b, 0);
 leppz=XX(:,2)';			MVUB=leppz/[xa' ya' za'];
 if dot(leppz, [0 0 1]) < 0
	leppz=-leppz;
 end

leppolar=acos(dot(leppz, [0 0 1]));
leplong=acos(dot(leppz, [1 0 0])/sin(leppolar));
if dot(leppz, [0 1 0]) < 0
	leplong=2*pi-leplong;
end

%%%%%%%%%	get the uncertainty line of Lepping axis
%[errpolar21, errlong21]=errlines([XX(:,2) XX(:,1) XX(:,3)], Lamerr(1,2));
%[errpolar23, errlong23]=errlines([XX(:,2) XX(:,3) XX(:,1)], Lamerr(2,2));
%%%%%%%%%

realy=realy*[xa' ya' za'];
if dot(realy, [0 0 1]) < 0
	realy=-realy;
end
ypolar=acos(dot(realy, [0 0 1]));
ylong=acos(dot(realy, [1 0 0])/sin(ypolar));
if dot(realy, [0 1 0]) < 0
	ylong=2*pi-ylong;
end

realz=realz*[xa' ya' za'];
if dot(realz, [0 0 1]) < 0
	realz=-realz;
end
zpolar=acos(dot(realz, [0 0 1]));
zlong=acos(dot(realz, [1 0 0])/sin(zpolar));
if dot(realz, [0 1 0]) < 0
	zlong=2*pi-zlong;
end

%*************************************************************

x_HT=vhtn*time';
Ba=b; %

for jf=1:3
[XBh,M13Bh]=BengtHu3(Ba, x_HT, jf); 	% to real case (0,L)

Bhuz=XBh(:,1);
Lamda3(jf)=M13Bh(2);
errbar(jf)=M13Bh(3)*pi/180;
if dot(Bhuz, [0 1]) < 0
	Bhuz=-Bhuz;
end

Bpolar(jf)=acos(dot(Bhuz, [0 1]));
Blong(jf)=0;
if dot(Bhuz, [1 0]) < 0
	Blong(jf)=pi-Blong(jf);
end
angSH=acos(dot([Bhuz(1) 0 Bhuz(2)], realz))*180/pi;
end

angLep=acos(dot(leppz, realz))*180/pi;


%%%%%%%%%%%%	Resampling  or Interpolating of B and p %%%%%%%%%%%%%%
if nx~=rn 
% 	bxi1=resample(b(:,1), nx, rn);
% 	byi1=resample(b(:,2), nx, rn);
% 	bzi1=resample(b(:,3), nx, rn);
% 	ps=resample(p, nx, rn);
	dts=[dt(1) ones(1,nx-1)*(sum(dt)/(nx-1))];
    intopt='spline';
    bxi1=interp1(cumsum(dt)', b(:,1), cumsum(dts)', intopt);
    byi1=interp1(cumsum(dt)', b(:,2), cumsum(dts)', intopt);
    bzi1=interp1(cumsum(dt)', b(:,3), cumsum(dts)', intopt);
    ps=interp1(cumsum(dt)', p, cumsum(dts)', intopt);
    
	clear b bb;
	b=[bxi1 byi1 bzi1];
	bb=sqrt(b(:,1).^2+b(:,2).^2+b(:,3).^2);
	nx=length(bxi1); %ceil((nx/rn)*length(Xa1));
else
	ps=p;
	dts=dt;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555


%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%			PSD
%%%if SR==1
	%spec=figure;subplot(211);
	%psd(bb, 512, 1/dts(2)); 
%%%end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% noise addition
%NL=0.05;
dBd=NL*mean(bb);
dB=randn(nx,3)*dBd;

b0=b+dB;
b=b0;
bb=sqrt(b(:,1).^2+b(:,2).^2+b(:,3).^2);

%figure;
%plot(cumsum(dts),b0,cumsum(dts),b0-dB,'--');

%%%%%%%%%%%%%%%%			PSD
%if NL~=0
%	figure(spec);subplot(212);
%	psd(bb, 512, 1/dts(2)); 
%end
%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%% INTERPOLATION OPTION FOR Pt vs A plot %%%%%%%%%%%%%%%%%%%

intopt='linear'; % 'linear' for 0109


%%%%%%%%%%%%%%%%
xal=[1 0 0]; % ;
yal=[0 1 0];
zal=[0 0 1];

%if nw==0
i1=[];
i2=nx;
%end



optres=[];
latmax=90;
latstep=latsteps; %5;
long0=0;
longstep=longsteps; %10;
long1=360;














%%%%%%%%%%%%%%%%%%%%%%%%%%% SCHLEIFE ÜBER ALLE RICHTUNGEN %%%%%%%%%%%%%

for alpha1=long0:longstep:long1
	for alpha2=0:latstep:latmax
 [xap, yap, zap]=rotatezo(xal, yal, zal, alpha1+90); %rotate around z by 0-360 deg
 [xs, ys, zs]=rotatexo(xap, yap, zap, alpha2); %rotate around x by 0-90 deg

 zs=real(zs);
 zs=zs/norm(zs);

 i1=[];	%choose symmetric interval for any z?
 i2=nx;


xsv=nvht-dot(nvht, zs)*zs; %xs,ys,zs vorläufiges Koord.-System
xs=xsv/norm(xsv);
ys=cross(zs, xs);

vhtxs=norm(xsv);
dxa=dts*vhtxs;
xi=cumsum(dxa);

 %A1=A(mid,:);
bxs=b*xs';  % B in xs,ys,zs projizieren
bys=b*ys';
bzs=b*zs';


% integrate A along the satellite trajectory

 dA(1)=0;
 for i=2:nx

    dA(i)=-(bys(i)+bys(i-1)).*dxa(i).*0.5;
 end
 A1=[cumsum(dA)];


 
 % 4 if Schleifen
 
if (alpha1==0 & alpha2==0) | isempty(i1)
    
%%%%% interval with equal A value at ends	
if bys(1) <= 0 | bys(2) < 0
imid0=find(A1(:)==max(A1(:)));
imid=imid0(1);
if A1(1)<A1(nx)
	i1=find(abs(A1(1:imid)-A1(nx))==min(abs(A1(1:imid)-A1(nx))));
	i2=nx;
else
	i1=1;
	i2=find(abs(A1(imid:nx)-A1(1))==min(abs(A1(imid:nx)-A1(1)))) ...
		+imid-1;
end
Aur=max([A1(i1) A1(i2)]);
else
imid0=find(A1(:)==min(A1(:)));
imid=imid0(1);
if A1(1)>A1(nx)
	i1=find(abs(A1(1:imid)-A1(nx))==min(abs(A1(1:imid)-A1(nx))));
	i2=nx;
else
	i1=1;
	i2=find(abs(A1(imid:nx)-A1(1))==min(abs(A1(imid:nx)-A1(1)))) ...
		+imid-1;
end
Aur=min([A1(i1) A1(i2)]);
end
% figure(3);hold off;plot(fliplr([i1:imid]),fliplr(A1(i1:imid)), 'o-');
%       hold on; plot([imid:i2],A1(imid:i2),'*-');
%       hold on; plot(A1); line([0 i2],[0 0]);
% %pause;
% 
end




if isempty(i1) | isempty(i2) | i1==nx | i2==1 | abs(i1-i2)<=nx/2 | ...
   ((alpha1==0 | alpha1==180 | alpha1==360) & alpha2==90) 
	residue=20;
else

%%%%%%%%%%%%%%%%%%%%
if adjp<=-1 
    Pt=bzs';
else
    Pt=bzs'.^2/(2*mu0)+ps';
end
%Pt=exp(-2*A(mid,:))/2;
%lnpt=log(Pt);

Au=linspace(A1(imid), Aur, round((i2-i1+1)/2));
[A1half, I1h]=sort(A1(i1:imid));
Pt1half=interp1(A1half, Pt(I1h+i1-1), Au, intopt);
[A2half, I2h]=sort(A1(imid:i2));
Pt2half=interp1(A2half, Pt(I2h+imid-1), Au, intopt); 

%pwsub=abs(Pt2half-Pt(i1:imid));
pwsub=abs(Pt2half-Pt1half);

%%%%%%[aptf, S]=polyfit(A1, lnpt, porder);
%%%%%%[aptf, S]=polyfit(A1(i1:i2), Pt', porder+2);
%%%%%%dapt=polyder(aptf);

%%%%%ptfit=exp(polyval(aptf, A1));
%%%%%ptfit=(polyval(aptf, A1(i1:i2)));

%%%%aptf=[a^2/2 0 0];
%%%%%dapt=polyder(aptf);

%***************************** RESIDUUM **************
	
	residue=norm(pwsub)/(max(Pt(i1:i2))-min(Pt(i1:i2))); %/sqrt(length(pwsub)); 
end	

	optres((alpha1-long0)/longstep+1, alpha2/latstep+1)=residue;
end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ENDE DER FOR SCHLEIFEN %%%%%%%%%%%%%%%%%5
















%%%%%%%%%%%%%%%%%%% FINDEN DER MINIMUM-RESIDUE DIRECTION %%%%%%%%%%%%%

Alpha1=[long0:longstep:long1]'*pi/180; 
Alpha2=[0:latstep:latmax]/latmax;
X=cos(Alpha1)*Alpha2;
Y=sin(Alpha1)*Alpha2;

[nlong, nlat]=size(X);
minres=min(min(optres(:, 1:nlat-1))) %-8 for survey
[J I]=find(optres==minres);
alpha1=Alpha1(J)*180/pi;
alpha2=Alpha2(I)*latmax;

min_direction=[alpha1 alpha2];


%optres=(optres-minres);

reslevel=minres;

Yppolar=[];
Yplong=[];
a1=alpha1(1)*pi/180;       %alpha1 to rad
a2=alpha2(1)*pi/180;       %alpha2 to rad 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minlong=[a1(1) leplong];
minlat=[a2(1) leppolar];
linc=[1 0 0; 0 0 1];
figure;











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FOR SCHLEIFE EIN %%%%%%%%%%%%%%%%%%

for kk=1:2 %% Für minimum-residue und lepping axis...
[xap, yap, zap]=rotatezo(xal, yal, zal, minlong(kk)*180/pi+90); %rotate around z by 0-360 deg
[xs, ys, zs]=rotatexo(xap, yap, zap, minlat(kk)*180/pi); %rotate around x by 0-30 deg

zs=real(zs);
zs=zs/norm(zs);

i1=[];	%choose symmetric interval for any z?
i2=nx;

xsv=nvht-dot(nvht, zs)*zs;
xs=xsv/norm(xsv);
ys=cross(zs, xs);

vhtxs=norm(xsv);
dxa=dts*vhtxs;

%%%%%%%%%%%%%%
% upper half plane.


 %A1=A(mid,:);
bxs=b*xs';
bys=b*ys';
bzs=b*zs';
% integrate A along the satellite trajectory

 dA(1)=0;
 for i=2:nx

    dA(i)=-(bys(i)+bys(i-1)).*dxa(i).*0.5;
 end
 A1=[cumsum(dA)];

%%%%%%%%%%%%
%n0=find(A1==max(A1));
%figure;
%plot(A1(1:n0), bzs(1:n0).^2/2, 'o-', ...
%A1(n0+1:nx), bzs(n0+1:nx).^2/2, '*-');
%figure;
%contour(x,y,lBl,30); axis([x(1) x(nr) y(1) y(ny)]);
%axis equal;
%quiver(x,y,dAdy,By); axis equal;
%%%%%%%%%%%%%%

if (alpha1==0 & alpha2==0) | isempty(i1)
 %%%%% interval with equal A value at ends	
 if bys(1) <= 0 | bys(2) < 0
 imid0=find(A1(:)==max(A1(:)));
 imid=imid0(1);
 if A1(1)<A1(nx)
	i1=find(abs(A1(1:imid)-A1(nx))==min(abs(A1(1:imid)-A1(nx))));
	i2=nx;
 else
	i1=1;
	i2=find(abs(A1(imid:nx)-A1(1))==min(abs(A1(imid:nx)-A1(1)))) ...
		+imid-1;
 end
 Aur=max([A1(i1) A1(i2)]);
 else
 imid0=find(A1(:)==min(A1(:)));
 imid=imid0(1);
 if A1(1)>A1(nx)
	i1=find(abs(A1(1:imid)-A1(nx))==min(abs(A1(1:imid)-A1(nx))));
	i2=nx;
 else
	i1=1;
	i2=find(abs(A1(imid:nx)-A1(1))==min(abs(A1(imid:nx)-A1(1)))) ...
		+imid-1;
 end
 Aur=min([A1(i1) A1(i2)]);
 end
end



% kk
% [i1 i2]

if isempty(i1) | isempty(i2) | i1==nx | i2==1 | abs(i1-i2)<=nx/2 | ...
   ((alpha1==0 | alpha1==180 | alpha1==360) & alpha2==90) 
	residue=20;
 else

 %%%%%%%%%%%%%%%%%%%%
 if adjp<=-1
    Pt=bzs';
 else
 Pt=bzs'.^2/(2*mu0)+ps';
 end
%Pt=exp(-2*A(mid,:))/2;
%lnpt=log(Pt);

Au=linspace(A1(imid), Aur, round((i2-i1+1)/2));
[A1half, I1h]=sort(A1(i1:imid));
Pt1half=interp1(A1half, Pt(I1h+i1-1), Au, intopt); %spline for 05/02
[A2half, I2h]=sort(A1(imid:i2));
Pt2half=interp1(A2half, Pt(I2h+imid-1), Au, intopt); 

%pwsub=abs(Pt2half-Pt(i1:imid));
pwsub=abs(Pt2half-Pt1half);

%%%%%%%%%%%%% Calculate uncertainty by error propagation 
%%%%%%%%%%%% as described in Hu, 3/2004

sig2_1st=var(Pt1half);
sig2_2nd=var(Pt2half);



%%%%%%[aptf, S]=polyfit(A1, lnpt, porder);
%%%%%%[aptf, S]=polyfit(A1(i1:i2), Pt', porder+2);
%%%%%%dapt=polyder(aptf);

%%%%%ptfit=exp(polyval(aptf, A1));
%%%%%ptfit=(polyval(aptf, A1(i1:i2)));

%%
n0=imid; %find(A1==min(A1));
hold on;
%plot(A1(i1:n0), bzs(i1:n0).^2/2, 'o-', ...
%A1(i1:n0), Pt2half, '*-', linc(kk));
hPtA=plot(Au, Pt1half, 'o-',...
Au, Pt2half, '*-');
set(hPtA, 'color', linc(kk,:));
end
%%%%aptf=[a^2/2 0 0];
%%%%%dapt=polyder(aptf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
 %%%%%%%%%%%% Gleichungen (A2) und (A4) in Hu, 3/2004
 
	min_residue(kk)=norm(pwsub)/(max(Pt(i1:i2))-min(Pt(i1:i2))); %/sqrt(length(pwsub)); 
    uncertainty(kk)=sqrt(sig2_1st+sig2_2nd)/(max(Pt(i1:i2))-min(Pt(i1:i2)));
end
%%%%%%%%%%%%%%%%%%%%%%%% FOR SCHLEIFE AUS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
















%figure(Hcloud);
%line([x(i1) x(i2)]/R0, [y(mid) y(mid)]/R0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%% calculate uncertainty in z axis orientation (Hu 3/2004, JGR)
J=J(1);
I=I(1);
Jm1=J-1; Jp1=J+1;
if J==length(Alpha1)
    Jp1=2;
end
if J==1
    Jm1=length(Alpha1)-1;
end

p_a1=polyfit([Alpha1(Jm1) Alpha1(J) Alpha1(Jp1)], [optres(Jm1,I) optres(J,I) optres(Jp1,I)], 2); 
p_a2=polyfit([Alpha2(I-1) Alpha2(I) Alpha2(I+1)]*latmax*pi/180, [optres(J,I-1) optres(J,I) optres(J,I+1)], 2); 
sig_a1=sqrt(uncertainty(1)/p_a1(1));
sig_a2=sqrt(uncertainty(1)/p_a2(1));

dZmin=[0 0 0];
dZmin(1)=sqrt((sin(a1)*sin(a2))^2*sig_a1^2+(cos(a2)*cos(a1))^2*sig_a2^2);
dZmin(2)=sqrt((sin(a2)*cos(a1))^2*sig_a1^2+(cos(a2)*sin(a1))^2*sig_a2^2);
dZmin(3)=sqrt(sin(a2)^2*sig_a2^2);

rtn=inv([xa' ya' za']);

dZzp=sqrt(dZmin.^2 * rtn.^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Zmin=[sin(a2)*cos(a1) sin(a2)*sin(a1) cos(a2)];
Zp(nw+1,:)=Zmin;
Coord(:,:,nw+1)=[xa' ya' za'];
angleZmin=acos(dot(Zmin, yal))*180/pi;
angZz(nw+1)=acos(dot(Zmin, realz))*180/pi;
Ymin=cross(Zmin, xal);
Ymin=Ymin/norm(Ymin);
Z00=cross(xal, Ymin);
for aa=0:5:180
	[xpp ypp zpp]=rotatey(xal, Ymin, Z00, -aa);
	if dot(xpp, [0 0 1]) < 0
	xpp=-xpp;
	end
	yppolar=acos(dot(xpp, [0 0 1]));
	if yppolar==0
		yplong=0;
	else
	yplong=acos(dot(xpp, [1 0 0])/sin(yppolar));
	end
	if dot(xpp, [0 1 0]) < 0
	yplong=2*pi-yplong;
	end
	Yppolar=[Yppolar yppolar];
	Yplong=[Yplong yplong];
end
Yppolar=[Yppolar fliplr(Yppolar)];
Yplong=[Yplong Yplong-pi];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







% ******************************* RESIDUE MAP ********************

%***************************** FIGURE ***********************************
hres=figure;
%pcolor(X,Y, optres);
plot(X,Y,'k.', 'markersize', 2);
line((a2/(pi/2))*cos(a1),(a2/(pi/2))*sin(a1), ...
'marker', '.', 'color', [0. 0. 0.], 'markersize', 12, 'linewidth', 1); 
text([1 0 -1.15 0], [0 1.05 0 -1.05], ['  0^\circ'; ' 90^\circ'; '180^\circ'; '270^\circ'], ...
	'fontname', 'times', 'fontsize', 8);
text([0 0 0 0], [Y(10,5) Y(10,9) Y(10,13) Y(10,17)], ['20^\circ'; ...
'40^\circ'; '60^\circ'; '80^\circ'], 'fontname', 'times', 'fontsize', 8, ...
'color', [0 0 0]);
%%view(90, 90);
%%colormap jet;
axis square; axis off;
%%hbar=colorbar; 
%%posbar=get(hbar, 'position'); 
%%set(hbar, 'position', [posbar(1) posbar(2)+posbar(4)/2 posbar(3)/2 posbar(4)/2]);
%%set(hbar, 'fontname', 'times', 'fontsize', 8);
% line((leppolar/(pi/2))*cos(leplong),(leppolar/(pi/2))*sin(leplong), ...
% 'marker', 'x', 'color', [0 0 0], 'markersize', 8, 'linewidth', 1); 
%line((ypolar/(pi/2))*cos(ylong),(ypolar/(pi/2))*sin(ylong), ...
%'marker', '*', 'color', [0 0 0], 'markersize', 6, 'linewidth', .5); 
%line((errpolar21/(pi/2)).*cos(errlong21),(errpolar21/(pi/2)).*sin(errlong21), ...
%'linestyle', '-', 'color', [1 0 0], 'markersize', 6, 'linewidth', .5); 
%line((errpolar23/(pi/2)).*cos(errlong23),(errpolar23/(pi/2)).*sin(errlong23), ...
%'linestyle', '-', 'color', [0 0 1], 'markersize', 6, 'linewidth', .5); 
% line([(Bpolar(1)/(pi/2))*cos(Blong(1))+errbar(1) ...
%  (Bpolar(1)/(pi/2))*cos(Blong(1))-errbar(1)] ...
% ,[(Bpolar(1)/(pi/2))*sin(Blong(1)) (Bpolar(1)/(pi/2))*sin(Blong(1))], ...
% 'marker', '+', 'color', [0 1 0], 'markersize', 5, 'linewidth', .5); 
% line((Bpolar(1)/(pi/2))*cos(Blong(1)),(Bpolar(1)/(pi/2))*sin(Blong(1)), ...
% 'marker', '+', 'color', [0 0 0], 'markersize', 8, 'linewidth', 1); 

%line((Bpolar(2)/(pi/2))*cos(Blong(2)),(Bpolar(2)/(pi/2))*sin(Blong(2)), ...
%'marker', 'x', 'color', [0 0 0], 'markersize', 5, 'linewidth', 1); 
%line((Bpolar(2)/(pi/2))*cos(Blong(2)),(Bpolar(2)/(pi/2))*sin(Blong(2)), ...
%'marker', 'o', 'color', [0 0 0], 'markersize', 5, 'linewidth', 1); 
%line([(Bpolar(2)/(pi/2))*cos(Blong(2))+errbar(2) ...
% (Bpolar(2)/(pi/2))*cos(Blong(2))-errbar(2)] ...
%,[(Bpolar(2)/(pi/2))*sin(Blong(2)) (Bpolar(2)/(pi/2))*sin(Blong(2))], ...
%'marker', '+', 'color', [0 1 0], 'markersize', 5, 'linewidth', .5); 
%line((Bpolar(3)/(pi/2))*cos(Blong(3)),(Bpolar(3)/(pi/2))*sin(Blong(3)), ...
%'marker', '+', 'color', [0 0 0], 'markersize', 5, 'linewidth', 1); 
%line((Bpolar(3)/(pi/2))*cos(Blong(3)),(Bpolar(3)/(pi/2))*sin(Blong(3)), ...
%'marker', 'o', 'color', [0 0 0], 'markersize', 5, 'linewidth', 1); 
%line([(Bpolar(3)/(pi/2))*cos(Blong(3))+errbar(3) ...
% (Bpolar(3)/(pi/2))*cos(Blong(3))-errbar(3)] ...
%,[(Bpolar(3)/(pi/2))*sin(Blong(3)) (Bpolar(3)/(pi/2))*sin(Blong(3))], ...
%%%'marker', '+', 'color', [0 0 1], 'markersize', 5, 'linewidth', .5); 
hold on;
[nlong, nlat]=size(X);


%%%%%%%%%%%% Konturlinien bei minres+uncertainty
reslevel=[0 uncertainty(1)]+reslevel;
[CS, H]=contour(X(:,1:nlat-1), Y(:,1:nlat-1), optres(:,1:nlat-1), reslevel, 'k-');
set(H, 'linewidth', .5, 'color', [.3 .3 .3]); %[.5 .5 .5]
%clabel(CS, H, 'fontsize', 8, 'fontname', 'times', 'color', 'k');
% line((zpolar/(pi/2))*cos(zlong),(zpolar/(pi/2))*sin(zlong), ...
% 'marker', 'o', 'color', [1 0 0], 'markersize', 8, 'linewidth', .5); 

xa=xal;
ya=Ymin;
za=Z00;
nw=nw+1;

relares=(min_residue(1)-min_residue(2))/(min_residue(1));
end


















%%% output figure to FIG file for now
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





AngZz=angZz(length(angZz));

%Zp
Zzp=Zp(size(Zp,1),:)/Coord(:,:,nw)
Bz=bgse*Zzp';
Iz=find(abs(Bz)==max(abs(Bz)));
Iz=Iz(1);
%Richtung so dass Bz>0 sein muss!!! 
Zzp=Bz(Iz)/abs(Bz(Iz))*Zzp;    


%%%%%%%%%%%%%%%%%%%%%%%%%%%% ERRORS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% calculate uncertainty in original and angular frame
%dZzp



%%%%%%%% ONLY for VOYAGER 2 60/1978
%if dirstr(1:3)== 'vo2' 
%    Zzp=-Zzp 
%end



clear r t n;
r=Zzp(1); t=Zzp(2); n=Zzp(3);

% FOR DELTA/LAMBDA representation like Hu 3/2004
%'delta' entspricht DELTA und 'phi' entspricht LAMBDA in Hu 3/2004
%bzw.
%'delta' entspricht delta und 'phi' entspricht lambda in Hu 8/2005

delta=acos(r)*180/pi;
phi=acos(t/sqrt(n*n+t*t));
if n < 0
    phi=2*pi-phi;
end
phi=phi*180/pi;
delta_phi=[delta phi];

%errors
ddelta2=dZzp(1)^2/(1-r*r);
dphi2=dZzp(2)^2*n*n/((n*n+t*t)^2)+dZzp(3)^2*t*t/((n*n+t*t)^2);

ddelta=sqrt(ddelta2)*180/pi;
dphi=sqrt(dphi2)*180/pi;
d_delta_phi=[ddelta dphi];



% FOR GSE representation: latitude-> theta and longitude-> phi2

theta=asin(n)*180/pi;
phi2=acos(r/sqrt(r.^2+t.^2))*180/pi;
if t < 0
	phi2=360-phi2;
else phi2;
end

theta_phi2=[theta phi2];

%Koordsys drehen, snst Formeln gleich wie oben 
%x->y, y->z, z-> x;

dtheta2=dZzp(2)^2/(1-t*t);
dphi4=dZzp(3)^2*r*r/((r*r+n*n)^2)+dZzp(1)^2*n*n/((r*r+n*n)^2);

dtheta=sqrt(dtheta2)*180/pi;
dphi3=sqrt(dphi4)*180/pi;
d_theta_phi2=[dtheta dphi3];

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% for multiple connected contours; examine alternative directions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%button 1 - linke Maustaste
%button 3 - rechte Maustaste
%button 2 - mittlere Maustaste oder Scrollrad
button=0;

while button ~= 3
    figure(hres);
[xg,yg, button]=ginput(1);
disg=sqrt((X-xg).^2+(Y-yg).^2);
[Jg Ig]=find(disg==min(disg(:)));
a1g=Alpha1(Jg)*180/pi;
a2g=Alpha2(Ig)*latmax;
alternative_z=[a1g a2g];

a1g=a1g(1)*pi/180;
a2g=a2g(1)*pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%disp('button!!!!')
%button
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

minlong=[a1(1) a1g];
minlat=[a2(1) a2g];
linc=[1 0 0; 0 0 1];
figure;



for kk=1:2 
 [xap, yap, zap]=rotatezo(xal, yal, zal, minlong(kk)*180/pi+90); %rotate around z by 0-360 deg
 [xs, ys, zs]=rotatexo(xap, yap, zap, minlat(kk)*180/pi); %rotate around x by 0-30 deg

 zs=real(zs);
 zs=zs/norm(zs);

 i1=[];	%choose symmetric interval for any z?
 i2=nx;


xsv=nvht-dot(nvht, zs)*zs;
xs=xsv/norm(xsv);
ys=cross(zs, xs);

vhtxs=norm(xsv);
dxa=dts*vhtxs;

%%%%%%%%%%%%%%
% upper half plane.


 %A1=A(mid,:);
bxs=b*xs';
bys=b*ys';
bzs=b*zs';
% integrate A along the satellite trajectory

 dA(1)=0;
 for i=2:nx

    dA(i)=-(bys(i)+bys(i-1)).*dxa(i).*0.5;
 end
 A1=[cumsum(dA)];

%%%%%%%%%%%%
%n0=find(A1==max(A1));
%figure;
%plot(A1(1:n0), bzs(1:n0).^2/2, 'o-', ...
%A1(n0+1:nx), bzs(n0+1:nx).^2/2, '*-');
%figure;
%contour(x,y,lBl,30); axis([x(1) x(nr) y(1) y(ny)]);
%axis equal;
%quiver(x,y,dAdy,By); axis equal;
%%%%%%%%%%%%%%

if (alpha1==0 & alpha2==0) | isempty(i1)
%%%%% interval with equal A value at ends	
if bys(1) <= 0 | bys(2) < 0
imid0=find(A1(:)==max(A1(:)));
imid=imid0(1);
if A1(1)<A1(nx)
	i1=find(abs(A1(1:imid)-A1(nx))==min(abs(A1(1:imid)-A1(nx))));
	i2=nx;
else
	i1=1;
	i2=find(abs(A1(imid:nx)-A1(1))==min(abs(A1(imid:nx)-A1(1)))) ...
		+imid-1;
end
Aur=max([A1(i1) A1(i2)]);
else
imid0=find(A1(:)==min(A1(:)));
imid=imid0(1);
if A1(1)>A1(nx)
	i1=find(abs(A1(1:imid)-A1(nx))==min(abs(A1(1:imid)-A1(nx))));
	i2=nx;
else
	i1=1;
	i2=find(abs(A1(imid:nx)-A1(1))==min(abs(A1(imid:nx)-A1(1)))) ...
		+imid-1;
end
Aur=min([A1(i1) A1(i2)]);
end
end




zs_1_red_zsg_2_blue=[kk i1 i2]

if isempty(i1) | isempty(i2) | i1==nx | i2==1 | abs(i1-i2)<=nx/2 | ...
   ((alpha1==0 | alpha1==180 | alpha1==360) & alpha2==90) 
	residue=20;
else

%%%%%%%%%%%%%%%%%%%%
if adjp<=-1
    Pt=bzs';
else
Pt=bzs'.^2/(2*mu0)+ps';
end
%Pt=exp(-2*A(mid,:))/2;
%lnpt=log(Pt);

Au=linspace(A1(imid), Aur, round((i2-i1+1)/2));
[A1half, I1h]=sort(A1(i1:imid));
Pt1half=interp1(A1half, Pt(I1h+i1-1), Au, intopt); %spline for 05/02
[A2half, I2h]=sort(A1(imid:i2));
Pt2half=interp1(A2half, Pt(I2h+imid-1), Au, intopt); 

%pwsub=abs(Pt2half-Pt(i1:imid));
pwsub=abs(Pt2half-Pt1half);

%%%%%%%%%%%%%%%%% Calculate uncertainty by error propagation
sig2_1st=var(Pt1half);
sig2_2nd=var(Pt2half);

%%%%%%[aptf, S]=polyfit(A1, lnpt, porder);
%%%%%%[aptf, S]=polyfit(A1(i1:i2), Pt', porder+2);
%%%%%%dapt=polyder(aptf);

%%%%%ptfit=exp(polyval(aptf, A1));
%%%%%ptfit=(polyval(aptf, A1(i1:i2)));

%%
n0=imid; %find(A1==min(A1));
hold on;
%plot(A1(i1:n0), bzs(i1:n0).^2/2, 'o-', ...
%A1(i1:n0), Pt2half, '*-', linc(kk));
hPtA=plot(Au, Pt1half, 'o-',...
Au, Pt2half, '*-');
set(hPtA, 'color', linc(kk,:));
end
%%%%aptf=[a^2/2 0 0];
%%%%%dapt=polyder(aptf);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	min_residue(kk)=norm(pwsub)/(max(Pt(i1:i2))-min(Pt(i1:i2))); %/sqrt(length(pwsub)); 
    uncertainty(kk)=sqrt(sig2_1st+sig2_2nd)/(max(Pt(i1:i2))-min(Pt(i1:i2)));
end


%alternative direction
Zmin=[sin(a2g)*cos(a1g) sin(a2g)*sin(a1g) cos(a2g)];
Zzpg=Zmin/Coord(:,:,nw);  %Coord => xa ya za
%save zsg.dat Zzpg -ascii;

 Bz=bgse*Zzpg';
 Iz=find(abs(Bz)==max(abs(Bz)));
 Iz=Iz(1);
 %Richtung so dass Bz>0 sein muss!!! 
 Zzpg=Bz(Iz)/abs(Bz(Iz))*Zzpg;


disp(' ------------------ 1st: alternative ---------------------')
if Zzp~=Zzpg 
    a=hu12n_opt(adjp, Zzpg);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



Bz=bgse*Zzp';



disp('----------------- 2nd: minimum residue ------------------')

%disp('Zzp bei Zeile 1158 changed!')
%Zzp=Zzpg;
savefiles=hu12n_opt(adjp, Zzp);


delta_phi=double(delta_phi);
d_delta_phi=double(d_delta_phi);
Zzp=double(Zzp);
dZzp=double(dZzp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%55
%chris 24.10.2006 ebenso 1150






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5





%%%%%%%% ONLY for VOYAGER 2 60/1978
%if dirstr(1:3)== 'vo2' 
%    Zzp=-Zzp 
%end

if savefiles=='y'
    save zs.dat Zzp -ascii;
    save dzs.dat dZzp -ascii;
    save zs_ang.dat delta_phi d_delta_phi -ascii;
    %type zs_ang.dat;
    disp('Invarianz-Achse in delta-lambda (GSE) + error')
    disp(delta_phi)
    disp(d_delta_phi)
    disp('Invarianz-Achse in theta-phi (GSE) + error')
    disp(theta_phi2)
    disp(d_theta_phi2)
end

disp('Bz>0? Mean [nT]:')
disp(mean(bgse*Zzp'))
disp('Max Bz [nT]')
disp(max(bgse*Zzp'))
disp('Invarianz Richtung in GSE:')
Zzp
disp('links minimumresidue, rechts alternative Stelle:')
min_residue
uncertainty
disp('Number of original data points:')
rn
disp('Number of grid points:')
nx

%save 'optcloud_var.mat'

%hu34n_2;



cd('..')

disp('------------------------- END optcloud.m ------------------')
diary off;
