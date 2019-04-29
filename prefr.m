% prefr.m is a program that carries out the minimum variance analysis and 
% deHoffmann Teller analysis. 
%
% written by Qiang Hu, edited by Christian Möstl 10/2006
% 
%
% a nice summary on MVA from Huttunen et al., 
% Annales Geophysicae (2005) 23: 1–17:
%
%"The MVA method can be applied satisfyingly to the directional changes of the
%magnetic field vector exceeding 30 degrees. The large ratio of the
%intermediate eigenvalue Lambda 2 to the minimum eigenvalue Lambda 3 indicates
%that the eigenvectors are well defined. We required
%that Lambda 2/ Lambda 3 is greater than 2, based on the analysis of Lepping
%and Behannon (1980). BX' , BY' , and BZ' correspond to
%the magnetic field components in the directions of maximum,
%intermediate and minimum variance. TheMVA analysis provides
%us with the estimation of the orientation of the MC axis
%(theta, phi). theta and phi are the latitudinal and longitudinal angels
%of the magnetic field vector in solar ecliptic coordinates;
%theta=90 degree is defined northward and phi=90 degree is defined eastward.
%The MC axis orientation corresponds to the direction of the
%intermediate variance that is seen from Eq. (1) as the axial
%component is zero at the boundaries of the MC. The radial
%component corresponds to the minimum variance direction
%and the azimuthal component corresponds to the maximum
%variance direction. The boundaries of MCs were determined
%by solar wind signatures (start of the smooth rotation of the
%magnetic field vector, drop in plasma beta, and plasma and
%field discontinuities) and by the eigenvalue ratio. In those
%cases where the boundaries defined by the different signatures
%disagreed we used the magnetic field rotation."
%
%other good introductions on MVA are 
%Bothmer and Schwenn (1998), Annales %Geophysicae, and
% the book "Analysis Methods for Multi-Spacecraft Data" 
% Editors Paschmann and Daly
% 2 chapters: Minimum variance and deHoffmanTeller analysis





%
% INPUT: ALLCASE.par (oberste Zeilen), 
%        CASE1.par (tells you the data filename, the beginning and the end of the flux rope interval and adjp)  
%
% OUTPUT: Figures: Daten, Hodogramm 
%       bxh.dat, byh.dat,bzh.dat, dht.dat
%       X - Eigenvektoren aus MVAB (MVUB)  
%       normal - Normale mit constraint

% Deutsche Anmerkungen: Christian Möstl
% X(:,1) - Maximum variance
% X(:,2) - Intermediate variance % auch Leppings first estimate axis (mit MVUB!)
% X(:,3) - minimum variance, ist estimator für Normalkomponente auf die TD
% normal - Normale aus MVAB (o. MVUB) mit constraint <B_n>=0
% VHT - HT Frame Geschwindigkeit in km/s, GSE
% vmV - Geschwindigkeiten im HT-Frame in km/s, GSE
% Va - lokale Alfvén Geschwindigkeit, km/s GSE
% pw - Linear Fit to Walen plot
% walen_slope =1 für RD; =0 für TD


%adjp=1.
%clear;
%close all;
format compact;
format short;

%vorzeichen=1;   %für normale
%vor=-1;    % fuer y-Achse umdrehen

f0=fopen('ALLCASE.par', 'rt');
%dirstr=fscanf(f0, '%s\n')
dirstr=deblank(fgets(f0)); %liest aus dem Anfang von Allcase Dir-Namen für das Event
fclose(f0);

cd(dirstr); % In Bearbeitungs-Dir für das jeweilige Event gehen  
%end

%diary('output_prefr.txt') % in Output stehen dann alle Parameter!

disp('******************************************************************')
disp('------------------------- START prefr.m ------------------')


fpar=fopen('CASE1.par', 'rt');
str=fscanf(fpar, '%c', 21); % str= Dateinamen der Daten, z.B.:..\data\winall291_95.mat
while feof(fpar) ~= 1
i12=fscanf(fpar,'%i %i %i\n', 3); 
end
fclose(fpar);
i1=i12(1); % Anfangspunkt der Daten
i2=i12(2); % Endpunkt der Daten
cd ..

load(str);
cd(dirstr)
CASE1=clouddata; % Daten werden geladen

i1
i2
disp(str);


%Alternativ: i1, i2 einlesen vom Benutzer

%i1=input('data interval starts at i1 =');
%i2=input('data interval ends at i2 =');

%diary(str(9:(length(str)-4)));
%%%%%%%%%%%%%diary('all_rho_aHT.out');

%%%% load the magnetic field only
%load(str2);
%MAG=eval(str2(9:(length(str2)-4)));
%[rM,cM]=size(MAG);

%mi1=input('interval to determine the normal start at mi1=');
%mi2=input('interval to determine the normal end at mi2=');

%Bm(:,1)=MAG(mi1:mi2,cM-2);
%Bm(:,2)=MAG(mi1:mi2,cM-1);
%Bm(:,3)=MAG(mi1:mi2,cM);

%%%%%%%%%%  Daten werden in einzelne Arrays geschrieben (zwecks Bearbeitung) %%%%%%%%%%%%%

Dtime=CASE1(:,7);
ltime=cumsum(Dtime);

%[ltime,lT]=issireform(str);

time=ltime(i1:i2)-ltime(i1);   
B(:,1)=CASE1(i1:i2,8);    % B in nT
B(:,2)=CASE1(i1:i2,9);
B(:,3)=CASE1(i1:i2,10);
v(:,1)=CASE1(i1:i2,13);   % v in km/s   
v(:,2)=CASE1(i1:i2,14);
v(:,3)=CASE1(i1:i2,15);
Np=CASE1(i1:i2,11);        %N in ccm^-3
Tp=CASE1(i1:i2,12);        %T in 10^6 K
if size(CASE1,2)>15
 Ne=CASE1(i1:i2,16);
 Te=CASE1(i1:i2,17);
end
%T=lT(i1:i2);

DD=CASE1(i1:i2,2);
hh=CASE1(i1:i2,4);
mmm=CASE1(i1:i2,5);
ssm=CASE1(i1:i2,6);
%timem=mmm*60+ssm;
%timei=mmi*60+ssi;
%display the original data
Y=CASE1(i1:i2,3); M=CASE1(i1:i2,1);  
Datem=datenum(double(Y),double(M),double(DD),double(hh),double(mmm),double(ssm));   %Datum in datenumber konvertieren
%Datei=datenum(Y,M,D,hh,mmi,ssi);

%Zeit ausgeben:

zeitraum=[dirstr(1:3),' ',num2str(CASE1(i1,1)),'/',num2str(CASE1(i1,2)),'/', ...
  num2str(CASE1(i1,3)),' ',num2str(CASE1(i1,4)),':', ...
 num2str(CASE1(i1,5)),':',num2str(round(CASE1(i1,6))),'-',num2str(CASE1(i2,1)),'/',...
num2str(CASE1(i2,2)),'/',num2str(CASE1(i2,3)),' ',num2str(CASE1(i2,4)), ...
  ':',num2str(CASE1(i2,5)),':',num2str(round(CASE1(i2,6)))]

%****************************** Daten darstellen *******************************


Bmag=sqrt(B(:,1).^2+B(:,2).^2+B(:,3).^2);

%%%%% MAG 
figure;       
subplot(211);
plot(Datem,Bmag,Datem,B(:,1),Datem,B(:,2),Datem,B(:,3),'-');ylabel('B [nT]');legend('|B|','B_x','B_y','B_z');
datetick('x',15, 'keeplimits');grid on;
if dirstr(1:2)=='cl'
    datetick('x',13,'keeplimits');grid;
end   



%%%%%% PLASMA
subplot(212);
[AX, h1, h2]=plotyy(Datem,Np,Datem,Tp);%Ne,'--',Datem,Te,'--');
datetick('x',15,'keeplimits');grid;
if dirstr(1:2)=='cl'
    datetick('x',13,'keeplimits');grid;
end   
set(AX(2), 'XTickLabel', []);
legend([h1 h2], 'proton density [ccm^-3]','proton temp [MK]');%, 'e^- density', 'e^- Temp.');

if dirstr(1:2)=='cl'
    legend([h1 h2], 'ion density [ccm^-3]','ion temp [MK]');%, 'e^- density', 'e^- Temp.');
end 
grid on;

title([dirstr(1:3),' ',num2str(CASE1(i1,1)),'/',num2str(CASE1(i1,2)),'/', ...
  num2str(CASE1(i1,3)),' ',num2str(CASE1(i1,4)),':', ...
 num2str(CASE1(i1,5)),':',num2str(round(CASE1(i1,6))),'-',num2str(CASE1(i2,1)),'/',...
num2str(CASE1(i2,2)),'/',num2str(CASE1(i2,3)),' ',num2str(CASE1(i2,4)), ...
  ':',num2str(CASE1(i2,5)),':',num2str(round(CASE1(i2,6)))]);




%%%%%%%%%%%% CHANGE THIS FOR MVAB OR MVUB

%%%%%%%%%%%%% normalize the B vectors into unit vectors (falls MVUB)
lBl=sqrt(B(:,1).^2+B(:,2).^2+B(:,3).^2);
Bhat(:,1)=B(:,1)./lBl;
Bhat(:,2)=B(:,2)./lBl;
Bhat(:,3)=B(:,3)./lBl;
%%%%%%%%%%%%%%%%%%%%%%
%Bhat=B;    %für MVAB


disp('number of magnetic field data points for MVUB:')
ndat_mag=length(lBl)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%*************************  MVAB (oder MVUB) Analyse *********************
clear M;

for i=1:3
	for j=1:3      %%%%%%%%%%%%%%%% Magnetic variance matrix %%%%%%%%%%%%%%%%%%%%%%
%          M(i,j)=mean((Bhat(:,i).*Bhat(:,j)))-mean(Bhat(:,i))*mean(Bhat(:,j)); 
		M(i,j)=mean((B(:,i).*B(:,j)))-mean(B(:,i))*mean(B(:,j));
	end
end

[X,D]=eig(M);            %Eigenwerte berechnen; D -> Eigenwert-Matrix

format long;

[limda,In]=sort([D(1,1) D(2,2) D(3,3)]) ; %Sortieren der Eigenwerte nach Grösse (In->Index)
X=[X(:,In(3)) X(:,In(2)) X(:,In(1))];   % Eigenvektoren sortieren damit Reihenfolge wie oben         


limda=(fliplr(limda))'
X                  % Eigenvektoren jeweils in Spalten

aBxi=[mean(B*X(:,1)) mean(B*X(:,2)) mean(B*X(:,3))]
%aBxi=[mean(Bm*X(:,1)) mean(Bm*X(:,2)) mean(Bm*X(:,3))]






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MVAB plus constraint <B_n> = 0 %-> normalerweise fuer tangentiale
% Diskontinuitäten
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Bi=B*X;                %%% Bi.. Komponenten im MVA Koord.Sys
%Bi=Bm*X;
a=(mean(Bi,1)).^2;

C(1)=sum(a);
C(2)=-a(1)*(limda(2)+limda(3))-a(2)*(limda(1)+limda(3))-a(3)*(limda(1)+limda(2));
C(3)=a(1)*limda(2)*limda(3)+a(2)*limda(1)*limda(3)+a(3)*limda(1)*limda(2);

Lim=sort(roots(C));
Limmin=Lim(1);

Gama=1/sqrt(sum(a'./((limda-Limmin).^2)));
ni=Gama*mean(Bi,1)'./(limda-Limmin);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
normal=X*ni;	% normal obtained with constraint
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

aBni=mean(B*normal) %sollte ~= 0 sein! da constraint
%aBni=mean(Bm*normal)








%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%% FIGURE
%Hodogram pair determined by MVAB analysis


[roBi coBi]=size(Bi);
figure;
subplot('Position',[0.1 0.1 .45 .8]);
plot(Bi(:,2),Bi(:,1));axis('equal');
line(Bi(1,2),Bi(1,1),'marker','d','markerfacecolor','r','markersize',12);
line(Bi(roBi,2),Bi(roBi,1),'marker','x','markersize',12,'markeredgecolor','black');
grid;xlabel('B_2');ylabel('B_1');
title(['Hodogram pair for data interval ',dirstr(1:3),' ',num2str(CASE1(i1,4)),':', ...
  num2str(CASE1(i1,5)),':',num2str(round(CASE1(i1,6))),'-',num2str(CASE1(i2,4)), ...
  ':',num2str(CASE1(i2,5)),':',num2str(round(CASE1(i2,6)))]);


subplot('Position',[.65 0.1 .3 .8]);
plot(Bi(:,3),Bi(:,1));
line(Bi(1,3),Bi(1,1),'marker','d','markerfacecolor','r','markersize',12);
line(Bi(roBi,3),Bi(roBi,1),'marker','x','markersize',12,'markeredgecolor','black');
axis('equal');grid;xlabel('B_3');ylabel('B_1');





%%%%%%%%%%%%%%%%%%%%%%%%%%%%% deHoffmann Teller Analyse %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Now determine HT velocity for both VHT=Const and
% the accelerating HT frame: VHT=VHT0+aHT*time 

r=length(time);
for m=1:r
	BB=B(m,1)^2+B(m,2)^2+B(m,3)^2;
	for i=1:3
		for j=1:3
			if i==j
				Kijm(i,j,m)=BB*(1-B(m,i)*B(m,j)/BB);
			else
				Kijm(i,j,m)=-B(m,i)*B(m,j);
			end
		end
	end
	K1ijm(:,:,m)=Kijm(:,:,m)*time(m);
	K2ijm(:,:,m)=Kijm(:,:,m)*(time(m)^2);
end

K0=sum(Kijm,3)/r;
K1=sum(K1ijm,3)/r;
K2=sum(K2ijm,3)/r;

for m=1:r
	Rb(:,m)=Kijm(:,:,m)*([v(m,1) v(m,2) v(m,3)]');
	Rb2(:,m)=Rb(:,m)*time(m);
end

VHT=inv(K0)*(sum(Rb, 2)/r)        %%% Lösung von VHT in Beobachtungs-Koordinaten
VHTi=VHT'*X	                      %%% Lösung von VHT in MVAB-Koordinaten (z invariant) ???????
VHTn=VHT'*normal                  %%% VHT in Richtung Normalkomponente 

RHS1=sum(Rb,2)/r;
RHS2=sum(Rb2,2)/r;

VHT0aHT=([K0 K1;K1 K2])\([RHS1;RHS2]);

VHT0=VHT0aHT(1:3)
aHT=VHT0aHT(4:6)

VHT0n=VHT0'*normal
aHTn=aHT'*normal

for m=1:r
	aVHT(m,:)=VHT0'+aHT'*time(m);
end
save aHT aHT;                  %************ DHT Analysis mit konstantem a! ****
save aVHT aVHT;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

avmV=v-aVHT;
for m=1:r
	vmV(m,:)=v(m,:)-VHT';
end
%vcB=sum((cross((vmV),B,2)).^2,2);
DV=mean(sum((cross((vmV),B,2)).^2, 2))
aDV=mean(sum((cross(avmV,B,2)).^2, 2))

E=-cross(v, B, 2);
for m=1:r
	EHT(m,:)=-cross(VHT', B(m,:));
	aEHT(m,:)=-cross(aVHT(m,:), B(m,:));
end
	
%%%%%%%%%%%%%%%%%%%%%%%% Scatter plot E=vxB usw.%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%figure;hold on;
%for m=1:r % Test: bis 50 oder so
%	h1=scatter(EHT(m,:)/1000,E(m,:)/1000);
%	h2=scatter(aEHT(m,:)/1000,E(m,:)/1000,'r*');
% end
% axis('square');axis([-20 20 -20 20]);
% title(['E_{HT} vs E_c for data interval ',num2str(CASE1(i1,4)),':', ...
%   num2str(CASE1(i1,5)),':',num2str(CASE1(i1,6)),'-',num2str(CASE1(i2,4)), ...
%   ':',num2str(CASE1(i2,5)),':',num2str(CASE1(i2,6))]);
% plot([-20 20],[-20 20]);
% xlabel('E_{HT} in mV/m');ylabel('E_c in mV/m');
%legend([h1 h2],'Const V_{HT}','Accel. V_{HT}');
%% 

E1=[E(1,:)' EHT(1,:)'];aE1=[E(1,:)' aEHT(1,:)'];
for m=2:r
	E1=cat(1,E1,[E(m,:)' EHT(m,:)']);aE1=cat(1,aE1,[E(m,:)' aEHT(m,:)']);
end

S=corrcoef(E1);aS=corrcoef(aE1);
CCoef=S(1,2)                %******* Korrelationskoeffizient dHT Analyse
aCCoef=aS(1,2)              %******* Korrelationskoeffizient dHT Analyse mit const. Verzoegerung

save CCoef;







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%*************************** Richtungen nach MVAB wichtig!!!! ************
X3n=dot(normal, X(:,3));  
% X(:,1) - Maximum variance
% X(:,2) - Intermediate variance % auch Leppings axis (mit MVUB!)
% X(:,3) - minimum variance, ist estimator für Normalkomponente auf die TD
% normal - Normale aus MVAB (o. MVUB) mit constraint <B_n>=0


%%*************************** Richtungen nach MVAB wichtig!!!! ************
%%if abs(X3n) < .9
%if limda(2)/limda(3) < 10.
%	XX(:,3)=normal;
%	XX2=[-(normal(2)*X(2,2)+normal(3)*X(3,2))/normal(1) X(2,2) X(3,2)]';
%	XX(:,2)=XX2/norm(XX2);
%	XX(:,1)=cross(XX(:,3)', XX(:,2)')';
%else	
%%%%%%% output the results into files (using X/-X as initial invariance axes)
%if normal(1) > 0 & X3n < 0	% assuming normal is pointing toward the sun%
%	XX=-X;
	%XX=X;
 %else 
%	XX=X;
% end
%end

%if X3n < 0.
% XX=-XX
%end


%%%%%%%%%%%%%%%%%%%%%%%% X Eigenvektoren -> RICHTUNGEN FESTLEGEN %%%%%%
XX=X;




%TEST 1 Haendigkeit der Eigenvektoren
zdir=cross(XX(:,1),XX(:,2));
if zdir ~= XX(:,3)    % wenn left-handed
    XX(:,3)=-XX(:,3); % machs right-handed
end

%TEST 2 Normale richtig
X3n=dot(normal, XX(:,3)) %Normale soll in z-Richtung schauen?
if X3n < 0.              %wenn nicht    
    XX(:,3)=-XX(:,3)     % flip z; damit ist z  festgelegt
   
    
    ydir=cross(XX(:,3),XX(:,1));
    if ydir ~= XX(:,2)    % wenn left-handed
     XX(:,2)=-XX(:,2); % machs right-handed
    end
 
end


XX2=X;

LepZlatGSE2=asin(XX2(3,2))*180/3.14159265
LepZlongGSE2=acos(XX2(1,2)/sqrt(XX2(1,2)^2+XX2(2,2)^2))*180/3.14159265;
if XX2(2,2) < 0
	LepZlongGSE2=360-LepZlongGSE2
else LepZlongGSE2
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%zdir=cross(XX(:,3),XX(:,1));
%if dot(zdir,XX(:,2)) < 0.
%	XX(:,2)=XX(:,2);
%end

Lambda1_zu_Lambda2=limda(1)/limda(2)
Lambda2_zu_Lambda3=limda(2)/limda(3)

LepZlatGSE=asin(XX(3,2))*180/3.14159265
LepZlongGSE=acos(XX(1,2)/sqrt(XX(1,2)^2+XX(2,2)^2))*180/3.14159265;
if XX(2,2) < 0
	LepZlongGSE=360-LepZlongGSE
else LepZlongGSE
end

	
fx=fopen('bxh.dat','w');
fprintf(fx, '%e %e %e', XX(:,3)');
fy=fopen('byh.dat','w');
fprintf(fy, '%e %e %e', XX(:,1)');
fz=fopen('bzh.dat','w');
fprintf(fz, '%e %e %e', XX(:,2)');

fht=fopen('dht.dat','w');
fprintf(fht, '%e %e %e', VHT');

fno=fopen('normal.dat','w');
fprintf(fno, '%e %e %e', normal');

fyps=fopen('mvaby.dat','w');
fprintf(fyps, '%e %e %e',XX(:,2)');

fclose('all');

save zeitraum;
save d_in_AU;

%break;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%************************************** Walen plot ******************
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
mu0=4.0*pi*1e-7; %SI Units
rho=Np*1.673*1e-27*1e6; 	% kg/m^3


%lokale Alfvengeschwindigkeit
for m=1:r
	Va(m,:)=(B(m,:)*1e-9/sqrt(mu0*rho(m)))/1000;   % km/s GSE
end

%FIGURE
 %figure;hold on;
 %for m=1:r
% 	scatter(Va(m,:),vmV(m,:));
% 	scatter(Va(m,:),avmV(m,:),'r*');
% end
% axis('equal');axis([-400 400 -400 400]);xlabel('V_A in km/s');
%&% ylabel('V - V_{HT} in km/s');
%FIGURE END


V1=[Va(1,:)' vmV(1,:)'];aV1=[Va(1,:)' avmV(1,:)'];
for m=2:r
	V1=cat(1,V1,[Va(m,:)' vmV(m,:)']);
	aV1=cat(1,aV1,[Va(m,:)' avmV(m,:)']);
end
pw=polyfit(V1(:,1), V1(:,2), 1);

%%%%%%%%%%%%% WALEN SLOPE 
walen_slope=pw(1)

save walen walen_slope;

%FIGURE %apw=polyfit(aV1(:,1), aV1(:,2), 1);
% plot([-400 400],[polyval(pw, -400) polyval(pw, 400)],[-400 400],[polyval(apw,-400) polyval(apw, 400)],'r--');
% grid;title(['Walen plot for data interval ',num2str(CASE1(i1,4)),':', ...
%   num2str(CASE1(i1,5)),':',num2str(CASE1(i1,6)),'-',num2str(CASE1(i2,4)), ...
%   ':',num2str(CASE1(i2,5)),':',num2str(CASE1(i2,6))]);
%FIGURE END





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%****************************************** PRESSURE ****************
figure;
%adjp=1.6; 	% IRM
adjp=1.;
ec=1.6*1e-19;
k=1.38*1e-23;
Pp=adjp*Np.*Tp.*k*1e6*1e6;		% Temperature is in 10^6 K

%PB=(sum(B.^2,2)*(1e-18)/(2*mu0));
Bmag=sqrt(B(:,1).^2+B(:,2).^2+B(:,3).^2);

PB=(Bmag.*1e-9).^2/(2*mu0);

%if exist('Pe','var')
if size(CASE1,2)>15
    Pe=adjp*Ne.*Te.*k*1e6*1e6;
    Ptotal=Pp+PB+Pe;
else
    Ptotal=Pp+PB;
end    

subplot(411)
plot(Pp);ylabel('Pp');grid;
if size(CASE1,2)>15
 subplot(412)
 plot(Pe);ylabel('Pe');grid;
end
subplot(413)
plot(PB);ylabel('B^2/2u_0');grid;
subplot(414)
plot(Ptotal);ylabel('P_{tot}');grid;title(['adjustment factor =',num2str(adjp)]);

%%%%%%%%%%%%%%%%%%%%%%% Pressure overplots %%%%%%%%%%%%%%%%%%%%%%%%55


%d=Datem-datenum(Y(1),0,0,0,0,0);
%figure;
%plot(d,Pe,'b',d,Pp,'r',d,Pe+Pp,'k',d,PB,'g');
%ylabel('Pressure [Pa]');
%xlabel('Day of Year');
%grid;
%title('ACE pressures');
%legend('Pe','Pp','Pe+Pp','PB');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%mean(Pp./Pe)
%mean(PB./Pp)
%mean(PB./Pe)
%mean(PB./(Pe+Pp))


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ?????????????????? 
B02=max(BB)*1e-18;
rho0=mean(rho);
a0=norm(aHT)*1e3;
L=1e6*1e3;

P0=max(Pp);

rie=rho0*a0*mu0*L/B02;
rip=rho0*a0*L/P0;

pw_rie_rip=[pw(1) rie rip]

%diary off;

cd('..');
pwd

% %******************************** Vektoren visualisieren*******************
% 		 
% %x=zeros(4,1);
% %y=zeros(4,1);
% %z=zeros(4,1);
% x=zeros(2,1);
% y=zeros(2,1);
% z=zeros(2,1);
% 
% %u=XX(1,[1 3]);
% %v=XX(2,[1 3]);
% %w=XX(3,[1 3]);
% 
% u=XX(1,[1]);
% v=XX(2,[1]);
% w=XX(3,[1]);
% 
% %u(1,4)=normal(1);
% %v(1,4)=normal(2);
% %w(1,4)=normal(3);
% 
% figure;
% h1=quiver3(x,y,z,u,v,w,'b');
% hold on;
% %h2=quiver3(0,0,0,normal(1),normal(2),normal(3));
% h3=quiver3(0,0,0,XX(1,2),XX(2,2),XX(3,2));
% h4=quiver3(0,0,0,XX(1,3),XX(2,3),XX(3,3),'y');
% h5=quiver3(0,0,0,1,0,0,'k');
% 
% set([h1 h3],'LineWidth',1.8);      
% set(h5,'LineWidth',3);      
% if dirstr(1:3)=='vo2' 
%     xlabel('R´');ylabel('T´');zlabel('N');
% else    
%  xlabel('X');ylabel('Y');zlabel('Z');
% end
% %legend([h1 h2 h3 h4 h5],'x -> max Var.','normal from constraint ~ z','y -> interm Var.','z -> min Var.')
% legend([h1 h3 h4 h5],'x -> max Var.','y -> interm Var.','z -> min Var.', 'Sunward direction')
% title('Vectors from MVUB in GSE')
% axis([-1 1 -1 1 -1 1]);
% axis vis3d;
% 


disp('-------------------------MAIN RESULTS ------------------')

VHT
norm_VHT=norm(VHT)
CCoef

%Ausdehnung d entlang Trajektorie
AU=149.6*1e6;
delta_t=(Datem(length(Datem))-Datem(1))*24*60*60;
d_in_AU=(delta_t*norm_VHT)/AU

%figure;
%plot(vmV)  %V im VHT-Frame



disp('axis from MVUB in GSE')

XX(:,2)

disp('Orientation angles: Longitude Phi and Latitude Theta:')
[phi,theta]=vectodeg(XX(:,2));
phi_round=round(phi)
theta_round=round(theta)

save la23 Lambda2_zu_Lambda3
if Lambda2_zu_Lambda3 > 2 
    disp('Eigenvector condition satisfied! L2/L3 greater 2') 
end

if Lambda2_zu_Lambda3 < 2 
    disp('Eigenvector condition NOT satisfied! L2/L3 less than 2') 
end

disp('------------------------- END prefr.m ------------------')
%diary off;
