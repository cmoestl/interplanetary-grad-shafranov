
function getwind(year,doybeg,doyend,write_file)
% READ WIND DATA for gui
% getwind(2010,44,55,1)

%from NSSDC: 1-min MAG ASCII and 92-sec SWE ASCII
% Christian Möstl IWF March 2008

%makes a .MAT File in Folder GS-code\DATA

close all;
%clear
format compact;
format short;

disp('*********************** START ************************')






%%use to get DOY from date
%DOY=datenum(year,1,21,0,0,0)-datenum(year,0,0,0,0,0)
%%%%%%%%%%%%%%
%year=1996;
%doybeg=272;  %73,74
%doyend=273;
%eg. DOY 140-144 covers 140,141,142,143
secres=120; %Common time resolution for the finally interpolated data
%write_file=0;  % =1 if data is to be saved           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DOYbeg=datenum(year,0,doybeg,0,0,0);
%datestr(DOYbeg);
datebeg=datevec(DOYbeg);

DOYend=datenum(year,0,doyend,0,0,0);
%datestr(DOYend);
dateend=datevec(DOYend);

disp(['FULL TIME interval is:      ', datestr(DOYbeg),'    ',datestr(DOYend)])

%%%%%%%%%%%%%%%%%%%%%%%%%%% MAG %%%%%%%%%%%%%%%%%%%%%%%%%%%%
scstring='win';

ndays=doyend-doybeg;
yearstr=num2str(year);


a=datenum(year, 0, doybeg, 0, 0, 0); 
b=datevec(a);
if b(2) < 10    
  monthstr=['0',num2str(b(2))];
else
monthstr=num2str(b(2));
end

twomonths=0;
c=datenum(year, 0, doyend, 0, 0, 0); 
d=datevec(c);
if d(2) ~= b(2)
   twomonths=1
   monthstr2=num2str(d(2),'%02i');
   strmonth2=strcat(yearstr,monthstr2,'_wind_mag_1min.asc')
end




if doybeg < 100    
  doybegstr=['0',num2str(doybeg)];
else
 doybegstr=num2str(doybeg);
end

if doyend < 100    
  doyendstr=['0',num2str(doyend)];
 else
   doyendstr=num2str(doyend);
end


%monthly mag files 1-min
strwin1=strcat(yearstr,monthstr,'_wind_mag_1min.asc')

%yearly swe files key parameters 2-min
strwin2=strcat('wind_kp_unspike',yearstr,'.txt')

%time=clock;
%time(4:6)

if twomonths==0 

    %%%load mag from the days of interest
    cd('WINDDATA/mag/1min_ascii')
    mfidata = load(strwin1);
    cd ..
    cd ..

    sm=length(mfidata);

    mday=mfidata(:,3);
    
    
    

    indm1=find(mday ==datebeg(3));
    indm2=find(mday ==dateend(3));

    anfang=min(indm1);
    ende=max(indm2);

    myear=mfidata(anfang:ende,1);
    mmonth=mfidata(anfang:ende,2);
    mday=mfidata(anfang:ende,3);
    mhour=mfidata(anfang:ende,4);
    mmin=mfidata(anfang:ende,5);
    msec=mfidata(anfang:ende,6);


    mdata=mfidata(anfang:ende,:); %%%%%%%%%!!!!!!!!!!!!!!!!!!!!


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if twomonths==1

    %%%load mag from the days of interest
    cd('WINDDATA/mag/1min_ascii')
    mfidata1 = load(strwin1);
    mfidata2 = load(strmonth2);
    cd ..
    cd ..
whos mfidata1
whos mfidata2
    sm1=length(mfidata1)
    sm2=length(mfidata2)
%    mdata=zeros(sm1+sm2,22);
    mfidata=zeros(sm1+sm2,22);
    mfidata(1:sm1,:)=mfidata1;
    mfidata(sm1+1:sm1+sm2,:)=mfidata2;
    
    mday=zeros(sm1+sm2,1);

    %convert to datenumbers    
    mday(1:sm1)=datenum(mfidata1(:,1),mfidata1(:,2),mfidata1(:,3),mfidata1(:,4),mfidata1(:,5),mfidata1(:,6))    ;
    mday(sm1+1:sm1+sm2)=datenum(mfidata2(:,1),mfidata2(:,2),mfidata2(:,3),mfidata2(:,4),mfidata2(:,5),mfidata2(:,6))    ;

    indm1=find(mday < DOYbeg);
    indm2=find(mday < DOYend);

    anfang=max(indm1)+1;
    ende=max(indm2);

myear=mfidata(anfang:ende,1);
mmonth=mfidata(anfang:ende,2);
mday=mfidata(anfang:ende,3);
mhour=mfidata(anfang:ende,4);
mmin=mfidata(anfang:ende,5);
msec=mfidata(anfang:ende,6);


mdata=mfidata(anfang:ende,:); %%%%%%%%%!!!!!!!!!!!!!!!!!!!!
%%%%%%%%%%%%%%%%%%%%%


end







%%%load swe from the days of interest

cd('swe_kp_unspike')
swedata = load(strwin2);


sp=length(swedata);
pDOY=swedata(:,2);

%%% DOYs auslesen
 indp1=find(pDOY > doybeg);
 indp2=find(pDOY < doyend);

pdata=swedata(min(indp1):max(indp2),:); %%%%%%%%%!!!!!!!!!!!!!!!!!!!!
pDOY=swedata(min(indp1):max(indp2),2);


disp('Einlesen der Dateien fertig!')
time=clock;
time(4:6)


%%%%%%%%%%%%%%%%%%

cd ..
cd ..


disp('********************************************************')

%%%%%%%%%%%%%%% Physical constants [SI] %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=1.38066*1e-23;
Mp=1.672623*1e-27;
mu0=4*pi*1e-7;
TtoeV=1.29*1e-4;    % aus kT=E, umrechnen T[k] -> E in [eV]
Re=6378;
%%%%%%%%%%%%%%%%%%%% Bad data values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%badmag=1000; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%mag
%1    Year            I4
% 2    Month           I3      01-12
% 3    Day             I3      01-31
% 4    Hour            I3      00-23
% 5    Min             I3      00-59
% 6    Sec             I3      at midpoint of average; typically 30
% 7    Bx,GSE          F10.3   nT  (same as Bx,GSM)
% 8    By,GSE          F10.3   nT
% 9    Bz,GSE          F10.3   nT
%fill values E10.2 and have the value -1.00E+31 
%18    Rx,GSE          F9.3    Re, GSE X component of Wind position vector
%19    Ry,GSE          F9.3    Re, Re = Earth radii
%20    Rz,GSE          F9.3    Re
%fill values E9.1 and have the value -1.0E+31 

Bdata(:,1)=mdata(:,7);%Bx usw. nT
Bdata(:,2)=mdata(:,8);%Bx usw. nT
Bdata(:,3)=mdata(:,9);%Bx usw. nT

xgse=mdata(:,18);%Re GSE
ygse=mdata(:,19);
zgse=mdata(:,20);

%swe
%Wd Format  Fill values      Parameters
% 1  I5     9999             Year
% 2  F12.6  99999.9999999    Decimal Day (day.fractionofday)
% 3  F8.1   99999.9          Flow Speed, km/s
% 4  F8.1   99999.9          Vx, GSE, km/s, flow velocity component
% 5  F8.1   99999.9          Vy, GSE, km/s
% 6  F8.1   99999.9          Vz, GSE, km/s
% 7  F9.0   9999999.         Proton Temperature, K.
% 8  F7.2   999.99           Proton Density, n/cc

vdata(:,1)=pdata(:,4); %vx km/s
vdata(:,2)=pdata(:,5); %vy km/s
vdata(:,3)=pdata(:,6); %vz km/s

%protons
T=pdata(:,7);    %K
N=pdata(:,8);        % 1/ccm
T=T/1e6; 	%T in 10^6 K;


%convert to Matlab datenumber
mDOY=datenum(myear,mmonth,mday,mhour,mmin,msec);
pDOY=pDOY+datenum(year,0,0,0,0,0);

badmag=-1.00E+31; %mag and pos from mag
badvel=99999.9; %swe
badT=9999999;   
badN=999.99;

% interpolate fill values


Bmag=sqrt(Bdata(:,1).^2+Bdata(:,2).^2+Bdata(:,3).^2);% nT, GSE    







[Ib,Jb]=find(Bdata(:,3)~=badmag & Bdata(:,2)~=badmag & Bdata(:,1)~=badmag & abs(Bmag) < 100);
if length(Ib) < length(mDOY)     %Wenn Fehler drinnen -> interpolieren
    Bdata=interp1(mDOY(Ib),Bdata(Ib,:),mDOY,'linear');
else Ib=0;
end
%disp(['bad mag data ',num2str(Ib)])

[Ip,Jp]=find(pdata(:,3)~=badvel);
if length(Ip) < length(pDOY)
   vdata=interp1(pDOY(Ip),vdata(Ip,:),pDOY,'linear');
else Ip=0;
end
   
[Ip,Jp]=find(T~=badT &  N~=badN);
if length(Ip) < length(pDOY)
   N=interp1(pDOY(Ip),N(Ip),pDOY,'linear');
   T=interp1(pDOY(Ip),T(Ip),pDOY,'linear');
else Ip=0;
end
%disp(['bad swe data ',num2str(Ip)])

%disp('1st')
%size(T)
%size(N)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% linear interpolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%e - 92 sec.
%p, B - 3 sec
%1 DOY Vektor mit gewünschtem Abstand -> darauf alles interpolieren!!


DOYvec=[DOYbeg:(secres/(60*60*24)):DOYend]';

%disp('stop')
%%%%%%%%%%% STOP CONDITIONS
if isempty(mDOY) 
    disp('Error: mDOY is zero!')
    return;
   end 
if isempty(Bdata) return;    
    disp('Error: Bdata is zero!')
    return;
end 
   
if isempty(vdata) return;    
    disp('Error: vdata is zero!')
    return;
   end 
if isempty(T) return;    
    disp('Error: T is zero!')
    return;
end 

 if isempty(N) return;    
    disp('Error: N is zero!')
    return;
  end 



Bdata=interp1(mDOY,Bdata,DOYvec,'linear');
vdata=interp1(pDOY,vdata,DOYvec,'linear');
N=interp1(pDOY,N,DOYvec,'linear'); 
T=interp1(pDOY,T,DOYvec,'linear'); 

disp('Interpolieren fertig!')
time=clock;
time(4:6)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



vx=vdata(:,1);               % km/s GSE 
vy=vdata(:,2);               % km/s GSE    
vz=vdata(:,3);               % km/s GSE
V=sqrt(vx.^2+vy.^2+vz.^2);   % km/s GSE 


bx=Bdata(:,1);              % nT, GSE    
by=Bdata(:,2);              % nT, GSE    
bz=Bdata(:,3);              % nT, GSE    
Bmag=sqrt(bx.^2+by.^2+bz.^2);% nT, GSE    

Dt=[1:length(Bmag)]';
Dt(:)=secres;

[YY,MM,DD,hh,mm,sec]=datevec(DOYvec);
%zu Testzwecken
%clouddata=[round(MM) round(DD) round(YY) round(hh) round(mm) round(sec) Dt bx by bz Nnew Tnew vx vy vz Nenew Tenew];
clouddata=[double(MM) double(DD) double(YY) double(hh) double(mm) double(sec) double(Dt) double(bx) double(by) double(bz) double(N) double(T) double(vx) double(vy) double(vz)];% double(Ne) double(Te)];

% 1. und letzter Datenpunkt ist NaN (durch Interpolieren)
%clouddata=clouddata(2:1440,:);


if year > 1999
 yearstrout=num2str(year-2000,'%02i') 
else
 yearstrout=num2str(year-1900);
end

%e.g. winall291_95.mat
filestr=strcat(scstring,'all',doybegstr,'_',yearstrout);
%%%%%%%%%%%%%%%%%%%% DATEN RAUSSCHREIBEN %%%%%%%%%%%%%%%%%%
if write_file == 1
 %outstr=strcat('c:\chris\Matlab\GS-code\DATA\winall', str(5:7),str(12:14),'_',num2str(secres),'.mat')
 outstr=strcat('DATA/',filestr,'.mat')
 save(outstr, 'clouddata');
end;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%T in [K] and  N in [m-3]
N1=N*1e6;    %[SI]
T1=T*1e6;   

pressure=N1*k.*T1;%+Ne1*k.*Te1; %[Pa]
    
%N [cm-3] T [K]
%pressure2=N*1.381.*T*1e-8+Ne*1.381.*Te*1e-8;


Bpre=(Bmag*1e-9).^(2)/(2*mu0);
beta=pressure./Bpre;
 
Bmax=ceil(max(Bmag)/2)*2;   %Plot-Obergrenzen für B,N,T,V
Nmax=ceil(max(N)/10)*10;
Tmax=ceil(max(T)/.1)*.1;
Vmax=ceil(max(V)/100)*100;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DATEN DARSTELLEN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%






figure;
subplot(511)


i=DOYvec;
day1=min(DOYvec);
day2=max(DOYvec);
xtickvec=datenum(year,0,doybeg,0:24:round(day2-day1+1)*24,00,00);
dk=6;

plot(i,Bmag,i,bx,i,by,i,bz);
axis([day1 day2 -Bmax Bmax]);
grid on;
disp('Mean position in Earth Radii (Re), GSE system:')
pos=['x: ',num2str(mean(xgse)),'  y: ',num2str(mean(ygse)), ' z: ',num2str(mean(zgse))]
xlabel(['WIND at mean position (GSE): ',pos]);
ylabel(['B [nT]' ]);
lh=legend('|B|', 'B_x', 'B_y', 'B_z');%, 'Location','Eastoutside');
%set(lh, 'FontSize',6)

%title([' ACE plasma bulk parameters from ',num2str(datebeg(1)),' ',num2str(datebeg(2)),' ',num2str(datebeg(3)),' -  ',num2str(dateend(1)),' ',num2str(dateend(2)),' ',num2str(dateend(3))])
title([' WIND plasma bulk parameters from ',datestr(datebeg),' to ',datestr(dateend)])

set(gca,'XTick',xtickvec,'TickDir','in');
datetick('x',dk,'keeplimits','keepticks');
set(gca,'XTickLabel',[]);



% subplot(612)
% plot(N);axis([1 length(N) 0 Nmax]);

subplot(512)
%axis([1 length(beta) 10^(-3) 10^2]);
semilogy(i, beta)
axis([day1 day2 0.01 10]);
ylabel(['plasma beta ']);
grid on;
set(gca,'XTick',xtickvec,'TickDir','in');
datetick('x',dk,'keeplimits','keepticks');
set(gca,'XTickLabel',[]);


subplot(513)
plot(i,N);
%axis([1 length(Bmag) 0 Vmax]);
axis([day1 day2 0 max(N)+10]);
ylabel('plasma density [#/ccm]');
grid on;
set(gca,'XTick',xtickvec,'TickDir','in');
datetick('x',dk,'keeplimits','keepticks');
set(gca,'XTickLabel',[]);


subplot(514)
plot(i,V);
%axis([1 length(Bmag) 0 Vmax]);
axis([day1 day2 0 Vmax]);
ylabel('plasma bulk V [km/s]');
grid on;

if (day2-day1) < 5
 xtickvec2=datenum(year,0,doybeg,0:6:round(day2-day1+1)*24,00,00);
end
if (day2-day1) >= 5
 xtickvec2=datenum(year,0,doybeg,0:12:round(day2-day1+1)*24,00,00);
end 
if (day2-day1) >= 10
 xtickvec2=0;
end 


set(gca,'XTick',xtickvec2,'TickDir','in');
dk2=15;
datetick('x',dk2,'keeplimits','keepticks');
%set(gca,'XTickLabel',[]);


subplot(515)
plot(i,T);
%axis([1 length(Bmag) 0 Vmax]);
axis([day1 day2 0 Tmax]);
ylabel('proton T [MK]');
set(gca,'XTick',xtickvec,'TickDir','in');
datetick('x',dk,'keeplimits','keepticks');
%set(gca,'XTickLabel',[]);
grid on;




























adjp=0;

%Grenzen einlesen, stehen dann in CASE1.par 

if write_file==1
    
[X,Y]=ginput(2);
%axes(axis3); Y=[1400; -980];
hbound=line([X'; X'], [Y Y]); set(hbound, 'Clipping', 'off','color',[0 0 0]);


i1=min(find(DOYvec>X(1)));
i2=max(find(DOYvec<X(2)));

X=[i1 i2 adjp]        %adjp = 0 initially      


%for directory
str=strcat(scstring,'_', doybegstr,'_',doyendstr,'_',yearstrout)

%for data file: e.g. ACEall316_04.mat
datstr=strcat(scstring,'all',doybegstr,'_',yearstrout,'.mat')



disp('lala***********************')
datstr
str
disp('lala***********************')


if exist(str, 'dir') == 0
    [suc, mess, messid]=mkdir(str);
    suc
    mess
end
cd(str);
disp('-----------------------')
pwd
disp('-----------------------')

%    if exist('CASE1.par', 'file') == 0
        parid=fopen('CASE1.par', 'wt');
        disp('*************filenew************')
        fprintf(parid, '%s\n', ['DATA/',datstr]);
        fclose(parid);
%    end
disp('*************file sowieso************')    
parid=fopen('CASE1.par', 'at');
fprintf(parid, '%i %i %i\n', X);
fclose(parid);

disp('mean position of WIND in the chosen interval (Re GSE):')
meanposx=mean(xgse(X(1):X(2)))
meanposy=mean(ygse(X(1):X(2)))
meanposz=mean(zgse(X(1):X(2)))

parid=fopen('winpos.txt', 'wt');
fprintf(parid, '%f %f %f\n', meanposx,meanposy,meanposz);
fclose(parid);


cd ..


disp('***********************')
str



end
whos clouddata

%%%%%%%%%%%%%%% MAKE IN SITU OUTPUT FOR IDL

%save('insitudata.dat', 'clouddata', '-ascii')

pwd



