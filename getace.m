
function getace(year,doybeg,doyend,write_file)
% This program takes merged data from ACE for a given time interval 
% and stores them in  a .MAT file in the correct format in folder \DATA,
% also creates an event folder (e.g. "ACE_345_355_08")
% 
%
%usage:
%getace(year, DOY begin, DOY end, write file yes/no=1/0)
%example: getace(2003, 323, 326, 0)
%
% Christian M�stl, Space Research Institute Graz
% last update November 2010
% based on a code by Qiang Hu



%clear;
format compact;
format short;

%year=2004;
%doybeg=305;
%doyend=310;
%write_file=0             %Falls Datei rausgeschrieben werden soll!!!!!


%use to get DOY from date
datebeg=datenum(year,0,doybeg,0,0,0);
datebeg=datevec(datebeg);

dateend=datenum(year,0,doyend-1,0,0,0);
dateend=datevec(dateend);


disp(['FULL TIME interval is:      ', datestr(datebeg),'    ',datestr(dateend)])

%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd('ACEDATA/4_min_merged_mag_plasma')

scstring='ACE';

a=datenum(year, 0, doybeg, 0, 0, 0); 
b=datevec(a);
if b(2) < 10    
  monthstr=['0',num2str(b(2))];
else
monthstr=num2str(b(2));
end

ndays=doyend-doybeg+1;
yearstr=num2str(year);
%monthstr='11';


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

 filestr=strcat('ace_m',yearstr,'.dat');
 
 f0=fopen(filestr, 'r');
  dat=fscanf(f0, '%f', [20 inf]);
 fclose(f0);
 
 dat=dat';
 
 %%% DOYs auslesen
 DOY=dat(:,2);
 ind1=find(DOY ==doybeg);
 ind2=find(DOY ==doyend);
 
 adata=dat(min(ind1):max(ind2),:); %%%%%%%%%!!!!!!!!!!!!!!!!!!!!
 s=size(adata);

 cd ..
 cd ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
format compact;

%str=strhel(36:49)

adjp=1;


disp('********************************************************')

%%%%%%%%%%%%%%% Physical constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k=1.381*1e-23;
Mp=1.673*1e-27;
mu0=4*pi*1e-7;


Y=adata(:,1);          %1- 1.Jan, bis 365
D=adata(:,2);
H=adata(:,3);          
M=adata(:,4);          
S=zeros(s(1),1);


xgse=adata(:,5)/6378;
ygse=adata(:,6)/6378;
zgse=adata(:,7)/6378;

Bdata=adata(:,9:11); %nT
N=adata(:,14);        % proton density 1/ccm
T=adata(:,15);        % proton temperature  K     
A=adata(:,16);        % alpha ratio   
V=adata(:,17);     %proton speed km/s
vdata=adata(:,18:20); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

last=length(D);
secres=240;


DOYbeg=datenum(Y(1),0,D(1),H(1),M(1),S(1));
DOYend=datenum(Y(1),0,D(last),H(last),M(last),S(last));

DOYvec=[DOYbeg:(secres/(60*60*24)):DOYend]';

DOY=datenum(Y, 0, D, H, M, S); 
% Weil Monat nicht gegeben, 0 setzen -> Days �ber 30 werden auf Monate �bergew�lzt!!

%%%%%%%%%%%%%%%%%%%% Bad data values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
badmag=9999.999;
badpos=badmag;
badv=9999.99;
badN=999.99;
badT=9999999;
bada=9.9999;
%%%%%%%%%%%%% INTERPOLIEREN UND RESAMPLEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bx=Bdata(:,1);

%Interpolieren bei schlechten Daten!!
[Ib,Jb]=find(bx~=badmag);
if length(Ib) < length(DOY)     %Wenn Fehler drinnen -> interpolieren
    bx=interp1(DOY(Ib),bx(Ib,:),DOYvec,'linear');
else 
bx=interp1(DOY,bx,DOYvec,'linear');

end

by=Bdata(:,2);

%Interpolieren bei schlechten Daten!!
[Ib,Jb]=find(by~=badmag);
if length(Ib) < length(DOY)     %Wenn Fehler drinnen -> interpolieren
    by=interp1(DOY(Ib),by(Ib,:),DOYvec,'linear');
else 
by=interp1(DOY,by,DOYvec,'linear');

end

bz=Bdata(:,3);

%Interpolieren bei schlechten Daten!!
[Ib,Jb]=find(bz~=badmag);
if length(Ib) < length(DOY)     %Wenn Fehler drinnen -> interpolieren
    bz=interp1(DOY(Ib),bz(Ib,:),DOYvec,'linear');
else 
bz=interp1(DOY,bz,DOYvec,'linear');

end

vx=vdata(:,1); 
vy=vdata(:,2); 
vz=vdata(:,3); 

[Ip,Jp]=find(V~=badv);
%evt
%[Ip,Jp]=find(V(Ip)~=0);
if length(Ip) < length(DOY)
    V=interp1(DOY(Ip),V(Ip),DOYvec,'linear');
    vx=interp1(DOY(Ip),vx(Ip),DOYvec,'linear');
    vy=interp1(DOY(Ip),vy(Ip),DOYvec,'linear');
    vz=interp1(DOY(Ip),vz(Ip),DOYvec,'linear');
else
   V=interp1(DOY,V,DOYvec,'linear');
   vx=interp1(DOY,vx,DOYvec,'linear');
   vy=interp1(DOY,vy,DOYvec,'linear'); 
   vz=interp1(DOY,vz,DOYvec,'linear'); 
end 


[In,Jn]=find(N~=badN);
if length(In) < length(DOY)
 if length(In) > 0    
     N=interp1(DOY(In),N(In),DOYvec,'linear');
 end
else    
    N=interp1(DOY,N,DOYvec,'linear');
end

[Ia,Ja]=find(A<1);
if length(Ia) < length(DOY)
    A=interp1(DOY(Ia),A(Ia),DOYvec,'linear');
else    
    A=interp1(DOY,A,DOYvec,'linear');
end

[It,Jt]=find(T~=badT);
if length(It) < length(DOY)
    T=interp1(DOY(It),T(It),DOYvec,'linear');
else    
    T=interp1(DOY,T,DOYvec,'linear');
end


T=T/1e6; 	%T in 10^6 K

Bmag=sqrt(bx.^2+by.^2+bz.^2);            % nT, SE    


for i=1:length(DOYvec)-1
Dt(i)=(DOYvec(i+1)-DOYvec(i))*24*3600; %% Delta t zwischen einzelnen Messungen in sec  
end

Dt(length(DOYvec))=round(mean(Dt(1:length(Dt)-1)));
Dt=Dt';
[YY,MM,DD,hh,mm,sec]=datevec(DOYvec);
clouddata=[round(MM) round(DD) round(YY) round(hh) round(mm) sec Dt bx by bz N T vx vy vz]; %Ne Te/1e6];

d=datenum(YY,MM,DD,hh,mm,sec)-datenum(YY,1,1,0,0,0)+1;

%%%%%%%%%% check for uniform sample interval
disp('Mean Delta t in seconds:')
mean(Dt)
%iwrong=(Dt~=240);
%disp('Wrong Delta t:');disp(size(iwrong));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if year < 2000 
    yearstrout=num2str(year-1900);
end
if year > 1999
    yearstrout=strcat('0',num2str(year-2000));
end    
if write_file == 1
 outstr=strcat('DATA/',scstring,'all', doybegstr,'_',yearstrout,'.mat')
 save(outstr, 'clouddata');
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pressure=adjp*N*1e6*k.*T*1e6; %[Pa]
Bpre=(Bmag*1e-9).^2/(2*mu0);
beta=pressure./Bpre;
 
day1=D(1);
day2=D(length(D));
Bmax=ceil(max(Bmag)/2)*2;   %Plot-Obergrenzen f�r B,N,T,V
Nmax=ceil(max(N)/10)*10;
Tmax=ceil(max(T)/.1)*.1;
Vmax=ceil(max(V)/100)*100;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5








%%%%%%%%% DATEN DARSTELLEN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%% B, T, N, V
%figure;
%subplot(311)
%plot(d,Bmag);% clouddata(:,8:10)); 
%axis([day1 day2 0 Bmax]);
%grid on;
%xlabel('Magnetic field |B|')
%ylabel('|B| nT')
%axis([day1 day2 0 Bmax]);
%grid on;
%xlabel('Magnetic field |B|')


%legend('|B|', 'B_X', 'B_Y', 'B_Z');

%subplot(312)
%[ax,h1,h2]=plotyy(d, N, d, T, @semilogy);
%axis(ax(1), [day1 day2 10^(-2) Nmax]);
%axis(ax(2), [day1 day2 10^(-3) Tmax]);
%legend([h1 h2], 'proton density','proton Temp.');
%xlabel('Np & Tp')



%subplot(313)

%plot(d,V,d, clouddata(:,13:15));
%hold on; %plot(V(:,2), '--');
%grid on;
%axis([day1 day2 -50 Vmax]);
%legend('V', 'v_x', 'v_y', 'v_z');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%2nd plot
%%%%%%%%%%%%%%%%%%%%5

figure;
subplot(511)


i=DOY;
day1=min(DOY);
day2=max(DOY);
xtickvec=datenum(year,0,doybeg,0:24:round(day2-day1+1)*24,00,00);
dk=6;

plot(i,Bmag,i,bx,i,by,i,bz);
axis([day1 day2 -Bmax Bmax]);
grid on;
disp('Mean position in Earth Radii (Re), GSE system:')
pos=['x: ',num2str(mean(xgse)),'  y: ',num2str(mean(ygse)), ' z: ',num2str(mean(zgse))]
xlabel(['ACE at mean position (GSE): ',pos]);
ylabel(['B [nT]' ]);
legend('|B|', 'B_X', 'B_Y', 'B_Z');

%title([' ACE plasma bulk parameters from ',num2str(datebeg(1)),' ',num2str(datebeg(2)),' ',num2str(datebeg(3)),' -  ',num2str(dateend(1)),' ',num2str(dateend(2)),' ',num2str(dateend(3))])
title([' ACE plasma bulk parameters from ',datestr(datebeg),' to ',datestr(dateend)])

set(gca,'XTick',xtickvec,'TickDir','in');
datetick('x',dk,'keeplimits','keepticks');
set(gca,'XTickLabel',[]);



% subplot(612)
% plot(N);axis([1 length(N) 0 Nmax]);

subplot(512)
%axis([1 length(beta) 10^(-3) 10^2]);
semilogy(i, beta)
axis([day1 day2 0 3]);
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
plot(i,A);
%axis([1 length(Bmag) 0 Vmax]);
axis([day1 day2 0 max(A)+0.05]);
ylabel('He/H ratio');
%set(gca,'XTick',xtickvec,'TickDir','in');
set(gca,'XTick',xtickvec,'TickDir','in');
datetick('x',dk,'keeplimits','keepticks');
%set(gca,'XTickLabel',[]);
grid on;




%axis3=subplot(513);
%plot(Lambda); axis([1 length(Lambda) 0 360]);

%subplot(512)
%plot(Delta); axis([1 length(Delta) -90 90]);

if write_file==1
    
    
%Grenzen einlesen, stehen dann in CASE1.par 
[X,Y]=ginput(2);
%axes(axis3); Y=[1400; -980];
hbound=line([X'; X'], [-1e8 +1e8]); set(hbound, 'Clipping', 'off','color',[0 0 0]);
%X=[round(X); 0]        %adjp = 0 initially      



[i,x1]=max(DOY.*(DOY < X(1)));
[i,x2]=max(DOY.*(DOY < X(2)));

X=[x1; x2; 0]        %adjp = 0 initially      

disp(['SELECTED (FLUX ROPE) TIME interval is:      ', datestr(DOY(x1)),'    ',datestr(DOY(x2))])


%for directory
str=strcat(scstring,'_', doybegstr,'_',doyendstr,'_',yearstrout)

%for data file: e.g. ACEall316_04.mat
datstr=strcat(scstring,'all',doybegstr,'_',yearstrout,'.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if exist(str, 'dir') == 0
    [suc, mess, messid]=mkdir(str);
    suc;
    mess;
end
    cd(str);
    pwd
    %if exist('CASE1.par', 'file') == 0
        parid=fopen('CASE1.par', 'wt');
        fprintf(parid, '%s\n', ['DATA/',datstr]);
        fclose(parid);
   % end
parid=fopen('CASE1.par', 'at');
fprintf(parid, '%i %i %i\n', X);
fclose(parid);

%disp('Mean position in Earth Radii (Re), GSE system:')
meanposx=mean(xgse(X(1):X(2)));
meanposy=mean(ygse(X(1):X(2)));
meanposz=mean(zgse(X(1):X(2)));

parid=fopen('acepos.txt', 'wt');
fprintf(parid, '%f %f %f\n', meanposx,meanposy,meanposz);
fclose(parid);


%generate file ALLCASE.par in event directory 


cd('..');
end




