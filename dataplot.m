% Plot data from ACE and WIND
% uses first line in ALLCASE.par


disp('--------------- START dataplot.m ------------------')

clear;
format short;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Parameters:

%how many data points are plotted before and after the i1 i2 interval in
%percent
%i1 i2 from CASE1.par

percbefore=40%60  %WIND Nov 20 2003    
percafter=40%90  %----// ------
%voff=-200;  % lower v =vmax-voff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f0=fopen('ALLCASE.par', 'rt');
%dirstr=fscanf(f0, '%s\n')
dirstr=deblank(fgets(f0));
fclose(f0);

%if strcmp(deblank(pwd), 'C:\Qiang\matlab\2001-2002')
cd(dirstr);
%end


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
CASE1=double(clouddata); % Daten werden geladen


s=size(clouddata)

before=i1-floor(percbefore*i1/100);
after=i2+floor(percafter*i2/100);
if before < 1
    before=1;
end
if after > s(1)
    after=s(1);
end   

i10=before;	 %i3=425;	%400 460 425 for 10/30 
i20=after; %5600;


%i3=1257+i10;

%for i=1:length(i12)/3
%i1(i)=i12(1+(i-1)*3);
%i2(i)=i12(2+(i-1)*3);
%out(i)=i12(3+(i-1)*3);
%end
%indout=find(out==1);
%i1=i1(indout); i2=i2(indout);
%i12=[i1' i2'];
%multi=size(i12,1);        % number of intervals to be plotted onto one plot
%i1=i1'; i2=i2';

indout=1;
out=1;
multi=2;

dk=15;

format long;

% adjp=1.3; %adjust p (only for oct 19, 1984)

% adjp=1.6;	% based on the assumption of the total pressure equalibrium across
		% a TD	-- IRM
 adjp=1;	%  0 for force-free and fit Bz; 1.5 for additional electron pressure
 adjbz=1.;	% Pt=p;
 order=1;	% order of the polynomial -- ptol = exp ( poly(A) )

 nden=1e6;  
 nv=1e3; 
 nb=1e-9; 
 k=1.38*1e-23; 
 miu=4.0*pi*1e-7;
 ec=1.6*1e-19;


%[ltime,lT]=issireform(str);
Dtime=CASE1(:,7);
ltime=cumsum(Dtime);


time=ltime(i10:i20)-ltime(i10);
bx=CASE1(i10:i20,8);   %nT
by=CASE1(i10:i20,9);
bz=CASE1(i10:i20,10);

vx=CASE1(i10:i20,13);
vy=CASE1(i10:i20,14);
vz=CASE1(i10:i20,15);
V=sqrt(vx.^2+vy.^2+vz.^2);

if size(CASE1,2)>15 
 Ne=CASE1(i10:i20, 16)*nden;
 Te=CASE1(i10:i20, 17);
pe=adjp*Ne.*Te*k*1e6; 
else
    Te=0;
    pe=0;
end

den=CASE1(i10:i20,11);
%T=CASE1(i1:i2,6);
%t=lT(i1:i2);
t=CASE1(i10:i20,12); %


n=den.*nden;
pp=adjp*n.*t.*k*1e6; 	%temperature is in 10^6 K
%%%%%%%%%%%%%%%%
p=pp+pe;
%%%%%%%%%%%%%%%%
bgse=[bx by bz];
b=bgse.*nb;
bb=sqrt(b(:,1).^2+b(:,2).^2+b(:,3).^2);
	PB=bb.^2/(2*miu);
	beta=p./PB;
pTotal=p+PB;


%%%%
Y=CASE1(:,3); M=CASE1(:,1); D=CASE1(:,2); hh=CASE1(:,4);
mmm=CASE1(:,5);
ssm=CASE1(:,6);
Datem0=datenum(Y,M,D,hh,mmm,ssm);
Datem=Datem0(i10:i20);
%Datei=datenum(Y,M,D,hh,mmi,ssi);




%*******************************************************************

hout=figure;
fpos=get(gcf, 'Position');
set(gcf, 'Position', [fpos(1) fpos(2) fpos(3) fpos(3)]);
left=.1;
bottom=.24;
width=.8;
hei=.15;
nsub=4;
space=0.01;
Date1=Datem(1);
Date2=max(Datem);


Bmax=ceil(max(bb/nb)/10)*10;
%Bmax=15
%Nmax=ceil(max(N)/10)*10;
%Tmax=ceil(max(t+Te)/.1)*.1;
[Ts, Is]=sort(t+Te);    % assuming bad values were marked by negative numbers
Is=find(Ts>0);
Tmin=floor(Ts(Is(1))*1e6/10)*10;
Tmax=ceil(max((t+Te)*1e6)/10)*30;
[Ns, Is]=sort(den);    
Is=find(Ns>0);
Nmin=10^floor(log10(Ns(Is(1))));
Nmax=ceil(max(den)/10)*10+10;
Vmax=ceil(max(V)/100)*100;
Vmin=floor(min(V)/100)*100;
[betas, Is]=sort(beta);    
Is=find(betas>0);
betam=10^floor(log10(betas(Is(1))));
betax=ceil(max(beta)/10)*10;



spacecraft=('Spacecraft?');

coord='GSE';


if dirstr(1:3) == ('sim')
    spacecraft='SIM';
    coord='RTN';
end


if dirstr(1:3) == ('ACE')
    spacecraft='ACE';
    coord='GSE';
end

if dirstr(1:3) == ('ace')
    spacecraft='ACE';
    coord='GSE';
end


if dirstr(1:3) == ('stb') 
    spacecraft='STEREO B';
   coord='RTN';
end

if dirstr(1:3) == ('sta') 
    spacecraft='STEREO A';
    coord='GSE';
   coord='RTN';
end



if dirstr(1:3) == ('win') 
    spacecraft='WIND';
end

if dirstr(1:3) == ('he1') 
    spacecraft='Helios 1';
end

if dirstr(1:3) == ('he2') 
    spacecraft='Helios 2';
end


if dirstr(1:3) == ('vo1') 
    spacecraft='Voyager 1';
        coord='RTN';
end


if dirstr(1:3) == ('vo2') 
    spacecraft='Voyager 2';
        coord='RTN';
end

if dirstr(1:3) == ('IMP') 
    spacecraft='IMP8';
end

if dirstr(1:3) == ('Pi1') 
    spacecraft='Pioneer 11';
    Bmax=1.2;
        coord='RTN';
end

if dirstr(1:3) == ('Uly') 
    spacecraft='Ulysses';
    coord='RTN';
end

if dirstr(1:3) == ('cl1') 
    spacecraft='Cluster 1';
    coord='GSE';
end

if dirstr(1:3) == ('cl3') 
    spacecraft='Cluster 3';
    coord='GSE';
end



plottitel=[spacecraft,' ',num2str(CASE1(i1,1)),'/',num2str(CASE1(i1,2)),'/', ...
  num2str(CASE1(i1,3)),' ',num2str(CASE1(i1,4)),':', ...
 num2str(CASE1(i1,5)),':',num2str(round(CASE1(i1,6))),'-',num2str(CASE1(i2,1)),'/',...
num2str(CASE1(i2,2)),'/',num2str(CASE1(i2,3)),' ',num2str(CASE1(i2,4)), ...
  ':',num2str(CASE1(i2,5)),':',num2str(round(CASE1(i2,6))),'  \Delta t= ',num2str(mean(Dtime)/60),' min'];

if (dirstr(1:2) == ('cl'))
    plottitel=[spacecraft,' ',num2str(CASE1(i1,1)),'/',num2str(CASE1(i1,2)),'/', ...
  num2str(CASE1(i1,3)),' ',num2str(CASE1(i1,4)),':', ...
 num2str(CASE1(i1,5)),':',num2str(round(CASE1(i1,6))),'-',num2str(CASE1(i2,1)),'/',...
num2str(CASE1(i2,2)),'/',num2str(CASE1(i2,3)),' ',num2str(CASE1(i2,4)), ...
  ':',num2str(CASE1(i2,5)),':',num2str(CASE1(i2,6)),'  \Delta t= ',num2str(mean(Dtime)),' sec'];

end





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% B plot

subplot('position', [left bottom+hei*(nsub-1) width hei]);

h1=plot(Datem,bx,'r-',Datem,by,'g-',Datem,bz,'b-', Datem, bb/nb, 'k-');
title(plottitel)
%axis([Date1 Date2 -inf inf]);
axis([Datem(1) Datem(length(Datem)) -Bmax Bmax]);
set(gca, 'fontname', 'times', 'fontangle', 'italic');
if coord=='RTN'
    hleg=legend('R','T','N');
else  
hleg=legend('X','Y','Z');
end
set(h1, 'linewidth', 1);
set(gca, 'fontname', 'times', 'fontangle', 'normal');
ylim=get(gca, 'YLim');
coloreye=[0 0 0; eye(3)];
for i=1:multi
line([Datem0(i12(i,:))'; Datem0(i12(i,:))'], [ylim' ylim'], ...
 'color', coloreye(i,:));
end



%line([Datem0(i12(:,1)'; Datem0(i12(:,2)'], [ylim(1)*ones(size(i2)); ylim(2)*ones(size(i2))], ...
% 'color', [0 0 0]);

% line([Datem0(i1p) Datem0(i1p)], [ylim(1) ylim(2)], ...
%  'color', [0 0 0]);
% line([Datem0(i2p) Datem0(i2p)], [ylim(1) ylim(2)], ...
%  'color', [0 0 0]);
if exist('i3')
line([Datem0(i3) Datem0(i3)], [ylim(1) ylim(2)], ...
 'color', [0 0 0], 'linestyle', '--');
end
ylabel('B (nT)', 'fontname', 'times', 'fontsize', 10);
datetick('x',dk,'keeplimits','keepticks');
grid;
set(gca, 'XTickLabel', []);
%%%%
subplot('position', [left bottom+hei*(nsub-2)-space width hei]);



%%%%%%%%%%%%%%%%%%%%%%%%%%% V PLOT

h1=plot(Datem,vx,'r-',Datem,vy,'g-',Datem,vz,'b-');
%h1=plot(Datem,V,'k-');

%axis([Date1 Date2 -inf inf]);
if coord=='RTN'
   axis([Datem(1) Datem(length(Datem)) Vmin Vmax]);
else  
  axis([Datem(1) Datem(length(Datem)) -Vmax 100]);
end
if dirstr(1:2) == ('cl') 
 axis([Datem(1) Datem(length(Datem)) -Vmax Vmax]);
end

%axis([Datem(1) Datem(length(Datem)) -1200 200]);

set(gca, 'fontname', 'times', 'fontangle', 'italic');
%hleg=legend('V_X','V_Y','V_Z');
set(h1, 'linewidth', 1);
set(gca, 'fontname', 'times', 'fontangle', 'normal');
ylim=get(gca, 'YLim');
for i=1:multi
line([Datem0(i12(i,:))'; Datem0(i12(i,:))'], [ylim' ylim'], ...
 'color', coloreye(i,:));
end
% line([Datem0(i1p) Datem0(i1p)], [ylim(1) ylim(2)], ...
%  'color', [0 0 0]);
% line([Datem0(i2p) Datem0(i2p)], [ylim(1) ylim(2)], ...
%  'color', [0 0 0]);
if exist('i3')
line([Datem0(i3) Datem0(i3)], [ylim(1) ylim(2)], ...
 'color', [0 0 0], 'linestyle', '--');
end
ylabel('Vp (km/s)', 'fontname', 'times', 'fontsize', 10);
datetick('x',dk,'keeplimits','keepticks');
grid;
set(gca, 'XTickLabel', []);
%Ytick=get(gca, 'YTick');
%c=length(Ytick);
set(gca,'YTick',[-1000:100:1000]);
set(gca,'YTickLabel',[-1000:100:1000]);

%set(gca, 'YTick', Ytick(1:c-1));
%%%%


subplot('position', [left bottom+hei*(nsub-3)-2*space width hei]);
%h1=plot(Datem, p/1e-9, 'k-');
%%axis([Date1 Date2 -inf inf]);
%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% T Plot
[AX, h1,h2]=plotyy(Datem, den, Datem, (t+Te)*1e6, @semilogy);
axis(AX(1),[Datem(1) Datem(length(Datem)) Nmin Nmax]);
axis(AX(2),[Datem(1) Datem(length(Datem)) Tmin Tmax]);
set(h1, 'linewidth', 1);
set(h2, 'linewidth', 1);
set(AX, 'fontname', 'times', 'fontangle', 'normal');

set(AX(1),'YTick',[10^-1;5;10;30]);
set(AX(1),'YTickLabel',[10^-1;5;10;30]);
set(AX(1), 'YminorTick', 'on');

%set(AX(1),'YGrid','on');
%set(AX,'YMinorGrid','off');

%set(AX(2),'YTick',[10^5;10^6]);
%set(AX(2),'YTickLabel',[10^5;10^6]);

ylim=get(gca, 'YLim');
line([Datem(1) Datem(length(Datem))], [ylim(2) ylim(2)]);
for i=1:multi
line([Datem0(i12(i,:))'; Datem0(i12(i,:))'], [ylim' ylim'], ...
 'color', coloreye(i,:));
end
% line([Datem0(i1p) Datem0(i1p)], [ylim(1) ylim(2)], ...
%  'color', [0 0 0]);
% line([Datem0(i2p) Datem0(i2p)], [ylim(1) ylim(2)], ...
%  'color', [0 0 0]);
if exist('i3')
line([Datem0(i3) Datem0(i3)], [ylim(1) ylim(2)], ...
 'color', [0 0 0], 'linestyle', '--');
end
ylabel('Np (cm^{-3})', 'fontname', 'times', 'fontsize', 10);
datetick('x',dk,'keeplimits','keepticks');
grid off
;
set(gca, 'XTickLabel', []);
Ytick=get(gca, 'YTick')
Yticklabel=get(gca, 'YTicklabel')
c=length(Ytick);
%set(gca, 'YminorTick', 'on');
%set(gca, 'YTick', [Nmin Ytick Nmax]);
%Yticklabel=num2str(Ytick);
%set(gca, 'YTickLabel', {num2str(Nmin); Yticklabel; num2str(Nmax)});
axes(AX(2)); grid off;
ylabel('Tp (K)', 'fontname', 'times', 'fontsize', 10);
datetick('x',dk,'keeplimits','keepticks');
set(gca, 'XTickLabel', []);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Te Plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Na/Np Plot

%subplot('position', [left bottom+hei*(nsub-5)-4*space width hei]);
%h1=semilogy(Datem, beta, 'k-');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Beta Plot %%%%%%%%%%%%%%%%%%%%



Datem0plot=datenum(Y,M,D,hh,mmm,ssm)-datenum(Y,1,1,0,0,0)+1;
Datemplot=Datem0plot(i10:i20);


subplot('position', [left bottom+hei*(nsub-4)-3*space width hei]);
h1=semilogy(Datemplot, beta, 'k-');

axis([Datemplot(1) Datemplot(length(Datemplot)) betam betax]);
set(h1, 'linewidth', 1);
set(gca, 'fontname', 'times', 'fontangle', 'normal');
ylim=get(gca, 'YLim'); 
for i=1:multi
line([Datem0plot(i12(i,:))'; Datem0plot(i12(i,:))'], [ylim' ylim'], ...
 'color', coloreye(i,:));
end
 %line([Datem0plot(i1p) Datem0plot(i1p)], [ylim(1) ylim(2)], ...
%  'color', [0 0 0]);
 %line([Datem0plot(i2p) Datem0plot(i2p)], [ylim(1) ylim(2)], ...
 % 'color', [0 0 0]);
 
if exist('i3')
line([Datem0plot(i3) Datem0plot(i3)], [ylim(1) ylim(2)], ...
 'color', [0 0 0], 'linestyle', '--');
end
ylabel('\beta', 'fontname', 'times', 'fontsize', 10);
%xlab=['Day of ',num2str(Y(1))];
%xlabel(xlab, 'fontname', 'times', 'fontsize', 10)
%Xtick=get(gca, 'XTick')
%set(gca, 'XTick', [323.2,323.4,324]);

%hxl=xlabel([datestr(Datem0plot(i1(1))) ' - ' datestr(Datem0plot(i2(1)))], ...
%    'fontname', 'times', 'fontsize', 10, 'color', coloreye(1,:));
%    interval=[time2DOY(datevec(Datem0plot(i1(1)))) time2DOY(datevec(Datem0plot(i2(1))))]
%if i>1
%    set(hxl, 'units', 'points');
%    xyl=get(hxl, 'position');
       %set(hxl, 'units', 'data');
%for i=2:multi
%    ht=text(xyl(1), xyl(2)-11*(i-1), [datestr(Datem0plot(i1(i))) ' - ' datestr(Datem0plot(i2(i)))],'units', 'points');
%    set(ht, 'fontname', 'times', 'fontsize', 10, 'color', coloreye(i,:), ...
%        'HorizontalAlignment', 'center', 'VerticalAlignment', 'cap');
%    interval=[time2DOY(datevec(Datem0plot(i1(i)))) time2DOY(datevec(Datem0plot(i2(i))))]
%   end
%end
%datetick('x',dk,'keeplimits','keepticks');grid;
%Ytick=get(gca, 'YTick')
%Yticklabel=get(gca, 'YTicklabel')
%c=length(Ytick);
set(gca, 'YminorTick', 'off');

set(gca,'YTick',[0.01;0.1;1;10]);
set(gca,'YTickLabel',[0.01;0.1;1;10]);

set(gca,'YGrid','on');
set(gca,'YMinorGrid','off');
%set(gca, 'YTick', [betam Ytick betax]);
%set(gca, 'YTickLabel', {num2str(log10(betam)); num2str(log10(Ytick)); num2str(log10(betax))});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Total pressure %%%%%%%%%%%%55

pTotal=pTotal*1e9;
PB=PB*1e9;
p=p*1e9;


%ptotmin=min(pTotal)-(min(pTotal)/10)
ptotmin=0
ptotmax=max(pTotal)+(max(pTotal)/10)


subplot('position', [left bottom+hei*(nsub-5)-4*space width hei]);
h1=plot(Datemplot, pTotal, 'k-', Datemplot, PB,'r-',Datemplot,p,'b-');

axis([Datemplot(1) Datemplot(length(Datemplot)) ptotmin ptotmax]);
set(h1, 'linewidth', 1);
set(gca, 'fontname', 'times', 'fontangle', 'normal');
ylim=get(gca, 'YLim'); 
for i=1:multi
line([Datem0plot(i12(i,:))'; Datem0plot(i12(i,:))'], [ylim' ylim'], ...
 'color', coloreye(i,:));
end
 %line([Datem0plot(i1p) Datem0plot(i1p)], [ylim(1) ylim(2)], ...
%  'color', [0 0 0]);
 %line([Datem0plot(i2p) Datem0plot(i2p)], [ylim(1) ylim(2)], ...
 % 'color', [0 0 0]);
 
hleg=legend('P_{tot}','P_B','P_p');
 
if exist('i3')
line([Datem0plot(i3) Datem0plot(i3)], [ylim(1) ylim(2)], ...
 'color', [0 0 0], 'linestyle', '--');
end
ylabel(' P_{tot} (nPa)', 'fontname', 'times', 'fontsize', 10);
xlab=['Day of ',num2str(Y(1))];
xlabel(xlab, 'fontname', 'times', 'fontsize', 10)



disp('boundaries at DOY:')
Datem0plot(i1)
Datem0plot(i2)





outgraphics=['data' str(6:8) str(12:17)];
set(hout, 'PaperPositionMode', 'manual', 'PaperPosition', [1 2 6 6.8]);
print(hout, '-depsc2', '-cmyk', outgraphics);

cd ..