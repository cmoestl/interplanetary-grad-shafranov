%function hu34n;
%
%nach h12n_opt -> hu34n.m

%load 'optcloud_var.mat' ;
%load 'hu12n_opt_var.mat' ;



%hu3.m 

%Normalisiert Gr�ssen und fittet Pt(A)


disp('--------------- START hu34n.m ------------------')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('hu3.m')
Re=6378;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% normalization
 A0=max(abs(A1));  %Max A 
 b0=bmax;          %Max B
 L0=A0./b0;
 p0=b0.*b0./miu;   %max B^2/mu0 

dxa1=dxa/L0;     %alle Gr�ssen normieren 
xL1=xL/L0;
Xa1=Xa./L0;
A11=A1./A0;
load Ab.dat;
Ab0=Ab/A0;
bxs1=bxs./b0;    % -s ist Invarianz System
bys1=bys./b0;
bzs1=bzs./b0;
bb1=bb./b0;     %   B^2 normiert
bz12=(bzs1.*bzs1)/(2.*1);
p1=p./p0;%
ptol1=ptol./p0;
Pt1=p1+bb1.*bb1./2.;




%%%%%%%%%%%%%%%%%%%% interpolation

Xi1=Xa1(1):(xL1./(nx-1)):Xa1(rn);
%bxi1=interp1(Xa1,bxs1,Xi1','spline');
%byi1=interp1(Xa1,bys1,Xi1','spline');
%bzi1=interp1(Xa1,bzs1,Xi1','spline');
%ptoli1=interp1(Xa1,ptol1,Xi1','spline');
%pi1=interp1(Xa1,p1,Xi1','spline');
%A2=interp1(Xa1, A11, Xi1', 'spline'); 
%%%%%%%%%%%%%%%%%%%% resampling
bxi1=resample(bxs1, nx, rn);
byi1=resample(bys1, nx, rn);
bzi1=resample(bzs1, nx, rn);
ptoli1=resample(ptol1, nx, rn);
pi1=resample(p1, nx, rn);
A2=resample(A11, nx, rn);

Vxsi=resample(Vxs, nx, rn)';
Vysi=resample(Vys, nx, rn)';
%Vxsi=interp1(Xa1, Vxs, Xi1, 'spline');
%Vysi=interp1(Xa1, Vys, Xi1, 'spline');

nx=length(bxi1); %ceil((nx/rn)*length(Xa1));
%Xi1=Xa1(1):(xL1./(nx-1)):Xa1(rn);

%%%%%%%%%%%%%%%%%%%%

bzi12=(bzi1.*bzi1)/(2.*1);
 ptoli1t=pi1+bzi12;

bbi1=sqrt(bxi1.^2.+byi1.^2.+bzi1.^2.);
Pti1=pi1+bbi1.*bbi1./2.;




%%%%%%%%%%%%%%%%%%%
 dxi1(1)=0;
for i=2:nx
 dxi1(i)=(Xi1(i)-Xi1(i-1));
end
dxi1=dxi1';
%%%%%%%%%%%%%%%%%%%


%integrate A along the satellite trajectory

 dAns(1)=0;
for i=2:(nx)
   dAns(i)=-(byi1(i)+byi1(i-1)).*dxi1(i).*0.5;
end
A1ns=[cumsum(dAns)];
%A2=A1ns';

%%%%%%% sort A2 in ascending order
%[A2 IN]=sort(A2);
%ptoli1=ptoli1(IN);
%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% fit ptol(A) ---- exponential or 2nd order poly.
%% theta=0
% theta=0
%porder=1;


[P,S]=polyfit(A11,ptol1,porder);

i0=1;
Ai0=A2(i0);
i02=length(A2);
Amax=max(A2);
Amin=min(A2);


	 lnptoli1=log(ptoli1);
     
%%CHANGE for fit! next2
%[aptf,S2]=polyfit(A2, lnptoli1, 1);


[aptf,S2]=polyfit([A2(i0:i02)], ptoli1(i0:i02), porder);


residue=sqrt((S2.normr)^2/(i02-i0+1))/(max(ptoli1)-min(ptoli1))
	%aptf=P;
%	if alpha1==0 & alpha2==0 
%		save optres.dat alpha1 alpha2 residue -ascii;
%	end


%%CHANGE for fit! next
ptfit=polyval(aptf,A2);
%	 ptfit=exp(polyval(aptf,A2));

dapt=polyder(aptf);%=3*a3,2*a2,a1
%        dptfit=polyval(dapt,A2);%=dapt(1)*A^2+dapt(2)*A+dapt(3)

%u0=ptfit(i0); u0p=polyval(dapt,A2(i0));
%il0=i0+di0; %length(A2);		%il0=6;
%Al0=A2(il0);
Al0=min(A2)+dAl0;          %.1;   %.5; 
u0=polyval(aptf, Al0);
u0p=polyval(dapt, Al0); %-.1; %.02;
co1=u0p/u0;
co2=log(u0)-co1*Al0;

exptail=exp(co1*[min(A2)-.4:.1:Al0]+co2);

if dirstr(1:3)=='cl5'

    if adjp==0
      clcoeff1=0.01;
      clcoeff2=2.3245;
    else
      clcoeff1=0.1;
      clcoeff2=1.5;  
    end
   exptail=clcoeff1*co1*[min(A2)-.4:.1:Al0]+co2+clcoeff2;
    
end    

Ar0=max(A2)-dAr0; % -.15 
ptr0=polyval(aptf, Ar0);
dptr0=polyval(dapt, Ar0); %.1; %-.02;  % -.1 for another solution
co1r=dptr0/ptr0;
co2r=log(ptr0)-co1r*Ar0;
rexptail=exp(co1r*[Ar0:.1:max(A2)+.4]+co2r);
%%%%%%%%%%%%%%%%%%%%%%%%%%
if adjp == 0
	aptfz=aptf; Al0z=Al0; co1z=co1; co2z=co2; Ar0z=Ar0; co1rz=co1r; co2rz=co2r;
    if dptr0>0
	save bz2fitup aptfz Al0z co1z co2z Ar0z co1rz co2rz;
    else
    save bz2fitdn aptfz Al0z co1z co2z Ar0z co1rz co2rz;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%% diaplay
%[P2,S2]=polyfit(A2,ptoli1,3);
%ptfit=polyval(P2, A2);
%figure;
%subplot(211);
%plot([Xi1]*L0,A2,Xa,A11,'ro');xlabel('x');ylabel('A');title(['\theta=',num2str(theta)]);
%subplot(212);plot([Xi1]*L0,ptoli1,Xa,ptol1,'ro');xlabel('x');ylabel('P_T');
hpta0=figure;
%plot(A11,ptol1,'ro-',A2(1:i0),exptail(1:i0),A2(i0:length(A2)),ptfit(i0:length(A2)));
plot(A2,ptoli1,'ro-', A2(i0:i02), ptfit(i0:i02), [min(A2)-.4:0.1:Al0], exptail, ...
    [Ar0:.1:max(A2)+.4], rexptail);
xlabel('A/A_0');ylabel('P_t /p_0');
%title(['WIND', str]);

%%%%%%%%%%%%%%%%%%%%%%%%
%II=input('Is the optimal invariance direction found? y/n/s [n]:', 's');
%if isempty(II)
%	II = 'n';
%end
%if II=='n'
%	clear;
%	hu12n;
%elseif II=='s'
%	break;
%else
disp('Continue to reconstruct the map... Please wait...');

%end































%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GS Solver und output in files

%hu4.m
%hu4.m

 hx=Xi1(2)-Xi1(1);;%uniform grids
 py=1.0*0.1/1; 
 hy=py*hx;

%nACE=round(abs(ACExyz(2))*Re*1000/L0/hy);

 %ny=nx;
 
 %ny=131; %141	131

 %mid=round(ny/2)-16; %+15 -6; %-16
 mid=round(ny/2)-dmid;

 %ny=mid+round(mid/2);
% mid=101;
%   mid=41;
%  mid=21;

 x=Xi1;

 for j=1:ny,
     y(j)=(j-mid).*hy;
 end


u=zeros(ny, nx);
u(mid,:)=A2';

udy=zeros(ny,nx);

	%%%udy(mid,:)=bxi1';
	udy(mid,:)=bxi1';
	udx(mid, :)=-byi1';

% upper half plane.

nl=1; nr=nx;

 for j=mid:ny-1,

    for i=nl:nr,

if u(j,i)<Al0
	dpt(j,i)=co1*exp(co1*u(j,i)+co2);
    
        %    for Cluster 
%    if dirstr(1:3)=='cl5'
%      dpt(j,i)=clcoeff1*co1*exp(co1*u(j,i)+co2+clcoeff2);
 %   end    
%elseif u(j,i)<Ab0   % SWITCH left right for different events
%    dpt(j,i)=0; %(horizontal)
        
elseif u(j,i)>Ar0
    dpt(j,i)=co1r*exp(co1r*u(j,i)+co2r);
else
       dpt(j, i)=polyval(dapt,u(j,i)); % fit p(A)
%		dpt(j,i)=exp(polyval(aptf,u(j,i)))*polyval(dapt,u(j,i));
end
%          if(dpt(j,i) > 0)
%	if(dpt(j,i) < 0)		%
%		if(u(j,i)>.08)
%          dpt(j,i)=0;
%          end

                rhs(j, i)=-dpt(j, i);

    end

% smoothing
    kk=1.00;
%	yfactor=(y(j)-y(mid))/(y(ny)-y(mid));
    yfactor=min([0.7 (y(j)-y(mid))/(y(ny)-y(mid))]);
%   yfactor=1;% everywhere is the same
    k1=kk*yfactor;k2=3-2*kk*yfactor;k3=kk*yfactor;
%    k1=.5*kk*yfactor;k2=3-kk*yfactor;k3=.5*kk*yfactor;
    for i=nl:nr,
    u1(j,i)=u(j,i);
    udy1(j,i)=udy(j,i);
    if ((i > nl) & (i < nr)) 
    u1(j,i)=(k1*u(j,i+1)+k2*u(j,i)+k3*u(j,i-1))/3;
    udy1(j,i)=(k1*udy(j,i+1)+k2*udy(j,i)+k3*udy(j,i-1))/3;
    end
    end
%%%%smoothing of the end points 
 
%%%    P1=polyfit([x(nr-2) x(nr-1) x(nr)],[u(j,nr-2) u(j,nr-1) u(j,nr)],2);
%%%    P2=polyfit([x(nl+2) x(nl+1) x(nl)],[u(j,nl+2) u(j,nl+1) u(j,nl)],2);
%%%    ul1=polyval(P2,x(nl)-hx);
%%%    ur1=polyval(P1,x(nr)+hx);
%%%
%%%	u1(j,nl)=(k1*u(j,nl+1)+(3-2*k1)*u(j,nl)+k3*ul1)/3;
%%%
%%%	u1(j,nr)=(k1*u(j,nr-1)+(3-2*k1)*u(j,nr)+k3*ur1)/3;
%%%   
%%%%%%%%%%%  

    u(j,:)=u1(j,:);
    udy(j,:)=udy1(j,:);

    for i=nl:nr,

  if (i == nl) 
   duxx=(1./(hx.*hx))*(-5*u(j,i+1)+2.0.*u(j,i)+4*u(j,i+2)-u(j,i+3));
  end 
  if (i == nr)
   duxx=(1./(hx.*hx))*(-5*u(j,i-1)+2.0.*u(j,i)+4*u(j,i-2)-u(j,i-3));
  end 
  if ((i > nl) & (i < nr)) 
        duxx=(1./(hx.*hx))*(u(j,i+1)-2.0.*u(j,i)+u(j,i-1));
  end

        aa(i)=rhs(j,i)-duxx;
   
        aa1(i)=0.5*aa(i).*((y(j+1)-y(j)).^2);
        u(j+1,i)=u(j,i)+udy(j,i).*(y(j+1)-y(j))+aa1(i);
        udy(j+1,i)=udy(j,i)+aa(i)*(y(j+1)-y(j));%Taylor expansion

    end

if j~=mid
    for i=nl:nr
     if i==nl
	    dux=(1./hx)*(u(j,i+1)-u(j,i));
     elseif  i==nr  
	    dux=(1./hx)*(u(j,i)-u(j,i-1));
     else
 	  dux=(1./(2.*hx))*(u(j,i+1)-u(j,i-1));
     end
	
	udx(j,i)=dux;
    end
end

end
    for i=nl:nr,
    u1(ny,i)=u(ny,i);
    udy1(ny,i)=udy(ny,i);
    if ((i > nl) & (i < nr)) 
    u1(ny,i)=(k1*u(ny,i+1)+k2*u(ny,i)+k3*u(ny,i-1))/3;
    udy1(ny,i)=(k1*udy(ny,i+1)+k2*udy(ny,i)+k3*udy(ny,i-1))/3;
    end
	end
	u(ny,:)=u1(ny,:);
	udy(ny,:)=udy1(ny,:);


    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% lower half plane.

nl=1; nr=nx;



for j=mid: -1: 2,

  for i=nl:nr,

    if u(j,i)<Al0
    %%%%%%% take exponential tails if A value is outside left       
    	dpt(j,i)=co1*exp(co1*u(j,i)+co2);

        
      %for Cluster 
      %if dirstr(1:3)=='cl5'
      %  dpt(j,i)=clcoeff1*co1*exp(co1*u(j,i)+co2+clcoeff2);
      %end    
      
       %elseif u(j,i)<Ab0   % SWITCH left right for different events
       %dpt(j,i)=0;            %horizontal dp/da outside 
       elseif u(j,i)>Ar0
       %%%%%%% take exponential tails if A value is outside right       
       dpt(j,i)=co1r*exp(co1r*u(j,i)+co2r);
      
      else
      %%%%%%% take polynomail if A value is inside       
       dpt(j, i)=polyval(dapt,u(j,i)); % fit p(A)
       % dpt(j,i)=exp(polyval(aptf,u(j,i)))*polyval(dapt,u(j,i));
    end

           %   if (dpt(j,i) > 0)
           %	if (dpt(j,i) < 0)		%
           %	if(u(j,i)>.08)
           %          dpt(j,i)=0;
            %         end
    %set right hand side
     rhs(j, i)=-dpt(j, i);
  end
   
    kk=1.00;
%	yfactor=(y(j)-y(mid))/(y(1)-y(mid));
    yfactor=min([0.7 (y(j)-y(mid))/(y(1)-y(mid))]);
%   yfactor=1;
   k1=kk*yfactor;k2=3-2*kk*yfactor;k3=kk*yfactor;
%     k1=.5*kk*yfactor;k2=3-kk*yfactor;k3=.5*kk*yfactor;

   for i=nl:nr,
    u1(j,i)=u(j,i);
    udy1(j,i)=udy(j,i);
     if ((i > nl) & (i < nr)) 
      u1(j,i)=(k1*u(j,i+1)+k2*u(j,i)+k3*u(j,i-1))/3;
      udy1(j,i)=(k1*udy(j,i+1)+k2*udy(j,i)+k3*udy(j,i-1))/3;
     end
    end
%%%%smoothing of the end points 
 
%%   P1=polyfit([x(nr-2) x(nr-1) x(nr)],[u(j,nr-2) u(j,nr-1) u(j,nr)],2);
%%    P2=polyfit([x(nl+2) x(nl+1) x(nl)],[u(j,nl+2) u(j,nl+1) u(j,nl)],2);
%%    ul1=polyval(P2,x(nl)-hx);
%%    ur1=polyval(P1,x(nr)+hx);
%%
%%	u1(j,nl)=(k1*u(j,nl+1)+(3-2*k1)*u(j,nl)+k3*ul1)/3;
%%
%%%%	u1(j,nr)=(k1*u(j,nr-1)+(3-2*k1)*u(j,nr)+k3*ur1)/3;
   
%%%%%%%%%%%  

    u(j,:)=u1(j,:);
    udy(j,:)=udy1(j,:);

  for i=nl:nr,
    %%%%%% end points special treatment
     if (i == nl) 
       duxx=(1./(hx.*hx))*(-5*u(j,i+1)+2.0.*u(j,i)+4*u(j,i+2)-u(j,i+3));
     end
     if (i == nr) 
       duxx=(1./(hx.*hx))*(-5*u(j,i-1)+2.0.*u(j,i)+4*u(j,i-2)-u(j,i-3));
     end
     %%% non end points
     if ((i > nl) & (i < nr)) 
        duxx=(1./(hx.*hx))*(u(j,i+1)-2.0.*u(j,i)+u(j,i-1));
     end

        aa(i)=rhs(j,i)-duxx;

        aa1(i)=0.5*aa(i).*((y(j-1)-y(j)).^2);
        u(j-1,i)=u(j,i)+udy(j,i).*(y(j-1)-y(j))+aa1(i);
        udy(j-1,i)=udy(j,i)+aa(i)*(y(j-1)-y(j));%Taylor expansion

  end 

      
   if j~=mid
     for i=nl:nr
      if i==nl
	     dux=(1./hx)*(u(j,i+1)-u(j,i));
      elseif  i==nr  
	     dux=(1./hx)*(u(j,i)-u(j,i-1));
      else
 	     dux=(1./(2.*hx))*(u(j,i+1)-u(j,i-1));
      end
	
   	 udx(j,i)=dux;
     end
   end

end

 for i=nl:nr,
    u1(1,i)=u(1,i);
    udy1(1,i)=udy(1,i);
  if ((i > nl) & (i < nr)) 
    u1(1,i)=(k1*u(1,i+1)+k2*u(1,i)+k3*u(1,i-1))/3;
    udy1(1,i)=(k1*udy(1,i+1)+k2*udy(1,i)+k3*udy(1,i-1))/3;
  end
end
u(1,:)=u1(1,:);
udy(1,:)=udy1(1,:);


    
    
dpt(ny,:)=dpt(ny-1,:);
dpt(1,:)=dpt(2,:);
    
%%%%%%%%%%%%%%%%%%%%%%
AU=1.49e11; %m
Xi=Xi1.*L0/AU;
bxi=bxi1.*b0;
byi=byi1.*b0;
bzi=bzi1.*b0; 
bbi=bbi1.*b0; 
pi=pi1.*p0; 
xe=x'.*L0/AU;
ye=y'.*L0/AU;
Ae=u.*A0;
Bx=udy.*b0/1e-9;
By=-udx*b0/1e-9;
yL=L0*(y(ny)-y(1))/AU;
xL=xL/AU;
%%%%%%%%%%%%%%%%%%%%%%
Xa=Xa/AU;
umin=min(u(:));
umax=max(u(:));

%%%%%%%%% total pressure adjustment (Hu and Sonnerup 2003)
%p0adj=PTol(1);
%load Ab.dat;
%chi=1;
%for i=2:ny-1,
 %   for j=nl:nr,
       %Bx -> udy, A-> u
        
  %    if  Ae(i,j) < Ab   %<> je nach Event checken
  %        padj=((Bx(i,j)*1e-9)^2+(By(i,j)*1e-9)^2)/(2*miu);  %no Bz for now
  %        deltap=abs(padj - p0adj);
  %        
  %       if deltap > chi*p0adj
  %          vorz=sign(p0adj - padj)*sign(Bx(i,j));
  %          Bx(i,j)=(Bx(i,j)*1e-9+(sqrt(2*chi*p0adj*miu))*vorz*abs(j-mid)/(ny/3) )*1e9;
  %          By(i,j)=(By(i,j)*1e-9+(sqrt(2*chi*p0adj*miu))*vorz*abs(j-mid)/(ny/3) )*1e9;
  %          Ae(i,j)=Ae(i,j)+(sqrt(2*chi*p0adj*miu))*vorz*abs(j-mid)/(ny/3);                      
  %       end   
  %   end 
 
%    end
%end

    










if dirstr(1:2)=='cl'
 xe=xe.*AU*1e-3;
 ye=ye.*AU*1e-3;
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%% fit Bz
if dptr0 > 0
load bz2fitup;
else
    load bz2fitdn;
end
for j=1:ny
	for i=1:nx
		if u(j,i) < Al0z
			Bz2(j,i)=exp(co1z*u(j,i)+co2z);
            
            if dirstr(1:3)=='cl5'
                    %    for Cluster 
                  Bz2(j,i)=clcoeff1*co1*exp(co1*u(j,i)+co2+clcoeff2);
           end  
             
        elseif u(j,i) > Ar0z
            Bz2(j,i)=exp(co1rz*u(j,i)+co2rz);
        else
			Bz2(j,i)=polyval(aptfz, u(j,i));
		end
	end
end
Bz=sqrt(Bz2*2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% determine the boundary (e-fold of the max.Bz)
%if abs(min(Bz))> abs(max(Bz) ....   %dann zeigt Bz in die Ebene
Bzmax=max(Bz(:));
[J0,I0]=find(Bz==Bzmax);
X0=xe(I0)
Y0=ye(J0)
% Bzbound=Bzmax/exp(1);
% Bzb2=Bzbound^2/2;
% 
% u1exp=linspace(umin, Al0z, 200);
% u2poly=linspace(Al0z, umax, 200);
% 
% Bz21exp=exp(co1z*u1exp+co2z);
% Bz22poly=polyval(aptfz, u2poly);
% 
% Jb1=find(abs(Bz21exp-Bzb2)<0.0005)
% Jb2=find(abs(Bz22poly-Bzb2)<0.0005)
% Jb3=find(abs(Bz22poly-(2.18*Bzbound)^2/2)<0.0007)
% Jb4=find(abs(Bz22poly-(2.5*Bzbound)^2/2)<0.00055)
% 
% if isempty(Jb1) 
% 	uB=u2poly(Jb2);
% else
% 	uB=u1exp(Jb1);
% end
% uB2=[u2poly(Jb3) u2poly(Jb4)];
% Ab2=uB2*A0;	
% 
%
if get_Ab==0
    adjp
    'Select Ab...'
figure(hpta0);
[uB, dum]=ginput(1);
Ab=uB*A0;
save Ab.dat Ab -ascii;
    'Ab.dat saved!'
end
load Ab.dat;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[maxbz,inz]=max(abs(bzs));
Bz=(bzs(inz)/maxbz)*Bz*b0/1e-9; % nT

Bz=double(Bz);
Bzmax=Bz(J0, I0)

Re=6378;
J00=b0/(miu*L0);
Jz=dpt*J00;        
%max_Jz=max(abs(Jz(:)))
%[J2, I2]=find(abs(Jz)==max_Jz)
Jz00=Jz(J0, I0)

















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAP %%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hout(1)=figure;

fpos=get(gcf, 'Position');
set(gcf, 'Position', [fpos(1) fpos(2) fpos(3) fpos(3)]);
%%%%subplot(311);
%subplot('position', [.13 .31 .66 .7])
%contourf(xe,ye,Bz,[10,20,30,40,50]); 
%shading flat;

%xLim=max(Bz(:))+max(Bz(:))/10;
%xe=-xe;
%ye=-ye;


contourf(xe,ye,Bz,50); 
shading flat;
colormap jet;


spacecraft=('Spacecraft?');
if dirstr(1:3) == ('ace') | ('ACE')
    spacecraft='ACE';
end

if dirstr(1:3) == ('stb') 
    spacecraft='STEREO B';
end

if dirstr(1:3) == ('sta') 
    spacecraft='STEREO A';
end

if dirstr(1:3) == ('win') 
    spacecraft='WIND';
end
if dirstr(1:3) == ('he1') 
    spacecraft='HELIOS 1';
end

if dirstr(1:3) == ('he2') 
    spacecraft='HELIOS 2';
end

if dirstr(1:3) == ('vo1') 
    spacecraft='Voyager 1';
end


if dirstr(1:3) == ('vo2') 
    spacecraft='Voyager 2';
end

if dirstr(1:3) == ('IMP') 
    spacecraft='IMP8';
end

if dirstr(1:3) == ('Pi1') 
    spacecraft='Pioneer 11';
end

if dirstr(1:3) == ('Uly') 
    spacecraft='Ulysses';
end

if dirstr(1:3) == ('cl1') 
    spacecraft='Cluster 1';
end

if dirstr(1:3) == ('cl2') 
    spacecraft='Cluster 2';
end

if dirstr(1:3) == ('cl3') 
    spacecraft='Cluster 3';
end

if dirstr(1:3) == ('cl4') 
    spacecraft='Cluster 4';
end


plottitel=[spacecraft,' ',num2str(CASE1(i1,1)),'/',num2str(CASE1(i1,2)),'/', ...
  num2str(CASE1(i1,3))];

if dirstr == ('cl3_234_234_01') 
    plottitel=('C3 22/08/2001 10:08:33-45 UT');
end

if dirstr == ('sta_141_144_07') 
    plottitel=('MC STEREO A   23 May 2007');
end

title(plottitel, 'fontsize',16);


hold on;

contour(xe,ye,Ae,40,'k'); 
axis('equal'); 
line(X0, Y0, 'marker','.','markersize',12,'color',[1 1 1]);
Lim=ceil(max(max(Bz))/5)*5;

set(gca,'CLim',[0 Lim])
hbar = colorbar('east',... 
  'Box','on',...
  'CLim',[0 Lim],...
  'YAxisLocation','right',...
  'FontSize',12,...
  'XLim',[-0.5 1.5],...
  'Location','manual');

%max(max(Bz))
%break;

axis([xe(1) xe(length(xe)) ye(1) ye(length(ye))]); 
set(gca, 'fontsize', 12);


% Ab=A2(length(A2))*A0;
hold on;
[cBound, hBound]=contour(xe,ye,Ae,[Ab Ab],'w-');
set(hBound, 'linewidth', 4);

set(gca, 'TickDir', 'out', 'XMinorTick', 'on', 'YMinorTick', 'on');



xlabel('x� (AU)', 'fontsize', 16);
ylabel('y� (AU)', 'fontsize', 16);
hcon=get(gcf, 'CurrentAxes');


hold on;
%arrows along spacecraft trajectory
hbt=quiver([xe' xe(5)], [zeros(size(xe')) ye(length(ye)-10)],[bxi1' 10e-9/b0],[byi1' 0],0.7);   
text(xe(2), ye(length(ye)-10), '10 nT','fontsize', 8, 'color', [1 1 0]);
set(hbt, 'linewidth', 1.5, 'color', [1 1 0]);

hvt=quiver([xe' xe(5)], [zeros(size(xe')) ye(length(ye)-15)], [Vxsi/1000 50], [Vysi/1000 0], 1);
set(hvt, 'linewidth', 1.5, 'color', [0 0.8 0]);   
text(xe(2), ye(length(ye)-15), '50 km/s','fontsize', 8, 'color', [0 0.8 0]); 

cd(dirstr)
save('variables3D', 'Bz', 'Ae', 'xe', 'ye');   
cd ..
%hvt=quiver([xe' xe(length(xe)-2)], [zeros(size(xe')) ye(length(ye)-10)], [bxi/1000 100], [byi/1000 0], .4);
%hvt=quiver([xe' xe(length(xe)-2)], [zeros(size(xe')) ye(length(ye)-10)], [Vxsi/1000 100], [Vysi/1000 0], .4);
%if dirstr(1:2)~='cl'
% hvt=quiver([xe' xe(length(xe)-2)], [zeros(size(xe')) ye(7)], [Vxsi/1000 100], [Vysi/1000 0], .1);
% set(hvt, 'linewidth', 1.0, 'color', [0 1 0]);
%end 

if dirstr(1:2)=='cl'
 xlabel('x� (km)', 'fontsize', 10);
 ylabel('y� (km)', 'fontsize', 10);
 hcon=get(gcf, 'CurrentAxes');
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%COLORBAR
pos=get(hcon, 'Position')
aspr=abs((ye(length(ye))-ye(1))/(xe(length(xe))-xe(1)));
barph=aspr*pos(3);
barpy=pos(2)+(pos(4)-aspr*pos(3))/2;
set(hbar, 'Position', [pos(1)+pos(3)+.02 barpy+0.05 pos(3)/40 barph-0.15], 'unit', 'normalized');


tx=(pos(3)+.02+pos(3)/30)/pos(3)*(xe(length(xe))-xe(1));
tx2=(pos(3)+.02)/pos(3)*(xe(length(xe))-xe(1));
text(tx2, ye(1)-.05*(ye(length(ye))-ye(1)), 'B_z (nT)', 'fontsize', 8, 'fontangle','o');



%hold on;
%contour(xe,ye,Ae,40,'k');

% Ab=A2(length(A2))*A0;

% 
% Ax=A2(10)*A0+1;                   % X-point
% [cBound2, hBound2]=contour(xe,ye,Ae,[Ax Ax],'m-');
% set(hBound2, 'linewidth', 1.5);
% xloc=xe(10)/xe(length(xe));
% time_x=xe(10)*AU/vxs/(24*3600)+Datem(1);
% date_x=datevec(time_x)


%%%%%%%%%%%%%%%%%%%%%%% Projecting second spacecraft example %%%%%%%%%%%%5

if dirstr=='ACE_076_082_01'
 cd('C:\chris\Matlab\GS-code\win_076_082_01')
 load winpos.txt;
 cd('C:\chris\Matlab\GS-code\ACE_076_082_01')
 load acepos.txt;
 cd .. 
 winposAU=winpos*Re/149.5e6;
 aceposAU=acepos*Re/149.5e6;
   sepvec=(winposAU-aceposAU)*[xs' ys' zs']   %zeigt zu Wind wenn so definiert!
   yWIN=sepvec(2)
   line([xe(1) xe(length(xe))], [yWIN yWIN], 'color', [1 1 0]);
    %%%% sepvec projeziert in reconstr- coord. zeigt alle Abst�nde.... 
    
end

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    



%%%%%%%%%%% Now calculate the content of the flux rope bounded by the box
Bz=abs(Bz);
flux=0.0;
current=0.0;
area=0;
garea=hx*L0*hy*L0;
for j=1:length(ye)
    for i=1:length(xe)
        if get_Ab==1
          if Ae(j,i) > Ab
            flux=flux+Bz(j,i)*1e-9*garea;
            current=current+Jz(j,i)*garea;
            area=area+garea;
          end
      elseif get_Ab==-1
          if Ae(j,i) < Ab
            flux=flux+Bz(j,i)*1e-9*garea;
            current=current+Jz(j,i)*garea;
          end
      else
            flux=flux+Bz(j,i)*1e-9*garea;
            current=current+Jz(j,i)*garea;
        end
    end
end
flux=abs(flux);
current=abs(current);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% FOR MAGNETIC CLOUDS

if dirstr(1:2)~='cl'
    
disp('poloidal flux in 10^21 Mx:')

  
garea;
%hx
%hy
%L0
hxinAU=L0*hx/1.496e11;
hyinAU=L0*hy/1.496e11; 

L_in_AU=1
L=L_in_AU*AU;
Am=Ae(J0,I0)
Ab
axflux=flux*1e8/1e21
polflux=abs(Am-Ab)*L*1e8/1e21
disp(['complete axial current in MA:'])
  current*1e-6
 

end


%%%%UNITS FOR CLUSTER
if dirstr(1:2)=='cl'

  Am=Ae(J0,I0);
  disp(['axial flux in Mx:'])
  axflux=flux*1e8       

  L=Re*1e3;
  Am=Ae(J0,I0);
  disp(['poloidal flux in Mx/Re'])
  polflux=abs(Am-Ab)*L*1e8   % per Re

  disp(['poloidal flux in T m = Wb/m'])
  polflux_unitlength=abs(Am-Ab)   % per m
  % Slavin 15-19 nA /m^2 "current density at the central axis"
  disp(['complete axial current in MA:'])
  current*1e-6
  disp(['axial current in nA/m^2:'])
  mean_current_density=current/(area)*1e9
  axial_current_density=abs(Jz(J0,I0))*1e9
    
end    

%%%%%%%%%%%%%%%%%%%%%%%%%% calculate the distance from Boundary to the center
%%%%% 2 segment boundary
%ns1=cBound(2,1);
%Boundpoints=[cBound(:,2:ns1+1) cBound(:,ns1+3:size(cBound,2))];
%radius=sqrt((Boundpoints(1,:)-X0).^2+(Boundpoints(2,:)-Y0).^2);

%aradi=mean(radius)/1.49e8	% AU in km

%%%%%%%%%%%%%%%%%%%%%
uv=.017;
GSE=eye(3);
Local=[xs' ys' zs'];
GSE_xyz=GSE*Local;
x0=xe(4);
y0=ye(length(ye)-20);

%%%%%%%%%%%%%    output;

%  line([x0 x0+uv*GSE_xyz(3,1)], [y0 y0+uv*GSE_xyz(3,2)], 'marker', '.', ...
%  'color',[1 .7 0],'linewidth', 1.5, 'markersize', 6);
%  line([x0 x0+uv*GSE_xyz(2,1)], [y0 y0+uv*GSE_xyz(2,2)], 'marker', '.', ...
%  'color',[0 1 0],'linewidth', 1.5, 'markersize', 6);
%  line([x0 x0+uv*GSE_xyz(1,1)], [y0 y0+uv*GSE_xyz(1,2)], 'marker', '.', ...
%  'color',[1 0 0],'linewidth', 1.5, 'markersize', 6);

% load normal.dat;
% sn=normal(1,:)*[xs' ys' zs']
% ang_zn=acos(dot(normal(1,:),zs))*180/3.14159265
% hn=quiver([x0 x0], [y0 y0], [sn(1) sn(1)], [sn(2) sn(2)], .045);
% hn=quiver([x0 x0 x0 x0], [y0 y0 y0 y0], [sn(1) GSE_xyz(1,1) GSE_xyz(2,1) GSE_xyz(3,1)], ...
% [sn(2) GSE_xyz(1,2) GSE_xyz(2,2) GSE_xyz(3,2)], .05);
%set(hn, 'linewidth', 2, 'color', [1 1 0]);
% set(hn(2), 'linewidth', 1.5, 'color', [1 0 0]);
% set(hn(3), 'linewidth', 1.5, 'color', [0 1 0]);
% set(hn(4), 'linewidth', 1.5, 'color', [1 .7 0]);
% sn=normal(2,:)*[xs' ys' zs']
% ang_zn=acos(dot(normal(2,:),zs))*180/3.14159265
% hold on; hn=quiver([xe(length(xe)-4) xe(length(xe)-4)], [ye(20) ye(20)], [sn(1) sn(1)], [sn(2) sn(2)], .045);
% % hn=quiver([x0 x0 x0 x0], [y0 y0 y0 y0], [sn(1) GSE_xyz(1,1) GSE_xyz(2,1) GSE_xyz(3,1)], ...
% % [sn(2) GSE_xyz(1,2) GSE_xyz(2,2) GSE_xyz(3,2)], .05);
% set(hn, 'linewidth', 2, 'color', [1 1 0]);


% subplot('position',[pos(1) pos(2)-.2 pos(3) .2]);
% [ax,h1,h2]=plotyy(Datem, bb/1e-9, Datem, vv);set(gca,'Xticklabel',[]);
% grid on;  datetick('x',15);
% Bmax=ceil(max(bb/1e-9)/10)*10;
% Vmax=ceil(max(vv)/100)*100;
% axis(ax(1),[Datem(1) Datem(length(Datem)) 0 Bmax]);%set(ax(1),'Xticklabel',[]);
% axis(ax(2),[Datem(1) Datem(length(Datem)) Vmax/2 Vmax]);set(ax(2),'Xticklabel',[]);
% xlabel(str);
% axes(ax(1)); ylabel('B (nT)');
% axes(ax(2));ylabel('V (km/s)');


%x_GSE=[1 0 0];
%hold on;
%hx_GSE=quiver([xe(length(xe)-3) 2*xe(length(xe))],[ye(length(ye)-10) 1*1e3], ...
%[x_GSE*xs' x_GSE*xs'],[x_GSE*ys' x_GSE*ys'],.1,'w-');
%set(hx_GSE, 'linewidth', 4);
%%%%%%%%%gtext('n');gtext('n');

%contourf(xe,ye,Ae);

%load byh.dat
%ya1=byh;

%quiver([.5*1e4 3*1e4],[0*1e2 1*1e3],[ya1*xs' ya*xs'],[ya1*ys' ya*ys'],.1,'k-');
%gtext('n');gtext('n');


%subplot(312);
%plot(Xa', bgse*ya1');ylabel('B_n (nT)');%axis([ixe1 ixe2 -inf inf]);
%axis([xe(1) xe(length(xe)) -inf inf]);
%grid;
%xlabel('x (km)');ylabel('B_n (nT)');

%subplot(313);
%plot(Xa', bxs/1e-9, ...
%'o-',Xa',bys/1e-9,'+-',Xa',bzs/1e-9,'x-',Xa',sqrt(bxs.^2+bys.^2+bzs.^2)/1e-9);
%axis([xe(1) xe(length(xe)) -inf inf]);
%grid;%axis([ixe1 ixe2 -inf inf]);
%legend('B_x','B_y','B_z','B');ylabel('nT');

%subplot(313);
%quiver3([0 0.5],[0 0.5],[0 0.5],[xs(1) ys(1)],[xs(2),ys(2)],[xs(3) ys(3)]);axis([-1 1 -1 1 -1 1]);axis('equal');
%view([1 0 0]);
%quiver([0],[0],[xs(2)],[xs(3)]);axis('equal');axis([-1 1 -1 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
ul=linspace(umin, Al0,20);
ur=linspace(Ar0, umax, 20);
uall=linspace(Al0, Ar0);
ptu=polyval(aptf, uall);
ptl=exp(co1*ul+co2);
if dirstr(1:3)=='cl5'
      ptl=clcoeff1*co1*ul+co2+clcoeff2;
end      
 

ptr=exp(co1r*ur+co2r);

%%%ul=linspace(umin, Ai0);
%%%ur=linspace(Ai0, umax);
%%%ptl=polyval(aptf, ul);
%%%ptr=exp(co1*ur+co2);























%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Pt(A) plot
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

Im=find(abs(A2)==max(abs(A2)));
hout(2)=figure;
hdata=plot(A2(1:Im)*A0, ptoli1(1:Im)*p0/1e-9, 'o--');
set(hdata, 'markersize', 8); hold on;
hdata=plot(A2(Im+1:length(A2))*A0, ptoli1(Im+1:length(A2))*p0/1e-9, '*--');
set(hdata, 'markersize', 8); 
set(gca, 'fontsize', 12);


hold on;plot(uall*A0, ptu*p0/1e-9,'k-', 'linewidth', 2);
hold on;plot(ul*A0, ptl*p0/1e-9, 'k--', 'linewidth', 2);
ulh=ul(1)*A0:0.001:Ab;
ph=zeros(length(ulh),1);
ph(:)=ptu(11)*p0/1e-9;
%hold on;plot(ulh, ph, 'k-', 'linewidth', 2);
hold on;plot(ur*A0, ptr*p0/1e-9, 'k--', 'linewidth', 2);
xlabel('A (T \cdot m)', 'fontsize', 16);
ylabel('P_t (nPa)', 'fontsize', 16);


if dirstr=='cl3_234_234_01'
 plottitel=('C3 P_t(A)') ;
end    

if dirstr=='sta_141_144_07'
   plottitel=('MC STEREO A   23 May 2007')
end


title(plottitel,'fontsize',16);




%%%%%%%%%%%%%%%%%%%%%%


  hRf=text(((umax-umin)*0.27+umin)*A0, max(ptu)*p0/1e-9, ['R_f=' num2str(residue, '%4.2f')]);
  set(hRf, 'fontsize', 18);


if dirstr=='stb_141_144_07'
    
  hRf=text(((umax-umin)*0.27+umin)*A0, max(ptu)*p0/1e-9, ['R_f=' num2str(residue, '%4.2f')]);
  set(hRf, 'fontsize', 18);
end  

if dirstr=='stb_142_144_07'
    
  hRf=text(((umax-umin)*0.27+umin)*A0, max(ptu)*p0/1e-9, ['R_f=' num2str(residue, '%4.2f')]);
  set(hRf, 'fontsize', 18);
end  


if dirstr=='sta_141_144_07'
    
  hRf=text(((umax-umin)*0.27+umin)*A0, max(ptu)*p0/1e-9, ['R_f=' num2str(residue, '%4.2f')]);
  set(hRf, 'fontsize', 18);
end  



%%%title(['\theta =', num2str(theta)], 'fontsize', 18);
%%axis([-1 .5 .025 .05]);
 grid on;

 %output;

ylim=get(gca, 'Ylim');
line([Ab Ab], [ylim(1) ylim(2)], 'color', [0 0 0]);
text(Ab, ylim(2)-(ylim(2)-ylim(1))*.2, 'A_b');


%%%%
%%% figure;
%%% fpos=get(gcf, 'Position');
%%% set(gcf, 'Position', [fpos(1) fpos(2) fpos(3) fpos(3)]);
%%% %subplot(211);
%%% hbt=quiver([xe' xe(2) 2e7*ones(1, length(xe))], [zeros(size(xe')) 6e6 6e6*ones(1, length(xe))],[bxi1' 5e-9/b0 byi1'],[bzi1' 0 bzi1'],0.6);
%%% %axis([xe(1) xe(length(xe)) -1e6 8e6]);
%%% axis equal; %grid on;
% set(hbt, 'linewidth', .5);
% line([1.6e7 2.4e7], [6e6 6e6], 'linewidth', 0.5);
% line([2e7 2e7], [3e6 10e6], 'linewidth', 0.5);
% %line([xe(1) xe(length(xe))], [0 0], 'linewidth', 0.5);
% line([xe(1) xe(1)], [-3e6 3e6], 'linewidth', 0.5);
% axis off;
% %subplot(212);
% hold on;
% [Ax, H1, H2]=plotyy(xe, zeros(size(xe)), xe, atan(sqrt(bxi1.^2+byi1.^2)./bzi1)*180/3.14159265);
% %H2=plot(xe, atan(sqrt(bxi1.^2+byi1.^2)./bzi1)*180/3.14159265);
% set(H2, 'marker', 'o', 'linewidth', 1, 'linestyle', '--', 'color', [.5 .5 .5]);
% %axis(Ax(2));
%xlabel('x (km)'); ylabel('pitch angle (^\circ)');
% %axis equal; 
% grid on;
% %axis([xe(1) xe(length(xe)) 0 90]);

cd(dirstr)

if adjp==1 
    parid=fopen('CASE1_good.par', 'at');
    fprintf(parid, '%s\n', '   ');
    fprintf(parid, '%s\n', str);
    fprintf(parid, '%s\n', [int2str(i12(1:2)') '  ' int2str(adjp) '  %hu12n good']);
    fprintf(parid, '%f %f %f\n', zs');
    fprintf(parid, '%s\n', 'ave_beta, Mta, Mta2, Rf, Y0, Bzmax, Jzmax, flux, current are:');
    fprintf(parid, '%f %f %f %f %f %f %e %e %e\n', ave_beta, Mta, Mta2, residue, Y0, Bzmax, Jz00, flux, current); 
    fclose(parid);
end






outgraphics=['map' str(6:8) str(12:17); 'pta' str(6:8), str(12:17)];
for i=1:length(hout)
    set(hout(i), 'PaperPositionMode', 'manual', 'PaperPosition', [1 3 4.5 4.5]);
    print(hout(i), '-depsc2', '-cmyk', outgraphics(i,:));
end

save('mapvar.mat','xe','ye', 'Ae', 'Bz');

residuef=residue;
save solvervariables ave_beta Mta Mta2 residuef Y0 Bzmax Jz00 axflux current zs polflux;


fluxstr=strcat(num2str(axflux,'%3.3f'),';', num2str(polflux,'%3.3f'));

H='R';
if Jz(J0,I0)*Bzmax < 0
    H='L';
end    
    

  intbegin=datestr(Datem(1));
  duration=datestr(Datem(length(Datem))-Datem(1));
  
  parid=fopen('00_results_for_table.txt', 'at');
  
  fprintf(parid, '    Begin            Duration     VHT (km/s) Bzmax (nT) (theta,phi) H  size (AU) Y0 (AU)   ax/polflux (1e21 Mx)  Iz (MA)\n');  
  fprintf(parid, ' %s (%s) & %i & %f &    (%i,%i) & %s & %s & %s &  %s &  %s\n',intbegin, duration , abs(round(norm(vht)/1000)), Bzmax, round(ZlatGSE), round(ZlongGSE), H, num2str(max(xe),'%1.3f'),num2str(Y0,'%1.3f'), fluxstr, num2str(round(current/1e6)));
  fclose(parid);

cd ..


zs
ys
xs

dht



r=zs(1); t=zs(2); n=zs(3);

% FOR DELTA/LAMBDA representation like Hu 3/2004
%'delta' entspricht DELTA und 'phi' entspricht LAMBDA in Hu 3/2004
%bzw.
%'delta' entspricht delta und 'phi' entspricht lambda in Hu 8/2005

delta=acos(r)*180/3.1415;
phi=acos(t/sqrt(n*n+t*t));
if n < 0
    phi=2*3.1415-phi;
end
phi=phi*180/pi;
delta_phi=[delta phi];


cd(dirstr)
 type 00_results_for_table.txt
cd ..
disp('--------------- END hu34n.m ------------------')


