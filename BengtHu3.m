function [Xt, M]=BengtHu3(Ba, x, j)
%Bengt-Hu approach in axis determination 2D MVAB

x=x-x(length(x))/2;
Lp=x(length(x));
x=x/Lp;
switch j
	case 1
	fac=(x)/((x(length(x))));
Bt(:,1)=(Ba(:,1)./Ba(:,2)).*(fac');
%Bt(:,2)=zeros(size(Ba(:,2)));
Bt(:,2)=(Ba(:,3)./Ba(:,2)).*(fac');
[Xt, L12, M]=mvab2D(Bt, 1);

	case 2
for k=1:800
kapa=pi/2*(k/800);
Bt(:,1)=(Ba(:,1)./Ba(:,2)).*sin(x'*sin(kapa))/(sin(kapa));
%Bt(:,2)=zeros(size(Ba(:,2)));
Bt(:,2)=(Ba(:,3)./Ba(:,2)).*sin(x'*sin(kapa))/(sin(kapa));

[Xt, L12, M]=mvab2D(Bt, 1);
LL0(k)=M(2);

end
K0=find(LL0==min(LL0));
kapa=pi/2*(K0(1)/800);
Bt(:,1)=(Ba(:,1)./Ba(:,2)).*sin(x'*sin(kapa))/(sin(kapa));
%Bt(:,2)=zeros(size(Ba(:,2)));
Bt(:,2)=(Ba(:,3)./Ba(:,2)).*sin(x'*sin(kapa))/(sin(kapa));

[Xt, L12, M]=mvab2D(Bt, 1);

	
	case 3
for k=1:400
kapa=pi/2*(k/400);
Bt(:,1)=(Ba(:,1)./Ba(:,2)).*sinh((x')*sin(kapa))/(sin(kapa));
%Bt(:,2)=zeros(size(Ba(:,2)));
Bt(:,2)=(Ba(:,3)./Ba(:,2)).*sinh(x'*sin(kapa))/(sin(kapa));

[Xt, L12, M]=mvab2D(Bt, 1);
LL0(k)=M(2);

end
K0=find(LL0==min(LL0));
kapa=pi/2*(K0(1)/400);
Bt(:,1)=(Ba(:,1)./Ba(:,2)).*sinh((x')*sin(kapa))/(sin(kapa));
%Bt(:,2)=zeros(size(Ba(:,2)));
Bt(:,2)=(Ba(:,3)./Ba(:,2)).*sinh(x'*sin(kapa))/(sin(kapa));

[Xt, L12, M]=mvab2D(Bt, 1);



	otherwise
	disp('Undefined normalization function.');
end

M=[M; L12*180/pi];

