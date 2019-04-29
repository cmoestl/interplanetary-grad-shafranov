function [X, L12, L0]=mvab2D(B, flag)
% if flag=1, original magnetic field vectors Bengt-Hu
% else if flag=0, unit normalized magnetic field vectors


if flag==0
%%%%%%%%%%%%% normalize the B vectors into unit vectors
lBl=sqrt(B(:,1).^2+B(:,2).^2);
Bhat(:,1)=B(:,1)./lBl;
Bhat(:,2)=B(:,2)./lBl;
%Bhat(:,3)=B(:,3)./lBl;
%Bhat=B;
%%%%%%%%%%%%%


elseif flag==1
	Bhat=B;
else
end

for i=1:2
	for j=1:2
		M(i,j)=mean((Bhat(:,i).*Bhat(:,j)))-mean(Bhat(:,i))*mean(Bhat(:,j));
%		M(i,j)=mean((Bm(:,i).*Bm(:,j)))-mean(Bm(:,i))*mean(Bm(:,j));
	end
end

[X,D]=eig(M);
%L0=abs(det(M));

format long;

[limda,In]=sort([D(1,1) D(2,2)]) ;
X=[X(:,In(2)) X(:,In(1))];

limda=(fliplr(limda))';

aBxi=[mean(Bhat*X(:,1)) mean(Bhat*X(:,2))];
%aBxi=[mean(Bm*X(:,1)) mean(Bm*X(:,2)) mean(Bm*X(:,3))]


L12=sqrt(limda(1)*limda(2)/((size(B,1)-1)*(limda(1)-limda(2)).^2)); %uncertainty in Radians
L0=limda;
