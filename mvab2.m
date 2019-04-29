function [X, M13]=mvab2(B, flag)
% if flag=1, original magnetic field vectors
% else if flag=0, unit normalized magnetic field vectors


if flag==0
%%%%%%%%%%%%% normalize the B vectors into unit vectors
lBl=sqrt(B(:,1).^2+B(:,2).^2+B(:,3).^2);
Bhat(:,1)=B(:,1)./lBl;
Bhat(:,2)=B(:,2)./lBl;
Bhat(:,3)=B(:,3)./lBl;
%Bhat=B;
%%%%%%%%%%%%%
elseif flag==3
	lBlt=sqrt(B(:,1).^2+B(:,3).^2);
	Bhat(:,1)=B(:,1)./lBlt;
	Bhat(:,2)=B(:,2);
	Bhat(:,3)=B(:,3)./lBlt;
elseif flag==2
	lBl=sqrt((B(:,2)).^2);
	Bhat(:,1)=B(:,1)./lBl;
	Bhat(:,2)=B(:,2);
	Bhat(:,3)=B(:,3)./lBl;

else
	Bhat=B;
end

for i=1:3
	for j=1:3
		M(i,j)=mean((Bhat(:,i).*Bhat(:,j)))-mean(Bhat(:,i))*mean(Bhat(:,j));
%		M(i,j)=mean((Bm(:,i).*Bm(:,j)))-mean(Bm(:,i))*mean(Bm(:,j));
	end
end


[X,D]=eig(M);


format long;

[limda,In]=sort([D(1,1) D(2,2) D(3,3)]) ;
X=[X(:,In(3)) X(:,In(2)) X(:,In(1))];

limda=(fliplr(limda))';

aBxi=[mean(Bhat*X(:,1)) mean(Bhat*X(:,2)) mean(Bhat*X(:,3))];
%aBxi=[mean(Bm*X(:,1)) mean(Bm*X(:,2)) mean(Bm*X(:,3))]

%if abs(limda(3)) <= 1e-6
%disp('max. to intermediate');
%M13=limda(1)/limda(2);
%else
%M13=limda(1)/limda(3);
%end

r=size(B,1)-1;
D21=sqrt(limda(3)*(limda(2)+limda(1)-limda(3))/(r*(limda(1)-limda(2))^2));
D23=sqrt(limda(3)*(limda(2)+limda(3)-limda(3))/(r*(limda(2)-limda(3))^2));
D13=sqrt(limda(3)*(limda(1)+limda(3)-limda(3))/(r*(limda(1)-limda(3))^2));

M13=[limda [D21; D23; D13]];
