function X=mvab_perp(B, flag, ev)
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
else
	Bhat=B;
end

delta=eye(length(ev));

for i=1:3
	for j=1:3
		M(i,j)=mean((Bhat(:,i).*Bhat(:,j)))-mean(Bhat(:,i))*mean(Bhat(:,j));
		P(i,j)=delta(i,j)-ev(i)*ev(j);
	end
end

[X,D]=eig(P*M*P);


format long;

[limda,In]=sort([D(1,1) D(2,2) D(3,3)]) ;
X=[X(:,In(3)) X(:,In(2)) X(:,In(1))];

limda=(fliplr(limda))';

aBxi=[mean(Bhat*X(:,1)) mean(Bhat*X(:,2)) mean(Bhat*X(:,3))];
%aBxi=[mean(Bm*X(:,1)) mean(Bm*X(:,2)) mean(Bm*X(:,3))]




X=[X limda]; 
