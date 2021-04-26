% % read the data file
% M6 = xlsread('41467_2017_2571_MOESM6_ESM.xlsx','B2:JE28599');
% 
% %%
% 
% % averge the 3 samples per day
% K=kron(eye(size(M6,2)/3),ones(3,1))/3;
% X=(M6*K)';
% 
% % retain only time series with positive numbers each day
% s=min(X);
% a=[find(isnan(s)),find(s==0)];
% X(:,a)=[];
% 
% %%
% 
% % arbitrary scale factor to make the spectrum O(1) (will remove this later)
% scale=100;
% m=mean(X*scale);
% 
% % subtract mean and compute spectrum of remainder
% Y=X*scale-(ones(88,1)*m);

load planktonY.mat

[P,W]=pcov(Y,8);

W(1)=[];
P(1,:)=[];

p=mean(P')';

s=@(v,b,cs2,g,x) -((v.*(-2.*b.^3+sqrt(2).*b.^2.*sqrt(b.^2-4.*cs2.*g-x.^2+sqrt(b.^4+(4.*cs2.*g+x.^2).^2+b.^2.*(-8.*cs2.*g+2.*x.^2)))-2.*b.*(x.^2+sqrt(b.^4+(4.*cs2.*g+x.^2).^2+b.^2.*(-8.*cs2.*g+2.*x.^2)))+sqrt(2).*sqrt(b.^2-4.*cs2.*g-x.^2+sqrt(b.^4+(4.*cs2.*g+x.^2).^2+b.^2.*(-8.*cs2.*g+2.*x.^2))).*(4.*cs2.*g+x.^2+sqrt(b.^4+(4.*cs2.*g+x.^2).^2+b.^2.*(-8.*cs2.*g+2.*x.^2)))))./(cs2.*(-2.*b.^3+sqrt(2).*b.^2.*sqrt(b.^2-4.*cs2.*g-x.^2+sqrt(b.^4+(4.*cs2.*g+x.^2).^2+b.^2.*(-8.*cs2.*g+2.*x.^2)))+sqrt(2).*sqrt(b.^2-4.*cs2.*g-x.^2+sqrt(b.^4+(4.*cs2.*g+x.^2).^2+b.^2.*(-8.*cs2.*g+2.*x.^2))).*(4.*cs2.*g+x.^2+sqrt(b.^4+(4.*cs2.*g+x.^2).^2+b.^2.*(-8.*cs2.*g+2.*x.^2)))-2.*b.*(-4.*cs2.*g.^2+x.^2+sqrt(b.^4+(4.*cs2.*g+x.^2).^2+b.^2.*(-8.*cs2.*g+2.*x.^2))))))

%%

options = fitoptions('Method', 'NonLinearLeastSquares');
options.Lower=[0,0,0,-1];
options.Upper=[Inf,Inf,Inf,1];
options.StartPoint=[0.3,1,0.3,0.1];
pfit=fit(W,p,s,options);
cv=coeffvalues(fit(W,p,s,options))

cv(2)-sqrt(cv(3))*(1+cv(4))


options.Lower=[0,0,0,-1];
options.Upper=[Inf,Inf,Inf,0];
options.StartPoint=[0.3,1,0.3,-0.1];
pfit=fit(W,p,s,options);
cvn=coeffvalues(fit(W,p,s,options))
%%

scale=100;
cv(1)=cv(1)/scale^2;
cvn(1)=cvn(1)/scale^2;


w=10.^(-2:0.1:1);
subplot(2,2,[1,3]);
loglog(W,p/scale^2,'o',w,s(cv(1),cv(2),cv(3),cv(4),w),'-k',w,s(cvn(1),cvn(2),cvn(3),cvn(4),w),'--k')
xlabel('\omega')
ylabel('\phi(\omega)')
legend('Mean power spectrum from data','Model best fit','Model best fit with \gamma<0')
subplot(2,2,2);
plot(W,p/scale^2,'o',w,s(cv(1),cv(2),cv(3),cv(4),w),'-k',w,s(cvn(1),cvn(2),cvn(3),cvn(4),w),'--k')
xlabel('\omega')
ylabel('\phi(\omega)')
axis([0,2,0,0.0003])
%%
N=500;
c=100;
sigma2=cv(3)/c;
b=cv(2);
gamma=cv(4);
[A,B]=makeAB(N,c,sigma2,gamma,b);
e=eig(A); 
subplot(2,2,4)
plot(e,'.k');
ellipse(sqrt(cv(3))*(1+cv(4)),sqrt(cv(3))*(1-cv(4)),0,-b,0,'k');
xlabel('Re\lambda')
ylabel('Im\lambda')
%%

