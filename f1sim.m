N=1000;
c=50;
gamma=-1;
sigma2=1/4/c;
mu=sqrt(2*sigma2/pi);
b=0.05;
%b=0.4+sqrt(c*sigma2)*(1+gamma);

[A,B]=makeAB(N,c,sigma2,gamma,b);
%%
T=2^10;
H=2^17;
t=(0:H-1)*T/H;

G=chol(B);
X1=OU_timeseries(A,G,t);
%%
w1=(0:0.01:1.5);
[Pe1,we1]=empirical_spectrum(X1(1:N,:),t,max(w));

%%
Nx=200;
Ny=800;
cx=20;
cy=5;
a=10;
b2=1;
d=1;
[A,B]=makeAB_bi(Nx,Ny,cx,cy,a,b2,d);

%e=eig(A);
%hist(imag(e),100);
%%
T=2^10;
H=2^17;
t=(0:H-1)*T/H;

G=chol(B);
X=OU_timeseries(A,G,t);

%%

w=(0:0.01:1.5);
[Pe,we]=empirical_spectrum(X(1:Nx,:),t,max(w));

%%
xs=(cx*a*b2-d)/(cx*cy*a^2+1);
ys=(cy*a*d+b2)/(cx*cy*a^2+1);

Px=[0,0;0,0];
Py=[0,0;0,0];
L=length(w);
qx=zeros(1,L);
qy=zeros(1,L);

for n=1:L
 
    XX=[0,xs-i*w(L+1-n);xs+i*w(L+1-n),-2*xs];
    YY=[0,ys-i*w(L+1-n);ys+i*w(L+1-n),-2*ys*(1+ys)];
    Z=[0,-ys;xs,xs*ys];
    
    for q=1:100
        
        Px=inv(XX-a^2*cx*Z*Py*Z');
        Py=inv(YY-a^2*cy*Z'*Px*Z);
        
    end
    
    qx(L+1-n)=real(Px(1,1));
    qy(L+1-n)=real(Py(1,1));
    
end
%%
subplot(2,2,1);
plot(X1(1,:));
xlim([0,10000]);
xlabel('$t$','interpreter','latex');
ylim([-20,20]);
ylabel('$\xi$','interpreter','latex','rotation',0);
set(gca,'xaxisLocation','top')
subplot(2,2,3);
plot(we1,mean(Pe1),'o',w1,(2*b+c*mu)/2/c/sigma2/b*real(sqrt(4*c*sigma2-w1.^2)),'k');
ylim([0,200]);
xlim([0,1.5]);
xlabel('$\omega$','interpreter','latex');
ylabel('$\phi(\omega)$','interpreter','latex','rotation',0,'position',[-0.2,110]);

subplot(2,2,2);
plot(X(1,:));
xlim([0,10000]);
xlabel('$t$','interpreter','latex');
ylim([-2,2]);
ylabel('$\xi$','interpreter','latex','rotation',0);
set(gca,'xaxisLocation','top')
subplot(2,2,4);
plot(we,mean(Pe),'o',w,qx,'k')
xlim([0,1.5]);
ylim([0,20]);
xlabel('$\omega$','interpreter','latex');
ylabel('$\phi(\omega)$','interpreter','latex','rotation',0,'position',[-0.2,11]);
