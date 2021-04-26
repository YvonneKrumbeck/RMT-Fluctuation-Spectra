% f3sim

N=500;
c=20;
gamma=-1;
sigma2=1/4/c;
mu=sqrt(2*sigma2/pi);
b=0.2;

[A,B]=makeAB(N,c,sigma2,gamma,b);

w=(0:0.01:1.5);

%%
%T=2^10;
%H=2^17;
%t=(0:H-1)*T/H;

%G=chol(B);
%X=OU_timeseries(A,G,t);
%%

%[Pe,we]=empirical_spectrum(X(1:N,:),t,max(w));

%%
S=sqrt((b^2+4*c*sigma2)^2+2*(b^2-4*c*sigma2)*w.^2+w.^4);
Q=sqrt(b^2+4*c*sigma2+S-w.^2);

Pm=(2*b+c*mu)*(sqrt(2)*Q-2*b)/4/b/c/sigma2;
rmx=(-2*b+sqrt(2)*sqrt(b^2+4*c*sigma2-w.^2+sqrt(16*b^2*c*sigma2+(b^2-4*c*sigma2+w.^2).^2)))/4/c/sigma2;
rmy=(-sqrt(2)*b^2*Q-4*b*w.^2+sqrt(2)*Q.*(-4*c*sigma2+S+w.^2))/8/b/c/sigma2./w;
rmy(1)=rmy(2); % fudge to fix w=0

Ps=zeros(N,length(w));
for n=1:N
    A2(n)=sum(A(n,:).^2)-A(n,n)^2;
    AB(n)=A(n,:)*B(n,:)'-A(n,n)*B(n,n);
    AA(n)=A(n,:)*A(:,n)-A(n,n)^2;
    Ps(n,:)=(Pm*A2(n)+2*rmx*AB(n)+B(n,n))./abs(A(n,n)+i*w+(rmx-i*rmy)*AA(n)).^2;
    %Ps(n,:)=(Pm*c*sigma2+2*b+c*mu)./abs(-b+i*w+(rmx-i*rmy)*(-c*sigma2)).^2;
end

%%

wd=(0:0.01:1.5);
Pd=OU_spectrum(A,B,wd);

%%

[~,n1]=max(Pd(:,1));
[~,n2]=min(Pd(:,1));

figure();
hold on;
for n=1:250
    p=plot(w,Ps(n,:),'LineWidth',1.5);
    p.Color=[0,0.45,0.74,0.1];
end
plot(wd,Pd(n1,:),'k','LineWidth',2);
plot(wd,Pd(n2,:),'k','LineWidth',2);
plot(wd,mean(Pd),'--k','LineWidth',2);
hold off;
xlim([0,1.5]);
xlabel('$\omega$','interpreter','latex');
ylabel('$\phi(\omega)$','interpreter','latex','rotation',0);

