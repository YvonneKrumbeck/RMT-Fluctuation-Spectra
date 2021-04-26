N=1000;
c=50;
sigma2=0.5;
mu=sqrt(2*sigma2/pi);
b0=2+sqrt(c*sigma2);
b1=2+2*sqrt(c*sigma2);

[A0,B0]=makeAB(N,c,sigma2,0,b0);
[A1,B1]=makeAB(N,c,sigma2,1,b1);

[A0i,B0i,b0i]=makeAB_inst(N,c,sigma2,0);
[A1i,B1i,b1i]=makeAB_inst(N,c,sigma2,1);

%%
T=2^8;
H=2^15;
t=(0:H-1)*T/H;

Ti=2^8;
Hi=2^16;
ti=(0:Hi-1)*Ti/Hi;

G0=chol(B0);
G1=chol(B1);
G0i=chol(B0i);
G1i=chol(B1i);

X0=OU_timeseries(A0,G0,t);
X1=OU_timeseries(A1,G1,t);
X0i=OU_timeseries(A0i,G0i,ti);
X1i=OU_timeseries(A1i,G1i,ti);
%%

w=(0:0.01:10);
[Pe0,we]=empirical_spectrum(X0,t,max(w));
[Pe1,we]=empirical_spectrum(X1,t,max(w));
wi=(0.01:0.01:100);
[Pe0i,wei]=empirical_spectrum(X0i,ti,max(wi));
[Pe1i,wei]=empirical_spectrum(X1i,ti,max(wi));
%%


S=sqrt((b1^2-4*c*sigma2)^2+2*(b1^2+4*c*sigma2)*w.^2+w.^4);
Q=sqrt(b1^2-4*c*sigma2+S-w.^2);


subplot(2,2,[1,3]);
plot(we,mean(Pe0),'o',we,mean(Pe1),'o',w,(2*b1+c*mu)*(-sqrt(2)*b1^2*Q-4*b1*w.^2+sqrt(2)*Q.*(4*c*sigma2+S+w.^2))/(8*b1*c*sigma2).*w.^(-2),w,(2*b0+c*mu)./(b0^2-c*sigma2+w.^2))
ylim([0,3])
xlabel('$\omega$','interpreter','latex');
ylabel('$\phi(\omega)$','interpreter','latex','rotation',0);

subplot(2,2,2);
plot(wei,mean(Pe0i),'o',wei,mean(Pe1i),'o',wi,(4*sqrt(c*sigma2)+c*mu)*sqrt(2).*(-wi.^2+sqrt(16*c*sigma2*wi.^2+wi.^4)).^(-1/2)/sqrt(c*sigma2)-(4*sqrt(c*sigma2)+c*mu)/2/(c*sigma2),wi,(2*b0i+c*mu)./((b0i^2-c*sigma2)+wi.^2))
ylim([0,5]);
xlim([0,10]);
xlabel('$\omega$','interpreter','latex');
ylabel('$\phi(\omega)$','interpreter','latex','rotation',0);

subplot(2,2,4);
loglog(wei,mean(Pe0i),'o',wei,mean(Pe1i),'o',wi,(4*sqrt(c*sigma2)+c*mu)*sqrt(2).*(-wi.^2+sqrt(16*c*sigma2*wi.^2+wi.^4)).^(-1/2)/sqrt(c*sigma2)-(4*sqrt(c*sigma2)+c*mu)/2/(c*sigma2),wi,(2*b0i+c*mu)./(wi.^2))
xlabel('$\omega$','interpreter','latex');
ylabel('$\phi(\omega)$','interpreter','latex','rotation',0);
