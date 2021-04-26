function [A,B,b]=makeAB_inst(N,c,sigma2,gamma)

mu=sqrt(2*sigma2/pi);

A=zeros(N);
B=zeros(N);
for n=1:N
    for m=n+1:N
        if rand*N<c
            A(n,m)=sqrt(sigma2)*randn;
            if 2*rand<1-gamma
                A(m,n)=-A(n,m);
                B(n,m)=-abs(A(n,m));
                B(m,n)=-abs(A(n,m));
            else
                A(m,n)=A(n,m);
            end          
        end
    end
end

b=max(real(eig(A)));

A=A-b*eye(N);
B=B+(2*b+c*mu)*eye(N);

end