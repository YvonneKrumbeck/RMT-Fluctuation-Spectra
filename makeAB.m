function [A,B]=makeAB(N,c,sigma2,gamma,b)

mu=sqrt(2*sigma2/pi);

A=zeros(N);
B=zeros(N);
for n=1:N
    for m=n+1:N
        if rand*N<c
            %r=mvnrnd([0,0],[1,gamma;gamma,1]);
            %A(n,m)=r(1);
            %A(m,n)=r(2);
            %if prod(r)<0
            %    B(n,m)=-abs(r(1));
            %    B(m,n)=-abs(r(2));
            %end
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


A=A-b*eye(N);
B=B+(2*b+c*mu)*eye(N);

end