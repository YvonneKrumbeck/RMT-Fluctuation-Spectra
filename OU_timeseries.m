function Xi=OU_timeseries(A,G,t)

    dt=t(2)-t(1);
    H=length(t);
    Adt=A*dt;
    Gsdt=G*sqrt(dt);
    
    N=length(A);
    M=size(G,2);
        
    xi=zeros(N,1);
    for h=1:H/4
        xi=xi+Adt*xi+Gsdt*randn(M,1);
    end
    Xi=zeros(N,H);

    wb=waitbar(0,'timeseries');
    
    for h=1:H
        Xi(:,h)=xi;
        xi=xi+Adt*xi+Gsdt*randn(M,1);
        
        if mod(h,64)==0
            waitbar(h/H);
        end
    end
    
    close(wb);
    
end
