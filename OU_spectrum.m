function S=OU_spectrum(A,B,w)

    A=full(A);
    B=full(B);

    N=length(A);
    L=length(w);

    [U,e]=eig(A); a=diag(e); UBU=U'*B*U; 

    S=zeros(N,L); 

    wb=waitbar(0,'spectrum');
    
    for s=1:L
        PS=U*diag(1./(a-i*w(s)))*UBU*diag(1./(a'+i*w(s)))*U'; 
        S(:,s)=real(diag(PS)); 
        
        waitbar(s/L);
    end
        
    close(wb);
end