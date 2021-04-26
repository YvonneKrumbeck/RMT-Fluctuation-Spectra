function [Se,we]=empirical_spectrum(X,t,wmax)

    H=length(t);
    T=max(t);

    %[Se,we]=pwelch(X'-mean(X'),H/4);
    [Se,we]=pwelch(X'-1,H/4);
    we=we*H/T;
    Se=Se'*T/H*pi;      
    
    q=find(we>wmax,1);
    we([1,q:end])=[];
    Se(:,[1,q:end])=[];

end