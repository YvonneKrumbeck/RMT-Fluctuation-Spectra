function [A,B]=makeAB_bi(Nx,Ny,cx,cy,a,b,d)

if Nx*cx~=Ny*cy
   error('degrees/dimensions mismatch');
end

xs=(cx*a*b-d)/(cx*cy*a^2+1);
ys=(cy*a*d+b)/(cx*cy*a^2+1);

success=0;
while success==0
    tries=0;
    R=zeros(Nx,Ny);
dx=cx*ones(1,Nx);
dy=cy*ones(1,Ny);
maxtries=10*Nx*cx;
while (sum(dx)>0)&&(tries<maxtries)
    tries=tries+1;
    x=randsample(Nx,1,true,dx);
    y=randsample(Ny,1,true,dy);
    if R(x,y)==0
        R(x,y)=a;
        dx(x)=dx(x)-1;
        dy(y)=dy(y)-1;
    end
end
if sum(dx)==0
    success=1;
end
    
end

A=[-xs*eye(Nx),xs*R;-ys*R',-ys*eye(Ny)];
B=[2*xs*(xs+d)*eye(Nx),-xs*ys*R;-xs*ys*R',2*ys*b*eye(Ny)];

end