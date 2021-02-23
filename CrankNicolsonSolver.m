function Psi=CrankNicolsonSolver(psi,Vfunc,x1,x2,NumX,tau,NumT,m,hbar,fileID)
    deltaT=tau/NumT;
    Tlist=linspace(0,tau,NumT+1);
    Xlist=linspace(x1,x2,NumX+1);
    Xlist=Xlist(1:NumX);
    Psi=psi(Xlist)';
    writematrix(Psi',fileID);
    for t=Tlist(1:end-1)
        U=TimeEvo(Vfunc,x1,x2,m,hbar,NumX,deltaT,t);
        Psi=U*Psi;
        writematrix(Psi',fileID,'WriteMode','append');
    end
end

function U=TimeEvo(Vfunc,x1,x2,m,hbar,NumX,deltaT,t)
    Vf=@(x) Vfunc(x,t);
    H=discreteHmt(Vf,m,hbar,x1,x2,NumX);
    s=1i*deltaT/2/hbar;
    E=eye(NumX);
    N=E-s*H;
    D=E+s*H;
    U=D\N;
end

function H=discreteHmt(Vfunc,m,hbar,x1,x2,NumX)
    deltaX=(x2-x1)/NumX;
    X=linspace(x1,x2,NumX+1);
    X=X(1:NumX);
    V=Vfunc(X);
    A=hbar^2/m/deltaX^2;
    B=-A/2;
    NumX=length(X);
    H=diag(A+V);
    OffD=diag(B.*ones(1,NumX-1),1);
    OffD(1,NumX)=B;
    H=H+OffD+OffD';
end