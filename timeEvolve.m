m=1;
hbar=1;
Vfunc=@(x,t) x.^2/2;
psi=@(x) (1/pi)^0.25*exp(-x.^2/2);
Psi=CrankNicolsonSolver(psi,Vfunc,-5,5,500,4*pi,500,m,hbar,'test.txt');