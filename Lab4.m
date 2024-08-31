clc; clear; close all;
h=6.626*1e-34;
hbar=h/2/pi;
eV=1.602*1e-19;
m0=9.1*1e-31;
meff=0.067*m0;
U0=0.280*eV;
w=5:5:100;
w=w*1e-9;
Ep=zeros(length(w),length(w));
IDX=1;
for idx=1:length(w)
xs=-1.5*w(idx); xd=1.5*w(idx); Nx=1000; x=linspace(xs,xd,Nx); dx=x(2)-x(1);
a=2*meff*dx^2/hbar^2;
U=zeros(1,Nx); U(x<-w(idx)/2)=U0; U(x>w(idx)/2)=U0;
phi=zeros(1,Nx);
phi(2)=1e-6;
Emin=0; Emax=U0; NE=85000; E=U0*logspace(-2,0,NE);
phid=zeros(1,NE);
for iE=1:NE
    for ix=2:Nx-1
        phi(ix+1)=2*phi(ix)-phi(ix-1)+a*(U(ix)-E(iE))*phi(ix);
    end
    phid(iE)=phi(Nx);
end
logphi=log(abs(phid));
contor=1;
for iE=2:NE-1
    if (logphi(iE-1)>logphi(iE))&&(logphi(iE)<logphi(iE+1))
        Ep(IDX,contor)=E(iE); contor=contor+1; %O matrice in care stochez valorile necesare
    end
end
IDX=IDX+1;
end
Ep=Ep/eV*1000;
h=figure(1);
hold on;
xlabel('Largimea gropii (w/nm)'); ylabel('Nivele de energie (Ep/meV)'); grid;
title('Dependenta nivelelor de energie (Ep(w))');
for iEp=1:length(w)+1
   plot(w,Ep(:,iEp),'color',[rand rand rand]);
end
set(h,'Units','Inches');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[8.2677 11.6929]);
print(h,'Laborator4_grafic','-dpdf','-r0','-bestfit');
%*Comentarii*
%Curbele nivelelor de energie 3,4,5,6..21 incep la valori mai mari ale largimii gropii
%Am trasat niste linii pentru a scoate in evidenta aceste puncte.
%Se pot observa anumite momente in care Ep creste putin.