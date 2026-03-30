% Program Axial_Bar_2
%
%Basic Input
%
format long %number of signigicant digits in output screen
N = input('Enter the number of elements in the beam ')
a = 0.1; %m
b = 0.5; %m
T = 10; %N/m
Q = 3; %Pa
wa = 0.01; %m
wb = 0.02; %m
ks = 85; %Pa/m

%
%INITALIZATIONS
%

numnp=N+1;  %number of nodes
numel=N;	%number of elements
numimp=2;	%number of imposed dof
numeq=numnp-numimp; %number of equations
idof = zeros(numnp,1); %vector with equation numbers
idof(1) = -1; %negative for first and last nodes
idof(numnp) = -2;
for i = 2:numnp-1
	idof(i)=i-1;
end
%intialize global stiffness matricies kff, kpf, kpp
%and global load vectors Rp, Rf

Kff=sparse(zeros(numeq,numeq));
Kpf=sparse(zeros(numimp,numeq));
Kpp=sparse(zeros(numimp,numimp));
Rf=zeros(numeq,1);
Rp=zeros(numimp,1);

%Initialize Global DOF vectors Uf and Up
Uf = zeros(numeq,1);
Up = zeros(numimp,1);
Up(1) = wa;  %fixed at the wall
Up(2)= wb; %applied end displacement


%create coordinate vector x(N+1) and connectivity table lm(2,N)

x=zeros(N+1,1);
x(1)=a;
for i = 1:N
	x(i+1) = x(i)+ (b-a)/N;
end

lm=zeros(2,N);
for i=1:N
	lm(1,i)=i;
	lm(2,i)=i+1;
end
%
%FINITE ELEMENT SOLUTION
%
%loop over the N elements to get local and global K & R
for i=1:N
	k=zeros(2,2);
	r=zeros(2,1);
	%get absolute length
	le=abs(x(i+1)-x(i));
	%get two node numbers for element i
	n1=lm(1,i);
	n2=lm(2,i);
	%get load on the element (=value of p(x) at the center of the element)

    r1 = x(i);
    r2 = x(i+1);
    h = r2 - r1;
    rm = (r1 + r2)/2;

    k(1,1) = ((pi*T*(r2+r1))/h) + 2*(pi*ks*h)/3;
    k(1,2) = -1*((pi*T*(r2+r1))/h) + (pi*ks*h)/3;
    k(2,1) = k(1,2);
    k(2,2) = ((pi*T*(r2+r1))/h) + 2*(pi*ks*h)/3;

    r(1) = 2*pi * Q*h*(2*r1 + r2)/6;
    r(2) = 2*pi * Q*h*(r1 + 2*r2)/6;
	%assemble into Kff, Kpf, Kpp and Rf and Rp
	iaux=zeros(2,1);
	iaux(1) = idof(n1);
	iaux(2) = idof(n2);
	for jj=1:2
		nj=iaux(jj); 
			for kk=1:2
				nk=iaux(kk);
				if nj>0 && nk>0
					Kff(nj,nk)=Kff(nj,nk)+k(jj,kk); 
				elseif nj<0 && nk>0
					Kpf(-nj,nk)=Kpf(-nj,nk)+k(jj,kk); 
				elseif nj<0 && nk<0
					Kpp(-nj,-nk)=Kpp(-nj,-nk)+k(jj,kk);
				end 
			end
			if nj>0 
				Rf(nj)=Rf(nj)+r(jj);
			else
				Rp(-nj)=Rp(-nj)+r(jj);
			end 
		end
end %end loop over elements
%account for displacement boundary conditions
Rf=Rf-(Kpf)'*Up;
%solve linear system
Uf=Kff\Rf;

%POSTPROCESSING
%
%
%Rebuild, display and plot the global displacement vector
u=zeros(numnp,1);
for i=1:numnp
	if idof(i)>0
		u(i)=Uf(idof(i));
	else
		u(i)=Up(-idof(i));
	end
end

%I did use clause to modify the plotting issues as I never understand the
%plotting functions if I didn't write them myself
%print and plot solutions
sprintf('%s','nodal displacements')
iii=1:1:numnp;
aux=[iii;x';u'];
sprintf('node %d  x=%8.3f  u=%15.5f\n', aux)

xx = a:(b-a)/500:b;

plot(x, u, 'ro-', 'linewidth', 2)
title('Axi-symmetric Membrane Deflection - YourName')
xlabel('r (m)')
ylabel('w(r) (m)')
legend('numerical','location','best')

%compute largest differences between exact and FE solutions
%not available for this problem, remove error calculation

%Obtain and print reactions at supports
sprintf('%s','support reactions')
Rp=Rp-(Kpf*Uf+Kpp*Up);
for i=1:numnp
    if idof(i)<0
        sprintf('node # %d  Reaction= %15.5f Newtons', i, Rp(-idof(i)))
    end
end




