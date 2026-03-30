% Program Axial_Bar_2
%
%Basic Input
%
format long %number of signigicant digits in output screen
N = input('Enter the number of elements in the beam ')
L=2
A=0.1
E=10000
po=500
uimposed=0.01

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
Up(1) = 0;  %fixed at the wall
Up(2)=uimposed; %applied end displacement


%create coordinate vector x(N+1) and connectivity table lm(2,N)

x=zeros(N+1,1);
x(1)=0;
for i = 1:N
	x(i+1) = x(i)+L/N;
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
	xm=(x(i)+x(i+1))/2; %location of element center
	pm=po*xm*(L-xm)/L^2;
	%compute local sm and local lv
	k(1,1)=E*A/le;
	k(1,2)=-k(1,1);
	k(2,1)=k(1,2);
	k(2,2)=k(1,1);
	r(1)=pm*le/2;
	r(2)=pm*le/2;
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
%compute exact solution
uex1=zeros(N+1,1);  %used for pointwise comparison
xx=zeros(501,1);	%used for plotting exact solution
xx=0:L/500:L;
uex2=zeros(501,1);  %used for plotting exact solution
eta=E*A*uimposed/(po*L^2);
uex1=po*L^2/(E*A)*(((x/L).^4)/12-((x/L).^3)/6+(1/12+eta)*(x/L));
uex2=po*L^2/(E*A)*(((xx/L).^4)/12-((xx/L).^3)/6+(1/12+eta)*(xx/L));
%print and plot solutions
sprintf('%s','nodal displacements')
iii=1:1:numnp;
%build display matrix the ' means transpose
aux=[iii;x';u';uex1'];
sprintf('node %d  x=%8.3f  u=%15.5f  uex=%15.5f\n', aux)
subplot(1,2,1),plot(x/L,E*A*u/(po*L^2),'ro-',xx/L,E*A*uex2/(po*L^2),'b-','linewidth',2)
title('Displacement')
xlabel('$\frac{x}{L}$','Interpreter','latex')
ylabel('$\frac{EAu(x)}{p_0L^2}$','Interpreter','latex')
legend('numerical','analytical','location', 'south')

%compute largest differences between exact and FE solutions
%at the nodes(*100/uimposed)
error=100*max(abs(u-uex1))/uimposed;
sprintf('maximum error on nodal displacement = %17.10f percent\n',error)

%Compute and plot stress distribution
stress = zeros(numel,1);
xc = zeros(numel,1);
for i=1:numel
    n1=lm(1,i);
    n2=lm(2,i);
    le=abs(x(n2)-x(n1));
    strain=(u(n2)-u(n1))/le;
    stress(i)=E*strain;
    xc(i)=(x(n2)-x(n1))/2;
end
sprintf('%s','axial stress distribution')
iii=1:1:numel;
aux=[iii;xc';stress'];
sprintf('element %d  xc=%8.4f  stress=%15.8f Pascal \n', aux)
xx=zeros(2*numel,1); %used to plot the FE stress distribution
stress2=zeros(2*numel,1); %used t plot the FE stress distirbution
for i = 1:numel
    xx(2*i-1)=x(lm(1,i));
    xx(2*i)=x(lm(2,i));
    stress2(2*i-1)=stress(i);
    stress2(2*i)=stress(i);
end
xxx=zeros(501,1); %used to plot the eact stress distribution
xxx=0:L/500:L;
sigma_ex=po*L/A*(((xxx/L).^3)/3-((xxx/L).^2)/2+(1/12+eta));
subplot(1,2,2),plot(xxx/L,A*sigma_ex/(po*L),'b-',xx/L,A*stress2/(po*L),'r.:','linewidth',2)
title('axial stress')
xlabel('$\frac{x}{L}$','Interpreter','latex')
ylabel('$\frac{A\sigma}{p_0L}$','Interpreter','latex')
legend('numerical','analytical','location', 'southwest')

%Obtain and print reactions at supports
sprintf('%s','support reactions')
Rp=Rp-(Kpf*Uf+Kpp*Up)
for i=1:numnp
    if idof(i)<0
        sprintf('node # %d  Reaction= %15.5f Newtons', i,Rp(-idof(i)))
    end
end




