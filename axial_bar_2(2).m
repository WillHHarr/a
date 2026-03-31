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

    k(1,1) = ((pi*T*(r2+r1))/h) + ks*pi*h*(r1/3 + h/12)*2;
    k(1,2) = -1*((pi*T*(r2+r1))/h) + (pi*ks*h)/3;
    k(2,1) = k(1,2);
    k(2,2) = ((pi*T*(r2+r1))/h) + ks*pi*h*(r1/3 + 7*h/12)*2;

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

%--- I did use AI for the post-processing... I will admit I am useless
%modifying matlab plotting code.
%POSTPROCESSING
%
% Rebuild the global displacement vector
u = zeros(numnp,1);
for i = 1:numnp
    if idof(i) > 0
        u(i) = Uf(idof(i));
    else
        u(i) = Up(-idof(i));
    end
end

% Print nodal displacements
fprintf('\nNodal Displacements:\n')
fprintf('Node    r (m)        w (m)\n')
for i = 1:numnp
    fprintf('  %3d   %8.4f   %15.8f\n', i, x(i), u(i))
end

% ----------------------------------------------------------------
% PLOT 1: Deflection for N = 1, 2, 3, 4 on the same graph
% ----------------------------------------------------------------
figure(1)
clf
hold on

colors  = {'r','b','g','m'};
markers = {'o','s','^','d'};

for Nplot = 1:4
    % Re-run the full FE solve for this N
    [x_plt, u_plt] = solve_membrane(Nplot, a, b, T, Q, ks, wa, wb);
    plot(x_plt, u_plt, ...
        [colors{Nplot} markers{Nplot} '-'], ...
        'LineWidth', 2, ...
        'MarkerSize', 7, ...
        'DisplayName', sprintf('N = %d', Nplot))
end

hold off
title('Axi-Symmetric Membrane Deflection - YourName')
xlabel('r (m)')
ylabel('w(r) (m)')
legend('show', 'Location', 'best')
grid on

% ----------------------------------------------------------------
% PLOT 2: Effect of foundation stiffness ks, N = 20
% ----------------------------------------------------------------
figure(2)
clf
hold on

ks_vals  = [0, 50, 500, 5000];
colors2  = {'b','r','g','k'};

for ii = 1:length(ks_vals)
    [x_plt, u_plt] = solve_membrane(20, a, b, T, Q, ks_vals(ii), wa, wb);
    plot(x_plt, u_plt, ...
        [colors2{ii} '-'], ...
        'LineWidth', 2, ...
        'DisplayName', sprintf('k_s = %g Pa/m', ks_vals(ii)))
end

hold off
title('Effect of Foundation Stiffness k_s - YourName')
xlabel('r (m)')
ylabel('w(r) (m)')
legend('show', 'Location', 'best')
grid on

% ----------------------------------------------------------------
% Print reactions at supports
% ----------------------------------------------------------------
fprintf('\nSupport Reactions:\n')
Rp = Rp - (Kpf*Uf + Kpp*Up);
for i = 1:numnp
    if idof(i) < 0
        fprintf('Node %d   Reaction = %15.5f N\n', i, Rp(-idof(i)))
    end
end


% ================================================================
% LOCAL FUNCTION: re-solves the membrane problem for a given N and ks
% so that Plot 1 and Plot 2 can loop over different parameters
% without re-running the whole script.
% ================================================================
function [x_out, w_out] = solve_membrane(Nval, a, b, T, Q, ks, wa, wb)

    numnp = Nval + 1;
    numimp = 2;
    numeq  = numnp - numimp;

    % idof vector
    idof = zeros(numnp,1);
    idof(1)     = -1;
    idof(numnp) = -2;
    for i = 2:numnp-1
        idof(i) = i-1;
    end

    % Global matrices
    Kff = sparse(zeros(numeq,numeq));
    Kpf = sparse(zeros(numimp,numeq));
    Rf  = zeros(numeq,1);

    % Imposed displacements
    Up    = zeros(numimp,1);
    Up(1) = wa;
    Up(2) = wb;

    % Node coordinates
    x = zeros(numnp,1);
    x(1) = a;
    for i = 1:Nval
        x(i+1) = x(i) + (b-a)/Nval;
    end

    % Connectivity
    lm = zeros(2,Nval);
    for i = 1:Nval
        lm(1,i) = i;
        lm(2,i) = i+1;
    end

    % Assembly
    for i = 1:Nval
        k = zeros(2,2);
        r = zeros(2,1);

        n1 = lm(1,i);
        n2 = lm(2,i);
        r1 = x(n1);
        r2 = x(n2);
        h  = r2 - r1;
        rm = (r1 + r2)/2;

        k(1,1) = 2*pi * (T*rm/h + ks*h*(r1/3 + h/12));
        k(1,2) = 2*pi * (-T*rm/h + ks*h*(r1/6 + h/12));
        k(2,1) = k(1,2);
        k(2,2) = 2*pi * (T*rm/h + ks*h*(r1/3 + 7*h/12));

        r(1) = 2*pi * Q*h*(2*r1 + r2)/6;
        r(2) = 2*pi * Q*h*(r1 + 2*r2)/6;

        iaux = [idof(n1); idof(n2)];
        for jj = 1:2
            nj = iaux(jj);
            for kk = 1:2
                nk = iaux(kk);
                if nj > 0 && nk > 0
                    Kff(nj,nk) = Kff(nj,nk) + k(jj,kk);
                elseif nj < 0 && nk > 0
                    Kpf(-nj,nk) = Kpf(-nj,nk) + k(jj,kk);
                end
            end
            if nj > 0
                Rf(nj) = Rf(nj) + r(jj);
            end
        end
    end

    % Apply BCs and solve
    Rf  = Rf - Kpf'*Up;
    Uf  = Kff\Rf;

    % Rebuild full solution
    w_out = zeros(numnp,1);
    for i = 1:numnp
        if idof(i) > 0
            w_out(i) = Uf(idof(i));
        else
            w_out(i) = Up(-idof(i));
        end
    end

    x_out = x;
end
