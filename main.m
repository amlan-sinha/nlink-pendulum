close all
clear all
clc

set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','Latex')

% N-Link Pendulum
% Simulate the dynamics of an n-link pendulum using:
% 1. Minimum Coordinates
% 2. Lagrange's Equations
% 3. Differential Algebraic Equations (DAEs)

% Amlan Sinha
% Cornell University

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           System parameters                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of links
n = 3;
% Initial conditions
th0 = ones(1,n)*pi/2;
dth0 = zeros(1,n);
% System parameters
p = struct('n',n,'g',1,'m',ones(1,n),'a',0.5*ones(1,n),'L',ones(1,n));
% Time span
time = [0 10];
% Flags
rederive = 1;
simulate = 1;
% Add the path to the folders
addpath('./DAEs');
addpath('./Lagrange');
addpath('./MinimumCoordinates');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                Solution                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T1,solAMB] = NLinkPend_AMB(n,rederive,simulate,p,[th0 dth0],time);
[T2,solLAG] = NLinkPend_LAG(n,rederive,simulate,p,[th0 dth0],time);
[T3,solDAE] = NLinkPend_DAE(n,rederive,simulate,p,[th0 dth0],time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                Graphics                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                                  Theta                                  %

figure
hold on
s = cell(1,n);
for i = 1:0.5*length(solAMB(1,:))
    plot(T1,solAMB(:,i),'LineWidth',3)
    s{i} = sprintf('Link %2.0f',i);
end
title('Solutions: AMB');
xlabel('$$time$$ (s)');
ylabel('$$\theta$$ (rad)');
legend(s,'Location','northeast');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
s = cell(1,n);
for i = 1:0.5*length(solLAG(1,:))
    plot(T2,solLAG(:,i),'LineWidth',3)
    s{i} = sprintf('Link %2.0f',i);
end
title('Solutions: Lagrange');
xlabel('$$time$$ (s)');
ylabel('$$\theta$$ (rad)');
legend(s,'Location','northeast');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
s = cell(1,n);
for i = 1:0.5*length(solDAE(1,:))
    plot(T3,solDAE(:,i),'LineWidth',3)
    s{i} = sprintf('Link %2.0f',i);
end
title('Solutions: DAE');
xlabel('$$time$$ (s)');
ylabel('$$\theta$$ (rad)');
legend(s,'Location','northeast');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                                  Error                                  %
figure
hold on
s = cell(1,n);
for i = 1:0.5*length(solAMB(1,:))
    plot(T1,solAMB(:,i)-solLAG(:,i),'LineWidth',3)
    s{i} = sprintf('Link %2.0f',i);
end
title('Comparing AMB with Lagrange');
xlabel('$$time$$ (s)');
ylabel('$$\theta$$ (rad)');
legend(s,'Location','northeast');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
s = cell(1,n);
for i = 1:0.5*length(solLAG(1,:))
    plot(T2,solLAG(:,i)-solDAE(:,i),'LineWidth',3)
    s{i} = sprintf('Link %2.0f',i);
end
title('Comparing Lagrange with DAE');
xlabel('$$time$$ (s)');
ylabel('$$\theta$$ (rad)');
legend(s,'Location','northeast');
hold off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
hold on
s = cell(1,n);
for i = 1:0.5*length(solDAE(1,:))
    plot(T3,solDAE(:,i)-solAMB(:,i),'LineWidth',3)
    s{i} = sprintf('Link %2.0f',i);
end
title('Comparing DAE with AMB');
xlabel('$$time$$ (s)');
ylabel('$$\theta$$ (rad)');
legend(s,'Location','northeast');
hold off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                          Conservation of energy                       %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[a,b] = size(solAMB);
b = 0.5*b;

AMB_y = zeros(a,b);  % y direction position for AMB
LAG_y = zeros(a,b);  % y direction position for LAG
DAE_y = zeros(a,b);  % y direction position for DAE

AMB_vx = zeros(a,b); % x direction velocity for AMB
AMB_vy = zeros(a,b); % y direction velocity for AMB

LAG_vx = zeros(a,b); % x direction velocity for LAG
LAG_vy = zeros(a,b); % y direction velocity for LAG

DAE_vx = zeros(a,b); % x direction velocity for DAE
DAE_vy = zeros(a,b); % y direction velocity for DAE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inertial frame velocity
for i = 1:b
	% Add the velocities for each link
    for j = 1:i-1
        AMB_vx(:,i) = AMB_vx(:,i) + p.L(j)*cos(solAMB(:,j)).*solAMB(:,n+j);
        AMB_vy(:,i) = AMB_vy(:,i) + p.L(j)*cos(solAMB(:,j)).*solAMB(:,n+j);
        
        LAG_vx(:,i) = LAG_vx(:,i) + p.L(j)*cos(solLAG(:,j)).*solLAG(:,n+j);
        LAG_vy(:,i) = LAG_vy(:,i) + p.L(j)*cos(solLAG(:,j)).*solLAG(:,n+j);
        
        DAE_vx(:,i) = DAE_vx(:,i) + p.L(j)*cos(solDAE(:,j)).*solDAE(:,n+j);
        DAE_vy(:,i) = DAE_vy(:,i) + p.L(j)*cos(solDAE(:,j)).*solDAE(:,n+j);        
        
        AMB_y(:,i) = AMB_y(:,i)  + p.L(j)*cos(solAMB(:,j));
        LAG_y(:,i) = LAG_y(:,i)  + p.L(j)*cos(solLAG(:,j));
        DAE_y(:,i) = DAE_y(:,i)  + p.L(j)*cos(solDAE(:,j));
    end
    % Add the velocities of the centers of mass with respect to the previous link
    AMB_vx(:,i) = AMB_vx(:,i) + p.a(i)*cos(solAMB(:,i)).*solAMB(:,n+i);
    AMB_vy(:,i) = AMB_vy(:,i) + p.a(i)*sin(solAMB(:,i)).*solAMB(:,n+i);

    LAG_vx(:,i) = LAG_vx(:,i) + p.a(i)*cos(solLAG(:,i)).*solLAG(:,n+i);
    LAG_vy(:,i) = LAG_vy(:,i) + p.a(i)*sin(solLAG(:,i)).*solLAG(:,n+i);

    DAE_vx(:,i) = DAE_vx(:,i) + p.a(i)*cos(solDAE(:,i)).*solDAE(:,n+i);
    DAE_vy(:,i) = DAE_vy(:,i) + p.a(i)*sin(solDAE(:,i)).*solDAE(:,n+i);

	AMB_y(:,i) = AMB_y(:,i) + p.a(i)*cos(solAMB(:,i));
    
    LAG_y(:,i) = LAG_y(:,i) + p.a(i)*cos(solLAG(:,i));
    
    DAE_y(:,i) = DAE_y(:,i) + p.a(i)*cos(solDAE(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Kinetic Energy
T_AMB = zeros(a,1);
T_LAG = zeros(a,1);
T_DAE = zeros(a,1);

for i = 1:b
    p.I(i) = 1/12*p.m(i)*p.L(i)^2;
    T_AMB = T_AMB + 0.5*abs(p.I(i).*solAMB(:,n+i).^2) + 0.5*abs(p.m(i)*(AMB_vx(:,i).*AMB_vy(:,i)));
    T_LAG = T_LAG + 0.5*abs(p.I(i).*solLAG(:,n+i).^2) + 0.5*abs(p.m(i)*(LAG_vx(:,i).*LAG_vy(:,i)));    
	T_DAE = T_DAE + 0.5*abs(p.I(i).*solDAE(:,n+i).^2) + 0.5*abs(p.m(i)*(DAE_vx(:,i).*DAE_vy(:,i)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Potential energy
V_AMB = zeros(a,1);
V_LAG = zeros(a,1);
V_DAE = zeros(a,1);

for i = 1:b
    V_AMB = V_AMB + p.m(i)*p.g*(sum(p.L)-AMB_y(:,i));
    V_LAG = V_LAG + p.m(i)*p.g*(sum(p.L)-LAG_y(:,i));    
    V_DAE = V_DAE + p.m(i)*p.g*(sum(p.L)-DAE_y(:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Total energy
E_AMB = V_AMB + T_AMB;
E_LAG = V_LAG + T_LAG;
E_DAE = V_DAE + T_DAE;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error in total energy
figure
clf
plot(T1,E_DAE-E_AMB,T1,E_DAE-E_LAG,T1,E_LAG-E_AMB,'LineWidth',3)
title('Error in Total Energies')
legend({'E_{DAE}-E_{AMB}','E_{DAE}-E_{LAG}','E_{LAG}-E_{AMB}'})
ylabel('difference')
xlabel('$$time$$ (s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error in potential energy
figure
clf
plot(T1,V_DAE-V_AMB,T1,V_DAE-V_LAG,T1,V_LAG-V_AMB,'LineWidth',3)
title('Error in Total Potential Energies')
legend({'U_{DAE}-U_{AMB}','U_{DAE}-U_{LAG}','U_{LAG}-U_{AMB}'})
ylabel('difference')
xlabel('$$time$$ (s)');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error in kinetic energy
figure
clf
plot(T1,T_DAE-T_AMB,T1,T_DAE-T_LAG,T1,T_LAG-T_AMB,'LineWidth',3)
title('Error in Total Kinetic Energies')
ylabel('difference')
legend({'K_{DAE}-K_{AMB}','K_{DAE}-K_{LAG}','K_{LAG}-K_{AMB}'})
xlabel('$$time$$ (s)');