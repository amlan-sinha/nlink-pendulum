close all
clear all
clc

set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultTextInterpreter','Latex')

% Amlan Sinha
% Cornell University

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           System parameters                           %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of links
n = 3;
% Initial conditions
th0  = ones(1,n)*pi/2;
dth0 = zeros(1,n);
% System parameters
p = struct('n',n,'g',1,'m',ones(1,n),'a',0.5*ones(1,n),'L',ones(1,n));
% Time span
time = [0 10];
% Flags
rederive = 1;
simulate = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                Solution                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T,solAMB] = NLinkPend_AMB(n,rederive,simulate,p,[th0 dth0],time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                Graphics                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                                  Theta                                  %

figure
plot(T,solAMB(:,1:3));
title('Solutions: AMB');
xlabel('$$time$$ (s)');
ylabel('$$\theta$$ (rad)');
legend(s,'Location','northeast');
hold off