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
n = 20;
% Initial conditions
th0  = linspace(pi/2+0.25,pi/2-0.25,n);
dth0 = zeros(1,n);
% System parameters
p = struct('n',n,'g',1,'m',ones(1,n),'a',0.5*ones(1,n),'L',ones(1,n));
% Time span
time = [0 10];
% Flags
rederive  = 1;
simulate  = 1;
constrain = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                Solution                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[T1,solDAE] = NBarLink_DAE(n,rederive,constrain,simulate,p,[th0 dth0],time);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                                Graphics                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%                                  Theta                                  %

figure
hold on
s = cell(1,n);
for i = 1:0.5*length(solDAE(1,:))
    plot(T1,solDAE(:,i),'LineWidth',3)
    s{i} = sprintf('Link %2.0f',i);
end
plot(T1,solDAE(:,1:3));
title('Solutions: DAE');
xlabel('$$time$$ (s)');
ylabel('$$\theta$$ (rad)');
legend(s,'Location','northeast');