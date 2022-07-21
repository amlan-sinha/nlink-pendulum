% Amlan Sinha
% Cornell University

function [T,sol] = NLinkPend_DAE(n,rederive,simulate,p,Y0,time)

set(0,'DefaultFigureWindowStyle','docked')

if rederive == 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                           System parameters                           %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    g       = sym('g',           'real');
	a       = sym('a',    [1 n], 'positive');
    L       = sym('L',    [1 n], 'positive');
    m       = sym('m',    [1 n], 'positive');
    
    t       = sym('th',   [1 n], 'real');
    dt      = sym('dth',  [1 n], 'real');
    ddt     = sym('ddth', [1 n], 'real');
    
    ddx     = sym('ddx',  [1 n], 'real');
    ddy     = sym('ddy',  [1 n], 'real');
    
    I = (1/12)*m.*L.^2;

    % Constraint Forces
    Fx      = sym('Fx',[1 n], 'real');
    Fy      = sym('Fy',[1 n], 'real');
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                                 Frames                                %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Intertial Frame
    e1 = [1 0 0]';
    e2 = [0 1 0]';
    e3 = [0 0 1]';
    
    % Body Frames
    ea1 = cell(1,n);
    ea2 = cell(1,n);
    ea3 = cell(1,n);
    
    for i = 1:n
        ea1{i} = +cos(t(i))*e1+sin(t(i))*e2;
        ea2{i} = -sin(t(i))*e1+cos(t(i))*e2;
        ea3{i} = e3;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                               Kinematics                              %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Position Vector
    rG     = cell(1,n);
    rE     = cell(1,n);
    rErelO = cell(1,n);
    rGrelO = cell(1,n);
    
    rG{1}     = a(1)*L(1)*ea1{1};
    rE{1}     =      L(1)*ea1{1};
    rGrelO{1} = rG{1};
    rErelO{1} = rE{1};
    
    for i = 2:n
        rG{i}     = a(i)*L(i)*ea1{i};
        rE{i}     =      L(i)*ea1{i};
        rGrelO{i} = rErelO{i-1}+rG(i);
        rErelO{i} = rErelO{i-1}+rE{i};
    end
    
    % Velocity Vector
    vG     = cell(1,n);
    vE     = cell(1,n);
    vErelO = cell(1,n);
    vGrelO = cell(1,n);
    
    vG{1}     = cross(dt(1)*ea3{1},rG{1});
    vE{1}     = cross(dt(1)*ea3{1},rE{1});
    vGrelO{1} = vG{1};
    vErelO{1} = vE{1};
    
    for i = 2:n
        vG{i}     = cross(dt(i)*ea3{i},rG{i});
        vE{i}     = cross(dt(i)*ea3{i},rE{i});
        vGrelO{i} = vErelO{i-1}+vG(i);
        vErelO{i} = vErelO{i-1}+vE{i};
    end
    
    % Acceleration Vector
    aG     = cell(1,n);
    aE     = cell(1,n);
    aErelO = cell(1,n);
    aGrelO = cell(1,n);
    aGrelOCartesian = cell(1,n); 
    
    aG{1}     = cross(ddt(1)*ea3{1},rG{1}) + cross(dt(1)*ea3{1},cross(dt(1)*ea3{1},rG{1}));
    aE{1}     = cross(ddt(1)*ea3{1},rE{1}) + cross(dt(1)*ea3{1},cross(dt(1)*ea3{1},rE{1}));
    aGrelO{1} = aG{1};
    aErelO{1} = aE{1};
    
    aGrelOCartesian{1} = ddx(1)*e1 + ddy(1)*e2;
    
    for i = 2:n
        aG{i}     = cross(ddt(i)*ea3{i},rG{i}) + cross(dt(i)*ea3{i},cross(dt(i)*ea3{i},rG{i}));
        aE{i}     = cross(ddt(i)*ea3{i},rE{i}) + cross(dt(i)*ea3{i},cross(dt(i)*ea3{i},rE{i}));
        aGrelO{i} = aErelO{i-1}+aG(i);
        aErelO{i} = aErelO{i-1}+aE{i};
        aGrelOCartesian{i} = ddx(i)*e1 + ddy(i)*e2;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%                                Dynamics                               %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    F  = cell(1,n);
    Fg = cell(1,n);
    
    for i = 1:n
        F{i}  = Fx(i)*e1 + Fy(i)*e2;
        Fg{i} = m(i)*g*e1;
    end
    
    eq  = cell(1,5*n);
    
	% Linear Momentum Balance    
    lmb = cell(1,n);
    
    lmb{1} = m(1)*aGrelOCartesian{1}-(F{1}+F{2})-Fg{1};
    eq{1}  = dot(lmb{1},e1);
    eq{2}  = dot(lmb{1},e2);
    for i = 2:n-1
        lmb{i} = m(i)*aGrelOCartesian{i}-(-F{i}+F{i+1})-Fg{i};
        eq{2*i-1} = dot(lmb{i},e1);
        eq{2*i}   = dot(lmb{i},e2);
    end
    lmb{n} = m(n)*aGrelOCartesian{n}-(-F{n})-Fg{n};
    eq{2*n-1}  = dot(lmb{n},e1);
    eq{2*n}    = dot(lmb{n},e2);
    
    % Angular Momentum Balance
    amb       = cell(1,n);
    
    amb{1}    = I(1)*ddt(1) - (cross(-rGrelO{1},F{1})+cross(rErelO{1}-rGrelO{1},F{2}));
    eq{2*n+1} = dot(amb{1},e3);
    for i = 2:n-1
        amb{i} = I(i)*ddt(i) - (cross(rErelO{i-1}-rGrelO{i},- F{i})+cross(rErelO{i}-rGrelO{i},F{i+1}));
        eq{2*n+i} = dot(amb{i},e3);
    end
    amb{n}  = I(n)*ddt(n) - (cross(rErelO{n-1}-rGrelO{n},-F{n}));
    eq{3*n} = dot(amb{n},e3);
    
    % Constraint Equation
    con = cell(1,n);
    for i = 1:n
        con{i} = aGrelO{i}-aGrelOCartesian{i};
        eq{3*n+2*i-1} = dot(con{i},e1);
        eq{3*n+2*i}   = dot(con{i},e2);
    end
    
    [A,b] = equationsToMatrix(eq,[ddt,ddx,ddy,Fx,Fy]);
    
    matlabFunction(A,'file','NLinkPend_A_DAE','Vars',{g,a,L,m,t,dt});
    matlabFunction(b,'file','NLinkPend_b_DAE','Vars',{g,a,L,m,t,dt});

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                           Equation of motion                          %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function dydt = odeDAE(~,y,p)
        
        th   = y(1:p.n)';
        dth  = y(p.n+1:2*p.n)';
        P    = NLinkPend_A_DAE(p.g,p.a,p.L,p.m,th,dth);
        q    = NLinkPend_b_DAE(p.g,p.a,p.L,p.m,th,dth);
        x    = P\q;
        
        dydt = [dth, x(1:p.n)']';
    end

tVector = linspace(time(1),time(end),(time(end)-time(1))*100);

opt.relTol = 1e-10;
opt.absTol = 1e-10;

% Integration
[T,sol] = ode45(@odeDAE,tVector,Y0,opt,p);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                               Animation                               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if simulate == 1
    
    theta = sol(:,1:p.n);
    
    xE    = zeros(size(theta));
    yE    = zeros(size(theta));
    xG    = zeros(size(theta));
    yG    = zeros(size(theta));
    
    xE(:,1) =  p.L(1)*sin(theta(:,1));
    yE(:,1) = -p.L(1)*cos(theta(:,1));
    xG(:,1) =  p.a(1)*p.L(1)*sin(theta(:,1));
    yG(:,1) = -p.a(1)*p.L(1)*cos(theta(:,1));
    
    for i = 2:p.n
        xE(:,i) =  p.L(i)*sin(theta(:,i))+xE(:,i-1);
        yE(:,i) = -p.L(i)*cos(theta(:,i))+yE(:,i-1);
        xG(:,i) =  p.a(i)*p.L(i)*sin(theta(:,i))+xE(:,i-1);
        yG(:,i) = -p.a(i)*p.L(i)*cos(theta(:,i))+yE(:,i-1);
    end
    
    col1 = [1 1 0]';
    col2 = [1 0 0]';
    colr = zeros(3,p.n);
    
    for i = 1:p.n
        colr(:,i) = (col1*(p.n-i)+col2*i)/p.n;
    end
    
    idx  = 1;
    minx = 1.25*min(min(xE));
    maxx = 1.25*min(max(xE));
    miny = 1.25*min(min(yE));
    maxy = 1.25*min(max(yE));
    
    limx = max(max(abs(xE)));
    
    arm  = cell(1,p.n);
    traj = cell(1,p.n);
    leg    = cell(1,p.n);
    
    figure
    clf
    hold on
    
    for i = 1:length(xE(1,:))
        traj{i} = plot(xG(:,i),yG(:,i),'--','Color',colr(:,i),'LineWidth',1,'Visible','off');
        leg{i}    = sprintf('arm %2.0f',i);
    end
    
    plot([minx minx],[miny maxy],'w',[maxx maxx],[miny maxy],'w')
    quiver([0 2*limx],[0 0],[1,0],[0,0],'k');
    quiver([0 0],[0 -2*limx],[0,0],[-1,0],'k');
    
    arm{1} = plot([0 xE(idx,1)],[0 yE(idx,1)],'Color',colr(:,1)','LineWidth',3);
    
    for i = 2:p.n
        arm{i} = plot([xE(idx,i-1) xE(idx,i)],[yE(idx,i-1) yE(idx,i)],'Color',colr(:,i)','LineWidth',3);
    end
    
    axis equal
    box on
    set(gcf,'Color','w')
    hold off
    pause(.5)
    
    while idx < length(tVector)
        
        idx = idx+1;
        set(arm{1},'XData',[0 xE(idx,1)],'YData',[0 yE(idx,1)])
        
        for i = 2:p.n
            set(arm{i},'XData',[xE(idx,i-1) xE(idx,i)],'YData',[yE(idx,i-1) yE(idx,i)])
        end
        
        pause(eps)
        drawnow
        
    end
    
    for i = 1:length(xE(1,:))
        set(traj{i},'Visible','on')
    end
    
    title([num2str(n) '-Link Pendulum'])
    xlabel('x (m)');
    ylabel('y (m)');
    legend(leg,'Location','East')
    
end

end