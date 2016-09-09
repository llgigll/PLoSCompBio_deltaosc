 function [sol,end_val] = balcirc_std4(q,u,fpos,w,k,t1,t2,t_end,I,Io_begin,Io_end,decay_tol,stable_low)
%% Balanced Circuitry for Persistent Activity
%  This program implements the balanced circuitry from Lim and Goldman
%  (2013) where the excitatory inputs are split between NMDA and AMPA
%  receptors. STD is also implemented and applied to AMPA and NMDA
%  dominated synapses. This is a reparameterization of balcirc_std.m

time_const = 'diff';
Io_type = 'const';
end_val = 0; % Tells event solver when to end the computation

% net_type = 'positive'; % amp, positive, derivative or hybrid
% time_const = 'diff'; % same or diff
% Io_type = 'const'; % pulse or const or const1 or osc
Jee = w+fpos;
Jei = k*w;
Jie = w;
Jii = k*w;
Jo  = [I;0;0;0;0;0;0;0;1;1];
% q = 0.5;
qee_a = 1-q;
qee_n = q;
% q=0.5;
qie_a = 1-q;
qie_n = q;


if strcmp(time_const,'diff')
    te  = 0.02;
    ti  = 0.01;
    tee_a = 0.005;
    tee_n = 0.1;
    tie_a = 0.005;
    tie_n = 0.1;
    tei = 0.01;
    tii = 0.01;
elseif strcmp(time_const,'same')
    te  = 0.02;
    ti  = 0.02;
    tee_a = 0.01;
    tee_n = 0.1;
    tie_a = 0.01;
    tie_n = 0.1;
    tei = 0.02;
    tii = 0.02;
end

std_off = 0;
u1 = u(1);
u2 = u(2);


tau = [1/te;1/ti;1/tee_a;1/tee_n;1/tie_a;1/tie_n;1/tei;1/tii;1/t1;1/t2];
t0 = 0;
y0 = [0.00001;zeros(7,1);1;1]; 

A = [-1/te,0,qee_a*Jee/te,qee_n*Jee/te,0,0,-Jei/te,0,0,0;...
     0,-1/ti,0,0,qie_a*Jie/ti,qie_n*Jie/ti,0,-Jii/ti,0,0;...
     1/tee_a,0,-1/tee_a,0,0,0,0,0,0,0;...
     1/tee_n,0,0,-1/tee_n,0,0,0,0,0,0;...
     1/tie_a,0,0,0,-1/tie_a,0,0,0,0,0;...
     1/tie_n,0,0,0,0,-1/tie_n,0,0,0,0;...
     0,1/tei,0,0,0,0,-1/tei,0,0,0;...
     0,1/tii,0,0,0,0,0,-1/tii,0,0;...
     -u1,0,0,0,0,0,0,0,-1/t1,0;...
     -u2,0,0,0,0,0,0,0,0,-1/t2];

t_past = [];
tspan = [t0 t_end]; 
options = odeset('MaxStep',0.001);%,'NonNegative',1:10);%,'Events',@events);
sol = ode45(@odefun,tspan,y0,options);


function dy = odefun(t,y)
    if t > Io_end(end) && end_val == 0
        end_val = y(1); 
    end
    
    if strcmp(Io_type,'pulse')
        Io = input_pulse(t);
    elseif strcmp(Io_type,'const')
        Io = input_const(t);
    elseif strcmp(Io_type,'const1')
        Io = input_const1(t);
    elseif strcmp(Io_type,'osc')
        Io = input_osc(t);
    end
    
    if std_off
        y(end-1:end) = 1;
    end
    
    B = A;
%     B(3,1) = A(3,1) * y(9);
%     B(4,1) = A(4,1) * y(10);
%     B(5,1) = A(5,1) * (y(10)+y(9))/2;
%     B(6,1) = A(6,1) * (y(10)+y(9))/2;
    
    B(1,3) = A(1,3) * y(9);
    B(1,4) = A(1,4) * y(10);
    B(2,5) = A(2,5) * (y(10)+y(9))/2;
    B(2,6) = A(2,6) * (y(10)+y(9))/2;
    
    
%     B(2,5) = A(2,5) * y(10);
%     B(2,6) = A(2,6) * y(9);
    
    B(9,1) = A(9,1) * y(9);
    B(10,1) = A(10,1) * y(10);

    dy =  B * y + tau .* Jo .* Io ;
    t_past = [t_past,t];
end

function [value,isterminal,direction] = events(t,y)
    Io = input_const(t);
    drdt = 1/te * (-y(1)+u1*y(9)*(1-qee)*Jee*y(3)+u2*y(10)*qee*Jee*y(4)...
        -Jei*y(7)+Io(1)*Jo(1));
    if end_val == 0;
        r = 1;
    else
        r = y(1) - decay_tol * (end_val-stable_low);
    end
    value = [drdt;r];
    isterminal = [0;1];   % Stop at local minimum
    direction = [-1;-1];  % [local minimum, local maximum]
end 

function Io = input_pulse(t)
    x = @(s)1*(exp(-s/0.1)-exp(-s/0.02));

    if t >= Io_begin
        Io = [x(t-Io_begin);x(t-Io_begin);0;0;0;0;0;0;1;1];
    else
        Io = [zeros(8,1);1;1];
    end
end

function Io = input_const(t)
    Io  = zeros(10,1);
    Io(9:10) = 1;
    for p = 1:length(Io_end)
        if t >= Io_begin(p) && t <= Io_end(p)
            Io(1) = Io(1)+(sigmf(t,[50,Io_begin(p)])-0.5)*2;
        elseif t > Io_end(p)
            Io(1) = Io(1)+sigmf(t,[-50,Io_end(p)])*2;
        end
    end
end

function Io = input_const1(t)
    if t >= Io_begin && t <= Io_end
        Io = [sigmf(t,[200,Io_begin+0.045]);sigmf(t,[200,Io_begin+0.045]);0;0;0;0;0;0;1;1];
    elseif t > Io_end
        Io = [1.5*sigmf(t,[200,Io_end+0.045]);sigmf(t,[200,Io_end+0.045]);0;0;0;0;0;0;1;1];
    end
end

function Io = input_osc(t)
    if t >= Io_begin && t <= Io_end
        Io = [sigmf(t,[200,Io_begin+0.045])*sin(t*40)^2;sigmf(t,[200,Io_begin+0.045]);0;0;0;0;0;0;1;1];
    elseif t > Io_end
        Io = [sigmf(t,[-200,Io_end+0.045]);sigmf(t,[-200,Io_end+0.045]);0;0;0;0;0;0;1;1];
    end
end

end


















