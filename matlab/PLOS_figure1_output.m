function [] = balcirc_basic_natneuroexplanation3()
%% Balanced Circuitry for Persistent Activity
%  Plot simulation results for all three versions of the model, from 2nd
%  order through 8th order.  There are four plots for the parameter sets
%  associated with values of the damping coefficient in the 2nd order
%  equation. We do plots for Overdamped, Critically Damped, Underdamped and
%  Un-damped systems.

timer = 0:0.001:0.08;
tlow = 0.005;
thigh = 0.01;
vals = 1/thigh*exp(-timer/thigh)-1/tlow*exp(-timer/tlow);
figure('color','w')
plot(timer,-vals/(max(abs(vals))),'b','LineWidth',3)
hold on
plot(timer,vals/(max(abs(vals))),'r','LineWidth',3)
plot([0,0.08],[0,0],'--k')
set(gca,'fontsize',16)
legend('Fast Excitation','Slow Excitation')
xlabel('Time (s)','fontsize',30)
ylabel('Normalized Recurrent Input','fontsize',30)


I = 1;

q = 0.3;
w = 30.0;
k = 1.5;

te  = 0.02;
t_a = 0.005;
t_n = 0.1;

colors =  [0,0.4470,0.7410;0.8500,0.3250,0.0980;0.9290,0.6940,0.1250];

undamped = -(1+te/(t_n))/w;
critically_damped = (2*sqrt(te*t_n)-(te+t_n))/(t_n)/w;
underdamped = undamped*0.87;
overdamped = 0.1;
qdiff_second = [undamped,underdamped,critically_damped, overdamped]

fun = @(x) crit_point_third(x, w, te, t_a, t_n);
omega = 1/sqrt(t_a*te+t_n*te+t_a*t_n)*1/(2*pi);
undamped_third = (t_n*t_a*te/((t_n-t_a)*(t_a*te+t_n*te+t_a*t_n))-(t_n+t_a+te)/(t_n-t_a))/w;
underdamped_third = undamped_third*0.8;
crit_damped_third = fzero(fun, [undamped_third*1.001,0]);
overdamped_third = 0.125;

qdiff_third = [undamped_third,underdamped_third,crit_damped_third, overdamped_third]

[crit_damped_eigth,undamped_eigth,qdiff_high]=critical_unstable_points(q,w,k);
qdiff_eigth = [undamped_eigth,undamped_eigth*0.87,crit_damped_eigth, 0.05]

Io_begin = 0;
Io_end = 5.0;
t0 = 0;
t_end = 1.0;

second_order(I, w, te, t_n, qdiff_second, Io_begin, Io_end, t0, t_end)
third_order(I, q, w, te, t_a, t_n, qdiff_third, Io_begin, Io_end, t0, t_end)
eigth_order(I, q, w, k, te, t_a, t_n, qdiff_eigth, Io_begin, Io_end, t0, t_end)


function [] = second_order(I, w, te, t_n, qdiff, Io_begin, Io_end, t0, t_end)

    tau_2 = [1/te;0];
    Jo_2  = [I/t_n;0];

    y0 = [0.0;0.0];

    ye = cell(length(qdiff),1);
    yt = cell(length(qdiff),1);
    for i = 1:length(qdiff)
        A_2 = [-(te+(1+w*qdiff(i))*t_n)/(te*t_n),-1/(te*t_n);...
             1,0];
        tspan = [t0 t_end]; 
        options = odeset('MaxStep',0.001);%,'NonNegative',1:2);
        [T,Y] = ode45(@odefun,tspan,y0,options);

        ye{i} = Y(:,2);
        yt{i} = T;
    end


    for i = 1:length(ye)
        figure('color','w','name','second order system')
        plot(yt{i},ye{i},'linewidth',3.0,'color',colors(1,:))
        xlabel('Time (s)','fontsize',30)
        ylabel('Rate (Hz)','fontsize',30)
        ylim(1.05*[min(ye{1}),max(ye{1})])
        set(gca,'fontsize',20)
    end



    function dy = odefun(t,y)
        Io = input_const(t);
        dy =  A_2 * y + tau_2 .* Jo_2 .* Io ;
    end


    function Io = input_const(t)
        if t >= Io_begin && t <= Io_end
            Io = [(sigmf(t,[100,Io_begin])-0.5)*2;0];
        elseif t > Io_end
            Io = [sigmf(t,[-100,Io_end])*2;0];
        end
    end
    
end



function [] = third_order(I, q, w, te, t_a, t_n, qdiff, Io_begin, Io_end, t0, t_end)
    Jee = w;
    Jie = w;

    Jo  = [I;0;0];
    tau = [1/te;zeros(2,1)];

    y0 = [0.00001;zeros(2,1)];

    ye = cell(length(qdiff),1);
    yt = cell(length(qdiff),1);
    for i = 1:length(qdiff)
        qdiff1 = qdiff(i);
        A = [-1/te,((1-q)*(Jee-Jie)-qdiff1*Jee)/te,(q*(Jee-Jie)+qdiff1*Jee)/te;...
             1/t_a,-1/t_a,0;...
             1/t_n,0,-1/t_n];

        t_past = [];
        tspan = [t0 t_end]; 
        options = odeset('MaxStep',0.001);%,'NonNegative',1:2);
        [T,Y] = ode45(@odefun,tspan,y0,options);

        ye{i} = Y(:,1);
        yt{i} = T;
    end


    for i = 1:length(ye)
        figure('color','w','name','third order system')
        plot(yt{i},ye{i},'linewidth',3.0,'color',colors(2,:))
        xlabel('Time (s)','fontsize',30)
        ylabel('Rate (Hz)','fontsize',30)
        ylim(1.05*[min(ye{1}),max(ye{1})])
        set(gca,'fontsize',20)
    end


    function dy = odefun(t,y)
        Io = input_const(t);
        dy =  A * y + tau .* Jo .* Io ;
        t_past = [t_past,t];
    end


    function Io = input_const(t)
        if t >= Io_begin && t <= Io_end
            Io = [(sigmf(t,[100,Io_begin])-0.5)*2;0;0];
        elseif t > Io_end
            Io = [sigmf(t,[-100,Io_end])*2;0;0];
        end
    end
end

function [] = eigth_order(I, q, w, k, te, t_a, t_n, qdiff, Io_begin, Io_end, t0, t_end)
%% Balanced Circuitry for Persistent Activity
%  This program implements the balanced circuitry from Lim and Goldman
%  (2013) where the excitatory inputs are split between NMDA and AMPA
%  receptors. No STD is used in this implementation.


ti  = 0.01;
tee_a = t_a;
tee_n = t_n;
tie_a = t_a;
tie_n = t_n;
tei = 0.01;
tii = 0.01;


Jee = w;
Jei = k*w;
Jie = w;
Jii = k*w;

Jo  = [I;0;0;0;0;0;0;0];
tau = [1/te;zeros(7,1)];

y0 = [0.00001;zeros(7,1)];

ye = cell(length(qdiff),1);
yt = cell(length(qdiff),1);

for i = 1:length(qdiff)
    qdiff1 = qdiff(i);

    A = [-1/te,0,((1-q)-qdiff1)*(Jee)/te,(q+qdiff1)*(Jee)/te,0,0,-Jei/te,0;...
         0,-1/ti,0,0,(1-q)*Jie/ti,q*Jie/ti,0,-Jii/ti;...
         1/tee_a,0,-1/tee_a,0,0,0,0,0;...
         1/tee_n,0,0,-1/tee_n,0,0,0,0;...
         1/tie_a,0,0,0,-1/tie_a,0,0,0;...
         1/tie_n,0,0,0,0,-1/tie_n,0,0;...
         0,1/tei,0,0,0,0,-1/tei,0;...
         0,1/tii,0,0,0,0,0,-1/tii];


    tspan = [t0 t_end]; 
    options = odeset('MaxStep',0.001);
    [T,Y] = ode45(@odefun,tspan,y0,options);

    ye{i} = Y(:,1);
    yt{i} = T;
end


for i = 1:length(ye)
    figure('color','w','name','eigth order system')
    plot(yt{i},ye{i},'linewidth',3.0,'color',colors(3,:))
    ylim(1.05*[min(ye{1}),max(ye{1})])
    xlabel('Time (s)','fontsize',30)
    ylabel('Rate (Hz)','fontsize',30)
    set(gca,'fontsize',20)
end

function dy = odefun(t,y)

    Io = input_const(t);

    dy =  A * y + tau .* Jo .* Io ;
end

function Io = input_const(t)
    if t >= Io_begin && t <= Io_end
        Io = [(sigmf(t,[100,Io_begin])-0.5)*2;0;0;0;0;0;0;0];
    elseif t > Io_end
        Io = [sigmf(t,[-100,Io_end])*2;0;0;0;0;0;0;0];
    end
end
end

function num_poles = crit_point_third(qdiff1, w, te, t_a, t_n)
    A = [-1/te,-qdiff1*w/te,qdiff1*w/te;...
         1/t_a,-1/t_a,0;...
         1/t_n,0,-1/t_n];

    B = [1/te;0;0];
    C = [1,0,0];
    D = 0;
    
    sys1 = ss(A,B,C,D);
    poles_sys = pole(sys1);
    
    imagpoles = imag(poles_sys)/(2*pi);
    num_poles = sum(imagpoles~=0)-1;

end

end
















