function [] = PLOS_figure1_rate()
%% Balanced Circuitry for Persistent Activity
%  Test an idea from my 08FEB2016 thesis committee meeting.  Can we
%  develop an experiment where we can separate out the impact of excitatory
%  inhibitory balance from temporal balance?

Io_type = 'const';
end_val = 0; % Tells event solver when to end the computation

I = 5;

q = 0.3;
w = 30.0;
k = 1.2;

te  = 0.02;
ti  = 0.01;
tee_a = 0.005;
tee_n = 0.1;
tie_a = 0.005;
tie_n = 0.1;
tei = 0.01;
tii = 0.01;

qdiff = [-0.017,0.0,0.1];


Jee = w;
Jei = k*w;
Jie = w;
Jii = k*w;
pv_hypo = -0.000;

qampa = 0.00;
qnmda = 0.00;

Jo  = [I;0;0;0;0;0;0;0];
tau = [1/te;zeros(7,1)];
Io_begin = 0;
Io_end = 5;
t0 = 0;
t_end = 3;
y0 = [0.00001;zeros(7,1)];

ye = cell(length(qdiff),1);
yi = cell(length(qdiff),1);
yt = cell(length(qdiff),1);
y_leg = cell(length(qdiff),1);
for i = 1:length(qdiff)
    qdiff1 = qdiff(i);
    qdiff2 = 0.0;
    A = [-1/te,0,((1-q)-qdiff1+qampa)*(Jee)/te,(q+qdiff1+qnmda)*(Jee)/te,0,0,-Jei/te,0;...
         0,-1/ti,0,0,((1-q)+pv_hypo-qdiff2)*Jie/ti,(q+qdiff2)*Jie/ti,0,-Jii/ti;...
         1/tee_a,0,-1/tee_a,0,0,0,0,0;...
         1/tee_n,0,0,-1/tee_n,0,0,0,0;...
         1/tie_a,0,0,0,-1/tie_a,0,0,0;...
         1/tie_n,0,0,0,0,-1/tie_n,0,0;...
         0,1/tei,0,0,0,0,-1/tei,0;...
         0,1/tii,0,0,0,0,0,-1/tii];

    [U,T]=schur(A)
    cond(A)
    rank(A)
    sprank(A)
    t_past = [];
    tspan = [t0 t_end]; 
    options = odeset('MaxStep',0.001,'NonNegative',1:2);
    [T,Y] = ode45(@odefun,tspan,y0,options);

    ye{i} = Y(:,1);
    yi{i} = Y(:,2);
    yt{i} = T;
    y_leg{i} = ['{\Delta}q = ', num2str(qdiff(i))];
end


figure('color','w')
for i = 1:length(ye)
    plot(yt{i},ye{i},'LineWidth',3.0)
    hold on
%     equilib = I/(1-w+(k*w^2/(1+k*w)));
%     plot([T(1),T(end)],[equilib,equilib],'k--')
end
xlabel('Time (s)','FontSize',20)
ylabel('Rate (Hz)','FontSize',20)
set(gca,'fontsize',16)
legend(y_leg{:})

% figure
% plot(T,Y(:,3)*((1-q)-qdiff1)*(Jee)+Y(:,4)*(q+qdiff1)*(Jee)-Y(:,7)*Jei-3)

% input_const(0.1)
% input_const(0.05)

% figure
% plot(yt{2},ye{2})
% hold on
% plot(yt{2},yi{2})
% xlabel('Time (s)')
% ylabel('Rate (Hz)')
% legend('E','I')


function dy = odefun(t,y)
    if t > Io_end && end_val == 0
        end_val = y(1); 
    end
    
    if strcmp(Io_type,'pulse')
        Io = input_pulse(t);
    elseif strcmp(Io_type,'const')
        Io = input_const(t);
    end
    B=A;
%     if t>1
%         qie = 0.5*(1-p);
%         qee = 0.5*(1+p);
%         B(1,3)=x1*(1-qee)*Jee/(te);
%         B(1,4)=x2*qee*Jee/(te);
%         B(2,5)=x2*(1-qie)*Jie/(ti);
%         B(2,6)=x1*qie*Jie/(ti);
%     end
    dy =  B * y + tau .* Jo .* Io ;
    t_past = [t_past,t];
end


function Io = input_pulse(t)
    x = @(s)1*(exp(-s/0.1)-exp(-s/0.02));

    if t >= Io_begin
        Io = [x(t-Io_begin);x(t-Io_begin);0;0;0;0;0;0];
    else
        Io = [zeros(8,1)];
    end
end


function Io = input_const(t)
    if t >= Io_begin && t <= Io_end
        Io = [(sigmf(t,[100,Io_begin])-0.5)*2;0;0;0;0;0;0;0];
    elseif t > Io_end
        Io = [sigmf(t,[-100,Io_end])*2;0;0;0;0;0;0;0];
    end
end

end


















