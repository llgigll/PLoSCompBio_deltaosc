% Produces the plots for Figure 3

u_hold = {[0.17,0.23],[0.2,0.2],[0.275,0.125]};
t_hold = {[0.5,0.5],[0.5,0.5],[0.5,0.5]};
% u_hold = {[0.2,0.2],[0.2,0.2],[0.2,0.2]};
% t_hold = {[0.5,0.5],[0.6,0.4],[0.8,0.2]};
sol = cell(length(u_hold),1);

fpos = 0;
k = 1.2;
I = 10;
decay_tol = 0.1;
w = 50;
q = 0.5;
Io_begin = 0.0;
Io_end = 1.0;
t_end = 2.0;

for i = 1:length(u_hold)
    u = u_hold{i};
    t1 = t_hold{i}(1);
    t2 = t_hold{i}(2);

    [sol{i},end_val] = balcirc_std4(q,u,fpos,w,k,t1,t2,t_end,I,Io_begin,Io_end,decay_tol,0);
end

figure('Color','w')
for i = 1:length(u_hold)
    ye = sol{i}.y(1,:);
    T = sol{i}.x;
    plot(T,ye,'LineWidth',2)
    hold on
end
xlabel('Time(s)','FontSize',30)
ylabel('Activity(Hz)','FontSize',30)
legend('{\Delta}u = -0.03','{\Delta}u = 0.0','{\Delta}u = 0.075')

trajectories( k,w,q,[1,w],[-0.2,0.5],sol )