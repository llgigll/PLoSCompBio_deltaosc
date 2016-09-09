function [] = trajectories( k,w,q,we_span,diff_span, std_vals )
%Trajectories: Produce trajectories for a set of simulations
%   

we = linspace(we_span(1),we_span(2),100);
q_low = diff_span(1);
q_high = diff_span(2);

qdiff_low = zeros(numel(we),1);
qdiff_high = zeros(numel(we),1);
qdiff_100 = zeros(numel(we),1);
qdiff_1 = zeros(numel(we),1);

for i = 1:numel(we)
    new_k = k*w/we(i);
    new_w = we(i);
    fun = @(diff) comp_roots(diff,new_k,new_w,q);
    
    if fun(q_low) < 0
        qdiff_low(i) = NaN;
    else
        qdiff_low(i) = fzero(fun, [q_low,0]);
    end

    if fun(q_high) < 0
        qdiff_high(i) = NaN;
    else
        qdiff_high(i) = fzero(fun, [0,q_high]);
    end
    
    fun = @(diff) (rise_time(diff,new_k,new_w,q) - 0.100);
    if isnan(qdiff_low(i)) && (fun(q_low)<0)
        qdiff_100(i) = fzero(fun,[q_low,q_high-q_low*0.01]);
    elseif isnan(qdiff_high(i)) && (fun(qdiff_low(i)-qdiff_low(i)*0.01)<0)
        qdiff_100(i) = fzero(fun,[qdiff_low(i)-qdiff_low(i)*0.01,q_high+qdiff_low(i)*0.01]);
    elseif fun(qdiff_low(i)-qdiff_low(i)*0.01)<0
        qdiff_100(i) = fzero(fun,[qdiff_low(i)-qdiff_low(i)*0.01,qdiff_high(i)+qdiff_low(i)*0.01]);
    else
        qdiff_100(i) = NaN;
    end
    
    fun = @(diff) (rise_time(diff,new_k,new_w,q) - 1.0);
    if isnan(qdiff_high(i)) && (fun(q_high)>0)
        qdiff_1(i) = fzero(fun,[0,q_high]);
    elseif isnan(qdiff_high(i))
        qdiff_1(i) = NaN;
    elseif (fun(qdiff_high(i)-qdiff_high(i)*0.001)>0) && (fun(0) < 0)
        qdiff_1(i) = fzero(fun,[0,qdiff_high(i)-qdiff_high(i)*0.001]);
    else
        qdiff_1(i) = NaN;
    end
end

ind_low = find(~isnan(qdiff_low));
ind_high = find(~isnan(qdiff_high));
ind_100 = find(~isnan(qdiff_100));
ind_1 = find(~isnan(qdiff_1)&(qdiff_1<0.2));


traj_we = cell(length(std_vals),1);
traj_diff = cell(length(std_vals),1);
traj_rise = cell(length(std_vals),1);
for i = 1:length(std_vals) 
    x1 = std_vals{i}.y(9,:);
    x2 = std_vals{i}.y(10,:);
    traj_we{i} = x1*q*w+x2*(1-q)*w;
    traj_diff{i} = x2-x1;
    traj_rise{i} = zeros(length(x1),1);
   
    for j = 1:length(x1)
        new_k = k*w/traj_we{i}(j);
        new_w = traj_we{i}(j);
        new_diff = traj_diff{i}(j);
        traj_rise{i}(j) = rise_time(new_diff,new_k,new_w,q);
    end
end





figure('Color','w')
for i = 1:length(std_vals)
    plot(traj_diff{i},traj_we{i},'LineWidth',2)
    hold on
end
area(qdiff_low(ind_low),we(ind_low),50.0);
area(qdiff_high(ind_high),we(ind_high),50.0);
plot(qdiff_100(ind_100),we(ind_100),'k-.');
plot(qdiff_1(ind_1),we(ind_1),'k-.');
xlabel('{\Delta}q','FontSize',30)
ylabel('W','FontSize',30)
title('Trajectories')
legend('AMPA Instability','NMDA Instability','Rise Time 0.1 (s)','Rise Time 1.0 (s)','{\Delta}u = -0.03','{\Delta}u = 0.0','{\Delta}u = 0.075')

figure('Color','w')
for i = 1:length(std_vals)
   plot(std_vals{i}.x,traj_rise{i},'LineWidth',2)
   hold on
end
xlabel('Time (s)','FontSize',30)
ylabel('Rise Time (s)','FontSize',30)
legend('{\Delta}u = -0.03','{\Delta}u = 0.0','{\Delta}u = 0.075')


end

function [real_root] = comp_roots(qdiff,k,w,q)

te  = 0.02;
ti  = 0.01;
tee_a = 0.005;
tee_n = 0.1;
tie_a = 0.005;
tie_n = 0.1;
tei = 0.01;
tii = 0.01;

qdiff1 = qdiff;
qdiff2 = 0;

Jee = w;
Jei = k*w;
Jie = w;
Jii = k*w;

A = [-1/te,0,((1-q)-qdiff1)*(Jee)/te,(q+qdiff1)*(Jee)/te,0,0,-Jei/te,0;...
     0,-1/ti,0,0,((1-q)-qdiff2)*Jie/ti,(q+qdiff2)*Jie/ti,0,-Jii/ti;...
     1/tee_a,0,-1/tee_a,0,0,0,0,0;...
     1/tee_n,0,0,-1/tee_n,0,0,0,0;...
     1/tie_a,0,0,0,-1/tie_a,0,0,0;...
     1/tie_n,0,0,0,0,-1/tie_n,0,0;...
     0,1/tei,0,0,0,0,-1/tei,0;...
     0,1/tii,0,0,0,0,0,-1/tii];

B = [1;0;0;0;0;0;0;0];

C = [0,0,(q-qdiff1)*(Jee)/te,(q+qdiff1)*(Jee)/te,0,0,-Jei/te,0];

D = 0;

sys = ss(A,B,C,D);

real_root = max(real(pole(sys)));

end




function [r_time] = rise_time(qdiff,k,w,q)

te  = 0.02;
ti  = 0.01;
tee_a = 0.005;
tee_n = 0.1;
tie_a = 0.005;
tie_n = 0.1;
tei = 0.01;
tii = 0.01;

qdiff1 = qdiff;
qdiff2 = 0;

Jee = w;
Jei = k*w;
Jie = w;
Jii = k*w;

A = [-1/te,0,((1-q)-qdiff1)*(Jee)/te,(q+qdiff1)*(Jee)/te,0,0,-Jei/te,0;...
     0,-1/ti,0,0,((1-q)-qdiff2)*Jie/ti,(q+qdiff2)*Jie/ti,0,-Jii/ti;...
     1/tee_a,0,-1/tee_a,0,0,0,0,0;...
     1/tee_n,0,0,-1/tee_n,0,0,0,0;...
     1/tie_a,0,0,0,-1/tie_a,0,0,0;...
     1/tie_n,0,0,0,0,-1/tie_n,0,0;...
     0,1/tei,0,0,0,0,-1/tei,0;...
     0,1/tii,0,0,0,0,0,-1/tii];

B = [1;0;0;0;0;0;0;0];

C = [0,0,(q-qdiff1)*(Jee)/te,(q+qdiff1)*(Jee)/te,0,0,-Jei/te,0];

D = 0;

sys = ss(A,B,C,D);

S = stepinfo(sys);
r_time = S.RiseTime;

end









