function [ stability ] = stability_qdiff( q, k_hold, w_hold )
%stability_qdiff: find the range of q for which the system is stable
%   Detailed explanation goes here

stability = struct('qdiff_low',{},'qdiff_high',{},'rise_time',{});
qdiff_low = zeros(length(k_hold),length(w_hold));
qdiff_high = zeros(length(k_hold),length(w_hold));
rise = zeros(length(k_hold),length(w_hold));
min_rise = zeros(length(k_hold),length(w_hold));
max_rise = zeros(length(k_hold),length(w_hold));

% q must be 0.5 or less for this to work
q_high = (1-q);
q_low = -q;

for i = 1:length(k_hold)
    k = k_hold(i);
    for j = 1:length(w_hold)
        
        w = w_hold(j);
        fun = @(qdiff) comp_roots(qdiff,k,w,q);

        if fun(q_low) < 0
            qdiff_low(i,j) = 0;
        else
            qdiff_low(i,j) = fzero(fun, [0,q_low]);
        end
        
        if fun(q_high) < 0
            qdiff_high(i,j) = q_high;
        else
            qdiff_high(i,j) = fzero(fun, [0,q_high]);
        end
        
        fun2 = @(qdiff) rise_time(qdiff,k,w,q);
        fun3 = @(qdiff) (-1*(rise_time(qdiff,k,w,q)+0.1));
        fun4 = @(qdiff) settling(qdiff,k,w,q);
        diff_high = qdiff_high(i,j)-abs(qdiff_high(i,j)*0.0001);
%         
%         [x,fval] = fminbnd(fun3,0,diff_high);
%         max_rise(i,j) = -fval;
%         if (k == 1.5) & (w==50)
%             fval
%         end
%         
%         hold = qdiff_high(i,j)*0.95:0.0001:qdiff_high(i,j);
%         holder = zeros(length(hold),1);
%         for m = 1:length(hold)
%             holder(m) = fun4(hold(m));
%         end
%         holder(isnan(holder)) = 0;
%         max_rise(i,j) = max(holder);

        max_rise(i,j) = rise_time(0.075,k,w,q)+0.100;
        
        diff_low = qdiff_low(i,j)+abs(qdiff_low(i,j)*0.0001);
        min_rise(i,j) = fun2(diff_low);
        if fun2(diff_high) < 0 || fun2(diff_low) > 0
            rise(i,j) = -0.3;
        else
            rise(i,j) = fzero(fun2, [diff_low,diff_high]);
        end
    end
end

stability(1).qdiff_low = qdiff_low;
stability(1).qdiff_high = qdiff_high;
stability(1).rise = rise;
stability(1).min_rise = min_rise;
stability(1).max_rise = max_rise;

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
r_time = S.RiseTime-0.100; % Compute difference from 100 ms rise time

end



function [settling_time] = settling(qdiff,k,w,q)

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

S = stepinfo(sys,'SettlingTimeThreshold',0.05);
settling_time = S.SettlingTime; % Compute difference from 100 ms rise time

end