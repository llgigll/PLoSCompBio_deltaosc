function [ stability ] = stability_oscillations_tn(t_n_hold, q_hold, k, w )
%stability_qdiff: find the resonant frequency at instability
%   Allows changing of the NMDA time constant.

stability = struct([]);
qdiff_low = zeros(length(t_n_hold),length(q_hold));
qdiff_high = zeros(length(t_n_hold),length(q_hold));
osc_low = zeros(length(t_n_hold),length(q_hold));
qdiff_osc = zeros(length(t_n_hold),length(q_hold));
qdiff_rise = zeros(length(t_n_hold),length(q_hold));
max_rise = zeros(length(t_n_hold), length(q_hold));
slope = zeros(length(t_n_hold), length(q_hold));

% q must be 0.5 or less for this to work

for i = 1:length(t_n_hold)
    t_n = t_n_hold(i);
    for j = 1:length(q_hold)
        
        q = q_hold(j);
        q_low = -q; % q_low represents the minimum value of qdiff
        q_high = 1-q; % q_high represents the maximum value of qdiff
        fun1 = @(qdiff) comp_roots(t_n,qdiff,k,w,q);

        if fun1(q_low) < 0
            qdiff_low(i,j) = 0;
        else
            qdiff_low(i,j) = fzero(fun1, [0,q_low]);
        end
        
        osc_low(i,j) =  freq_low(t_n,qdiff_low(i,j),k,w,q);
        
        
        if fun1(q_high) < 0
            qdiff_high(i,j) = q_high;
        else
            qdiff_high(i,j) = fzero(fun1, [0,q_high]);

        end
        max_rise(i,j) = risetime(t_n,0.9*qdiff_high(i,j),k,w,q);
        rise_at_zero = risetime(t_n,0,k,w,q);
        slope(i,j) = (max_rise(i,j)-rise_at_zero)./(0.9*qdiff_high(i,j));
      
        
        fun2 = @(qdiff) num_poles(t_n,qdiff,q,w,k);

        if (qdiff_low(i,j)==0) && (fun2(q_low) < 0)
            qdiff_osc(i,j) = q_low;
                
        elseif (qdiff_low(i,j)==0) && (fun2(q_low) > 0)
            
            qdiff_osc(i,j) = fzero(fun2, [q_low,0]);
        else
            qdiff_osc(i,j) = fzero(fun2, [qdiff_low(i,j),0]);
        end
         
        
        qdiff_rise(i,j) = risetime(t_n,qdiff_osc(i,j),k,w,q);
        
    end
end

stability(1).qdiff_low = qdiff_low;
stability(1).osc_low = osc_low;
stability(1).qdiff_osc = qdiff_osc;
stability(1).qdiff_rise = qdiff_rise;
stability(1).slope = slope;

end

 function rise = risetime(t_n,diff,k,w,q)
    te  = 0.02;
    ti  = 0.01;
    tee_a = 0.005;
    tee_n = t_n;
    tie_a = 0.005;
    tie_n = t_n;
    tei = 0.01;
    tii = 0.01;
    qdiff1 = diff;
    qdiff2 = 0;
    Jee = w;
    Jei = k*w;
    Jie = w;
    Jii = k*w;
    fpos = 0;

    A = [-1/te,0,((1-q)-qdiff1)*(Jee+fpos)/te,(q+qdiff1)*(Jee+fpos)/te,0,0,-Jei/te,0;...
         0,-1/ti,0,0,((1-q)-qdiff2)*Jie/ti,(q+qdiff2)*Jie/ti,0,-Jii/ti;...
         1/tee_a,0,-1/tee_a,0,0,0,0,0;...
         1/tee_n,0,0,-1/tee_n,0,0,0,0;...
         1/tie_a,0,0,0,-1/tie_a,0,0,0;...
         1/tie_n,0,0,0,0,-1/tie_n,0,0;...
         0,1/tei,0,0,0,0,-1/tei,0;...
         0,1/tii,0,0,0,0,0,-1/tii];
    B = [1/te;0;0;0;0;0;0;0];
%     C = [0,0,(q-qdiff1)*(Jee+fpos)/te,(q+qdiff1)*(Jee+fpos)/te,0,0,-Jei/te,0];
    C = [1,0,0,0,0,0,0,0];
    D = 0;

    sys1 = ss(A,B,C,D);
    S = stepinfo(sys1);
    rise = S.RiseTime;
 end

function num = num_poles(t_n,diff,q,w,k)
te  = 0.02;
ti  = 0.01;
tee_a = 0.005;
tee_n = t_n;
tie_a = 0.005;
tie_n = t_n;
tei = 0.01;
tii = 0.01;
qdiff1 = diff;
qdiff2 = 0;
Jee = w;
Jei = k*w;
Jie = w;
Jii = k*w;
fpos = 0;

A = [-1/te,0,((1-q)-qdiff1)*(Jee+fpos)/te,(q+qdiff1)*(Jee+fpos)/te,0,0,-Jei/te,0;...
     0,-1/ti,0,0,((1-q)-qdiff2)*Jie/ti,(q+qdiff2)*Jie/ti,0,-Jii/ti;...
     1/tee_a,0,-1/tee_a,0,0,0,0,0;...
     1/tee_n,0,0,-1/tee_n,0,0,0,0;...
     1/tie_a,0,0,0,-1/tie_a,0,0,0;...
     1/tie_n,0,0,0,0,-1/tie_n,0,0;...
     0,1/tei,0,0,0,0,-1/tei,0;...
     0,1/tii,0,0,0,0,0,-1/tii];
B = [1/te;0;0;0;0;0;0;0];
% C = [0,0,(q-qdiff1)*(Jee+fpos)/te,(q+qdiff1)*(Jee+fpos)/te,0,0,-Jei/te,0];
C = [1,0,0,0,0,0,0,0];
D = 0;

sys1 = ss(A,B,C,D);
poles1 = pole(sys1);
% num = sum(imag(poles1)~=0)-2;
holder1 = abs(imag(poles1))>0;
holder2 = abs(imag(poles1))<60.00;
num = sum(holder1 & holder2)/2-0.5;
end


function [real_root] = comp_roots(t_n,qdiff,k,w,q)

te  = 0.02;
ti  = 0.01;
tee_a = 0.005;
tee_n = t_n;
tie_a = 0.005;
tie_n = t_n;
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

B = [1/te;0;0;0;0;0;0;0];

% C = [0,0,(q-qdiff1)*(Jee)/te,(q+qdiff1)*(Jee)/te,0,0,-Jei/te,0];
C = [1,0,0,0,0,0,0,0];

D = 0;

sys = ss(A,B,C,D);

real_root = max(real(pole(sys)));

end

function freq = freq_low(t_n,qdiff,k,w,q)

te  = 0.02;
ti  = 0.01;
tee_a = 0.005;
tee_n = t_n;
tie_a = 0.005;
tie_n = t_n;
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

B = [1/te;0;0;0;0;0;0;0];

% C = [0,0,(q-qdiff1)*(Jee)/te,(q+qdiff1)*(Jee)/te,0,0,-Jei/te,0];
C = [1,0,0,0,0,0,0,0];

D = 0;

sys = ss(A,B,C,D);
poles = pole(sys);

[m,ind] = max(real(poles));
freq = abs(imag(poles(ind)))/(2*pi);

end