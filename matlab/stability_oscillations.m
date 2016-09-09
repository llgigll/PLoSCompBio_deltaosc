function [ stability ] = stability_oscillations(q, k_hold, w_hold )
%stability_qdiff: find the resonant frequency at instability
%   Detailed explanation goes here

stability = struct([]);
qdiff_low = zeros(length(k_hold),length(w_hold));
osc_low = zeros(length(k_hold),length(w_hold));
qdiff_osc = zeros(length(k_hold),length(w_hold));
qdiff_rise = zeros(length(k_hold),length(w_hold));

% q must be 0.5 or less for this to work
q_low = -q;

for i = 1:length(k_hold)
    k = k_hold(i);
    for j = 1:length(w_hold)
        
        w = w_hold(j);
        fun = @(qdiff) comp_roots(qdiff,k,w,q);

        if fun(q_low) < 0
            qdiff_low(i,j) = 0; % If system is stable for smallest value of delta q
        else
            qdiff_low(i,j) = fzero(fun, [0,q_low]);
        end
        
        osc_low(i,j) =  freq_low(qdiff_low(i,j),k,w,q); % Resonant Frequency at instability
        
        
      
        
        fun = @(qdiff) num_poles(qdiff,q,w,k);
        holder = linspace(qdiff_low(i,j)+0.001*q,0,200);
        num_hold = zeros(length(holder),1);
        
        for ii = 1:length(holder)
           num_hold(ii) = fun(holder(ii)); 
        end
        f_high = holder(find(num_hold==0,1,'first'));
        f_low = qdiff_low(i,j)+0.001*q;
        if fun(f_low) <= 0
            qdiff_osc(i,j) = qdiff_low;
        else
            qdiff_osc(i,j) = fzero(fun, [f_low,f_high]);
        end
        
        qdiff_rise(i,j) = risetime(qdiff_osc(i,j),k,w,q);
        
    end
end

stability(1).qdiff_low = qdiff_low;
stability(1).osc_low = osc_low;
stability(1).qdiff_osc = qdiff_osc;
stability(1).qdiff_rise = qdiff_rise;

end

 function rise = risetime(diff,k,w,q)
    te  = 0.02;
    ti  = 0.01;
    tee_a = 0.005;
    tee_n = 0.1;
    tie_a = 0.005;
    tie_n = 0.1;
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
    B = [1;0;0;0;0;0;0;0];
    C = [0,0,(q-qdiff1)*(Jee+fpos)/te,(q+qdiff1)*(Jee+fpos)/te,0,0,-Jei/te,0];
    D = 0;

    sys1 = ss(A,B,C,D);
    S = stepinfo(sys1);
    rise = S.RiseTime;
 end

function num = num_poles(diff,q,w,k)
te  = 0.02;
ti  = 0.01;
tee_a = 0.005;
tee_n = 0.1;
tie_a = 0.005;
tie_n = 0.1;
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
B = [1;0;0;0;0;0;0;0];
C = [0,0,(q-qdiff1)*(Jee+fpos)/te,(q+qdiff1)*(Jee+fpos)/te,0,0,-Jei/te,0];
D = 0;

sys1 = ss(A,B,C,D);
poles1 = pole(sys1);
num = sum(imag(poles1)~=0)-2;
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

function freq = freq_low(qdiff,k,w,q)

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
poles = pole(sys);

[m,ind] = max(real(poles));
freq = abs(imag(poles(ind)))/(2*pi);

end