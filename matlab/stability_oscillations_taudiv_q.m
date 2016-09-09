function [ stability ] = stability_oscillations_taudiv_q(tau_div, q, k, w)
%stability_qdiff: find the resonant frequency at instability
%   Detailed explanation goes here

stability = struct('qdiff_low','osc_low');
qdiff_low = zeros(length(q),length(tau_div));
osc_low = zeros(length(q),length(tau_div));

% q must be 0.5 or less for this to work
for i = 1:length(q)
    q_low = -q(i);
    for j = 1:length(tau_div)
        
        fun = @(qdiff) comp_roots(qdiff,k,w,q(i),tau_div(j));

        if fun(q_low) < 0
            qdiff_low(i,j) = 0;
        else
            qdiff_low(i,j) = fzero(fun, [0,q_low]);
        end
        
        osc_low(i,j) =  freq_low(qdiff_low(i,j),k,w,q(i),tau_div(j));
    end
end

stability(1).qdiff_low = qdiff_low;
stability(1).osc_low = osc_low;


end

function [real_root] = comp_roots(qdiff,k,w,q,tau_div)

te  = 0.02*tau_div;
ti  = 0.01*tau_div;
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



function freq = freq_low(qdiff,k,w,q,tau_div)

te  = 0.02*tau_div;
ti  = 0.01*tau_div;
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