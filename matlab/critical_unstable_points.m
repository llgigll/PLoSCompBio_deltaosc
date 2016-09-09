function [qdiff_osc,qdiff_low,qdiff_high] = critical_unstable_points(q,w,k)
% Find the delta oscillatory bifurcation and the unstable points as a
% function of qdiff for the eigth order system.


q_high = q;
q_low = -q;

fun = @(qdiff) comp_roots(qdiff,k,w,q);

% Find the delta oscillatory instability
if fun(q_low) < 0
    qdiff_low = 0;
else
    qdiff_low = fzero(fun, [q_low,0]);
end

% Find the gamma oscillatory instability
if fun(q_high) < 0
    qdiff_high = q_high;
else
    qdiff_high = fzero(fun, [0,q_high]);
end

% Find bifurcation near the oscillatory instability
fun = @(diff) num_poles(diff,q,w,k);
diff_val_hold = linspace(qdiff_low+0.001*q,qdiff_high-0.001*q,100);

holder = linspace(qdiff_low+0.001*q,0,100);
num_hold = zeros(length(holder),1);
for i = 1:length(holder)
   num_hold(i) = fun(holder(i)); 
end
f_high = holder(find(num_hold==0,1,'first'));
f_low = diff_val_hold(1);
if fun(f_low) <= 0
    qdiff_osc = qdiff_low;
else
    qdiff_osc = fzero(fun, [f_low,f_high]);
end


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

function [real_root] = comp_roots(diff,k,w,q)
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



 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 