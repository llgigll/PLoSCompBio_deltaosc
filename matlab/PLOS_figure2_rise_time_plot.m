function [] = rise_time_plot()
gca
hold on
q = 0.45;
w = 30;
k = 1.2;
global div
div = 1;

q_high = q;
q_low = -q;

fun = @(qdiff) comp_roots(qdiff,k,w,q);

if fun(q_low) < 0
    qdiff_low = 0;
else
    qdiff_low = fzero(fun, [q_low,0]);
end

if fun(q_high) < 0
    qdiff_high = q_high;
else
    qdiff_high = fzero(fun, [0,q_high]);
end


diff_val_hold = linspace(qdiff_low+0.001*q,qdiff_high-0.001*q,100);


rise_time = zeros(length(diff_val_hold),1);
for i = 1:length(diff_val_hold)
    diff = diff_val_hold(i);
    rise_time(i) = risetime(diff,q,w,k);
end

fun = @(diff) (risetime(diff,q,w,k)-0.100);
f_low = diff_val_hold(find(rise_time<0.100,1,'last'));
f_high = diff_val_hold(find(rise_time>0.100,1,'first'));
qdiff_100 = fzero(fun, [f_low,f_high]);

fun = @(diff) num_poles(diff,q,w,k);

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
qdiff_osc_time = risetime(qdiff_osc,q,w,k);

figure('Color','w')
plot(diff_val_hold(1:end-1),rise_time(1:end-1),'LineWidth',3)
hold on
plot([qdiff_100,qdiff_100],[-max(rise_time)*0.025,max(rise_time)*1.05],'--','LineWidth',3.0)
plot(qdiff_high,rise_time(end-1)*1.01,'rs','MarkerSize',14,'LineWidth',3)
plot(qdiff_low,rise_time(1),'rs','MarkerSize',14,'LineWidth',3)
plot(qdiff_osc,qdiff_osc_time,'gs','MarkerSize',14,'LineWidth',3)
xlabel('\Delta q','FontSize',20)
ylabel('Rise Time (s)','FontSize',20)
set(gca,'fontsize',16)
xlim([qdiff_low-0.05*qdiff_high,qdiff_high+0.05*qdiff_high])
ylim([-max(rise_time)*0.025,max(rise_time)*1.05])


% normed_diff = diff_val_hold/((qdiff_high+0.05*qdiff_high)-(qdiff_low-0.05*qdiff_high));
% normed_rise = rise_time/(max(rise_time)*1.05);
% ind1 = find(normed_diff>0.55,1,'first');
% ind2 = find(normed_diff<0.45,1,'last');
% 
% arrow_len = ceil(length(normed_diff)/10)
% height = -0.05;
% 
% max(rise_time)*1.05
% 
% normed_diff(ind1)
% normed_diff(ind1+arrow_len)
% normed_rise(ind1)
% normed_rise(ind1+arrow_len)

% h=annotation('textarrow',[normed_diff(ind1),normed_diff(ind1+arrow_len)]...
%     ,[normed_rise(ind1)+height,normed_rise(ind1+arrow_len)+height]);
% 
% xdiff = (diff_val_hold(ind1+arrow_len)-diff_val_hold(ind1));
% ydiff = rise_time(ind1+arrow_len)-rise_time(ind1);
% h = text(diff_val_hold(ind1),rise_time(ind1),...
%     'NMDA','VerticalAlignment','Bottom');
% set(h, 'rotation', 0)
% 
% 
% h=annotation('textarrow',[normed_diff(ind2),normed_diff(ind2-arrow_len)]...
%     ,[normed_rise(ind2)+height,normed_rise(ind2-arrow_len)+height],...
%     'String','AMPA','VerticalAlignment','Bottom','HorizontalAlignment','Left');

% figure
% opt = stepDataOptions('InputOffset',0,'StepAmplitude',1);
% step(make_sys(qdiff_100,k,w,q),opt)
% figure
% opt = stepDataOptions('InputOffset',0,'StepAmplitude',1);
% step(make_sys(qdiff_osc,k,w,q),opt)
end
 
 
 function rise = risetime(diff,q,w,k)
    global div
    te  = 0.02/div;
    ti  = 0.01/div;
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
    global div
    te  = 0.02/div;
    ti  = 0.01/div;
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
    C = [0,0,(1-q-qdiff1)*(Jee+fpos)/te,(q+qdiff1)*(Jee+fpos)/te,0,0,-Jei/te,0];
    D = 0;

    sys1 = ss(A,B,C,D);
    poles1 = pole(sys1);
    num = sum(imag(poles1)~=0)-2;
  end

function [real_root] = comp_roots(diff,k,w,q)
    global div
    te  = 0.02/div;
    ti  = 0.01/div;
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
    C = [0,0,((1-q)-qdiff1)*(Jee)/te,(q+qdiff1)*(Jee)/te,0,0,-Jei/te,0];
    D = 0;

    sys = ss(A,B,C,D);

    real_root = max(real(pole(sys)));

end

function [sys] = make_sys(diff,k,w,q)
    global div
    te  = 0.02/div;
    ti  = 0.01/div;
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
    B = [1/te;0;0;0;0;0;0;0];
    C = [1,0,0,0,0,0,0,0];
    D = 0;

    sys = ss(A,B,C,D);

end
   

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 