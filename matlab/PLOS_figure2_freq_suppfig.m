%% Control Theory approach to negative derivative feedback
%  Find the transfer function of the linear negative derivative feedback
%  model. Use the transfer function to examine the root loci and the bode
%  plots of the system. Determine the phase and gain margins and examine
%  how they change for particular parameters in the model.


% Plots figure 2 panel E and supplementary figure panel H.


te  = 0.02;
ti  = 0.01;
tee_a = 0.005;
tee_n = 0.1;
tie_a = 0.005;
tie_n = 0.1;
tei = 0.01;
tii = 0.01;


opts = bodeoptions;
opts.FreqUnits = 'Hz';
opts.MagUnits = 'abs';
opts.PhaseWrapping = 'off';
opts.PhaseMatching = 'off';

q = 0.3;
w_val = 30;
fpos_hold = 0.0;

% Top set of values are used for the supplementary figure and the bottom
% set of values are used for Figure 2.
diff_val_hold = fliplr([-0.02,-0.0175,-0.015,-0.0125,-0.01,-0.0075,-0.005,-0.0025,0.0,0.005,0.01,0.025,0.05,0.1,0.14,0.15]);
% diff_val_hold = fliplr([-0.02,-0.01,0.0,0.05,0.14]);

set(0,'DefaultAxesColorOrder','factory')
legender = cell(length(diff_val_hold),1);
sys_plot = cell(2*length(diff_val_hold),1);
sys_hold = cell(length(diff_val_hold),1);
sys2_hold = cell(length(diff_val_hold),1);
all_poles = cell(length(diff_val_hold),1);
pole_holder = zeros(length(diff_val_hold),1);
imag_poles1 = zeros(length(diff_val_hold),1);
imag_poles2 = zeros(length(diff_val_hold),1);
for i = 1:length(diff_val_hold)
    diff_val = diff_val_hold(i);
    w = realp('w',w_val);
    qdiff1 = realp('qdiff1',diff_val);
    qdiff2 = realp('qdiff2',0.0);
    fpos = realp('fpos',fpos_hold);

    k = 1.2;
    Jee = w;
    Jei = k*w;
    Jie = w;
    Jii = k*w;

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
    C2 = [1,0,0,0,0,0,0,0];
    
    D = 0;

    sys1 = ss(A,B,C,D);
    sys2 = ss(A,B,C2,D);
    
    sys_hold{i} = sys1;
    sys2_hold{i} = sys2;
    
    sys_plot{2*i-1} = sys1;
    poles1 = pole(sys1);
    poles2 = pole(sys2);
    all_poles{i} = poles1;
    [m,ind] = max(real(poles1));
    if m > 0
        sys_plot{2*i} = '--';
    else
        sys_plot{2*i} = '-';
    end
    pole_holder(i) = poles1(ind);
    
    imag_poles1(i) = sum(imag(poles1)~=0);
    imag_poles2(i) = sum(imag(poles2)~=0);

    legender{i} = num2str(diff_val);

end




wout = logspace(-2.0,2.3010,250);
wout_rad = wout*2*pi;
rise_time = zeros(length(diff_val_hold),1);
stable = zeros(length(diff_val_hold),1);
all_mag = zeros(length(wout),length(diff_val_hold));
all_phase = zeros(length(wout),length(diff_val_hold));
for i = 1:length(diff_val_hold)
    [mag,phase] = bode(sys_hold{i},wout_rad);  
    all_mag(:,i) = squeeze(mag);
    all_phase(:,i) = squeeze(phase);
    S = stepinfo(sys_hold{i});

    rise_time(i) = S.RiseTime;
    stable(i) = max(real(pole(sys_hold{i})))<0;
end



colors = jet(length(diff_val_hold));
figure('Color','w')
for i = 1:length(diff_val_hold)
    poles_imag = imag(all_poles{i}(1:5))/(2*pi);
    poles_real = real(all_poles{i}(1:5));
    plot(poles_real,poles_imag,'o','color',colors(i,:),'linewidth',2.0,'markersize',10)
    hold on
end
plot([0,0],[-80,80],'k--','linewidth',2.0)

xlim([-80,10])
ylim([-80,80])
xlabel('Real Axis','FontSize',30)
ylabel('Imaginary Axis','FontSize',30)
% legend('0.15','0.14','0.1','0.05','0.025','0.01','0.005','0.0','-0.0025','-0.005','-0.0075','-0.01','-0.0125','-0.015','-0.0175','-0.02')


figure('Color','w')
for i = 1:length(diff_val_hold)
    poles_imag = imag(all_poles{i}(3:5))/(2*pi);
    poles_real = real(all_poles{i}(3:5));
    plot(poles_real,poles_imag,'o','color',colors(i,:),'linewidth',2.0,'markersize',10)
    hold on
end
plot([0,0],[-4,4],'k--','linewidth',2.0)
xlim([-12,1])
ylim([-4,4])
xlabel('Real Axis','FontSize',30)
ylabel('Imaginary Axis','FontSize',30)
% legend('0.15','0.14','0.1','0.05','0.025','0.01','0.005','0.0','-0.0025','-0.005','-0.0075','-0.01','-0.0125','-0.015','-0.0175','-0.02')
legend('{\Delta}q = -0.02','{\Delta}q = -0.01','{\Delta}q = 0.0','{\Delta}q = 0.05','{\Delta}q = 0.14')
 

figure('Color','w')
semilogx(wout,all_mag,'LineWidth',3)
set(gca,'Fontsize',20)
xlim([0.01,300])
xlabel('Frequency (Hz)','FontSize',30)
ylabel('Magnitude (abs)','FontSize',30)
% legend('0.15','0.14','0.1','0.05','0.025','0.01','0.005','0.0','-0.0025','-0.005','-0.0075','-0.01','-0.0125','-0.015','-0.0175','-0.02')

legend('{\Delta}q = -0.02','{\Delta}q = -0.01','{\Delta}q = 0.0','{\Delta}q = 0.05','{\Delta}q = 0.14')




 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 