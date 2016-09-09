%% STD Plot

load('~/Documents/Brian_Python/Goldman_Spiking_Model/PNAS_goldman_std_tfunc_p20')
delta_u = [-0.075,-0.07,-0.06,-0.05,-0.025,0,0.15];
leg_hold = cell(length(delta_u),1);
for i = 1:length(delta_u)
    leg_hold{i} = ['{\Delta}u = ',num2str(delta_u(i))];
end

step_size = freqs(2)-freqs(1);

inds = (freqs>0.75)&(freqs<=15.0);
figure('Color','w')
for i = 1:length(delta_u)
  plot(freqs(inds==1),t_func_total(i,inds==1),'linewidth',3.0)
  hold on
end
set(gca,'fontsize',20)
xlabel('Frequency (Hz)','fontsize',30)
ylabel('Amplitude (V)','fontsize',30)
legend(leg_hold{:},'fontsize',14)


figure('Color','w')
errorbar(delta_u,t_func_int_mean,t_func_int_std,'linewidth',2.0)
set(gca,'fontsize',20)
xlabel('{\Delta}u','fontsize',30)
%ylabel('Mean Distance','fontsize',30)

%% No STD Plot


load('~/Documents/Brian_Python/Goldman_Spiking_Model/PNAS_goldman_nostd_tfunc_p20')
q = [-0.04,-0.03,-0.02,-0.01,0.0,0.20];
leg_hold = cell(length(q),1);
for i = 1:length(q)
    leg_hold{i} = ['{\Delta}q = ',num2str(q(i))];
end

step_size = freqs(2)-freqs(1);


inds = (freqs>0.75)&(freqs<=15.0);
figure('Color','w')
for i = 1:length(q)
  plot(freqs(inds==1),t_func_total(i,inds==1),'linewidth',3.0)
  hold on
end
set(gca,'fontsize',20)
xlabel('Frequency (Hz)','fontsize',30)
ylabel('Amplitude (V)','fontsize',30)
legend(leg_hold{:},'fontsize',14)



figure('Color','w')
errorbar(q,t_func_int_mean,t_func_int_std,'linewidth',2.0)
set(gca,'fontsize',20)
xlabel('{\Delta}q','fontsize',30)
ylabel('Mean Area (V{\times}Hz)','fontsize',25)

%% Multisine Input Plot

% figure('Color','w')
% for i = 1:length(delta_u)
%   plot(freqs((freqs>=0)&(freqs<=15.0)),input_fft_total(i,(freqs>=0)&(freqs<=15.0)),'linewidth',3.0)
%   hold on
% end
% set(gca,'fontsize',20)
% ylabel('Magnitude (abs)','fontsize',30)
% xlabel('Frequency (Hz)','fontsize',30)
% ylim([0,1.05*max(input_fft_total(1,:))])

% figure('Color','w')
% for i = 1:length(q)
%   plot(freqs((freqs>=0)&(freqs<=15.0)),input_fft_total(i,(freqs>=0)&(freqs<=15.0)),'linewidth',3.0)
%   hold on
% end
% set(gca,'fontsize',20)
% ylabel('Magnitude (abs)','fontsize',30)
% xlabel('Frequency (Hz)','fontsize',30)
% ylim([0,1.05*max(input_fft_total(1,:))])

% load('~/Documents/Brian_Python/Goldman_Spiking_Model/PNAS_goldman_nostd_fft_p20_qee0.225_1.mat')
% t = linspace(0,5,10000);
% figure('Color','w')
% plot(t,Pe_input(1,:)*10^3,'linewidth',3.0)
% set(gca,'fontsize',20)
% xlabel('Time (s)','fontsize',30)
% ylabel('Modulation (mV)','fontsize',30)












