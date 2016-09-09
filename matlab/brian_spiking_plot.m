%% Load spiking files from Brian
%  Loads a set of .mat files produced in Brian and averages and plots them.
%  I did two separate sets of runs for the STD version. The first set was
%  one run per network instantiation and the second was ten runs for the
%  first network in each setup.
std = false;
nostd = true;
if std
    
figure('Color','w');
num_runs = 1;
num_reps = 1;
delta_hold = {-0.095};
delta_u = {'-0.044','0.0','0.15'};
for k = 1:length(delta_u) % delta_q below is a naming error on the files
    for i = num_runs:num_runs
        for j = num_reps:num_reps
            load(['~/Documents/Brian_Python/Goldman_Spiking_Model/',...
                'PNAS_goldman_std_rate_p20_delta_q_',delta_u{k},'_',num2str(i),'rep',num2str(j),'.mat'])
%              load(['~/Documents/Brian_Python/Goldman_Spiking_Model/',...
%                 'PNAS_goldman_std_rate_p20_w75_delta_q_',num2str(delta_u),'_',num2str(i),'rep',num2str(j),'.mat'])
%             if (i == 1) && (j == 1)
%                 rate_ave = zeros(size(Pe_rate));
%             end
            rate_ave = Pe_rate;
%             figure
%             plot(Pe_time,Pe_rate)
%             xlabel('Time (s)')
%             ylabel('Rate (Hz)')
%             title(['delta u=',num2str(delta_u),'  Run ',num2str(i)])
        end
    end
    rate_ave = rate_ave/num_reps;
    filter = fspecial('gaussian',[30 1],1.0); % gaussian kernel where s= size of contour

    rate_ave = conv(rate_ave, filter,'same');
    plot(Pe_time-0.45,rate_ave,'LineWidth',2)
    hold on

end
xlabel('Time (s)','FontSize',30)
ylabel('Rate (Hz)','FontSize',30)
set(gca,'fontsize',20)
xlim([0,2.0])
legend('{\Delta}u = -0.044','{\Delta}u = 0.0','{\Delta}u = 0.15')
end

if nostd

num_runs = 1;
num_reps = 1;
figure('Color','w');
for qee = {0.26,0.30,0.50} % delta_q below is a naming error on the files
    
    for i = num_runs:num_runs
        for j = num_reps:num_reps

            load(['~/Documents/Brian_Python/Goldman_Spiking_Model/',...
                'PNAS_goldman_nostd2_rate_p20_qee_',num2str(qee{1}),'_',num2str(i),'rep',num2str(j),'.mat'])
            if (i == 1) && (j == 1)
                rate_ave = zeros(size(Pe_rate));
            end
            rate_ave = rate_ave + Pe_rate;
%             figure
%             plot(Pe_time,Pe_rate)
%             xlabel('Time (s)')
%             ylabel('Rate (Hz)')
%             title(['qee',num2str(qee{1}),'  Run ',num2str(i)])
        end
    end
    
    rate_ave = rate_ave/num_reps;
    filter = fspecial('gaussian',[30 1],2.0); % gaussian kernel where s= size of contour

    %rate_ave = conv(rate_ave, filter,'same');
    plot(Pe_time(1:end-5)-0.5,rate_ave(1:end-5),'LineWidth',2)
    hold on

end
xlabel('Time (s)','FontSize',30)
ylabel('Rate (Hz)','FontSize',30)
set(gca,'fontsize',20)
xlim([0,3.0])
legend('{\Delta}q = -0.04','{\Delta}q = 0.0','{\Delta}q = 0.2')
end











