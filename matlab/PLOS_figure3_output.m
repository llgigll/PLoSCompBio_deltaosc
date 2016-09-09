% Computes all of the output for Figure 3 uses stability_oscillations_tn

q = [0.0:0.001:0.015,0.025:0.025:0.5];
% q = [0.0,0.01,0.1,0.2];
t_n = [0.075,0.1,0.15,0.4];
k = 1.2;
w = 30;

stability = stability_oscillations_tn(t_n,q,k,w);

% figure('Color','w');
% pcolor(q,t_n,stability(1).osc_low); colorbar;
% set(gca,'YDir','normal')
% xlabel('q'); ylabel('\tau_n'); title('Resonant Frequency at AMPA Instability')
% colormap(flipud(jet))
% shading interp
% h = get(gca,'XLabel'); set(h,'FontSize',20);
% h = get(gca,'YLabel'); set(h,'FontSize',20);
% h = get(gca,'Title'); set(h,'FontSize',16); 
% set(gca,'fontsize',16)

stab_ind = cell(length(t_n),1);
figure('Color','w');
for i = 1:length(t_n)
    stab_ind{i} = find(stability(1).osc_low(i,:)~=0);
%     stab_ind{i} = find(stability(1).osc_low(i,:)<100000);
    plot(q(stab_ind{i}),stability(1).osc_low(i,stab_ind{i}),'LineWidth',3);
    hold on
end
xlabel('q'); ylabel('Resonant Frequency (Hz)');
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
set(gca,'fontsize',16)

figure('Color','w');
for i = 1:length(t_n)
    stab_ind{i} = find(stability(1).qdiff_low(i,:)~=0);
    plot(q(stab_ind{i}),stability(1).qdiff_low(i,stab_ind{i}),'LineWidth',3);
    hold on
end
plot(q,-q,'k--','LineWidth',3)
ylim([min(min(stability(1).qdiff_low)),0])
% plot(q,stability(1).qdiff_low);
xlabel('q'); ylabel('\Delta q');
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
legend(['\tau^{nmda}=',num2str(t_n(1))],['\tau^{nmda}=',num2str(t_n(2))],...
    ['\tau^{nmda}=',num2str(t_n(3))],['\tau^{nmda}=',num2str(t_n(4))])
set(gca,'fontsize',16)

figure('Color','w');
plot(q,stability(1).qdiff_rise,'LineWidth',3);
xlabel('q'); ylabel('Rise time (s)');
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
set(gca,'fontsize',16)

figure('Color','w');
plot(q,stability(1).slope./100,'LineWidth',3);
xlabel('q'); ylabel('Rise time slope (s/\Deltaq)');
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
set(gca,'fontsize',16)


figure('Color','w');
plot(q,stability(1).qdiff_osc,'LineWidth',3);
xlabel('q'); ylabel('\Deltaq Oscillations');
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
set(gca,'fontsize',16)




