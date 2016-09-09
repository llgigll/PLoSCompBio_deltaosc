% Computes panel G on the supplementary figure.

q = 0.1:0.05:0.5;
tau_div = 0.2:0.1:1.0;
k = 1.2;
w = 30;

stability = stability_oscillations_taudiv_q(tau_div,q,k,w);

figure('Color','w');
pcolor(tau_div,q,stability(1).osc_low); colorbar;
set(gca,'YDir','normal')
xlabel('f_{\tau}'); ylabel('q'); title('Resonant Frequency at AMPA Instability')
colormap(flipud(jet))
shading interp
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
h = get(gca,'Title'); set(h,'FontSize',16); 
set(gca,'fontsize',16)









