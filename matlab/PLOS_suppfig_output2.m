% Computes panels B, C and F on the supplementary figure.

q = 0.30;
k_hold = 2.5:-0.02:1.1;
w_hold = 5:2.5:50;

stability = stability_oscillations(q,k_hold,w_hold);

figure('Color','w');
pcolor(w_hold,k_hold,-stability(1).qdiff_low); colorbar;
set(gca,'YDir','normal')
xlabel('w'); ylabel('k'); title('{\Delta}q to AMPA Instability')
colormap(flipud(jet))
shading interp
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
h = get(gca,'Title'); set(h,'FontSize',16); 
set(gca,'fontsize',16)

figure('Color','w');
pcolor(w_hold,k_hold,stability(1).osc_low); colorbar;
set(gca,'YDir','normal')
xlabel('w'); ylabel('k'); title('Resonant Frequency at AMPA Instability')
colormap(flipud(jet))
shading interp
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
h = get(gca,'Title'); set(h,'FontSize',16); 
set(gca,'fontsize',16)

figure('Color','w');
pcolor(w_hold,k_hold,-stability(1).qdiff_osc); colorbar;
set(gca,'YDir','normal')
xlabel('w'); ylabel('k'); title('{\Delta}q to Delta Oscillations')
colormap(flipud(jet))
shading interp
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
h = get(gca,'Title'); set(h,'FontSize',16); 
set(gca,'fontsize',16)

figure('Color','w');
pcolor(w_hold,k_hold,stability(1).qdiff_rise); colorbar;
set(gca,'YDir','normal')
xlabel('w'); ylabel('k'); title('Rise Time (s) at Bifurcation')
colormap(flipud(jet))
shading interp
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
h = get(gca,'Title'); set(h,'FontSize',16); 
set(gca,'fontsize',16)









