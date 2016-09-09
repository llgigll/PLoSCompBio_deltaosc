% Computes panels A, D and E on the supplementary figure.
q = 0.3;
k_hold = 2.5:-0.02:1.1;
w_hold = [5:2.5:50];

stability = stability_qdiff(q,k_hold,w_hold);

% figure;imagesc(w_hold,k_hold,-stability(1).qdiff_low); colorbar;
% set(gca,'YDir','normal')
% xlabel('w'); ylabel('k'); title('AMPA Dominated')
% 
% figure; imagesc(w_hold,k_hold,stability(1).qdiff_high); colorbar;
% set(gca,'YDir','normal')
% xlabel('w'); ylabel('k'); title('NMDA Dominated')
stability(1).rise(stability(1).rise==-0.3) = NaN;

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

stab = -stability(1).rise+stability(1).qdiff_high;
stab(stability(1).qdiff_high==max(max(stability(1).qdiff_high))) = max(max(stab));


figure('Color','w');pcolor(w_hold,k_hold,stab); colorbar;
set(gca,'YDir','normal')
xlabel('w'); ylabel('k'); title('{\Delta}q to NMDA Instability')
colormap(flipud(jet))
shading interp
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
h = get(gca,'Title'); set(h,'FontSize',16); 
set(gca,'fontsize',16)



figure('Color','w');pcolor(w_hold,k_hold,stability(1).max_rise); colorbar;
set(gca,'YDir','normal')
xlabel('w'); ylabel('k'); title('Rise Time (s)')
colormap jet
shading interp
h = get(gca,'XLabel'); set(h,'FontSize',20);
h = get(gca,'YLabel'); set(h,'FontSize',20);
h = get(gca,'Title'); set(h,'FontSize',16); 
set(gca,'fontsize',16)


% figure;imagesc(w_hold,k_hold,stability(1).min_rise+0.1); colorbar;
% set(gca,'YDir','normal')
% xlabel('w'); ylabel('k'); title('Min Rise Time')


% figure
% [w_mesh,k_mesh] = meshgrid(w_hold,k_hold);
% equilib = 1./(1-w_mesh+(k_mesh.*w_mesh.^2./(1+k_mesh.*w_mesh)));
% imagesc(w_hold,k_hold,equilib); colorbar;
% set(gca,'YDir','normal')
% xlabel('w'); ylabel('k'); title('Gain')