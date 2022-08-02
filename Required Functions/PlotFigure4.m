function PlotFigure4(ratios,Matrix_Sizes)

x = Matrix_Sizes./4; y = ratios;

% Read in Data
Files=dir('Figure4Data'); 
Files = Files(~[Files.isdir]);
for k=1:length(Files)
   load(['Figure4Data/Figure4Data_MatSiz',num2str(Matrix_Sizes(k))],'mean_Kappa_error','std_Kappa_error')
   Data(:,:,k) = squeeze(mean_Kappa_error);
   StdData(:,:,k) = squeeze(std_Kappa_error);
end

Fontsize = 12; Biasat = 120;
fig = figure('Color','w','Units','normalized','Position',[0.3,0.2,0.8292,0.573148148148148]);
tiledlayout(1,2,'Padding','none'); ax1=nexttile;
imagesc(x,y,squeeze(Data(Biasat,:,:)-Biasat),[-20 20]); xlabel('Phase Encoding Lines in TurboFLASH Readout','Fontsize',Fontsize); ylabel('\alpha:\beta Ratio','Fontsize',Fontsize); xticks([5:5:60]); xticklabels([5:5:60]); yticks([5:5:40]); yticklabels([5:5:40])
set(gca, 'YDir','normal'); colormap(ax1,bluewhitered); cb = colorbar; cb.Label.String = ['Bias at nominal \alpha = ',num2str(Biasat),char(176)]; cb.FontSize = Fontsize;
hold on; plot(9,10,'kx','MarkerSize',20,'Linewidth',3); plot(36,20,'kx','MarkerSize',20,'Linewidth',3);
text(9,11,'2D','HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',Fontsize);
text(36,21,'3D','HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',Fontsize);
ax2=nexttile;
imagesc(x,y,squeeze(StdData(Biasat,:,:)),[0 10]); xlabel('Phase Encoding Lines in TurboFLASH Readout','Fontsize',Fontsize); ylabel('\alpha:\beta Ratio','Fontsize',Fontsize); xticks([5:5:60]); xticklabels([5:5:60]); yticks([5:5:40]); yticklabels([5:5:40])
set(gca, 'YDir','normal'); colormap(ax2,bluewhitered); cb = colorbar; cb.Label.String = ['Standard Deviation of Bias at nominal \alpha = ',num2str(Biasat),char(176)]; cb.FontSize = Fontsize;
hold on; plot(9,10,'kx','MarkerSize',20,'Linewidth',3); plot(36,20,'kx','MarkerSize',20,'Linewidth',3);
text(9,11,'2D','HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',Fontsize);
text(36,21,'3D','HorizontalAlignment','center','VerticalAlignment','bottom','Fontsize',Fontsize);
saveas(fig,'CombinedTrainvsRatiowithStd.png');