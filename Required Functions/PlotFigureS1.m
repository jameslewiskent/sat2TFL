function PlotFigureS1(T1s)
% Plot Figure S1

load('Figure3Data/Figure3DataTopRow','PP_FAs','IT_FAs','Kappa_error_results');

Noise_n = 2;
% Read in satTFL
satTFLData = squeeze(Kappa_error_results{1,1});
Std_satTFLData = squeeze(Kappa_error_results{2,1});
% Read in Short TR
shortTRData = squeeze(Kappa_error_results{1,2});
Std_shortTRData = squeeze(Kappa_error_results{2,2});
% Read in Sandwich
sat2TFLData = squeeze(Kappa_error_results{1,3});
Std_sat2TFLData = squeeze(Kappa_error_results{2,3});

Values = meshgrid(PP_FAs,IT_FAs)';
cbdiffstring = ['Standard Deviation of Bias in \alpha (',char(176),')']; %'SD of \alpha Bias';

Font = 'Calibri';
tick_size = 12;
caxis2_minlim = 0; 
caxis2_maxlim = 10;
x_label_fontsize = 16;
y_label_fontsize = 16; 
title_fontsize = 16;
tick_size = 12;
cb_label_size = 12;
x_label_text = ['\alpha_{nom} (',char(176),')']; 
y_label_text = ['\beta_{nom} (',char(176),')'];
fig = figure('Units','normalized','Position',[0.1609375,0.061111111111111,0.771875,0.837037037037037]);set(fig,'color','w');
tiledlayout(2,3,'TileSpacing', 'Compact','Padding', 'Compact')

nexttile()
imagesc(PP_FAs,IT_FAs,Std_satTFLData(:,:,Noise_n))
set(gca,'XTick',[0,20,40,60,80,100,120,140,160,180],'YTick',[0,5,10,15,20],'Fontsize',tick_size,'FontName','Calibri','YDir','normal')
ylabel(y_label_text,'Fontsize',y_label_fontsize,'FontName',Font)
xlabel(x_label_text,'Fontsize',x_label_fontsize,'FontName',Font)
title('a) satTFL','Fontsize',title_fontsize,'FontName',Font);
caxis([caxis2_minlim caxis2_maxlim]);
xlim([0 180])
ylim([0 20])
colormap(gca,bluewhitered)
% Plot isoratios
line([200,0], [40,0], 'Color', 'k','Linestyle','-');
line([200,0], [20,0], 'Color', 'k','Linestyle','--');
line([400,0], [20,0], 'Color', 'k','Linestyle','-.');
line([600,0], [20,0], 'Color', 'k','Linestyle',':');
leg = legend('5:1','10:1','20:1','30:1','Location','NorthEast');
title(leg,'\alpha:\beta Ratio','FontName',Font)

nexttile()
imagesc(PP_FAs,IT_FAs,Std_shortTRData(:,:,Noise_n))
set(gca,'XTick',[0,20,40,60,80,100,120,140,160,180],'YTick',[0,5,10,15,20],'Fontsize',tick_size,'FontName','Calibri','YDir','normal')
ylabel(y_label_text,'Fontsize',y_label_fontsize,'FontName',Font)
xlabel(x_label_text,'Fontsize',x_label_fontsize,'FontName',Font)
title('b) Short TR satTFL','Fontsize',title_fontsize,'FontName',Font);
caxis([caxis2_minlim caxis2_maxlim]);
colormap(gca,bluewhitered)
xlim([0 180])
ylim([0 20])
% Plot isoratios
line([200,0], [40,0], 'Color', 'k','Linestyle','-');
line([200,0], [20,0], 'Color', 'k','Linestyle','--');
line([400,0], [20,0], 'Color', 'k','Linestyle','-.');
line([600,0], [20,0], 'Color', 'k','Linestyle',':');


nexttile()
imagesc(PP_FAs,IT_FAs,Std_sat2TFLData(:,:,Noise_n))
set(gca,'XTick',[0,20,40,60,80,100,120,140,160,180],'YTick',[0,5,10,15,20],'Fontsize',tick_size,'FontName','Calibri','YDir','normal')
ylabel(y_label_text,'Fontsize',y_label_fontsize,'FontName',Font) % ,'Position',[-6.677741106166408,10.000009536743166,1]
xlabel(x_label_text,'Fontsize',x_label_fontsize,'FontName',Font)
title('c) Sandwich','Fontsize',title_fontsize,'FontName',Font);
xlim([0 180])
ylim([0 20])
caxis([caxis2_minlim caxis2_maxlim]);
colormap(gca,bluewhitered)
% Plot isoratios
line([200,0], [40,0], 'Color', 'k','Linestyle','-');
line([200,0], [20,0], 'Color', 'k','Linestyle','--');
line([400,0], [20,0], 'Color', 'k','Linestyle','-.');
line([600,0], [20,0], 'Color', 'k','Linestyle',':');
    hcb2 = colorbar('Location','EastOutside','AxisLocation','out');
hcb2.Layout.Tile = 'east'; 
    set(hcb2.XLabel,{'String','Rotation','Fontsize'},{cbdiffstring,90,cb_label_size})
colormap(gca,bluewhitered)

% Plot isoratios
line([200,0], [40,0], 'Color', 'k','Linestyle','-');
line([200,0], [20,0], 'Color', 'k','Linestyle','--');
line([400,0], [20,0], 'Color', 'k','Linestyle','-.');
line([600,0], [20,0], 'Color', 'k','Linestyle',':');

% 

load('Figure3Data/Figure3DataBottomRow','Kappa_error_results');

% Read in satTFL
satTFLData = squeeze(Kappa_error_results{1,1});
Std_satTFLData = squeeze(Kappa_error_results{2,1});
% Read in Short TR
shortTRData = squeeze(Kappa_error_results{1,2});
Std_shortTRData = squeeze(Kappa_error_results{2,2});
% Read in Sandwich
sat2TFLData = squeeze(Kappa_error_results{1,3});
Std_sat2TFLData = squeeze(Kappa_error_results{2,3});

Yticks = 0.5:0.5:3;
y_label_text = 'T_1 (s)'; %['\beta_{nom},  degrees'];
nexttile()
imagesc(PP_FAs,T1s,Std_satTFLData(:,:,Noise_n)')
set(gca,'XTick',[0,20,40,60,80,100,120,140,160,180],'YTick',Yticks,'Fontsize',tick_size,'FontName','Calibri','YDir','normal')
ylabel(y_label_text,'Fontsize',y_label_fontsize,'FontName',Font) %,'Position',[-8.17275771746209,9.937704240793003,1]
xlabel(x_label_text,'Fontsize',x_label_fontsize,'FontName',Font)
caxis([caxis2_minlim caxis2_maxlim]);
xlim([0 180])
ylim([min(T1s) max(T1s)])
colormap(gca,bluewhitered)


nexttile()
imagesc(PP_FAs,T1s,Std_shortTRData(:,:,Noise_n)')
set(gca,'XTick',[0,20,40,60,80,100,120,140,160,180],'YTick',Yticks,'Fontsize',tick_size,'FontName','Calibri','YDir','normal')
ylabel(y_label_text,'Fontsize',y_label_fontsize,'FontName',Font)
xlabel(x_label_text,'Fontsize',x_label_fontsize,'FontName',Font)
caxis([caxis2_minlim caxis2_maxlim]);
colormap(gca,bluewhitered)
xlim([0 180])
ylim([min(T1s) max(T1s)])


nexttile()
imagesc(PP_FAs,T1s,Std_sat2TFLData(:,:,Noise_n)')
set(gca,'XTick',[0,20,40,60,80,100,120,140,160,180],'YTick',Yticks,'Fontsize',tick_size,'FontName','Calibri','YDir','normal')
ylabel(y_label_text,'Fontsize',y_label_fontsize,'FontName',Font) 
xlabel(x_label_text,'Fontsize',x_label_fontsize,'FontName',Font)
xlim([0 180])
ylim([min(T1s) max(T1s)])
caxis([caxis2_minlim caxis2_maxlim]);
colormap(gca,bluewhitered)


saveas(fig,'Figure3Data\FigureS1.png');
end
