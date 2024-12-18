%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mouse_LV = load('fitted_mouse_LV_8.mat');
porcine_LV = load('fitted_porcine_LV_8.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
h1 = subplot(2,2,1);
ax = gca;
errorbar(porcine_LV.pCa,porcine_LV.data_ktr(:,1),porcine_LV.data_ktr(:,2),...
         '.','MarkerSize',35,'MarkerEdgeColor','blue','Color','blue','LineWidth',2); hold on
errorbar(mouse_LV.pCa,mouse_LV.data_ktr(:,1),mouse_LV.data_ktr(:,2),...
         '.','MarkerSize',35,'MarkerEdgeColor','black','Color','black','LineWidth',2);
set(gca, 'XDir','reverse')
ax.FontSize = 16;
xlabel('pCa','FontSize',20)
ax.XAxis.LineWidth = 2;
ylabel('Rate of force redevelopment (s^{-1})','FontSize',20)
ax.YAxis.LineWidth = 2;
legend('Porcine LV','Mouse LV','FontSize',18,'Location','best')
ax.TitleHorizontalAlignment = 'left';
title('A','FontSize',24)
xlim([4.5 6.2])
%%%%%%%%%%%%%%%%%%%
h2 = subplot(2,2,2);
ax = gca;
errorbar(porcine_LV.pCa,porcine_LV.data_relktr(:,1),porcine_LV.data_relktr(:,2),...
         '.','MarkerSize',35,'MarkerEdgeColor','blue','Color','blue','LineWidth',2); hold on
errorbar(mouse_LV.pCa,mouse_LV.data_relktr(:,1),mouse_LV.data_relktr(:,2),...
         '.','MarkerSize',35,'MarkerEdgeColor','black','Color','black','LineWidth',2);
set(gca, 'XDir','reverse')
ax.FontSize = 16;
xlabel('pCa','FontSize',20)
ax.XAxis.LineWidth = 2;
ylabel('Relative ktr','FontSize',20)
ax.YAxis.LineWidth = 2;
legend('Porcine LV','Mouse LV','FontSize',18,'Location','best')
ax.TitleHorizontalAlignment = 'left';
title('B','FontSize',24)
xlim([4.5 6.2])
%%%%%%%%%%%%%%%%%%%
h3 = subplot(2,2,3);
ax = gca;
errorbar(porcine_LV.data_relf(:,1),porcine_LV.data_ktr(:,1),porcine_LV.data_ktr(:,2),...
         porcine_LV.data_ktr(:,2),porcine_LV.data_relf(:,2),porcine_LV.data_relf(:,2),...
         '.','MarkerSize',35,'MarkerEdgeColor','blue','Color','blue','LineWidth',2); hold on
errorbar(mouse_LV.data_relf(:,1),mouse_LV.data_ktr(:,1),mouse_LV.data_ktr(:,2),...
         mouse_LV.data_ktr(:,2),mouse_LV.data_relf(:,2),mouse_LV.data_relf(:,2),...
         '.','MarkerSize',35,'MarkerEdgeColor','black','Color','black','LineWidth',2);
ax.FontSize = 16;
xlabel('Relative force (P/P_0)','FontSize',20)
ax.XAxis.LineWidth = 2;
ylabel('Rate of force redevelopment (s^{-1})','FontSize',20)
ax.YAxis.LineWidth = 2;
legend('Porcine LV','Mouse LV','FontSize',18,'Location','best')
ax.TitleHorizontalAlignment = 'left';
title('C','FontSize',24)
%%%%%%%%%%%%%%%%%%%%
h4 = subplot(2,2,4);
ax = gca;
errorbar(porcine_LV.data_relf(:,1),porcine_LV.data_relktr(:,1),porcine_LV.data_relktr(:,2),...
         porcine_LV.data_relktr(:,2),porcine_LV.data_relf(:,2),porcine_LV.data_relf(:,2),...
         '.','MarkerSize',35,'MarkerEdgeColor','blue','Color','blue','LineWidth',2); hold on
errorbar(mouse_LV.data_relf(:,1),mouse_LV.data_relktr(:,1),mouse_LV.data_relktr(:,2),...
         mouse_LV.data_relktr(:,2),mouse_LV.data_relf(:,2),mouse_LV.data_relf(:,2),...
         '.','MarkerSize',35,'MarkerEdgeColor','black','Color','black','LineWidth',2);
ax.FontSize = 16;
xlabel('Relative force (P/P_0)','FontSize',20)
ax.XAxis.LineWidth = 2;
ylabel('Relative ktr','FontSize',20)
ax.YAxis.LineWidth = 2;
legend('Porcine LV','Mouse LV','FontSize',18,'Location','best')
ax.TitleHorizontalAlignment = 'left';
title('D','FontSize',24)
%%%%%%%%%%%%%%%%%%%
% set(0,'Units','normalized')
% set(h4,'position',[.57 .08 .38 .38])
% set(h3,'position',[.09 .08 .38 .38])
% set(h2,'position',[.57 .56 .38 .38])
% set(h1,'position',[.09 .56 .38 .38])
% set(gcf, 'PaperSize', [12 12], 'PaperPosition', [0 0 12 12])
% print('Figure1_Revised','-djpeg')
%%%%%%%%%%%%%%%%%%%%