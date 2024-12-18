close all
clc
%%%%%%%%%%%%%%%%%%%%%
% range of pCa values
n = 500;
pCa_start = 6.2;
pCa_end = 4.5;
pCa = linspace(pCa_start,pCa_end,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the fitted parameters for mouse LV data
mouse_LV = load('fitted_mouse_LV_2.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the fitted parameters for porcine LV data
porcine_LV = load('fitted_porcine_LV_2.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TF model (alpha,alpha_bar,beta,beta_bar,u1,u2,z1,z2,v,w)
% Compute rel ktr and P/P_0 using the fitted parameters
% in mouse LV data and porcine LV data
fitted_relktr_mouse_LV = rel_rate_force_redev(mouse_LV.relktr,pCa);
fitted_relktr_porcine_LV = rel_rate_force_redev(porcine_LV.relktr,pCa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure
ax = gca;
%%%%%%%%%
plot(pCa,fitted_relktr_mouse_LV,'k-',...
     'markersize',10,'LineWidth',4), hold on
errorbar(mouse_LV.pCa,mouse_LV.data_relktr(:,1),mouse_LV.data_relktr(:,2),...
         '.','MarkerSize',40,'MarkerEdgeColor','black',...
             'Color','black','LineWidth',2.5), hold on
%%%%%%%%%
plot(pCa,fitted_relktr_porcine_LV,'b-',...
    'markersize',10,'LineWidth',4), hold on
errorbar(porcine_LV.pCa,porcine_LV.data_relktr(:,1),porcine_LV.data_relktr(:,2),...
         '.','MarkerSize',40,'MarkerEdgeColor','blue',...
             'Color','blue','LineWidth',2.5); hold on
set(gca, 'XDir','reverse')
%%%%%%%%%
ax.FontSize = 20;
xlabel('pCa','FontName', 'Times','FontSize',24);
ax.XAxis.LineWidth = 1.5;
ylabel('Relative ktr','Interpreter','latex','FontName',... 
              'Times','FontSize',24,'FontWeight','bold');
ax.YAxis.LineWidth = 1.5;
legend('Mouse LV fitted solution','Mouse LV data',...
       'Porcine LV fitted solution','Porcine LV data',...
       'FontSize', 20,'LineWidth',1.5,'Location','best');
title('Parameter set 2: relktr vs pCa','FontSize',30)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(0,'Units','normalized')
% set(gcf, 'PaperSize', [10 8], 'PaperPosition', [0 0 10 8])
% print('Figure5B','-djpeg')