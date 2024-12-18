close all
clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% range of pCa values
n = 500;
pCa_start = 6.3;
pCa_end = 4.5;
pCa = linspace(pCa_start,pCa_end,n);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the fitted parameters for mouse LV data
mouse_LV = load('fitted_mouse_LV_8.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extract the fitted parameters for porcine LV data
porcine_LV = load('fitted_porcine_LV_8.mat');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute rel ktr and P/P_0 using the TF model with the fitted parameters
% in mouse LV data and porcine LV data with different u2
param_ktr_mouse_u2 = mouse_LV.ktr;
param_ktr_porcine_u2 = porcine_LV.ktr;
u2_mouse = [mouse_LV.ktr(16); 4; 7; 10];
u2_porcine = [porcine_LV.ktr(16); 16; 13; 10];
m1 = length(u2_mouse);
m2 = length(u2_porcine);
fitted_ktr_mouse_u2 = zeros(n,m1);
fitted_ktr_porcine_u2 = zeros(n,m2);
%%%%%%%%%%%%%
for i = 1:m1
    param_ktr_mouse_u2(16) = u2_mouse(i);
    fitted_ktr_mouse_u2(:,i) = rate_force_redev_8(param_ktr_mouse_u2,pCa);
end
%%%%%%%%%%%%
for i = 1:m2
    param_ktr_porcine_u2(16) = u2_porcine(i);
    fitted_ktr_porcine_u2(:,i) = rate_force_redev_8(param_ktr_porcine_u2,pCa);    
end
%%%%%%%%%%%%
fitted_relf_mouse = rel_force(mouse_LV.relf,pCa);
fitted_relf_porcine = rel_force(porcine_LV.relf,pCa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute rel ktr and P/P_0 using the TF model with the fitted parameters
% in mouse LV data and porcine LV data with different z2
param_ktr_mouse_z2 = mouse_LV.ktr;
param_ktr_porcine_z2 = porcine_LV.ktr;
z2_mouse = [mouse_LV.ktr(18); 1.5; 2];
z2_porcine = [porcine_LV.ktr(18); 1.5; 2];
m3 = length(z2_mouse);
m4 = length(z2_porcine);
fitted_ktr_mouse_z2 = zeros(n,m3);
fitted_ktr_porcine_z2 = zeros(n,m4);
%%%%%%%%%%%%%
for i = 1:m3
    param_ktr_mouse_z2(18) = z2_mouse(i);
    fitted_ktr_mouse_z2(:,i) = rate_force_redev_8(param_ktr_mouse_z2,pCa);
end
%%%%%%%%%%%%
for i = 1:m4
    param_ktr_porcine_z2(18) = z2_porcine(i);
    fitted_ktr_porcine_z2(:,i) = rate_force_redev_8(param_ktr_porcine_z2,pCa);    
end
%%%%%%%%%%%%
% fitted_relf_mouse_z2 = rel_force(mouse_LV.relf,pCa);
% fitted_relf_porcine_z2 = rel_force(porcine_LV.relf,pCa);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figures
colors = [1,0,0
0,1,0
0,0,1
0,1,1
1,0,1
1,1,0
0,0,0];
% colororder(colors);
set(gcf,'DefaultAxesColorOrder',colors)
%%%%%%%%%%%%%%%%%%%%
clf
h1 = subplot(2,2,1);
ax = gca;
L1 = cell(m1+1,1);
for i = 1:m1
    plot(fitted_relf_mouse,fitted_ktr_mouse_u2(:,i),...
                           'markersize',10,'LineWidth',3), hold on
    if i == 1
        L1{i} = strcat('fitted solution - u_2 = ', num2str(u2_mouse(i)));
    else
        L1{i} = strcat('solution - u_2 = ', num2str(u2_mouse(i)));
    end
end
errorbar(mouse_LV.data_relf(:,1),mouse_LV.data_ktr(:,1),...
         mouse_LV.data_ktr(:,2),mouse_LV.data_ktr(:,2),...
         mouse_LV.data_relf(:,2),mouse_LV.data_relf(:,2),...
         '.','MarkerSize',30,'MarkerEdgeColor','black',...
             'Color','black','LineWidth',1.5), hold on
L1{m1+1} = strcat('Mouse LV data');
ax.FontSize = 14;
xlabel('Relative force','FontName', 'Times','FontSize',18);
ax.XAxis.LineWidth = 1.5;
ylabel('Rate of force redevelopment ($s^{-1}$)','Interpreter','latex','FontName', 'Times',...
                                       'FontSize',18,'FontWeight','bold');
ax.YAxis.LineWidth = 1.5;
legend(L1,'FontSize',14,'Location','northwest')
ax.TitleHorizontalAlignment = 'left';
title('A','FontSize',20)
% set(gca, 'XDir','reverse')
%%%%%%%%%%%%%%%%%%%%
h2 = subplot(2,2,2);
ax = gca;
L2 = cell(m2+1,1);
for i = 1:m2
    plot(fitted_relf_porcine,fitted_ktr_porcine_u2(:,i),...
                               'markersize',10,'LineWidth',3), hold on
    if i == 1
        L2{i} = strcat('fitted solution - u_2 = ', num2str(u2_porcine(i)));
    else
        L2{i} = strcat('solution - u_2 = ', num2str(u2_porcine(i)));
    end
end
errorbar(porcine_LV.data_relf(:,1),porcine_LV.data_ktr(:,1),...
         porcine_LV.data_ktr(:,2),porcine_LV.data_ktr(:,2),...
         porcine_LV.data_relf(:,2),porcine_LV.data_relf(:,2),...
         '.','MarkerSize',30,'MarkerEdgeColor','black',...
             'Color','black','LineWidth',1.5); hold on
L2{m2+1} = strcat('Porcine LV data');
ax.FontSize = 14;
xlabel('Relative force','FontName', 'Times','FontSize',18);
ax.XAxis.LineWidth = 1.5;
ylabel('Rate of force redevelopment ($s^{-1}$)','Interpreter','latex','FontName', 'Times',...
                                   'FontSize',18,'FontWeight','bold');
ax.YAxis.LineWidth = 1.5;
legend(L2,'FontSize',14,'Location','north')
ax.TitleHorizontalAlignment = 'left';
title('B','FontSize',20)
axis([0 1 0 4])
% set(gca, 'XDir','reverse')
%%%%%%%%%%%%%%%%%%%%
h3 = subplot(2,2,3);
ax = gca;
L3 = cell(m3+1,1);
for i = 1:m3
    plot(fitted_relf_mouse,fitted_ktr_mouse_z2(:,i),...
                           'markersize',10,'LineWidth',3), hold on
    if i == 1
        L3{i} = strcat('fitted solution - z_2 = ', num2str(z2_mouse(i)));
    else
        L3{i} = strcat('solution - z_2 = ', num2str(z2_mouse(i)));
    end
end
errorbar(mouse_LV.data_relf(:,1),mouse_LV.data_ktr(:,1),...
         mouse_LV.data_ktr(:,2),mouse_LV.data_ktr(:,2),...
         mouse_LV.data_relf(:,2),mouse_LV.data_relf(:,2),...
         '.','MarkerSize',30,'MarkerEdgeColor','black',...
             'Color','black','LineWidth',1.5), hold on
L3{m3+1} = strcat('Mouse LV data');
ax.FontSize = 14;
xlabel('Relative force','FontName', 'Times','FontSize',18);
ax.XAxis.LineWidth = 1.5;
ylabel('Rate of force redevelopment ($s^{-1}$)','Interpreter','latex','FontName', 'Times',...
                                       'FontSize',18,'FontWeight','bold');
ax.YAxis.LineWidth = 1.5;
legend(L3,'FontSize',14,'Location','northwest')
ax.TitleHorizontalAlignment = 'left';
title('C','FontSize',20)
% set(gca, 'XDir','reverse')
%%%%%%%%%%%%%%%%%%%%
h4 = subplot(2,2,4);
ax = gca;
L4 = cell(m4+1,1);
for i = 1:m4
    plot(fitted_relf_porcine,fitted_ktr_porcine_z2(:,i),...
                               'markersize',10,'LineWidth',3), hold on
    if i == 1
        L4{i} = strcat('fitted solution - z_2 = ', num2str(z2_porcine(i)));
    else
        L4{i} = strcat('solution - z_2 = ', num2str(z2_porcine(i)));
    end
end
errorbar(porcine_LV.data_relf(:,1),porcine_LV.data_ktr(:,1),...
         porcine_LV.data_ktr(:,2),porcine_LV.data_ktr(:,2),...
         porcine_LV.data_relf(:,2),porcine_LV.data_relf(:,2),...
         '.','MarkerSize',30,'MarkerEdgeColor','black',...
             'Color','black','LineWidth',1.5); hold on
L4{m4+1} = strcat('Porcine LV data');
ax.FontSize = 14;
xlabel('Relative force','FontName', 'Times','FontSize',18);
ax.XAxis.LineWidth = 1.5;
ylabel('Rate of force redevelopment ($s^{-1}$)','Interpreter','latex','FontName', 'Times',...
                                   'FontSize',18,'FontWeight','bold');
ax.YAxis.LineWidth = 1.5;
legend(L4,'FontSize',14,'Location','north')
ax.TitleHorizontalAlignment = 'left';
title('D','FontSize',20)
axis([0 1 0 4])
% set(gca, 'XDir','reverse')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set(0,'Units','normalized')
% set(h4,'position',[.57 .08 .38 .38])
% set(h3,'position',[.09 .08 .38 .38])
% set(h2,'position',[.57 .56 .38 .38])
% set(h1,'position',[.09 .56 .38 .38])
% set(gcf, 'PaperSize', [12 12], 'PaperPosition', [0 0 12 12])
% print('Figure7_Revised','-djpeg')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%