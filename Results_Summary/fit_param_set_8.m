%% Options of execution - paramter set 8
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fit_ktr_porcine_LV_data = 0; 
Fit_ktr_mouse_LV_data = 0;
%
Fit_rel_SSF_porcine_LV_data = 0; 
Fit_rel_SSF_mouse_LV_data = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit model to ktr vs pCa data in porcine LV

if Fit_ktr_porcine_LV_data
    
    close all
    clc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mechanical data from porcine LV
    pCa_porcine_LV = [6.1; 6.0; 5.9; 5.8; 5.7; 5.6; 5.5; 5.4; 4.5];
    ktr_mean_porcine_LV = [2.43; 1.79; 1.47; 1.09; 1.21; 1.67; 2.22; 2.51; 3.32];
    ktr_SEM_porcine_LV = [0.15; 0.15; 0.14; 0.11; 0.09; 0.09; 0.11; 0.11; 0.17];
    %%%%%%%%%%%%
    % Error function - the norm of the difference between actual observations
    % (data) and predicted observations (the solution of the model)
    ftns_1 = @(para_fit) norm(ktr_mean_porcine_LV - rate_force_redev_8(para_fit,pCa_porcine_LV));
    Parms_1 = 20; % the number of fitting parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial values of the parameters to be fitted
    pCa_50 = 5.65;
    k0_BC = 24.35112;
    kCa_BC = 243.53604;
    k0_CB = 1592.75304;
    kCa_CB = 117.51833;
    f0_CM1 = 90.99207;
    f0_M1C = 691.79199;
    k_M1M2 = 2.67635;
    k_M2M1 = 2.55643;
    k_M2C = 0.69478;
    alpha = 1;
    alpha_bar = 0.83950;
    beta = 1;
    beta_bar = 1;
    u1 = 1.10367;
    u2 = 17.08321;
    z1 = 1.00074;
    z2 = 1.00000;
    v = 1.03745;
    w = 1.00000;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the vector of parameters to be fitted
    para_fit_1 = zeros(1,20);
    
    para_fit_1(1) = pCa_50;
    para_fit_1(2) = k0_BC;
    para_fit_1(3) = kCa_BC;
    para_fit_1(4) = k0_CB;
    para_fit_1(5) = kCa_CB;
    para_fit_1(6) = f0_CM1;
    para_fit_1(7) = f0_M1C;
    para_fit_1(8) = k_M1M2;
    para_fit_1(9) = k_M2M1;
    para_fit_1(10) = k_M2C;
    para_fit_1(11) = alpha;
    para_fit_1(12) = alpha_bar;
    para_fit_1(13) = beta;
    para_fit_1(14) = beta_bar;
    para_fit_1(15) = u1;
    para_fit_1(16) = u2;
    para_fit_1(17) = z1;
    para_fit_1(18) = z2;
    para_fit_1(19) = v;
    para_fit_1(20) = w;
    lower_bound_1 = [5.65;zeros(13,1);ones(6,1)];
    upper_bound_1 = [5.65;Inf(9,1);ones(4,1);Inf(6,1)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % t0 = clock;
    % fprintf('\nStart Time: %4d-%02d-%02d %02d:%02d:%07.4f\n',t0)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % opts = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    % % Fit the model to the mechanical data from porcine LV
    % [fitted_para_1,resnorm_1,residual_1,exitflag_1,output_1,lambda_1,jacobian_1] = lsqcurvefit( ...
    %     @rate_force_redev,para_fit_1,pCa_porcine_LV,ktr_mean_porcine_LV,lower_bound_1,...
    %     upper_bound_1,opts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimization options for ps
    opts = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
    % pattern-search algorithm
    [fitted_para_1,fval_1,exitflag_1,output_1]=patternsearch(ftns_1,para_fit_1,...
                                    [],[],[],[],lower_bound_1,upper_bound_1,[],opts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print out values of the new fitted parameters
    fprintf('fitted parameters for porcine LV data:\n')
    fprintf('pCa50 = %f \n',fitted_para_1(1))
    fprintf('k0_BC = %f, kCa_BC = %f, k0_CB = %f, kCa_CB = %f \n', fitted_para_1(2:5))
    fprintf('f0_CM1 = %f, f0_M1C = %f \n', fitted_para_1(6:7))
    fprintf('k_M1M2 = %f, k_M2M1 = %f, k_M2C = %f \n', fitted_para_1(8:10))
    fprintf('alpha = %f, alpha_bar = %f, beta = %f, beta_bar = %f \n', fitted_para_1(11:14))
    fprintf('u1 = %f, u2 = %f \n', fitted_para_1(15:16))
    fprintf('z1 = %f, z2 = %f \n', fitted_para_1(17:18))
    fprintf('v = %f, w = %f \n', fitted_para_1(19:20))
    fprintf('RMSE = %f', sqrt(resnorm_1))
    % fprintf('RMSE = %f', fval_1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the fitted rate of force redevelopment
    n = 500;
    pCa_start = 6.3;
    pCa_end = 4.5;
    p_Ca = linspace(pCa_start,pCa_end,n);
    fitted_ktr_1 = rate_force_redev(fitted_para_1,p_Ca);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % t1 = clock;
    % fprintf('\nStop Time: %4d-%02d-%02d %02d:%02d:%07.4f\n', t1)
    % 
    % GA_Time = etime(t1,t0)
    % 
    % % Show the time elapsed for the ga
    % DT_GA_Time = datetime([0 0 0 0 0 GA_Time], 'Format', 'HH:mm:ss.SSSS');
    % fprintf('\nElapsed Time: %23.15E\t\t%s\n\n', GA_Time, DT_GA_Time)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the mechanical data and the fitted solution in the same graph
    clf
    h = axes('Position',[0 0 1 1],'Visible','off');
    axes('Position',[.1 .1 .62 .8]);
    ax = gca;
    errorbar(pCa_porcine_LV,ktr_mean_porcine_LV,ktr_SEM_porcine_LV,'.','MarkerSize',25,...
         'MarkerEdgeColor','blue','Color','blue','LineWidth',1); hold on
    plot(p_Ca,fitted_ktr_1,'k-','LineWidth',2);
    set(gca, 'XDir','reverse');
    ax.FontSize = 16
    xlabel('pCa','FontSize',16);
    ax.XAxis.LineWidth = 1.5;
    ylabel('Rate of force redevelopment ($ktr$)','Interpreter','latex','FontSize',16);
    ax.YAxis.LineWidth = 1.5;
    legend('Porcine LV data','fitted solution','Location','best');
    
    title('D', 'FontSize',20)
    ax.TitleHorizontalAlignment = 'left';
    str(1) = {'Fitted parameter set:'};
    str(2) = {[' $pCa_{50}$ = ', num2str(fitted_para_1(1))]};
    str(3) = {[' $k^0_{BC}$ = ', num2str(fitted_para_1(2))]};
    str(4) = {[' $k^{Ca}_{BC}$ = ', num2str(fitted_para_1(3))]};
    str(5) = {[' $k^0_{CB}$ = ', num2str(fitted_para_1(4))]};
    str(6) = {[' $k^{Ca}_{CB}$ = ', num2str(fitted_para_1(5))]};
    str(7) = {[' $f^0_{CM_1}$ = ', num2str(fitted_para_1(6))]};
    str(8) = {[' $f^0_{M_1C}$ = ', num2str(fitted_para_1(7))]};
    str(9) = {[' $k_{M_1M_2}$ = ', num2str(fitted_para_1(8))]};
    str(10) = {[' $k_{M_2M_1}$ = ', num2str(fitted_para_1(9))]};
    str(11) = {[' $k_{M_2C}$ = ', num2str(fitted_para_1(10))]};
    str(12) = {[' $\alpha$ = ', num2str(fitted_para_1(11))]};
    str(13) = {[' $\bar{\alpha}$ = ', num2str(fitted_para_1(12))]};
    str(14) = {[' $\beta$ = ', num2str(fitted_para_1(13))]};
    str(15) = {[' $\bar{\beta}$ = ', num2str(fitted_para_1(14))]};
    str(16) = {[' $u_1$ = ', num2str(fitted_para_1(15))]};
    str(17) = {[' $u_2$ = ', num2str(fitted_para_1(16))]};
    str(18) = {[' $z_1$ = ', num2str(fitted_para_1(17))]};
    str(19) = {[' $z_2$ = ', num2str(fitted_para_1(18))]};
    str(20) = {[' $v$ = ', num2str(fitted_para_1(19))]};
    str(21) = {[' $w$ = ', num2str(fitted_para_1(20))]};
    str(22) = {[' RMSE = ', num2str(sqrt(resnorm_1))]};
    % str(22) = {[' RMSE = ', num2str(fval_1)]};
    set(gcf,'CurrentAxes',h)
    text(.77,.5,str,'Interpreter','latex','FontSize',16)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set(0,'Units','normalized')
    % set(gcf, 'PaperSize', [10 8], 'PaperPosition', [0 0 10 8])
    % print('fig2D','-dpdf')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit model to ktr vs pCa data in mourse LV

if Fit_ktr_mouse_LV_data
    
    close all
    clc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mechanical data from mouse LV
    pCa_mouse_LV = [6.2; 6.1; 6.0; 5.9; 5.8; 5.7; 5.6; 5.5; 5.4; 4.5];
    ktr_mean_mouse_LV = [3.20; 3.16; 3.42; 4.35; 7.03; 11.04; 16.40; 20.70; 22.52; 36.23];
    ktr_SEM_mouse_LV = [0.28; 0.39; 0.48; 0.54; 0.79; 1.02; 1.19; 1.24; 1.20; 1.79];
    %%%%%%%%%%%%
    % Error function - the norm of the difference between actual observations
    % (data) and predicted observations (the solution of the model)
    ftns_3 = @(para_fit) norm(ktr_mean_mouse_LV - rate_force_redev_8(para_fit,pCa_mouse_LV));
    Parms_1 = 20; % the number of fitting parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial values of the parameters to be fitted
    pCa_50 = 5.64;
    k0_BC = 0.00000;
    kCa_BC = 1.9452;
    k0_CB = 2.3635;
    kCa_CB = 0.3963;
    f0_CM1 = 42.9689;
    f0_M1C = 50.2944;
    k_M1M2 = 42.7133;
    k_M2M1 = 44.7155;
    k_M2C = 47.8151;
    alpha = 0.9998;
    alpha_bar = 0.8758;
    beta = 0.5589;
    beta_bar = 0.3856;
    u1 = 1.0054;
    u2 = 2.2786;
    z1 = 1.0047;
    z2 = 1.0118;
    v = 1.0082;
    w = 1.1034;
    % RMSE = 1.04423
    % Create the vector of parameters to be fitted
    para_fit_3 = zeros(1,20);
    
    para_fit_3(1) = pCa_50;
    para_fit_3(2) = k0_BC;
    para_fit_3(3) = kCa_BC;
    para_fit_3(4) = k0_CB;
    para_fit_3(5) = kCa_CB;
    para_fit_3(6) = f0_CM1;
    para_fit_3(7) = f0_M1C;
    para_fit_3(8) = k_M1M2;
    para_fit_3(9) = k_M2M1;
    para_fit_3(10) = k_M2C;
    para_fit_3(11) = alpha;
    para_fit_3(12) = alpha_bar;
    para_fit_3(13) = beta;
    para_fit_3(14) = beta_bar;
    para_fit_3(15) = u1;
    para_fit_3(16) = u2;
    para_fit_3(17) = z1;
    para_fit_3(18) = z2;
    para_fit_3(19) = v;
    para_fit_3(20) = w;
    lower_bound_3 = [5.63;zeros(13,1);ones(6,1)];
    upper_bound_3 = [5.64;Inf(9,1);ones(4,1);Inf(6,1)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % t0 = clock;
    % fprintf('\nStart Time: %4d-%02d-%02d %02d:%02d:%07.4f\n',t0)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % opts = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    % % Fit the model to the mechanical data from porcine LV
    % [fitted_para_3,resnorm_3,residual_3,exitflag_3,output_3,lambda_3,jacobian_3] = lsqcurvefit( ...
    %     @rate_force_redev,para_fit_3,pCa_mouse_LV,ktr_mean_mouse_LV,lower_bound_3,...
    %     upper_bound_3,opts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimization options for ps
    opts = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
    % pattern-search algorithm
    [fitted_para_3,fval_3,exitflag_3,output_3]=patternsearch(ftns_3,para_fit_3,...
                                    [],[],[],[],lower_bound_3,upper_bound_3,[],opts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print out values of the new fitted parameters
    fprintf('fitted parameters for mouse LV data:\n')
    fprintf('pCa50 = %f \n',fitted_para_3(1))
    fprintf('k0_BC = %f, kCa_BC = %f, k0_CB = %f, kCa_CB = %f \n', fitted_para_3(2:5))
    fprintf('f0_CM1 = %f, f0_M1C = %f \n', fitted_para_3(6:7))
    fprintf('k_M1M2 = %f, k_M2M1 = %f, k_M2C = %f \n', fitted_para_3(8:10))
    fprintf('alpha = %f, alpha_bar = %f, beta = %f, beta_bar = %f \n', fitted_para_3(11:14))
    fprintf('u1 = %f, u2 = %f \n', fitted_para_3(15:16))
    fprintf('z1 = %f, z2 = %f \n', fitted_para_3(17:18))
    fprintf('v = %f, w = %f \n', fitted_para_3(19:20))
    % fprintf('RMSE = %f', sqrt(resnorm_3))
    fprintf('RMSE = %f', fval_3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the fitted rate of force redevelopment
    n = 500;
    pCa_start = 6.3;
    pCa_end = 4.5;
    p_Ca = linspace(pCa_start,pCa_end,n);
    fitted_ktr_3 = rate_force_redev(fitted_para_3,p_Ca);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the mechanical data and the fitted solution in the same graph
    clf
    h = axes('Position',[0 0 1 1],'Visible','off');
    axes('Position',[.1 .1 .62 .8]);
    ax = gca;
    errorbar(pCa_mouse_LV,ktr_mean_mouse_LV,ktr_SEM_mouse_LV,'.','MarkerSize',25,...
         'MarkerEdgeColor','blue','Color','blue','LineWidth',1); hold on
    plot(p_Ca,fitted_ktr_3,'k-','LineWidth',2);
    set(gca, 'XDir','reverse');
    ax.FontSize = 16
    xlabel('pCa','FontSize',16);
    ax.XAxis.LineWidth = 1.5;
    ylabel('Rate of force redevelopment ($ktr$)','Interpreter','latex','FontSize',16);
    ax.YAxis.LineWidth = 1.5;
    legend('Mouse LV data','fitted solution','Location','best');
    
    title('C', 'FontSize',20)
    ax.TitleHorizontalAlignment = 'left';
    str(1) = {'Fitted parameter set:'};
    str(2) = {[' $pCa_{50}$ = ', num2str(fitted_para_3(1))]};
    str(3) = {[' $k^0_{BC}$ = ', num2str(fitted_para_3(2))]};
    str(4) = {[' $k^{Ca}_{BC}$ = ', num2str(fitted_para_3(3))]};
    str(5) = {[' $k^0_{CB}$ = ', num2str(fitted_para_3(4))]};
    str(6) = {[' $k^{Ca}_{CB}$ = ', num2str(fitted_para_3(5))]};
    str(7) = {[' $f^0_{CM_1}$ = ', num2str(fitted_para_3(6))]};
    str(8) = {[' $f^0_{M_1C}$ = ', num2str(fitted_para_3(7))]};
    str(9) = {[' $k_{M_1M_2}$ = ', num2str(fitted_para_3(8))]};
    str(10) = {[' $k_{M_2M_1}$ = ', num2str(fitted_para_3(9))]};
    str(11) = {[' $k_{M_2C}$ = ', num2str(fitted_para_3(10))]};
    str(12) = {[' $\alpha$ = ', num2str(fitted_para_3(11))]};
    str(13) = {[' $\bar{\alpha}$ = ', num2str(fitted_para_3(12))]};
    str(14) = {[' $\beta$ = ', num2str(fitted_para_3(13))]};
    str(15) = {[' $\bar{\beta}$ = ', num2str(fitted_para_3(14))]};
    str(16) = {[' $u_1$ = ', num2str(fitted_para_3(15))]};
    str(17) = {[' $u_2$ = ', num2str(fitted_para_3(16))]};
    str(18) = {[' $z_1$ = ', num2str(fitted_para_3(17))]};
    str(19) = {[' $z_2$ = ', num2str(fitted_para_3(18))]};
    str(20) = {[' $v$ = ', num2str(fitted_para_3(19))]};
    str(21) = {[' $w$ = ', num2str(fitted_para_3(20))]};
    % str(22) = {[' RMSE = ', num2str(sqrt(resnorm_3))]};
    str(22) = {[' RMSE = ', num2str(fval_3)]};
    set(gcf,'CurrentAxes',h)
    text(.77,.5,str,'Interpreter','latex','FontSize',14)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set(0,'Units','normalized')
    % set(gcf, 'PaperSize', [10 8], 'PaperPosition', [0 0 10 8])
    % print('fig2C','-dpdf')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit model to rel force vs pCa data in porcine LV

if Fit_rel_SSF_porcine_LV_data
    
    close all
    clc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % relative force in porcine LV data
    pCa_porcine_LV = [6.1; 6.0; 5.9; 5.8; 5.7; 5.6; 5.5; 5.4; 4.5];
    rel_SSF_mean_porcine_LV = [0.021; 0.040; 0.081; 0.205; 0.470; 0.723; 0.847; 0.877; 1];
    rel_SSF_SEM_porcine_LV = [0.003; 0.005; 0.008; 0.015; 0.027; 0.023; 0.012; 0.008; 0];
    %%%%%%%%%%%%
    % Error function - the norm of the difference between actual observations
    % (data) and predicted observations (the solution of the model)
    ftns_4 = @(para_fit) norm(rel_SSF_mean_porcine_LV - rel_force(para_fit,pCa_porcine_LV));
    Parms_4 = 20; % the number of fitting parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial values of the parameters to be fitted
    pCa_50 = 5.68;
    k0_BC = 0.00000;
    kCa_BC = 1.25272;
    k0_CB = 9.41757;
    kCa_CB = 0.00000;
    f0_CM1 = 0.50694;
    f0_M1C = 2.94072;
    k_M1M2 = 3.74260;
    k_M2M1 = 5.66271;
    k_M2C = 0.12256;
    alpha = 1.00000;
    alpha_bar = 0.74030;
    beta = 0.91648;
    beta_bar = 0.94096;
    u1 = 1.04719;
    u2 = 14.43902;
    z1 = 1;
    z2 = 2.22580;
    v = 1.04729;
    w = 1;
    % RMSE = 0.03094
    % Create the vector of parameters to be fitted
    para_fit_4 = zeros(1,20);
    
    para_fit_4(1) = pCa_50;
    para_fit_4(2) = k0_BC;
    para_fit_4(3) = kCa_BC;
    para_fit_4(4) = k0_CB;
    para_fit_4(5) = kCa_CB;
    para_fit_4(6) = f0_CM1;
    para_fit_4(7) = f0_M1C;
    para_fit_4(8) = k_M1M2;
    para_fit_4(9) = k_M2M1;
    para_fit_4(10) = k_M2C;
    para_fit_4(11) = alpha;
    para_fit_4(12) = alpha_bar;
    para_fit_4(13) = beta;
    para_fit_4(14) = beta_bar;
    para_fit_4(15) = u1;
    para_fit_4(16) = u2;
    para_fit_4(17) = z1;
    para_fit_4(18) = z2;
    para_fit_4(19) = v;
    para_fit_4(20) = w;
    lower_bound_4 = [5.68;zeros(13,1);ones(6,1)];
    upper_bound_4 = [5.69;Inf(9,1);ones(4,1);Inf(6,1)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % t0 = clock;
    % fprintf('\nStart Time: %4d-%02d-%02d %02d:%02d:%07.4f\n',t0)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % opts = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    % % Fit the model to the mechanical data from porcine LV
    % [fitted_para_4,resnorm_4,residual_4,exitflag_4,output_4,lambda_4,jacobian_4] = lsqcurvefit( ...
    %     @rel_force,para_fit_4,pCa_porcine_LV,rel_SSF_mean_porcine_LV,lower_bound_4,...
    %     upper_bound_4,opts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimization options for ps
    opts = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
    % pattern-search algorithm
    [fitted_para_4,fval_4,exitflag_4,output_4]=patternsearch(ftns_4,para_fit_4,...
                                    [],[],[],[],lower_bound_4,upper_bound_4,[],opts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % t1 = clock;
    % fprintf('\nStop Time: %4d-%02d-%02d %02d:%02d:%07.4f\n', t1)
    % 
    % GA_Time = etime(t1,t0)
    % 
    % % Show the time elapsed for the ga
    % DT_GA_Time = datetime([0 0 0 0 0 GA_Time], 'Format', 'HH:mm:ss.SSSS');
    % fprintf('\nElapsed Time: %23.15E\t\t%s\n\n', GA_Time, DT_GA_Time)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print out values of the new fitted parameters
    fprintf('fitted parameters for porcine LV data:\n')
    fprintf('pCa50 = %f \n',fitted_para_4(1))
    fprintf('k0_BC = %f, kCa_BC = %f, k0_CB = %f, kCa_CB = %f \n', fitted_para_4(2:5))
    fprintf('f0_CM1 = %f, f0_M1C = %f \n', fitted_para_4(6:7))
    fprintf('k_M1M2 = %f, k_M2M1 = %f, k_M2C = %f \n', fitted_para_4(8:10))
    fprintf('alpha = %f, alpha_bar = %f, beta = %f, beta_bar = %f \n', fitted_para_4(11:14))
    fprintf('u1 = %f, u2 = %f \n', fitted_para_4(15:16))
    fprintf('z1 = %f, z2 = %f \n', fitted_para_4(17:18))
    fprintf('v = %f, w = %f \n', fitted_para_4(19:20))
    fprintf('RMSE = %f', sqrt(resnorm_4))
    % fprintf('RMSE = %f', fval_4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the fitted relative force
    n = 500;
    pCa_start = 6.3;
    pCa_end = 4.5;
    p_Ca = linspace(pCa_start,pCa_end,n);
    fitted_relf_4 = rel_force(fitted_para_4,p_Ca);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the mechanical data and the fitted solution in the same graph
    clf
    h = axes('Position',[0 0 1 1],'Visible','off');
    axes('Position',[.1 .1 .62 .8]);
    ax = gca;
    errorbar(pCa_porcine_LV,rel_SSF_mean_porcine_LV,rel_SSF_SEM_porcine_LV,'.','MarkerSize',25,...
         'MarkerEdgeColor','blue','Color','blue','LineWidth',1); hold on
    plot(p_Ca,fitted_relf_4,'k-','LineWidth',2);
    set(gca, 'XDir','reverse');
    ax.FontSize = 16
    xlabel('pCa','FontSize',16);
    ax.XAxis.LineWidth = 1.5;
    ylabel('Relative force','FontSize',16);
    ax.YAxis.LineWidth = 1.5;
    legend('Porcine LV data','fitted solution','Location','best');
    
    title('B', 'FontSize',20)
    ax.TitleHorizontalAlignment = 'left';
    str(1) = {'Fitted parameter set:'};
    str(2) = {[' $pCa_{50}$ = ', num2str(fitted_para_4(1))]};
    str(3) = {[' $k^0_{BC}$ = ', num2str(fitted_para_4(2))]};
    str(4) = {[' $k^{Ca}_{BC}$ = ', num2str(fitted_para_4(3))]};
    str(5) = {[' $k^0_{CB}$ = ', num2str(fitted_para_4(4))]};
    str(6) = {[' $k^{Ca}_{CB}$ = ', num2str(fitted_para_4(5))]};
    str(7) = {[' $f^0_{CM_1}$ = ', num2str(fitted_para_4(6))]};
    str(8) = {[' $f^0_{M_1C}$ = ', num2str(fitted_para_4(7))]};
    str(9) = {[' $k_{M_1M_2}$ = ', num2str(fitted_para_4(8))]};
    str(10) = {[' $k_{M_2M_1}$ = ', num2str(fitted_para_4(9))]};
    str(11) = {[' $k_{M_2C}$ = ', num2str(fitted_para_4(10))]};
    str(12) = {[' $\alpha$ = ', num2str(fitted_para_4(11))]};
    str(13) = {[' $\bar{\alpha}$ = ', num2str(fitted_para_4(12))]};
    str(14) = {[' $\beta$ = ', num2str(fitted_para_4(13))]};
    str(15) = {[' $\bar{\beta}$ = ', num2str(fitted_para_4(14))]};
    str(16) = {[' $u_1$ = ', num2str(fitted_para_4(15))]};
    str(17) = {[' $u_2$ = ', num2str(fitted_para_4(16))]};
    str(18) = {[' $z_1$ = ', num2str(fitted_para_4(17))]};
    str(19) = {[' $z_2$ = ', num2str(fitted_para_4(18))]};
    str(20) = {[' $v$ = ', num2str(fitted_para_4(19))]};
    str(21) = {[' $w$ = ', num2str(fitted_para_4(20))]};
    str(22) = {[' RMSE = ', num2str(sqrt(resnorm_4))]};
    % str(22) = {[' RMSE = ', num2str(fval_4)]};
    set(gcf,'CurrentAxes',h)
    text(.77,.5,str,'Interpreter','latex','FontSize',16)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set(0,'Units','normalized')
    % set(gcf, 'PaperSize', [10 8], 'PaperPosition', [0 0 10 8])
    % print('fig2B','-dpdf')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit model to rel force vs pCa data in mouse LV

if Fit_rel_SSF_mouse_LV_data
    
    close all
    clc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mechanical data from mouse LV
    pCa_mouse_LV = [6.2; 6.1; 6.0; 5.9; 5.8; 5.7; 5.6; 5.5; 5.4; 4.5];
    rel_SSF_mean_mouse_LV = [0.043; 0.093; 0.189; 0.348; 0.541; 0.705; 0.824; 0.896; 0.930; 1];
    rel_SSF_SEM_mouse_LV = [0.006; 0.015; 0.021; 0.024; 0.011; 0.010; 0.007; 0.011; 0.009; 0];
    %%%%%%%%%%%%
    % Error function - the norm of the difference between actual observations
    % (data) and predicted observations (the solution of the model)
    ftns_6 = @(para_fit) norm(rel_SSF_mean_mouse_LV - rel_force(para_fit,pCa_mouse_LV));
    Parms_6 = 20; % the number of fitting parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial values of the parameters to be fitted
    pCa_50 = 5.79380;
    k0_BC = 0.00000;
    kCa_BC = 7.91031;
    k0_CB = 16.52960;
    kCa_CB = 7.03212;
    f0_CM1 = 0.19514;
    f0_M1C = 2.53993;
    k_M1M2 = 5.27425;
    k_M2M1 = 11.32683;
    k_M2C = 0.29365;
    alpha = 0.94244;
    alpha_bar = 0.89144;
    beta = 0.92448;
    beta_bar = 0.99662;
    u1 = 1.00856;
    u2 = 17.56493;
    z1 = 1.00096;
    z2 = 4.51009;
    v = 1.00051;
    w = 1.00235;
    % RMSE = 0.00760
    % Create the vector of parameters to be fitted
    para_fit_6 = zeros(1,20);
    
    para_fit_6(1) = pCa_50;
    para_fit_6(2) = k0_BC;
    para_fit_6(3) = kCa_BC;
    para_fit_6(4) = k0_CB;
    para_fit_6(5) = kCa_CB;
    para_fit_6(6) = f0_CM1;
    para_fit_6(7) = f0_M1C;
    para_fit_6(8) = k_M1M2;
    para_fit_6(9) = k_M2M1;
    para_fit_6(10) = k_M2C;
    para_fit_6(11) = alpha;
    para_fit_6(12) = alpha_bar;
    para_fit_6(13) = beta;
    para_fit_6(14) = beta_bar;
    para_fit_6(15) = u1;
    para_fit_6(16) = u2;
    para_fit_6(17) = z1;
    para_fit_6(18) = z2;
    para_fit_6(19) = v;
    para_fit_6(20) = w;
    lower_bound_6 = [5.79;zeros(13,1);ones(6,1)];
    upper_bound_6 = [5.81;Inf(9,1);ones(4,1);Inf(6,1)];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % t0 = clock;
    % fprintf('\nStart Time: %4d-%02d-%02d %02d:%02d:%07.4f\n',t0)
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % opts = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    % % Fit the model to the mechanical data from porcine LV
    % [fitted_para_6,resnorm_6,residual_6,exitflag_6,output_6,lambda_6,jacobian_6] = lsqcurvefit( ...
    %     @rel_force,para_fit_6,pCa_mouse_LV,rel_SSF_mean_mouse_LV,lower_bound_6,...
    %     upper_bound_6,opts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimization options for ps
    opts = optimoptions('patternsearch','Display','iter','PlotFcn',@psplotbestf);
    % pattern-search algorithm
    [fitted_para_6,fval_6,exitflag_6,output_6]=patternsearch(ftns_6,para_fit_6,...
                                    [],[],[],[],lower_bound_6,upper_bound_6,[],opts);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % t1 = clock;
    % fprintf('\nStop Time: %4d-%02d-%02d %02d:%02d:%07.4f\n', t1)
    % 
    % GA_Time = etime(t1,t0)
    % 
    % % Show the time elapsed
    % DT_GA_Time = datetime([0 0 0 0 0 GA_Time], 'Format', 'HH:mm:ss.SSSS');
    % fprintf('\nElapsed Time: %23.15E\t\t%s\n\n', GA_Time, DT_GA_Time)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print out values of the new fitted parameters
    fprintf('fitted parameters for mouse LV data:\n')
    fprintf('pCa50 = %f \n',fitted_para_6(1))
    fprintf('k0_BC = %f, kCa_BC = %f, k0_CB = %f, kCa_CB = %f \n', fitted_para_6(2:5))
    fprintf('f0_CM1 = %f, f0_M1C = %f \n', fitted_para_6(6:7))
    fprintf('k_M1M2 = %f, k_M2M1 = %f, k_M2C = %f \n', fitted_para_6(8:10))
    fprintf('alpha = %f, alpha_bar = %f, beta = %f, beta_bar = %f \n', fitted_para_6(11:14))
    fprintf('u1 = %f, u2 = %f \n', fitted_para_6(15:16))
    fprintf('z1 = %f, z2 = %f \n', fitted_para_6(17:18))
    fprintf('v = %f, w = %f \n', fitted_para_6(19:20))
    fprintf('RMSE = %f', sqrt(resnorm_6))
    % fprintf('RMSE = %f', fval_6)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the fitted rate of force redevelopment
    n = 500;
    pCa_start = 6.3;
    pCa_end = 4.5;
    p_Ca = linspace(pCa_start,pCa_end,n);
    fitted_relf_6 = rel_force(fitted_para_6,p_Ca);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot the mechanical data and the fitted solution in the same graph
    clf
    h = axes('Position',[0 0 1 1],'Visible','off');
    axes('Position',[.1 .1 .62 .8]);
    ax = gca;
    errorbar(pCa_mouse_LV,rel_SSF_mean_mouse_LV,rel_SSF_SEM_mouse_LV,'.','MarkerSize',25,...
         'MarkerEdgeColor','blue','Color','blue','LineWidth',1); hold on
    plot(p_Ca,fitted_relf_6,'k-','LineWidth',2);
    set(gca, 'XDir','reverse');
    ax.FontSize = 16
    xlabel('pCa','FontSize',16);
    ax.XAxis.LineWidth = 1.5;
    ylabel('Relative force','FontSize',16);
    ax.YAxis.LineWidth = 1.5;
    legend('Mouse LV data','fitted solution','Location','best');
    
    title('A', 'FontSize',20)
    ax.TitleHorizontalAlignment = 'left';
    str(1) = {'Fitted parameter set:'};
    str(2) = {[' $pCa_{50}$ = ', num2str(fitted_para_6(1))]};
    str(3) = {[' $k^0_{BC}$ = ', num2str(fitted_para_6(2))]};
    str(4) = {[' $k^{Ca}_{BC}$ = ', num2str(fitted_para_6(3))]};
    str(5) = {[' $k^0_{CB}$ = ', num2str(fitted_para_6(4))]};
    str(6) = {[' $k^{Ca}_{CB}$ = ', num2str(fitted_para_6(5))]};
    str(7) = {[' $f^0_{CM_1}$ = ', num2str(fitted_para_6(6))]};
    str(8) = {[' $f^0_{M_1C}$ = ', num2str(fitted_para_6(7))]};
    str(9) = {[' $k_{M_1M_2}$ = ', num2str(fitted_para_6(8))]};
    str(10) = {[' $k_{M_2M_1}$ = ', num2str(fitted_para_6(9))]};
    str(11) = {[' $k_{M_2C}$ = ', num2str(fitted_para_6(10))]};
    str(12) = {[' $\alpha$ = ', num2str(fitted_para_6(11))]};
    str(13) = {[' $\bar{\alpha}$ = ', num2str(fitted_para_6(12))]};
    str(14) = {[' $\beta$ = ', num2str(fitted_para_6(13))]};
    str(15) = {[' $\bar{\beta}$ = ', num2str(fitted_para_6(14))]};
    str(16) = {[' $u_1$ = ', num2str(fitted_para_6(15))]};
    str(17) = {[' $u_2$ = ', num2str(fitted_para_6(16))]};
    str(18) = {[' $z_1$ = ', num2str(fitted_para_6(17))]};
    str(19) = {[' $z_2$ = ', num2str(fitted_para_6(18))]};
    str(20) = {[' $v$ = ', num2str(fitted_para_6(19))]};
    str(21) = {[' $w$ = ', num2str(fitted_para_6(20))]};
    str(22) = {[' RMSE = ', num2str(sqrt(resnorm_6))]};
    % str(22) = {[' RMSE = ', num2str(fval_6)]};
    set(gcf,'CurrentAxes',h)
    text(.77,.5,str,'Interpreter','latex','FontSize',16)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set(0,'Units','normalized')
    % set(gcf, 'PaperSize', [10 8], 'PaperPosition', [0 0 10 8])
    % print('fig2A','-dpdf')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end