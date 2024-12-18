%% Options of execution - paramter set 1
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mechanical data from porcine LV
    pCa_porcine_LV = [6.1; 6.0; 5.9; 5.8; 5.7; 5.6; 5.5; 5.4; 4.5];
    ktr_mean_porcine_LV = [2.43; 1.79; 1.47; 1.09; 1.21; 1.67; 2.22; 2.51; 3.32];
    ktr_SEM_porcine_LV = [0.15; 0.15; 0.14; 0.11; 0.09; 0.09; 0.11; 0.11; 0.17];
    %%%%%%%%%%%%
    % Error function - the norm of the difference between actual observations
    % (data) and predicted observations (the solution of the model)
    ftns_1 = @(para_fit) norm(ktr_mean_porcine_LV - rate_force_redev_1(para_fit,pCa_porcine_LV));
    Parms_1 = 20; % the number of fitting parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial values of the parameters to be fitted
    pCa_50 = 5.65;
    k0_BC = 0;
    kCa_BC = 208.4857;
    k0_CB = 5907.6526;
    kCa_CB = 246.3005;
    f0_CM1 = 103.8106;
    f0_M1C = 157.6488;
    k_M1M2 = 14.6009;
    k_M2M1 = 1.6750;
    k_M2C = 0.0367;
    alpha = 1;
    alpha_bar = 1;
    beta = 1;
    beta_bar = 1;
    u1 = 1;
    % u2 = 15.09571;
    u2 = 1;
    z1 = 1;
    % z2 = 1.13477;
    z2 = 1;
    v = 1;
    w = 1;
    % RMSE = 0.27265
    para_fit_1=[pCa_50;k0_BC;kCa_BC;k0_CB;kCa_CB;f0_CM1;f0_M1C;k_M1M2;k_M2M1;k_M2C;
                alpha;alpha_bar;beta;beta_bar;u1;u2;z1;z2;v;w];
    lower_bound_1 = [0 % pCa_50
                     0    % k0_BC
                     0    % kCa_BC
                     0    % k0_CB
                     0    % kCa_CB
                     0    % f0_CM1
                     0    % f0_M1C
                     0    % k_M1M2
                     0    % k_M2M1
                     0    % k_M2C
                     1    % alpha
                     1    % alpha_bar
                     1    % beta
                     1    % beta_bar
                     1    % u1
                     1    % u2
                     1    % z1
                     1    % z2
                     1    % v
                     1    % w
                     ];
    upper_bound_1 = [Inf % pCa_50
                     Inf    % k0_BC
                     Inf    % kCa_BC
                     Inf    % k0_CB
                     Inf    % kCa_CB
                     Inf    % f0_CM1
                     Inf    % f0_M1C
                     Inf    % k_M1M2
                     Inf    % k_M2M1
                     Inf    % k_M2C
                     1    % alpha
                     1    % alpha_bar
                     1    % beta
                     1    % beta_bar
                     1    % u1
                     1    % u2
                     1    % z1
                     1    % z2
                     1    % v
                     1    % w
                     ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimization options for ps
    opts = optimoptions('patternsearch', 'Display', 'iter',...
                        'Cache','off', 'CompletePoll','on',...
                        'InitialMeshSize',0.1, 'MaxIterations',2,...
                        'MaxFunEvals',1e5, 'ScaleMesh','off',...
                        'MeshTolerance', 1e-10,...
                       'PlotFcn',@psplotbestf);
    % pattern-search algorithm

    rmse = 10;
    while rmse > 8
        [fitted_para_1,fval_1,exitflag_1,output_1]=patternsearch(ftns_1,para_fit_1,...
                                    [],[],[],[],lower_bound_1,upper_bound_1,[],opts);
        para_fit_1 = fitted_para_1;
        rmse = fval_1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Genetic Algorithm
    % % func = @(para_fit) norm(Ctrl - total_host_cells(para_fit,Time));
    % Params = 10;
    % PopSz = 50;
    % optshf = optimoptions('ga','Display','iter','PopulationSize',PopSz,...
    %                       'InitialPopulationMatrix',...
    %                        randi(1E+4,PopSz,Params)*1E-3,...
    %                        'MaxGenerations',500,...
    %                        'FunctionTolerance',1E-10,...
    %                        'HybridFcn',@fmincon,...
    %                        'PlotFcn',@gaplotbestf);
    % [fitted_para_1,fval_1,exitflag_1,output_1,population_1,scores_1]=ga(ftns_1,Params,...
    %                        [],[],[],[],lower_bound_1,upper_bound_1,[],[],optshf);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % fprintf('RMSE = %f', sqrt(resnorm_1))
    fprintf('RMSE = %f', fval_1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fitted_result = [fitted_para_1;fval_1];
    % save('fitted_ktr_pCa_porcineLV_1.mat','fitted_result')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the fitted rate of force redevelopment
    n = 500;
    pCa_start = 6.1;
    pCa_end = 4.5;
    p_Ca = linspace(pCa_start,pCa_end,n);
    fitted_ktr_1 = rate_force_redev(fitted_para_1,p_Ca);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    ylabel('Rate of force redevelopment ($s^{-1}$)','Interpreter','latex','FontSize',16);
    ax.YAxis.LineWidth = 1.5;
    legend('Porcine LV data','fitted solution','Location','best');
    
    title('Parameter set 1 - ktr vs pCa Porcine LV', 'FontSize',20)
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
    % str(22) = {[' RMSE = ', num2str(sqrt(resnorm_1))]};
    str(22) = {[' RMSE = ', num2str(fval_1)]};
    set(gcf,'CurrentAxes',h)
    text(.77,.5,str,'Interpreter','latex','FontSize',16)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set(0,'Units','normalized')
    % set(gcf, 'PaperSize', [10 8], 'PaperPosition', [0 0 10 8])
    % print('Fig0D','-dpdf')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit model to ktr vs pCa data in mourse LV
if Fit_ktr_mouse_LV_data
    
    close all
    clc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mechanical data from mouse LV
    pCa_mouse_LV = [6.2; 6.1; 6.0; 5.9; 5.8; 5.7; 5.6; 5.5; 5.4; 4.5];
    ktr_mean_mouse_LV = [3.20; 3.16; 3.42; 4.35; 7.03; 11.04; 16.40; 20.70; 22.52; 36.23];
    ktr_SEM_mouse_LV = [0.28; 0.39; 0.48; 0.54; 0.79; 1.02; 1.19; 1.24; 1.20; 1.79];
    %%%%%%%%%%%%
    % Error function - the norm of the difference between actual observations
    % (data) and predicted observations (the solution of the model)
    ftns_3 = @(para_fit) norm(ktr_mean_mouse_LV - rate_force_redev_1(para_fit,pCa_mouse_LV));
    Parms_3 = 20; % the number of fitting parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial values of the parameters to be fitted
    pCa_50 = 5.63;
    k0_BC = 0.01663;
    kCa_BC = 0.4370;
    k0_CB = 2.3633;
    kCa_CB = 0.3979;
    f0_CM1 = 72.6072;
    f0_M1C = 44.2949;
    k_M1M2 = 50.6703;
    k_M2M1 = 53.5794;
    k_M2C = 55.7214;
    alpha = 1;
    alpha_bar = 1;
    beta = 1;
    beta_bar = 1;
    u1 = 1;
    u2 = 1;
    z1 = 1;
    z2 = 1;
    v = 1;
    w = 1;
    para_fit_3 = [pCa_50;k0_BC;kCa_BC;k0_CB;kCa_CB;f0_CM1;f0_M1C;k_M1M2;k_M2M1;k_M2C;...
                  alpha;alpha_bar;beta;beta_bar;u1;u2;z1;z2;v;w];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lower_bound_3 = [5.59 % pCa_50
                     0    % k0_BC
                     0    % kCa_BC
                     0    % k0_CB
                     0    % kCa_CB
                     0    % f0_CM1
                     0    % f0_M1C
                     0    % k_M1M2
                     0    % k_M2M1
                     0    % k_M2C
                     1    % alpha
                     1    % alpha_bar
                     1    % beta
                     1    % beta_bar
                     1    % u1
                     1    % u2
                     1    % z1
                     1    % z2
                     1    % v
                     1    % w
                     ];
    upper_bound_3 = [5.63 % pCa_50
                     Inf    % k0_BC
                     Inf    % kCa_BC
                     Inf    % k0_CB
                     Inf    % kCa_CB
                     Inf    % f0_CM1
                     Inf    % f0_M1C
                     Inf    % k_M1M2
                     Inf    % k_M2M1
                     Inf    % k_M2C
                     1    % alpha
                     1    % alpha_bar
                     1    % beta
                     1    % beta_bar
                     1    % u1
                     1    % u2
                     1    % z1
                     1    % z2
                     1    % v
                     1    % w
                     ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % optimization options for ps
    opts = optimoptions('patternsearch', 'Display', 'iter',...
                        'Cache','off', 'CompletePoll','on',...
                        'InitialMeshSize',0.1, 'MaxIterations',500,...
                        'MaxFunEvals',1e5, 'ScaleMesh','off',...
                        'MeshTolerance', 1e-10,...
                       'PlotFcn',@psplotbestf);
    % pattern-search algorithm
    
    rmse = 20;
    while rmse > 15
        [fitted_para_3,fval_3,exitflag_3,output_3]=patternsearch(ftns_3,para_fit_3,...
                                    [],[],[],[],lower_bound_3,upper_bound_3,[],opts);
        para_fit_3 = fitted_para_3;
        rmse = fval_3;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % Genetic Algorithm
    % % func = @(para_fit) norm(Ctrl - total_host_cells(para_fit,Time));
    % Params = 20;
    % PopSz = 50;
    % optshf = optimoptions('ga','Display','iter','PopulationSize',PopSz,...
    %                       'InitialPopulationMatrix',...
    %                        randi(1E+4,PopSz,Params)*1E-3,...
    %                        'MaxGenerations',500,...
    %                        'FunctionTolerance',1E-10,...
    %                        'HybridFcn',@fmincon,...
    %                        'PlotFcn',@gaplotbestf);
    % [fitted_para_3,fval_3,exitflag_3,output_3,population_3,scores_3]=ga(ftns_3,Params,...
    %                        [],[],[],[],lower_bound_3,upper_bound_3,[],[],optshf);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fitted_result = [fitted_para_3;fval_3];
    % save('fitted_ktr_pCa_mouseLV_1.mat','fitted_result')
    % End timing
    % elapsedTime = toc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the fitted rate of force redevelopment
    n = 500;
    pCa_start = 6.2;
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
    plot(p_Ca,fitted_ktr_3,'k-','LineWidth',2); hold on
    
    set(gca, 'XDir','reverse');
    ax.FontSize = 16
    xlabel('pCa','FontSize',16);
    ax.XAxis.LineWidth = 1.5;
    ylabel('Rate of force redevelopment ($s^{-1}$)','Interpreter','latex','FontSize',16);
    ax.YAxis.LineWidth = 1.5;
    legend('Mouse LV data','fitted solution','Location','best');

    title('Parameter set 1 - ktr vs pCa Mouse LV', 'FontSize',20)
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
    str(22) = {[' RMSE = ', num2str(fval_3)]};
    % str(22) = {[' RMSE = ', num2str(sqrt(resnorm_3))]};
    % str(12) = {[' RMSE = ', num2str(fval_3)]};
    set(gcf,'CurrentAxes',h)
    text(.77,.5,str,'Interpreter','latex','FontSize',16)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set(0,'Units','normalized')
    % set(gcf, 'PaperSize', [10 8], 'PaperPosition', [0 0 10 8])
    % print('Fig0C','-dpdf')
    % 
    % % Display elapsed time
    % fprintf('Elapsed time: %.2f seconds\n', elapsedTime);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fit model to rel force vs pCa data in porcine LV
if Fit_rel_SSF_porcine_LV_data
    
    close all
    clc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % relative force in porcine LV data
    pCa_porcine_LV = [6.1; 6.0; 5.9; 5.8; 5.7; 5.6; 5.5; 5.4; 4.5];
    rel_SSF_mean_porcine_LV = [0.021; 0.040; 0.081; 0.205; 0.470; 0.723; 0.847; 0.877; 1];
    rel_SSF_SEM_porcine_LV = [0.003; 0.005; 0.008; 0.015; 0.027; 0.023; 0.012; 0.008; 0];
    %%%%%%%%%%%%
    % Error function - the norm of the difference between actual observations
    % (data) and predicted observations (the solution of the model)
    ftns_4 = @(para_fit) norm(rel_SSF_mean_porcine_LV - rel_force(para_fit,pCa_porcine_LV));
    Parms_4 = 20; % the number of fitting parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial values of the parameters to be fitted
    pCa_50 = 5.71;
    k0_BC = 0.00000;
    kCa_BC = 1.34302;
    k0_CB = 9.42520;
    kCa_CB = 0.00000;
    f0_CM1 = 0.47649;
    f0_M1C = 2.98872;
    k_M1M2 = 3.73584;
    k_M2M1 = 5.67186;
    k_M2C = 0.18277;
    alpha = 1;
    alpha_bar = 1;
    beta = 1;
    beta_bar = 1;
    u1 = 1;
    u2 = 1;
    z1 = 1;
    z2 = 1;
    v = 1;
    w = 1;
    para_fit_4 = [pCa_50;k0_BC;kCa_BC;k0_CB;kCa_CB;f0_CM1;f0_M1C;k_M1M2;k_M2M1;k_M2C;
                  alpha;alpha_bar;beta;beta_bar;u1;u2;z1;z2;v;w];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lower_bound_4 = [5.69  % pCa_50   
                     0    % k0_BC
                     0    % kCa_BC
                     0    % k0_CB
                     0    % kCa_CB
                     0    % f0_CM1
                     0    % f0_M1C
                     0    % k_M1M2
                     0    % k_M2M1
                     0    % k_M2C
                     1    % alpha
                     1    % alpha_bar
                     1    % beta
                     1    % beta_bar
                     1    % u1
                     1    % u2
                     1    % z1
                     1    % z2
                     1    % v
                     1    % w
                     ];
    upper_bound_4 = [5.71 % pCa_50
                     Inf    % k0_BC
                     Inf    % kCa_BC
                     Inf    % k0_CB
                     Inf    % kCa_CB
                     Inf    % f0_CM1
                     Inf    % f0_M1C
                     Inf    % k_M1M2
                     Inf    % k_M2M1
                     Inf    % k_M2C
                     1    % alpha
                     1    % alpha_bar
                     1    % beta
                     1    % beta_bar
                     1    % u1
                     1    % u2
                     1    % z1
                     1    % z2
                     1    % v
                     1    % w
                     ];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Genetic Algorithm
    % func = @(para_fit) norm(Ctrl - total_host_cells(para_fit,Time));
    PopSz = 100;
    optshf = optimoptions('ga','Display','iter','PopulationSize',PopSz,...
                          'InitialPopulationMatrix',...
                           randi(1E+4,PopSz,Parms_4)*1E-3,...
                           'MaxGenerations',500,...
                           'FunctionTolerance',1E-10,...
                           'HybridFcn',@fmincon,...
                           'PlotFcn',@gaplotbestf);
    [fitted_para_4,fval_4,exitflag_4,output_4,population_4,scores_4]=ga(ftns_4,Parms_4,...
                           [],[],[],[],lower_bound_4,upper_bound_4,[],[],optshf);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % optimization options for ps
    % opts = optimoptions('patternsearch', 'Display', 'iter',...
    %                     'Cache','off', 'CompletePoll','on',...
    %                     'InitialMeshSize',0.1, 'MaxIterations',500,...
    %                     'MaxFunEvals',1e5, 'ScaleMesh','off',...
    %                     'MeshTolerance', 1e-10,...
    %                    'PlotFcn',@psplotbestf);
    % % pattern-search algorithm
    % rmse = 2;
    % while rmse > 0.619
    %     [fitted_para_4,fval_4,exitflag_4,output_4]=patternsearch(ftns_4,para_fit_4,...
    %                                 [],[],[],[],lower_bound_4,upper_bound_4,[],opts);
    %     para_fit_4 = fitted_para_4;
    %     rmse = fval_4;
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fitted_result = [fitted_para_4';fval_4];
    save('fitted_relf_pCa_porcineLV_1.mat','fitted_result')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % fprintf('RMSE = %f', sqrt(resnorm_4))
    fprintf('RMSE = %f', fval_4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the fitted relative force
    n = 500;
    pCa_start = 6.1;
    pCa_end = 4.5;
    p_Ca = linspace(pCa_start,pCa_end,n);
    fitted_relf_4 = rel_force(fitted_para_4,p_Ca);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    ylabel('Relative force','Interpreter','latex','FontSize',16);
    ax.YAxis.LineWidth = 1.5;
    legend('Porcine LV data','fitted solution','Location','best');
    
    title('Parameter set 1 - relf vs pCa Porcine LV', 'FontSize',20)
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
    % str(22) = {[' RMSE = ', num2str(sqrt(resnorm_4))]};
    str(22) = {[' RMSE = ', num2str(fval_4)]};
    set(gcf,'CurrentAxes',h)
    text(.77,.5,str,'Interpreter','latex','FontSize',16)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set(0,'Units','normalized')
    % set(gcf, 'PaperSize', [10 8], 'PaperPosition', [0 0 10 8])
    % print('Fig0B','-dpdf')
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
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % initial values of the parameters to be fitted
    pCa_50 = 5.79;
    k0_BC = 0.00000;
    kCa_BC = 7.8858;
    k0_CB = 16.6774;
    kCa_CB = 6.9937;
    f0_CM1 = 0.1832;
    f0_M1C = 2.5028;
    k_M1M2 = 7.7427;
    k_M2M1 = 13.7014;
    k_M2C = 0.2915;
    alpha = 1;
    alpha_bar = 1;
    beta = 1;
    beta_bar = 1;
    u1 = 1;
    u2 = 1;
    z1 = 1;
    z2 = 1;
    v = 1;
    w = 1;
    para_fit_6 = [pCa_50;k0_BC;kCa_BC;k0_CB;kCa_CB;f0_CM1;f0_M1C;k_M1M2;k_M2M1;k_M2C;
                  alpha;alpha_bar;beta;beta_bar;u1;u2;z1;z2;v;w];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    lower_bound_6 = [5.77  % pCa_50   
                     0    % k0_BC
                     0    % kCa_BC
                     0    % k0_CB
                     0    % kCa_CB
                     0    % f0_CM1
                     0    % f0_M1C
                     0    % k_M1M2
                     0    % k_M2M1
                     0    % k_M2C
                     1    % alpha
                     1    % alpha_bar
                     1    % beta
                     1    % beta_bar
                     1    % u1
                     1    % u2
                     1    % z1
                     1    % z2
                     1    % v
                     1    % w
                     ];
    upper_bound_6 = [5.79 % pCa_50
                     Inf    % k0_BC
                     Inf    % kCa_BC
                     Inf    % k0_CB
                     Inf    % kCa_CB
                     Inf    % f0_CM1
                     Inf    % f0_M1C
                     Inf    % k_M1M2
                     Inf    % k_M2M1
                     Inf    % k_M2C
                     1    % alpha
                     1    % alpha_bar
                     1    % beta
                     1    % beta_bar
                     1    % u1
                     1    % u2
                     1    % z1
                     1    % z2
                     1    % v
                     1    % w
                     ];
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Genetic Algorithm
    % func = @(para_fit) norm(Ctrl - total_host_cells(para_fit,Time));
    PopSz = 50;
    optshf = optimoptions('ga','Display','iter','PopulationSize',PopSz,...
                          'InitialPopulationMatrix',...
                           randi(1E+4,PopSz,Parms_6)*1E-3,...
                           'MaxGenerations',500,...
                           'FunctionTolerance',1E-10,...
                           'HybridFcn',@fmincon,...
                           'PlotFcn',@gaplotbestf);
    [fitted_para_6,fval_6,exitflag_6,output_6,population_6,scores_6]=ga(ftns_6,Parms_6,...
                           [],[],[],[],lower_bound_6,upper_bound_6,[],[],optshf);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % optimization options for ps
    % opts = optimoptions('patternsearch', 'Display', 'iter',...
    %                     'Cache','off', 'CompletePoll','on',...
    %                     'InitialMeshSize',0.1, 'MaxIterations',500,...
    %                     'MaxFunEvals',1e5, 'ScaleMesh','off',...
    %                     'MeshTolerance', 1e-10,...
    %                    'PlotFcn',@psplotbestf);
    % % pattern-search algorithm
    % rmse = 2;
    % while rmse > 0.554
    %     [fitted_para_6,fval_6,exitflag_6,output_6]=patternsearch(ftns_6,para_fit_6,...
    %                                 [],[],[],[],lower_bound_6,upper_bound_6,[],opts);
    %     para_fit_6 = fitted_para_6;
    %     rmse = fval_6;
    % end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fitted_result = [fitted_para_6';fval_6];
    save('fitted_relf_pCa_mouseLV_1.mat','fitted_result')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    % fprintf('RMSE = %f', sqrt(resnorm_6))
    fprintf('RMSE = %f', fval_6)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute the fitted rate of force redevelopment
    n = 500;
    pCa_start = 6.2;
    pCa_end = 4.5;
    p_Ca = linspace(pCa_start,pCa_end,n);
    fitted_relf_6 = rel_force(fitted_para_6,p_Ca);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    
    title('Parameter set 1 - relf vs pCa mouse LV', 'FontSize',20)
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
    % str(22) = {[' RMSE = ', num2str(sqrt(resnorm_6))]};
    str(22) = {[' RMSE = ', num2str(fval_6)]};
    set(gcf,'CurrentAxes',h)
    text(.77,.5,str,'Interpreter','latex','FontSize',16)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % set(0,'Units','normalized')
    % set(gcf, 'PaperSize', [10 8], 'PaperPosition', [0 0 10 8])
    % print('Fig0A','-dpdf')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
