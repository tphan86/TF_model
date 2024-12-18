%% Resolution Matrix for ktr vs pCa in porcine LV - parameter set 2

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
lambda = 1e-2;
ftns_1 = @(para_fit) norm(ktr_mean_porcine_LV ...
                  - rate_force_redev_2(para_fit,pCa_porcine_LV))...
                  + lambda*norm(para_fit);
Parms_1 = 20; % the number of fitting parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial values of the parameters to be fitted
pCa_50 = 5.65;
k0_BC = 0;
kCa_BC = 243.03571;
k0_CB = 1694.75258;
kCa_CB = 0.00045;
f0_CM1 = 91.00800;
f0_M1C = 691.67383;
k_M1M2 = 3.20091;
k_M2M1 = 2.50001;
k_M2C = 0.73673;
alpha = 1;
alpha_bar = 1;
beta = 1;
beta_bar = 1;
u1 = 1;
u2 = 15.09571;
z1 = 1;
z2 = 1.13477;
v = 1;
w = 1;
para_fit_1= [pCa_50;k0_BC;kCa_BC;k0_CB;kCa_CB;f0_CM1;f0_M1C;k_M1M2;k_M2M1;k_M2C;...
             u1;u2;z1;z2];
lower_bound_1 = [5.65 % pCa_50
                 0    % k0_BC
                 0    % kCa_BC
                 0    % k0_CB
                 0    % kCa_CB
                 0    % f0_CM1
                 0    % f0_M1C
                 0    % k_M1M2
                 0    % k_M2M1
                 0    % k_M2C
                 % 1    % alpha
                 % 1    % alpha_bar
                 % 1    % beta
                 % 1    % beta_bar
                 1    % u1
                 1    % u2
                 1    % z1
                 1    % z2
                 % 1    % v
                 % 1    % w
                 ];
upper_bound_1 = [5.65 % pCa_50
                 Inf    % k0_BC
                 Inf    % kCa_BC
                 Inf    % k0_CB
                 Inf    % kCa_CB
                 Inf    % f0_CM1
                 Inf    % f0_M1C
                 Inf    % k_M1M2
                 Inf    % k_M2M1
                 Inf    % k_M2C
                 % 1    % alpha
                 % 1    % alpha_bar
                 % 1    % beta
                 % 1    % beta_bar
                 Inf    % u1
                 Inf   % u2
                 Inf    % z1
                 Inf    % z2
                 % 1    % v
                 % 1    % w
                 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opts = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');

% Fit the model to the mechanical data from porcine LV
[fitted_para_1,resnorm_1,residual_1,exitflag_1,output_1,lambda_1,jacobian_1] = lsqcurvefit( ...
    @rate_force_redev_2,para_fit_1,pCa_porcine_LV,ktr_mean_porcine_LV,lower_bound_1,...
    upper_bound_1,opts);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% [fitted_para_1, resnorm_1, residual_1, exitflag_1, output_1, lambda_1, jacobian_1] = lsqcurvefit( ...
% @(para_fit, pCa_porcine_LV) lsqcurvefit_wrapper(para_fit, pCa_porcine_LV, ktr_mean_porcine_LV, lambda), ...
% para_fit_1, pCa_porcine_LV, ktr_mean_porcine_LV, lower_bound_1, upper_bound_1, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print out values of the new fitted parameters
fprintf('fitted parameters for porcine LV data:\n')
fprintf('pCa50 = %f \n',fitted_para_1(1))
fprintf('k0_BC = %f, kCa_BC = %f, k0_CB = %f, kCa_CB = %f \n', fitted_para_1(2:5))
fprintf('f0_CM1 = %f, f0_M1C = %f \n', fitted_para_1(6:7))
fprintf('k_M1M2 = %f, k_M2M1 = %f, k_M2C = %f \n', fitted_para_1(8:10))
% fprintf('alpha = %f, alpha_bar = %f, beta = %f, beta_bar = %f \n', fitted_para_1(11:14))
fprintf('u1 = %f, u2 = %f \n', fitted_para_1(11:12))
fprintf('z1 = %f, z2 = %f \n', fitted_para_1(13:14))
% fprintf('v = %f, w = %f \n', fitted_para_1(19:20))
fprintf('RMSE = %f', sqrt(resnorm_1))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Compute the data resolution matrix
% jacobian_inv = pinv(jacobian_1);
% data_resolution_matrix_1 = jacobian_1 * jacobian_inv;
% 
% % Plot the data resolution matrix as a heatmap
% figure;
% imagesc(data_resolution_matrix_1);
% colorbar;
% title({['Data Resolution Matrix'] 
%        ['ktr vs pCa in porcine LV']
%        ['Parameter set 2']});
% xlabel('Data Points');
% ylabel('Data Points');
% xticks(1:9);
% yticks(1:9);
% axis square;
% set(gca, 'FontSize', 14);
% 
% % Save the heatmap as a JPEG file
% saveas(gcf, 'data_res_matrix_ktr_porcine_2.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the model resolution matrix
J = jacobian_1;
model_resolution_matrix_1 = pinv(J'*J)*(J'*J);

% Plot the model resolution matrix as a heatmap
figure;
model_res_matrix_1_resized = imresize(model_resolution_matrix_1, [14, 14]);
imagesc(model_res_matrix_1_resized);
axis equal
axis tight
colorbar;
title({['Model Resolution Matrix'] 
       ['ktr vs pCa in porcine LV']
       ['Parameter set 2']});
xlabel('Parameters');
ylabel('Parameters');
xticks(1:14);
yticks(1:14);
axis square;
set(gca, 'FontSize', 14);

% % Save the heatmap as a JPEG file
% saveas(gcf, 'model_res_matrix_ktr_porcine_2.jpg');