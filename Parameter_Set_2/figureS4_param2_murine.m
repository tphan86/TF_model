close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mechanical data from mouse LV
pCa_mouse_LV = [6.2; 6.1; 6.0; 5.9; 5.8; 5.7; 5.6; 5.5; 5.4; 4.5];
ktr_mean_mouse_LV = [3.20; 3.16; 3.42; 4.35; 7.03; 11.04; 16.40; 20.70; 22.52; 36.23];
ktr_SEM_mouse_LV = [0.28; 0.39; 0.48; 0.54; 0.79; 1.02; 1.19; 1.24; 1.20; 1.79];

%%%%%%%%%%%%
% Error function - the norm of the difference between actual observations
% (data) and predicted observations (the solution of the model)
ftns_3 = @(para_fit) norm(ktr_mean_mouse_LV - rate_force_redev_2(para_fit,pCa_mouse_LV));
Parms_3 = 20; % the number of fitting parameters

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initial values of the parameters to be fitted
pCa_50 = 5.63;
k0_BC = 0.01663;
kCa_BC = 1.93703;
k0_CB = 2.36327;
kCa_CB = 0.39788;
f0_CM1 = 43.00724;
f0_M1C = 50.29489;
k_M1M2 = 39.57034;
k_M2M1 = 47.37937;
k_M2C = 47.72140;
alpha = 1;
alpha_bar = 1;
beta = 1;
beta_bar = 1;
u1 = 1;
u2 = 2.22077;
z1 = 1;
z2 = 1;
v = 1;
w = 1;
para_fit_3 = [pCa_50;k0_BC;kCa_BC;k0_CB;kCa_CB;f0_CM1;f0_M1C;k_M1M2;k_M2M1;k_M2C;...
              u1;u2;z1;z2];
lower_bound_3 = [5.61 % pCa_50
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
upper_bound_3 = [5.65 % pCa_50
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
[fitted_para_3,resnorm_3,residual_3,exitflag_3,output_3,lambda_3,jacobian_3] = lsqcurvefit( ...
    @rate_force_redev_2,para_fit_3,pCa_mouse_LV,ktr_mean_mouse_LV,lower_bound_3,...
    upper_bound_3,opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Print out values of the new fitted parameters
fprintf('fitted parameters for mouse LV data:\n')
fprintf('pCa50 = %f \n',fitted_para_3(1))
fprintf('k0_BC = %f, kCa_BC = %f, k0_CB = %f, kCa_CB = %f \n', fitted_para_3(2:5))
fprintf('f0_CM1 = %f, f0_M1C = %f \n', fitted_para_3(6:7))
fprintf('k_M1M2 = %f, k_M2M1 = %f, k_M2C = %f \n', fitted_para_3(8:10))
% fprintf('alpha = %f, alpha_bar = %f, beta = %f, beta_bar = %f \n', fitted_para_3(11:14))
fprintf('u1 = %f, u2 = %f \n', fitted_para_3(11:12))
fprintf('z1 = %f, z2 = %f \n', fitted_para_3(13:14))
% fprintf('v = %f, w = %f \n', fitted_para_3(19:20))
fprintf('RMSE = %f', sqrt(resnorm_3))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Compute the data resolution matrix
% jacobian_inv = pinv(jacobian_3);
% data_resolution_matrix_2 = jacobian_3 * jacobian_inv;
% 
% % Plot the data resolution matrix as a heatmap
% figure;
% imagesc(data_resolution_matrix_2);
% colorbar;
% title({['Data Resolution Matrix'] 
%        ['ktr vs pCa in mouse LV']
%        ['Parameter set 2']});
% xlabel('Data Points');
% ylabel('Data Points');
% xticks(1:10);
% yticks(1:10);
% axis square;
% set(gca, 'FontSize', 14);
% 
% % Save the heatmap as a JPEG file
% saveas(gcf, 'data_res_matrix_ktr_mouse_2.jpg');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the model resolution matrix
model_resolution_matrix_2 = jacobian_inv * jacobian_3;

% Plot the model resolution matrix as a heatmap
figure;
model_res_matrix_2_resized = imresize(model_resolution_matrix_2, [14, 14]);
imagesc(model_res_matrix_2_resized);
axis equal
axis tight
colorbar;
title({['Model Resolution Matrix'] 
       ['ktr vs pCa in mouse LV']
       ['Parameter set 2']});
xlabel('Parameters');
ylabel('Parameters');
xticks(1:14);
yticks(1:14);
axis square;
set(gca, 'FontSize', 14);

% % Save the heatmap as a JPEG file
% saveas(gcf, 'model_res_matrix_ktr_mouse_2.jpg');