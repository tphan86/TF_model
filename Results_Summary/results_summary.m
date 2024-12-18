%% Options of execution
results_summary_porcine_LV = 0;
results_summary_mouse_LV = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results Summary in Porcine LV - Parameter sets 1, 2, and 8

if results_summary_porcine_LV
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract the fitted parameters for porcine LV data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    porcine_LV_1 = load('fitted_porcine_LV_1.mat');
    porcine_LV_2 = load('fitted_porcine_LV_2.mat');
    porcine_LV_8 = load('fitted_porcine_LV_8.mat');
    param_set = zeros(3,size(porcine_LV_8.ktr, 1));
    param_set(1,:) = porcine_LV_1.ktr';
    param_set(2,:) = porcine_LV_2.ktr';
    param_set(3,:) = porcine_LV_8.ktr';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the residual fitted parameter set for porcine LV data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pCa_50 = param_set(:,1);
    ref_para = param_set(:,2:10);
    ref_para_resid = ref_para;
    ref_para_resid(:,6) = param_set(:,6) * (1E8 - 1);
    near_neighbor_para = param_set(:,11:20);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute (C,M1,M2,K,N,lambda_cyc,lambda_cyc_M2,ktr) at pCa = 4.5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    func = @filament_function;
    pCa = 4.5;
    B = zeros(3,1);
    C = zeros(3,1);
    M1 = zeros(3,1);
    M2 = zeros(3,1);
    K = zeros(3,1);
    N = zeros(3,1);
    lambda_cyc = zeros(3,1);
    lambda_cyc_M2 = zeros(3,1);
    ktr = zeros(3,1);
    fprintf('fitted data for porcine LV at maximal Ca2+ activation\n')
    for i = 1:3
        [C(i),M1(i),M2(i),K(i),N(i),lambda_cyc(i),lambda_cyc_M2(i),ktr(i)] =...
        TF_model(pCa,pCa_50(i),ref_para(i,:),ref_para_resid(i,:),near_neighbor_para(i,:),func);
        B(i) = 1 - C(i) - M1(i) - M2(i);
        if i ~= 3
            fprintf('For parameter set %d:\n', i)
            fprintf('B = %.4f, C = %.4f, M1 = %.4f, M2 = %.4f\n', B(i), C(i), M1(i), M2(i))
            fprintf('K = %f, N = %.4f\n', K(i), N(i))
            fprintf('lambda_cyc = %.4f, lambda_cyc_M2 = %.4f\n', lambda_cyc(i), lambda_cyc_M2(i))
            fprintf('ktr = %.4f', ktr(i))
        else
            fprintf('For parameter set 8:\n')
            fprintf('B = %.4f, C = %.4f, M1 = %.4f, M2 = %.4f\n', B(i), C(i), M1(i), M2(i))
            fprintf('K = %f, N = %.4f\n', K(i), N(i))
            fprintf('lambda_cyc = %.4f, lambda_cyc_M2 = %.4f\n', lambda_cyc(i), lambda_cyc_M2(i))
            fprintf('ktr = %.4f', ktr(i))
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Results Summary in Murine LV - Parameter sets 1, 2, and 8

if results_summary_mouse_LV
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Extract the fitted parameters for mouse LV data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    mouse_LV_1 = load('fitted_mouse_LV_1.mat');
    mouse_LV_2 = load('fitted_mouse_LV_2.mat');
    mouse_LV_8 = load('fitted_mouse_LV_8.mat');
    param_set = zeros(3,size(mouse_LV_8.ktr, 1));
    param_set(1,:) = mouse_LV_1.ktr';
    param_set(2,:) = mouse_LV_2.ktr';
    param_set(3,:) = mouse_LV_8.ktr';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create the residual fitted parameter set for mouse LV data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    pCa_50 = param_set(:,1);
    ref_para = param_set(:,2:10);
    ref_para_resid = ref_para;
    ref_para_resid(:,6) = param_set(:,6) * (1E8 - 1);
    near_neighbor_para = param_set(:,11:20);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Compute (C,M1,M2,K,N,lambda_cyc,lambda_cyc_M2,ktr) at pCa = 4.5
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    func = @filament_function;
    pCa = 4.5;
    B = zeros(3,1);
    C = zeros(3,1);
    M1 = zeros(3,1);
    M2 = zeros(3,1);
    K = zeros(3,1);
    N = zeros(3,1);
    lambda_cyc = zeros(3,1);
    lambda_cyc_M2 = zeros(3,1);
    ktr = zeros(3,1);
    fprintf('fitted data for mouse LV at maximal Ca2+ activation\n')
    for i = 1:3
        [C(i),M1(i),M2(i),K(i),N(i),lambda_cyc(i),lambda_cyc_M2(i),ktr(i)] =...
        TF_model(pCa,pCa_50(i),ref_para(i,:),ref_para_resid(i,:),near_neighbor_para(i,:),func);
        B(i) = 1 - C(i) - M1(i) - M2(i);
        if i ~= 3
            fprintf('For parameter set %d:\n', i)
            fprintf('B = %.4f, C = %.4f, M1 = %.4f, M2 = %.4f\n', B(i), C(i), M1(i), M2(i))
            fprintf('K = %f, N = %.4f\n', K(i), N(i))
            fprintf('lambda_cyc = %.4f, lambda_cyc_M2 = %.4f\n', lambda_cyc(i), lambda_cyc_M2(i))
            fprintf('ktr = %.4f', ktr(i))
        else
            fprintf('For parameter set 8:\n')
            fprintf('B = %.4f, C = %.4f, M1 = %.4f, M2 = %.4f\n', B(i), C(i), M1(i), M2(i))
            fprintf('K = %f, N = %.4f\n', K(i), N(i))
            fprintf('lambda_cyc = %.4f, lambda_cyc_M2 = %.4f\n', lambda_cyc(i), lambda_cyc_M2(i))
            fprintf('ktr = %.4f', ktr(i))
        end
    end
end