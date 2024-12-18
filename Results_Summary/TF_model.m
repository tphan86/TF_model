function [C,M1,M2,K,N,lambda_cyc,lambda_cyc_M2,ktr] = TF_model(pCa,pCa_50,ref_para,ref_para_resid,...
                                              neighbor_interaction_para,func)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Ca = 10.^(-pCa);
    Ca_Ca50 = 10.^(-pCa + pCa_50);
    initial = [0;0;0];
    tspan = [0,5];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    k0_BC = ref_para(1);
    kCa_BC = ref_para(2);
    k0_CB = ref_para(3);
    kCa_CB = ref_para(4);
    f0_CM1 = ref_para(5);
    f0_M1C = ref_para(6);
    alpha = neighbor_interaction_para(1);
    alpha_bar = neighbor_interaction_para(2);
    beta = neighbor_interaction_para(3);
    beta_bar = neighbor_interaction_para(4);
    u1 = neighbor_interaction_para(5);
    u2 = neighbor_interaction_para(6);
    z1 = neighbor_interaction_para(7);
    z2 = neighbor_interaction_para(8);
    v = neighbor_interaction_para(9);
    w = neighbor_interaction_para(10);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    p = [ref_para,Ca_Ca50,neighbor_interaction_para];
    p_resid = [ref_para_resid,Ca_Ca50,neighbor_interaction_para];
    [~,y] = ode45(@(r,y) func(r,y,p),tspan,initial);
    [~,z] = ode15s(@(s,z) func(s,z,p_resid),tspan,initial);
    C = y(end,1);
    M1 = y(end,2);
    M2 = y(end,3);
    lambda_cyc = C + M1 + M2;
    lambda_cyc_M2 = M2 / lambda_cyc;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    B = 1 - C - M1 - M2;
    k11_BC = k0_BC + (kCa_BC - k0_BC)*Ca/(Ca_Ca50 + Ca);
    k11_CB = k0_CB + (kCa_CB - k0_CB)*Ca/(Ca_Ca50 + Ca);
    k_BC = k11_BC*(alpha*(1 - B*(1 - u1^(-1)) ...
              + (M1 + M2)*(u2 - 1))^2 ...
              +(1-alpha)*(1 + M2*(exp(w-1)-1))^2);
    k_CB = k11_CB*(alpha_bar*(1 + B*(u1 - 1) ...
              - (M1 + M2)*(1 - u2^(-1)))^2 ...
              + (1-alpha_bar)*(1 + M2*(exp(-w+1)-1))^2);
    K = k_BC / k_CB;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    k_CM1 = f0_CM1*(beta*(1 - B*(1 - z1^(-1)) ...
               + (M1 + M2)*(z2 - 1))^2 ...
               + (1-beta)*(1 + M2*(exp(v-1)-1))^2);
    k_M1C = f0_M1C*(beta_bar*(1 + B*(z1 - 1) ...
               - (M1 + M2)*(1 - z2^(-1)))^2 ...
               + (1-beta_bar)*(1 + M2*(exp(-v+1)-1))^2);
    N = k_CM1 / k_M1C;

    %%%%%%%%%%%%%%%%%
    Yzero = z(end,:);
    tspan1 = [0,20];
    F_resid = z(end,end);
    threshold = F_resid + 0.5*(M2 - F_resid);
    opt = odeset('Events',@(t,u)  force_redevelop(t,u,threshold));
    [t,~] = ode45(@(t,x) func(t,x,p),tspan1,Yzero,opt);
    ktr = - log(1/2)/t(end);
    %%%%%%%%%%%%%%%%%%%%%%%%
end