function ktr = rate_force_redev_2(para_fit,pCa)

    pCa_50 = para_fit(1);
    k0_BC = para_fit(2);
    kCa_BC = para_fit(3);
    k0_CB = para_fit(4);
    kCa_CB = para_fit(5); 
    f0_CM1 = para_fit(6);
    f0_M1C = para_fit(7);
    k_M1M2 = para_fit(8);
    k_M2M1 = para_fit(9);
    k_M2C = para_fit(10);
    % alpha = para_fit(11);
    % alpha_bar = para_fit(12);
    % beta = para_fit(13);
    % beta_bar = para_fit(14);
    u1 = para_fit(11);
    u2 = para_fit(12);
    z1 = para_fit(13);
    z2 = para_fit(14);
    % v = para_fit(19);
    % w = para_fit(20);
    alpha = 1;
    alpha_bar = 1;
    beta = 1;
    beta_bar = 1;
    v = 1;
    w = 1;

    n = length(pCa);
    Ca_Ca50 = 10.^(-pCa + pCa_50);
    initial_value = [0;0;0];
    tspan = [0,5];
    neighbor_interaction_para = [alpha,alpha_bar,beta,beta_bar,u1,u2,z1,z2,v,w];
    ref_para = [k0_BC,kCa_BC,k0_CB,kCa_CB,f0_CM1,f0_M1C,k_M1M2,k_M2M1,k_M2C];
    ref_para_resid = [k0_BC,kCa_BC,k0_CB,kCa_CB,f0_CM1,(10^8-1)*f0_CM1,k_M1M2,k_M2M1,k_M2C];
    ktr = zeros(n,1);

    for i = 1:n
        p = [ref_para,Ca_Ca50(i),neighbor_interaction_para];
        p_resid = [ref_para_resid,Ca_Ca50(i),neighbor_interaction_para];
        func = @filament_function;
        [~,y] = ode15s(@(r,y) func(r,y,p),tspan,initial_value);
        [~,z] = ode15s(@(s,z) func(s,z,p_resid),tspan,initial_value);
        F_ss = y(end,end);
        Yzero = z(end,:);
        tspan1 = [0,20];
        F_resid = z(end,end);
        threshold = F_resid + 0.5*(F_ss - F_resid);
        opt = odeset('Events',@(t,u)  force_redevelop(t,u,threshold));
        [t,~] = ode15s(@(t,x) func(t,x,p),tspan1,Yzero,opt);
        ktr(i) = - log(1/2)/t(end);
    end
    
end