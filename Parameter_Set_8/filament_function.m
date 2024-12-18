function yprime = filament_function(t,y,p)
    % Note: y(1) represents C, the closed state.
    %       y(2) represents M1, the non-generating-force strongly bound state.
    %       y(3) represents M2, the generating force strongly bound state.
    % We assume that R_T = 1, the fixed total number of actin-myosin sites.
    % p = [k_BC^0,k_BC^Ca,k_CB^0,k_CB^Ca,f_CM1^0,f_M1C^0,k_M1M2,k_M2M1,k_M2C,
    % Ca/Ca50,alpha,alpha_bar,beta,beta_bar,u1,u2,z1,z2,v,w].
    k_BC = p(1) + (p(2) - p(1))*p(10)/(1+p(10));
    k_CB = p(3) + (p(4) - p(3))*p(10)/(1+p(10));
    k_CM1 = p(5);
    k_M1C = p(6);
    k_M1M2 = p(7);
    k_M2M1 = p(8);
    k_M2C = p(9);
    alpha = p(11);
    alpha_bar = p(12);
    beta = p(13);
    beta_bar = p(14);
    u1 = p(15);
    u2 = p(16);
    z1 = p(17);
    z2 = p(18);
    v = p(19);
    w = p(20);
    k_BC = k_BC*(alpha*(1-(1-y(1)-y(2)-y(3))*(1-u1^(-1)) + (y(2)+y(3))*(u2-1))^2 ...
                 + (1-alpha)*(1+y(3)*(exp(w-1)-1))^2);
    k_CB = k_CB*(alpha_bar*(1+(1-y(1)-y(2)-y(3))*(u1-1) - (y(2)+y(3))*(1-u2^(-1)))^2 ...
                 + (1-alpha_bar)*(1+y(3)*(exp(-w+1)-1))^2);
    k_CM1 = k_CM1*(beta*(1-(1-y(1)-y(2)-y(3))*(1-z1^(-1)) + (y(2)+y(3))*(z2-1))^2 ...
                 + (1-beta)*(1+y(3)*(exp(v-1)-1))^2);
    k_M1C = k_M1C*(beta_bar*(1+(1-y(1)-y(2)-y(3))*(z1-1) - (y(2)+y(3))*(1-z2^(-1)))^2 ...
                 + (1-beta_bar)*(1+y(3)*(exp(-v+1)-1))^2);
    yprime = zeros(3,1);
    yprime(1) = k_BC*(1-y(1)-y(2)-y(3)) + k_M2C*y(3) + k_M1C*y(2)... 
                - (k_CB + k_CM1)*y(1);
    yprime(2) = k_CM1*y(1) + k_M2M1*y(3) - (k_M1C + k_M1M2)*y(2);
    yprime(3) = k_M1M2*y(2) - (k_M2C + k_M2M1)*y(3);

end