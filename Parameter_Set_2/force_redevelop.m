function [value, isterminal, direction] = force_redevelop(t,u,threshold)
    % This function stops the integration as the third component
    % passes a threshold.
    value = u(3) - threshold;
    isterminal = 1;
    direction = 0;
end