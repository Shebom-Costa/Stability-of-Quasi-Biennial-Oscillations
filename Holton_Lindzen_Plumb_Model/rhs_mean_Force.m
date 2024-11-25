% RHS for HL model, which consists of the two most promiment waves types in
% the stratosphere

function du_dt = rhs_mean_Force(u, Re, dz, rdz2, nz)
    k1 = 1; k2 = 3;    % zonal wavenumbers for the two waves
    c1 = 1; c2 = -1;   % horizontal phase speed for the two waves  
    gwdamp = 1;     % GW damping rate
    ub0 = 0;

    % Extract the zonal wind profile
    ub = u';
    u0 = zeros(1,nz);
    % Kelvin wave
    gw1z = 1./(k1.*(1-ub).^2);
    gw10 = gwdamp./(k1.*(ub0-c1).^2);

    % Rossby-Gravity wave
    beta = 32;                      % Change beta as required
    amplit = 1;
    gw20 = amplit*gwdamp./(k2.*(ub0-c2).^2).*(1./(k2*(ub0-c2))-1);
    gpart = amplit*gwdamp./(k2*(1+ub).^2);
    rpart = beta./(k2^2*(1+ub))-1;
    gw2z = gpart.*rpart;

    % Force
    fgw1z = 1./(k1.*(1-u0).^2);
    fgw10 = gwdamp./(k1.*(ub0-c1).^2);

    fgw20 = amplit*gwdamp./(k2.*(ub0-c2).^2).*(1./(k2*(ub0-c2))-1);
    fgpart = amplit*gwdamp./(k2*(1+u0).^2);
    frpart = beta./(k2^2*(1+u0))-1;
    fgw2z = fgpart.*frpart;
    Force = zeros(1,nz);

    % Evaluate effect of gravity wave forcing
    F_gw1z = zeros(1,nz);                  % Vector for wave 1 forcing at time n
    F_gw2z = zeros(1,nz);                  % Vector for wave 2 forcing at time n 
    for j = 1:nz
            F_gw1z(j) = gw1z(j) * exp(-(0.5 * (gw10 + gw1z(j)) + sum(gw1z(1:j-2)))*dz);
            F_gw2z(j) = -gw2z(j) * exp(-(0.5 * (gw20 + gw2z(j)) + sum(gw2z(1:j-2)))*dz);
            Force(j) = fgw1z(j) * exp(-(0.5 * (fgw10 + fgw1z(j)) + sum(fgw1z(1:j-2)))*dz)-fgw2z(j) * exp(-(0.5 * (fgw20 + fgw2z(j)) + sum(fgw2z(1:j-2)))*dz);
    end

    % Adding wave plus mean flow diffusion terms to get total forcing
    du_dt = zeros(nz, 1);
    du_dt(2:end-1) = F_gw1z(2:end-1) + F_gw2z(2:end-1) - Force(2:end-1) + (1/Re) * rdz2 * (ub(3:end) - 2.*ub(2:end-1) + ub(1:end-2));
    du_dt(1) = 0;               % bottom boundary condition
    du_dt(nz) = du_dt(nz-1);    % top boundary condition
end