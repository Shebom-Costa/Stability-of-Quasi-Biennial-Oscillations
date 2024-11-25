% RHS_mean with alpha included, here alpha is the ratio between the viscous
% damping and the thermal damping

function du_dt = rhs_mean_a(u, Re, a, dz, rdz2, nz)
    k1 = 1; k2 = 1;    % zonal wavenumbers for the two waves
    c1 = 1; c2 = -1;   % horizontal phase speed for the two waves  
    ub0 = 0;

    % Extract the zonal wind profile
    ub = u';

    % rhs of mean wind equation
    gw1z_a = a./(1-ub).^4;
    gw1z = (1-a)./(1.*(1-ub).^2);
    gw2z_a = a./(1+ub).^4;
    gw2z = (1-a)./(1.*(1+ub).^2);
    gw10 = (1-a)./(k1.*(ub0-c1).^2);
    gw20 = (1-a)./(k2.*(ub0-c2).^2);
    gw10_a = a./(k1.*(ub0-c1).^4);
    gw20_a = a./(k2.*(ub0-c2).^4);

    % Evaluate effect of gravity wave forcing
    F_gw1z = zeros(1,nz);                  % Vector for wave 1 forcing at time n
    F_gw2z = zeros(1,nz);                  % Vector for wave 2 forcing at time n
    for j = 1:nz
            F_gw1z(j) = gw1z(j) * exp(-(0.5 * (gw10 + gw1z(j)) + sum(gw1z(1:j-2)))*dz) + gw1z_a(j) * exp(-(0.5 * (gw10_a + gw1z_a(j)) + sum(gw1z_a(1:j-2)))*dz);
            F_gw2z(j) = -gw2z(j) * exp(-(0.5 * (gw20 + gw2z(j)) + sum(gw2z(1:j-2)))*dz) - gw2z_a(j) * exp(-(0.5 * (gw20_a + gw2z_a(j)) + sum(gw2z_a(1:j-2)))*dz);
    end

    % Adding wave plus mean flow diffusion terms to get total forcing
    du_dt = zeros(nz, 1);
    du_dt(2:end-1) = F_gw1z(2:end-1) + F_gw2z(2:end-1) + (1/Re) * rdz2 * (ub(3:end) - 2.*ub(2:end-1) + ub(1:end-2));
    du_dt(1) = 0;               % bottom boundary condition
    du_dt(nz) = du_dt(nz-1);    % top boundary condition
end