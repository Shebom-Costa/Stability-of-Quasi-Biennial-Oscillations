% Define constants and parameters
% Non-dimensional physical parameters
Re = 4.7619 ;      % Reynolds number which controls everything
a = 0.6;         % alpha; ratio between viscous and newtonian damping   
ztop = 4  ;

% Numerical values, note that z=0 and z=ztop correspond to j=0 and j=nz-1
nz = 120;                               % number of vertical gridpoints
nzm1 = nz - 1;
dz = 1 * ztop / nz;                     % grid size
dz2 = dz^2;
rdz2 = 1 / dz2;
dt = 2.1e-04;                           % Timestep
time_stop = 80;                         % Final time
nsteps = time_stop / dt;
nplot1 = 2000000;                       % Plot u-movie every nplot1 
nplot2 = 25;                            % Plot contour every so many time steps
nplot2len = nsteps / nplot2;            % Size of the plot array
t0 = 0;                                 % Initial time
tspan = [t0, time_stop];                % Initial and final times

% Now initialize some of the basic arrays for later use
z = linspace(0, ztop, nz);              % Vertical prediction levels
zplot = linspace(0, ztop, nz);          % Plotting levels in vertical
ubplot = zeros(nz, round(nplot2len));   % Time series for mean wind
timeplot = zeros(1, round(nplot2len));

ub0 = 0;                                % Mean wind at z = 0
% ub =  0.1*sin(pi * z / 4) / 1;        % Initial mean wind profile, c1 = 1
% ub = zeros(1, nz);
ub = z * 0.1;                          % Shear flow 

options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'Refine', 1);
% ODE45 solving the numerical approximation
[t, u_sol] = ode45(@(t, u) rhs_mean(u, Re, dz, rdz2, nz), tspan, ub');

% Loop over time steps
idx = 1; 
for i = 1:nsteps
% Store mean wind profile at selected time steps
     if mod(i, nplot2) == 0 && idx <nplot2len
         i_acl = min(i, length(u_sol(:, 1)));
         ubplot(:, idx) = u_sol(i_acl, :);
         timeplot(idx) = t(i_acl);
         idx = idx + 1;  % Increment idx
     end
end

% Interpolate the solution at the specific z value over time
t_s = numel(t);
t_s2 = 2*t_s;
tsamp = linspace(t0,time_stop,t_s2);
u_zstar = interp1(t,u_sol(:,20),tsamp,"spline");
u_zstar2 = interp1(t,u_sol(:,80),tsamp,"spline");

% Plot period of u(z_star), u(z_star2), and u(z_star3)
figure;
plot(tsamp, u_zstar, 'LineWidth', 2, 'DisplayName', 'z*=0.5');
hold on;
plot(tsamp, u_zstar2, 'LineWidth', 1.5, 'Color', 'r', 'DisplayName', 'z*=2');
% yline(0);
hold off;
title('Zonal Period');
xlabel('Time');
xlim([0,60])
ylabel('u(z*)');
legend('show');
grid on;

% Fourier Analysis to Calculate Period

% Compute the Fast Fourier Transform

% Fs = 1 / (t(2) - t(1));  % Sampling frequency 
% Since dt varys with each step, we take an mean of the time-steps between
% points

Fs = 1/mean(diff(tsamp)); % Sampling frequency

L = length(u_zstar);  % Length of wave
% Increasing the resolution of the frequency
%L = 2^nextpow2(length(u_zstar));

Y = fft(u_zstar,L);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = (Fs*(0:(L/2))/L);   % Frequency

[~, max_idx] = max(P1); % Find the index of the maximum amplitude

% Extract the corresponding frequency and calculate the period and
% amplitude
dominant_frequency = f(max_idx);
period = 1 / dominant_frequency;
amplitude = P1(max_idx);
fprintf('Dominant Frequency: %.4f\n', dominant_frequency);
fprintf('Estimated Period: %.4f\n', period);
fprintf('Amplitude: %.4f\n', amplitude);

% Compute the Fast Fourier Transform u_zstar2
L1 = length(u_zstar2);      % Length of wave  
Y1 = fft(u_zstar2);
P2s = abs(Y1/L1);
P1s = P2s(1:L1/2+1);
P1s(2:end-1) = 2*P1s(2:end-1);
fs = (Fs*(0:(L1/2))/L1);    % Frequency


[~, max_idx_s] = max(P1s);   % Find the index of the maximum amplitude (dominant frequency)

% Extract the corresponding frequency and calculate the period and
% amplitude
dominant_frequency_s = fs(max_idx_s);
period2 = 1 / dominant_frequency_s;
amplitude = P1s(max_idx_s);
fprintf('Dominant Frequency: %.4f\n', dominant_frequency_s);
fprintf('Estimated Period: %.4f\n', period);
fprintf('Amplitude: %.4f\n', amplitude);

% Plotting the results
figure;
pcolor(timeplot,z, ubplot);
shading interp;
colormap jet;
colorbar;
xlabel('Time');
xlim([20,time_stop]);
ylabel('Height');
title('Zonal Mean Zonal Wind');

% Create a meshgrid for time and height
[timegrid, zgrid] = meshgrid(timeplot, z);

% % Plotting the results using contourf
figure;
contourf(timegrid, zgrid, ubplot, 20, 'LineStyle', 'none');
colormap jet;
colorbar;
xlabel('Time');
xlim([20,time_stop]);
ylabel('Height');
 title('Zonal Mean Zonal Wind');

% Plotting the results using surf
figure;
surf(timeplot, z, ubplot,'EdgeColor', 'none');
colormap jet;
colorbar;
xlabel('Time');
ylabel('Height');
zlabel('Zonal Mean Zonal Wind');
title('Zonal Mean Zonal Wind - 3D Surface Plot');

