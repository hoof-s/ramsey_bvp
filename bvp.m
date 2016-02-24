clear all; clc;

% parameters
a = 0.50;           % alpha
d = 0.05;           % delta
r = 0.05;           % rho
g = (d+r)/(a*d);    % gamma

% steady state
k_ss = (a/(d+r))^(1/(1-a));
c_ss = k_ss^a - d*k_ss;

% initial condition
k_0 = 0.5*k_ss;
c_0 = k_0^a - d*k_0;

% differential equations with y(1) = k and y(2) = c
dy = @(t,y) [y(1)^a - d*y(1) - y(2); ...        % dk/dt
             y(2)/g*(a*y(1)^(a-1) - d - r)];    % dc/dt
         
% boundary conditions
bc = @(y0,yT) [y0(1) - k_0;     % initial condition
               yT(2) - c_ss];   % terminal condition

% approximation
T = 400;                                            % time horizon
solinit = bvpinit(linspace(0,T,10), [k_0 c_0]);     % inititalization
sol = bvp4c(dy, bc, solinit);                       % call solver
t = linspace(0,T)';                                 % time axis
y = deval(sol,t)';                                  % evaluate solution
k_a = y(:,1);                                       % transform variables
c_a = y(:,2);                                       % transform variables

% true solution
k_t = (1/(d*g) + (k_0^(1-a) - 1/(d*g))*exp(d*(a-1).*t)).^(1/(1-a));
c_t = (1-1/g)*k_t.^a;

% error
norm_k = norm(k_a-k_t,Inf);
norm_c = norm(c_a-c_t,Inf);

% figures
subplot(211)
plot(t, k_a)
xlabel('t')
ylabel('k(t)')
subplot(212)
plot(t, c_a)
xlabel('t')
ylabel('c(t)')
suptitle('Capital k(t) and Consumption c(t) trajectories')






