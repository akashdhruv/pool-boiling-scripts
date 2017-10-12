function[l,t,u,rho_d,mu_d,cp_d,k_d,alpha_d,Re,Pr,Pe,St,Fr,We,Ts,Tb,Tw,Abar,Bbar,Cbar] = getScalingNumbers(sig,g,rho_l,rho_g,mu_l,mu_g,cp_l,cp_g,k_l,k_g,h_gl,Tsat,Tbulk,Twall,A,Rg)

% Non-dimensionalization script

% Author - A. V. Dhruv

% Assuming Bond number, Bo and Froude number, Fr equal to 1.0

l = sqrt((1.0*sig)/(g*(rho_l - rho_g))); % Length scale
u = sqrt(g*l);                           % Velocity scale
t = l/u;                                 % Time scale

dt = Twall - Tbulk;                      % Wall superheat

Tw = (Twall - Tbulk)/dt;                 % Wall temperature
Tb = (Tbulk - Tbulk)/dt;                 % Bulk temperature
Ts = (Tsat - Tbulk)/dt;                  % Saturation temperature

rho_d = rho_l/rho_g;                     % Density ratio
mu_d  = mu_l/mu_g;                       % Viscosity ratio
cp_d  = cp_l/cp_g;                       % Specific heat ratio
k_d   = k_l/k_g;                         % Thermal conductivity ratio
alpha_d = k_d/(rho_d*cp_d);              % Thermal diffusivity ratio

Re    = (rho_l*u*l)/mu_l;                % Reynolds number
Fr    = u/sqrt(g*l);                     % Froude number
Pr    = (mu_l*cp_l)/k_l;                 % Prandtl number
St    = (cp_l*dt)/h_gl;                  % Stefan (Jakob) number
We    = (rho_l*u*u*l)/sig;               % Weber number
Pe    = Re*Pr;                           % Peclet number
Abar  = A/(mu_l*u*l*l);                  % Dimensionless Hamakar constant
Bbar  = (Tsat*u*u)/(h_gl*dt);            % Dimensionless number #1       

Cbar  = (k_l*Tsat*sqrt(2*pi*Rg*Tsat))...
        /(2*l*h_gl*h_gl*rho_l);          % Dimensionless number #2
    
end