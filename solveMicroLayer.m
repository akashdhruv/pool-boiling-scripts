clear
close all

%% Introduction

% Microlayer and contact line model used to account for evaporative
% heat flux in nucleate boiling simulation. The model consits of a 4th
% order differential equation which is reduced to a set of 4 first oder
% differential equations (as shown below) and solved using an euler 
% explicit method.

% % % Equations

% z1 = theta; 
% z2 = dtheta/dr;  
% z3 = K(Re/We) + Abar/(theta^3);
% z4 = (RePr)/(3*St) * theta^3 * d/dr[K(Re/We) + Abar/(theta^3)];

% Where, 
% z1 = Film thickness
% z2 = Local slope of the interface
% z3 = Pc = Pl - Pv = capillary pressure
% z4 = Mass flow rate
% K  = Interface curvature = theta''/[1+(theta')^2]^3/2

% dz1/dr = z2;
% dz2/dr = {[z3 - Abar/(theta^3)]/(Re/We)}*(1+z2^2)^3/2
% dz3/dr = (3*St)/(Re*Pr) * z4/z1^3;
% dz4/dr = [Tsat - Twall + (Bbar/Re)*z3]/[z1 + Cbar/rho]

% Where,
% dz4/dr = Evaporative heat flux

% % % Boundary conditions

% Boundary conditions are vaguely defined across all literature. However 
% five solid boundary conditions are given by,

% % At lower limit

% theta     = theta0 (lower limit of film thickness ~ 1e-9 m)
% theta'    = 0
% theta'''  = 0
% Heat flux = 0

% which implies,

% z1        = theta0
% z2        = 0
% z4        = 0

% % At upper limit

% theta     = thetaMacro (upper limit of film thickness ~ dy/2)
% theta''   = 0
% Heat flux = 0

% which implies,

% z1        = thetaMacro (lower limit of film thickness ~ 1e-9 m)
% z3        = 0

% Since a complete set of boundary condition is missing on either sides,
% all the references have solved these equations with assumptions which has
% led to different profiles. Also, recently, Yazdani .et.al have come up 
% with an analytical solution using matched asymptotic theory

% Following is the soultion based on boundary conditions used by Dhir.et.al

% Author - A. V. Dhruv %

%% Parameters for micro and macro region

% % Fluid - Vapor parameters
sig    = 0.059;       % Suface tension
g      = 9.8;         % Acceleration due to gravity
rho_l  = 958.4;       % Density of liquid
rho_g  = 0.59;        % Density of vapor
mu_l   = 280.0e-6;    % Viscosity of liquid
mu_g   = 12.6e-6;     % Viscosity of vapor
cp_l   = 4216.0;      % Specific heat of liquid
cp_g   = 2030.0;      % Specific heat of vapor
k_l    = 0.679;       % Thermal conducitivity of liquid
k_g    = 0.025;       % Thermal conducitivity of vapor
h_gl   = 2260.0e3;    % Latent heat of vaporization
Ts     = 373.12;      % Saturation temperature
Tb     = 373.12;      % Bulk liquid temperature
Tw     = 379.3120;    % Wall temperature
A      = 8.5e-21;     % Hamakar constant
Rg     = 461.5;       % Gas constant for water in (J/KgK)

[lo,to,uo,rho_d,mu_d,cp_d,k_d,alpha_d,Re,Pr,Pe,St,Fr,We,Tsat,Tbulk,Twall,Abar,Bbar,Cbar] = getScalingNumbers(sig,g,rho_l,rho_g,mu_l,mu_g,cp_l,cp_g,k_l,k_g,h_gl,Ts,Tb,Tw,A,Rg);

rho = 1/rho_d;        % rho = rho_g/rho_l

% % Macro region
lx  = 1;              % Length of Domain - Dimensionless 
nx  = 100;            % Number of points on macro-scale
dx  = lx/nx;          % dx = dy on macro scale

% % Micro region
h   = dx/2;           % Maximum film thickness for micro scale
b   = 0.3*dx;         % Arbitary contact line for a macro cell
psi = 38*(pi/180);    % Contact angle
R   = h/tan(psi);     % Film radius - upper limit
R0  = 0;              % Film radius - lower limit

%% Solution

% % Integration steps
N      = 5000;
step   = R/N;

% % Arrays
z1     = zeros(N,1);  % z1
z2     = zeros(N,1);  % z2
z3     = zeros(N,1);  % z3
z4     = zeros(N,1);  % z4

q      = zeros(N,1);  % Heat flux
T_int  = zeros(N,1);  % Interface temperature
K      = zeros(N,1);  % Interface curvature

% % Initial conditions at upper boundary

z1(1) = h;
z2(1) = tan(psi);
z3(1) = Abar/(h^3);
z4(1) = 0.0;
step  = -step;        % change direction for time stepping
R     = R+R0;         % swap upper and lower limits
R0    = R-R0;         % swap upper and lower limits
R     = R-R0;         % swap upper and lower limits

% % Initial conditions at lower boundary

% z1(1) = ((Abar*Bbar)/(Re*(Twall-Tsat)))^(1/3);
% z2(1) = 0.0;
% z3(1) = ((Twall-Tsat)*Re)/Bbar;
% z4(1) = 0.0;

% % Solve using 1st order Euler method

for i = 2:N
     
    z1(i) = z1(i-1) + step*z2(i-1);
    z2(i) = z2(i-1) + step*((z3(i-1)-Abar/(z1(i-1)^3))/(Re/We))*((1+z2(i-1)^2)^(3/2));
    z3(i) = z3(i-1) + step*(3*St/Pe)*(z4(i-1)/(z1(i-1)^3));
    z4(i) = z4(i-1) + step*(Tsat - Twall + (Bbar/Re)*z3(i-1))/(z1(i-1) + (Cbar/rho));

end

%% Post Processing

% % Calculate required quantities

for i=1:N
   
    q(i) = (Twall - Tsat - (Bbar/Re)*z3(i))/(z1(i) + (Cbar/rho));
    T_int(i) = Twall - q(i)*z1(i);
    K(i) = (z3(i) - Abar/(z1(i)^3))/(Re/We);
    
end

x = linspace(R0,R,N);

qflux = sum(abs(step).*q);

% % Plots

% Curvature
figure
plot(x,-z3,'-b')
xlim([min(R0,R) max(R0,R)])
xlabel('r/lo')
ylabel('Pc')

% Film thickness
figure
plot(x,z1,'-b')
xlim([min(R0,R) max(R0,R)])
xlabel('r/lo')
ylabel('y/lo')

% Heat flux
figure
plot(x,q,'-r')
xlim([min(R0,R) max(R0,R)])
xlabel('r/lo')
ylabel('Nu')

% Interface temperature
figure
plot(x,T_int,'-k')
xlim([min(R0,R) max(R0,R)])
xlabel('r/lo')
ylabel('\Theta')

