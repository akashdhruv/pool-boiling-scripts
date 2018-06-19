clc
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
% sig    = 0.059;       % Suface tension
% g      = 9.8;         % Acceleration due to gravity
% rho_l  = 958.4;       % Density of liquid
% rho_g  = 0.59;        % Density of vapor
% mu_l   = 280.0e-6;    % Viscosity of liquid
% mu_g   = 12.6e-6;     % Viscosity of vapor
% cp_l   = 4216.0;      % Specific heat of liquid
% cp_g   = 2030.0;      % Specific heat of vapor
% k_l    = 0.679;       % Thermal conducitivity of liquid
% k_g    = 0.025;       % Thermal conducitivity of vapor
% h_gl   = 2260.0e3;    % Latent heat of vaporization
% Ts     = 373.12;      % Saturation temperature
% Tb     = 373.12;      % Bulk liquid temperature
% Tw     = 379.3120;    % Wall temperature
% A      = 8.5e-21;     % Hamakar constant
% Rg     = 461.5;       % Gas constant for water in (J/KgK)

sig   = 0.014707;
g     = 9.8;
rho_l = 1508.4;
rho_g = 7.4048;
mu_l  = (3.25e-7)*rho_l;
mu_g  = (1.39e-6)*rho_g;
cp_l  = 940.28;
cp_g  = 691.30;
k_l   = 0.063671;
k_g   = 0.0095023;
h_gl  = 144350;
Ts    = 47 + 273;
Tb    = 47 + 273;
Tw    = 72 + 273;
A     = 1.0e-20;
Rg    = 8.314/0.187376;


% sig    = 8.4e-3;
% g      = 9.8;
% rho_l  = 1621.2;
% rho_g  = 13.491; 
% mu_l   = 4.13e-4;
% mu_g   = 1.19e-5;
% cp_l   = 1106.7;
% cp_g   = 924.81;
% k_l    = 5.4165e-2;
% k_g    = 1.3778e-2;
% h_gl   = 83562;
% Tb     = 50 + 273;
% Tw     = 90 + 273;
% Ts     = 56 + 273;
% Rg = 8.314/0.33804;
% A = 1e-20;

[lo,to,uo,rho_d,mu_d,cp_d,k_d,alpha_d,Re,Pr,Pe,St,Fr,We,Tsat,Tbulk,Twall,Abar,Bbar,Cbar] = getScalingNumbers(sig,g,rho_l,rho_g,mu_l,mu_g,cp_l,cp_g,k_l,k_g,h_gl,Ts,Tb,Tw,A,Rg);

rho = 1/rho_d;        % rho = rho_g/rho_l

% % Macro region
lx  = 2;              % Length of Domain - Dimensionless 
nx  = 40; %100;       % Number of points on macro-scale
dx  = lx/nx;          % dx = dy on macro scale

% % Micro region
h   = dx/2;           % Maximum film thickness for micro scale
b   = 0.3*dx;         % Arbitrary contact line for a macro cell
psi = 50*(pi/180);    % Contact angle
R   = h/tan(psi);     % Film radius - upper limit
R0  = 0;              % Film radius - lower limit

%% Solution

% Integration steps
% N      = 2000;
% step   = R/N;

step   = 0.2d-4;
N      = floor(R/step);

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

% x     = (R-R0).*(1-cos(linspace(pi/2,0,N)));
x = linspace(R,R0,N);

z     = zeros(1,4);

% % Solve using 1st order Euler method

for i = 2:N
    
    step = x(i)-x(i-1);
    
    % Euler explicit
    z = [z1(i-1),z2(i-1),z3(i-1),z4(i-1)];
    a = step*getFun(z,Re,Pe,We,St,Abar,Bbar,Cbar,rho,Tsat,Twall);
    
    z = z + a;
    
%     % RK - 4
%     z = [z1(i-1),z2(i-1),z3(i-1),z4(i-1)];
%     a = step*getFun(z,Re,Pe,We,St,Abar,Bbar,Cbar,rho,Tsat,Twall);
%     b = step*getFun(z+a./2,Re,Pe,We,St,Abar,Bbar,Cbar,rho,Tsat,Twall);
%     c = step*getFun(z+b./2,Re,Pe,We,St,Abar,Bbar,Cbar,rho,Tsat,Twall);
%     d = step*getFun(z+c,Re,Pe,We,St,Abar,Bbar,Cbar,rho,Tsat,Twall);
%     
%     z = z + (a + 2.*b + 2.*c + d)./6;
    
    z1(i) = z(1);
    z2(i) = z(2);
    z3(i) = z(3);
    z4(i) = z(4);

end

%% Post Processing

% % Calculate required quantities

for i=1:N
   
    q(i) = (Twall - Tsat - (Bbar/Re)*z3(i))/(z1(i) + (Cbar/rho));
    T_int(i) = Twall - q(i)*z1(i);
    K(i) = (z3(i) - Abar/(z1(i)^3))/(Re/We);
    
end

St
qflux = sum(abs(step).*q)*dx
fflux = sum(step.*K)*dx
Tsat

% % Plots

% Curvature
figure
plot(x,-K,'-b')
xlim([min(R0,R) max(R0,R)])
xlabel('r/lo')
ylabel('K')

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


