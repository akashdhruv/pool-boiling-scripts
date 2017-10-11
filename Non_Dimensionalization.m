
clear
close all
clc

%% Ethanol 101.3 kPa

% sig    = 0.018;
% g      = 9.8;
% rho_l  = 757;
% rho_g  = 1.435;
% mu_l   = 4.29e-4;
% mu_g   = 1.04e-5;
% cp_l   = 3000.0;
% cp_g   = 1830.0;
% k_l    = 0.154;
% k_g    = 0.020;
% dt     = 1.8939;
% h_gl   = 9.63e5;
% 
% l   = sqrt(sig/(g*(rho_l - rho_g)));
% lam = (2*pi)*sqrt(3*sig/(g*(rho_l - rho_g)));
% u = sqrt(g*l);
% 
% rho_d = rho_l/rho_g;
% mu_d  = mu_l/mu_g;
% cp_d  = cp_l/cp_g;
% k_d   = k_l/k_g;
% alpha_d = k_d/(rho_d*cp_d);
% 
% 
% Re    = (rho_l*u*l)/mu_l;
% Fr    = u/sqrt(g*l);
% Pr    = (mu_l*cp_l)/k_l;
% St    = (cp_l*dt)/h_gl;
% We    = (rho_l*u*u*l)/sig;
% Pe    = Re*Pr;
% 
% dt_v = Re/max(1,rho_d/mu_d);
% dt_t = Re*Pr/max(1,1/alpha_d);
% 
% t = l/u;
% nx = 40;
% dx = (l*10)/nx;
% cfl = 0.03;
% 
% dt_c = 0.1*dx/u;
% dt_v = dt_v*cfl*dx*dx;
% dt_t = dt_t*cfl*dx*dx;

%% Water 101.3 kPa

% sig    = 0.059;
% g      = 9.8;
% rho_l  = 958.4;
% rho_g  = 0.59;
% mu_l   = 280.0e-6;
% mu_g   = 12.6e-6;
% cp_l   = 4216.0;
% cp_g   = 2030.0;
% k_l    = 0.679;
% k_g    = 0.025;
% dt     = 6.2;
% h_gl   = 2260.0e3;
% Tsat   = 373.12;
% A      = -8.5e-21;
% Rg     = 461.5;
% 
% l = sqrt((1.0*sig)/(g*(rho_l - rho_g)));
% u = sqrt(g*l);
% 
% rho_d = rho_l/rho_g;
% mu_d  = mu_l/mu_g;
% cp_d  = cp_l/cp_g;
% k_d   = k_l/k_g;
% alpha_d = k_d/(rho_d*cp_d);
% 
% Re    = (rho_l*u*l)/mu_l;
% Fr    = u/sqrt(g*l);
% Pr    = (mu_l*cp_l)/k_l;
% St    = (cp_l*dt)/h_gl;
% We    = (rho_l*u*u*l)/sig;
% Pe    = Re*Pr;
% Abar  = A/(mu_l*u*l*l);
% Bbar  = (Tsat*u*u)/(h_gl*dt);
% Cbar  = (k_l*Tsat*sqrt(2*pi*Rg*Tsat))/(2*l*h_gl*h_gl*rho_l);
% 
% % dt_v = Re/max(1,rho_d/mu_d);
% % dt_t = Re*Pr/max(1,1/alpha_d);
% 
% t = l/u;
% 
% % nx = 40;
% % dx = (l)/nx;
% % cfl = 0.03;
% % 
% % dt_c = 0.1*dx/u;
% % dt_v = dt_v*cfl*dx*dx;
% % dt_t = dt_t*cfl*dx*dx;
% 
% alpha_l = k_l/(rho_l*cp_l);
% nu_l    = mu_l/rho_l;
% 
% % beta    = 0.000752;
% beta    = 0.000214;
% 
% dtheta = 7.14*(((nu_l*alpha_l)/(g*beta*dt))^(1/3)) ;
% dtheta = dtheta/l;
% 
% l_d = l/ 0.0025;
% t_d = t/0.0160;

%% Water Pressure Ratio = 0.99

% sig    = 0.00007;
% g      = 9.8;
% rho_l  = 402.4;
% rho_g  = 242.7;
% mu_l   = 46.7e-6;
% mu_g   = 32.38e-6;
% cp_l   = 2180.0;
% cp_g   = 3520.0;
% k_l    = 0.545;
% k_g    = 0.538;
% dt     = 5.0;
% h_gl   = 276.4e3;
% 
% l   = sqrt(sig/(g*(rho_l - rho_g)));
% lam = (2*pi)*sqrt(3*sig/(g*(rho_l - rho_g)));
% u = sqrt(g*l);
% 
% rho_d = rho_l/rho_g;
% mu_d  = mu_l/mu_g;
% cp_d  = cp_l/cp_g;
% k_d   = k_l/k_g;
% alpha_d = k_d/(rho_d*cp_d);
% 
% Re    = (rho_l*u*l)/mu_l;
% Fr    = u/sqrt(g*l);
% Pr    = (mu_l*cp_l)/k_l;
% St    = (cp_l*dt)/h_gl;
% We    = (rho_l*u*u*l)/sig;
% Pe    = Re*Pr;
% 
% dt_v = Re/max(1,rho_d/mu_d);
% dt_t = Re*Pr/max(1,1/alpha_d);
% 
% t = l/u;
% nx = 40;
% dx = (l*10)/nx;
% cfl = 0.03;
% 
% dt_c = 0.1*dx/u; 
% dt_v = dt_v*cfl*dx*dx;
% dt_t = dt_t*cfl*dx*dx;

%% Virtual Fluid (Gibou)

% sig    = 0.1;
% g      = 9.8;
% rho_l  = 200.0;
% rho_g  = 5.0;
% mu_l   = 0.1;
% mu_g   = 0.005;
% cp_l   = 400.0;
% cp_g   = 200.0;
% k_l    = 40.0;
% k_g    = 1.0;
% dt     = 5.0;
% 
% h_gl   = 10^4;
% 
% l   = sqrt(sig/(g*(rho_l - rho_g)));
% lam = (2*pi)*sqrt(3*sig/(g*(rho_l - rho_g)));
% %l = lam;
% u = sqrt(g*l);
% 
% % l = 1;
% % u = 1;
% 
% rho_d = rho_l/rho_g;
% mu_d  = mu_l/mu_g;
% cp_d  = cp_l/cp_g;
% k_d   = k_l/k_g;
% alpha_d = k_d/(rho_d*cp_d);
% 
% 
% Re    = (rho_l*u*l)/mu_l;
% Fr    = u/sqrt(g*l);
% Pr    = (mu_l*cp_l)/k_l;
% St    = (cp_l*dt)/h_gl;
% We    = (rho_l*u*u*l)/sig;
% Pe    = Re*Pr;
% 
% dt_v = Re/max(1,rho_d/mu_d);
% dt_t = Re*Pr/max(1,1/alpha_d);
% 
% t = l/u;
% nx = 40;
% dx = (l*10)/nx;
% cfl = 0.03;
% 
% dt_c = 0.1*dx/u;
% dt_v = dt_v*cfl*dx*dx;
% dt_t = dt_t*cfl*dx*dx;
% 
% alpha_l = k_l/(rho_l*cp_l);
% nu_l    = mu_l/rho_l;
% 
% %beta    = 0.000214;
% %beta    = 0.000752;
% beta     = 1;
% 
% dtheta = 7.14*(((nu_l*alpha_l)/(g*beta*dt))^(1/3)) ;
% 
% t_ref = [0.355,0.405,0.455,0.505];
% t_num = t_ref/t;
% 
% ly = [0.016,0.0425,0.07110,0.09871,0.1260,0.1535]/l;
% ty = linspace(0.405,0.655,6);
% 
% tx = [13.1,17.6,19.6]*t;

%% Virtual Fluid (Needs Reference)

% sig    = 0.000182;
% g      = 9.8;
% rho_l  = 730.8;
% rho_g  = 310.9;
% mu_l   = 56.72e-6;
% mu_g   = 21.71e-6;
% cp_l   = 5128.0;
% cp_g   = 4445.0;
% k_l    = 0.052;
% k_g    = 0.0415;
% dt     = 30.0;
% h_gl   = 54.6e3;
% 
% l = sqrt(sig/(g*(rho_l - rho_g)));
% 
% lam = (2*pi)*sqrt(3*sig/(g*(rho_l - rho_g)));
% u = sqrt(g*l);
% 
% Re    = (rho_l*u*l)/mu_l;
% Fr    = u/sqrt(g*l);
% Pr    = (mu_l*cp_l)/k_l;
% St    = (cp_l*dt)/h_gl;
% We    = (rho_l*u*u*l)/sig;
% rho_d = rho_l/rho_g;
% mu_d  = mu_l/mu_g;
% cp_d  = cp_l/cp_g;
% k_d   = k_l/k_g;
% alpha_d = k_d/(rho_d*cp_d);
% 
% dt_v = Re/max(1,rho_d/mu_d);
% dt_t = Re*Pr/max(1,1/alpha_d);
% 
% t = l/u;
% nx = 40;
% dx = (l*10)/nx;
% cfl = 0.03;
% 
% dt_c = 0.1*dx/u; 
% dt_v = dt_v*cfl*dx*dx;
% dt_t = dt_t*cfl*dx*dx;

%% FC-72

% sig    = 0.0084;
% 
% g      = 9.8;
% 
% rho_l  = 1621.2;
% rho_g  = 13.491;
% 
% mu_l   = 4.13e-4;
% mu_g   = 1.19e-5;
% 
% cp_l   = 1106.7;
% cp_g   = 924.81;
% 
% k_l    = 0.054165;
% k_g    = 0.013778;
% 
% Tsat   = 57;
% Twall  = 76;
% Tblk  = 52;
% 
% dt = Twall-Tblk;
% 
% Tw = (Twall-Tblk)/(Twall-Tblk);
% Ts = (Tsat-Tblk)/(Twall-Tblk);
% Tb = (Tblk-Tblk)/(Twall-Tblk);
% 
% h_gl   = 83562;
% 
% l   = sqrt(sig/(g*(rho_l - rho_g)));
% lam = (2*pi)*sqrt(3*sig/(g*(rho_l - rho_g)));
% u = sqrt(g*l);
% t = l/u;
% 
% rho_d = rho_l/rho_g;
% mu_d  = mu_l/mu_g;
% cp_d  = cp_l/cp_g;
% 
% k_d   = k_l/k_g;
% alpha_d = k_d/(rho_d*cp_d);
% alpha_l = k_l/(rho_l*cp_l);
% alpha_g = k_g/(rho_g*cp_g);
% 
% 
% Re    = (rho_l*u*l)/mu_l;
% Fr    = u/sqrt(g*l);
% Pr    = (mu_l*cp_l)/k_l;
% St    = (cp_l*dt)/h_gl;
% We    = (rho_l*u*u*l)/sig;
% Pe    = Re*Pr;
% 
% D_eq_init = 0.1955*l*1000;
%% FC-72 properties at 89 kPa

% sig    = 9.27e-3;
% 
% g      = 9.8;
% gz     = 9.8*2.5*(10^-7);
% 
% rho_l  = 1641;
% rho_g  = 9.14;
% 
% mu_l   = 5.09e-4;
% mu_g   = 5.09e-4;
% 
% cp_l   = 1080.0;
% cp_g   = 1080.0;
% 
% k_l    = 5.35e-2;
% k_g    = 5.35e-2;
% 
% Tl   = 43.5;
% Tw   = 57;
% Tsat = 45;
% 
% dt     = Tw-Tl;
% 
% Ts     = (Tsat-Tl)/(Tw-Tl);
% 
% h_gl   = 87789;
% 
% l   = sqrt((1*sig)/(g*(rho_l - rho_g)));
% lam = (2*pi)*sqrt(3*sig/(g*(rho_l - rho_g)));
% u = sqrt(g*l);
% 
% rho_d = rho_l/rho_g;
% mu_d  = mu_l/mu_g;
% cp_d  = cp_l/cp_g;
% 
% k_d     = k_l/k_g;
% alpha_d = k_d/(rho_d*cp_d);
% alpha_l = k_l/(rho_l*cp_l);
% 
% R_initial     = 2.5581e-3/l;
% theta_initial = 7.44e-3/l;
% 
% 
% Re    = (rho_l*u*l)/mu_l;
% Fr    = u/sqrt(g*l);
% Pr    = (mu_l*cp_l)/k_l;
% St    = (cp_l*dt)/h_gl;
% We    = (rho_l*u*u*l)/sig;
% Pe    = Re*Pr;
% gs     = Fr*sqrt(g/gz);
% gs     = 1/(gs^2);
% 
% t = l/u;
% 
% Domain_x = 0.0190/l;
% Domain_y = 0.0761/l;
% 
% t_d = 0.0088/t;

%% R-113 properties 

sig    = 0.014707;

g      = 9.8;

rho_l  = 1508.4;
rho_g  = 7.4048;

mu_l = (3.25e-7)*rho_l;
mu_g = (1.39e-6)*rho_g;

cp_l   = 940.28;
cp_g   = 691.30;

k_l    = 0.063671;
k_g    = 0.0095023;

Tsat = 47.5;
Tw   = 61;

dt     = Tw-Tsat;

h_gl   = 144350;

Rg = 8.314/0.187376;
A = 8.5e-21;

l   = sqrt(sig/(g*(rho_l - rho_g)));
u   = sqrt(g*l);

rho_d = rho_l/rho_g;
mu_d  = mu_l/mu_g;
cp_d  = cp_l/cp_g;

k_d   = k_l/k_g;
alpha_d = k_d/(rho_d*cp_d);
alpha_l = k_l/(rho_l*cp_l);
CP_d = (rho_l*cp_l)/(rho_g*cp_g);

R_initial = 50e-6/l;

theta_initial = 352e-6/l;

bl_thickness = 1.35e-5/l;

Domain_x = 2.5e-3/l;
Domain_y = 4.0e-3/l;

Re    = (rho_l*u*l)/mu_l;
Fr    = u/sqrt(g*l);
Bo    = (g*(rho_l-rho_g)*l*l)/sig;
Pr    = (mu_l*cp_l)/k_l;
St    = (cp_l*dt)/h_gl;
We    = (rho_l*u*u*l)/sig;
Pe    = Re*Pr;
Abar  = A/(mu_l*u*l*l);
Bbar  = (Tsat*u*u)/(h_gl*dt);
Cbar  = (k_l*Tsat*sqrt(2*pi*Rg*Tsat))/(2*l*h_gl*h_gl*rho_l);

t = l/u;

nu_l = mu_l/rho_l;

beta = 0.0026;

dtheta = 7.14*(((nu_l*alpha_l)/(g*beta*dt))^(1/3)) ;

dtheta = dtheta/l;

% R_d = 4.0e-4;
% R_c = (sqrt(27)/2)*St*rho_d*alpha_l*sqrt((rho_l*R_d)/sig);
% t_c = (9/4)*St*rho_d*alpha_l*((rho_l*R_d)/sig);
% 
% u_c = R_c/t_c;
% 
% Re_c    = (rho_l*u_c*R_c)/mu_l;
% Fr_c    = u_c/sqrt(g*R_c);
% Pr_c    = (mu_l*cp_l)/k_l;
% St_c    = (cp_l*dt)/h_gl;
% We_c    = (rho_l*u_c*u_c*R_c)/sig;
% Pe_c    = Re_c*Pr_c;
% 
% R_initial_c = 50e-6/R_c;
% 
% theta_initial_c = 352e-6/R_c;

%% R-113 v2

% sig    = 0.014707;
% 
% g      = 9.8;
% 
% rho_l  = 1508.4;
% rho_g  = 7.4048;
% 
% mu_l = (3.25e-7)*rho_l;
% mu_g = (1.39e-6)*rho_g;
% 
% cp_l   = 940.28;
% cp_g   = 691.30;
% 
% k_l    = 0.063671;
% k_g    = 0.0095023;
% 
% Tb   = 32.0;
% Tsat = 47.5;
% Tw   = 61.0;
% 
% dt  = Tw-Tb;
% Ts  = (Tsat-Tb)/dt;
% 
% h_gl   = 144350;
% 
% rho_d = rho_l/rho_g;
% mu_d  = mu_l/mu_g;
% cp_d  = cp_l/cp_g;
% 
% k_d   = k_l/k_g;
% alpha_d = k_d/(rho_d*cp_d);
% alpha_l = k_l/(rho_l*cp_l);
% 
% Pr    = (mu_l*cp_l)/k_l;
% St    = (cp_l*dt)/h_gl;
% 
% % u = (sig/mu_l)*(1/10);
% % l = (10*mu_l)/(rho_l*u);
% 
% % R = 4.0e-4;
% % l = (sqrt(27)/2)*St*rho_d*alpha_l*sqrt((rho_l*R)/sig);
% % t = (9/4)*St*rho_d*alpha_l*((rho_l*R)/sig);
% 
% % l = 1e-3;
% % t = 1e-3;
% % u = l/t;
% 
% l   = sqrt(sig/(g*(rho_l - rho_g)));
% lam = (2*pi)*sqrt(3*sig/(g*(rho_l - rho_g)));
% u = sqrt(g*l);
% 
% R_initial = 50e-6/l;
% 
% theta_initial = 352e-6/l;
% 
% bl_thickness = 1.35e-5/l;
% 
% Domain_x = 2.5e-3/l;
% Domain_y = 4.0e-3/l;
% 
% Re    = (rho_l*u*l)/mu_l;
% Fr    = u/sqrt(g*l);
% Bo    = (g*(rho_l-rho_g)*l*l)/sig;
% We    = (rho_l*u*u*l)/sig;
% Pe    = Re*Pr;
% 
% t = l/u;

%% FC-72 FC-72 FC-72

% sig    = 8.4e-3;
% 
% g      = 9.8;
% gz     = 9.8*2.5*(10^-7);
% 
% rho_l  = 1621.2;
% rho_g  = 13.491;
% 
% mu_l   = 4.13e-4;
% mu_g   = 1.19e-5;
% 
% cp_l   = 1106.7;
% cp_g   = 924.81;
% 
% k_l    = 5.4165e-2;
% k_g    = 1.3778e-2;
% 
% Tl   = 43.5;
% Tw   = 57;
% Tsat = 50;
% 
% dt     = Tw-Tl;
% 
% Ts     = (Tsat-Tl)/(Tw-Tl);
% 
% h_gl   = 83562;
% 
% l   = sqrt(sig/(g*(rho_l - rho_g)));
% lam = (2*pi)*sqrt(3*sig/(g*(rho_l - rho_g)));
% u = sqrt(g*l);
% 
% rho_d = rho_l/rho_g;
% mu_d  = mu_l/mu_g;
% cp_d  = cp_l/cp_g;
% 
% k_d     = k_l/k_g;
% alpha_d = k_d/(rho_d*cp_d);
% alpha_l = k_l/(rho_l*cp_l);
% 
% R_initial     = 2.5581e-3/l;
% theta_initial = 7.44e-3/l;
% 
% 
% Re    = (rho_l*u*l)/mu_l;
% Fr    = u/sqrt(g*l);
% Pr    = (mu_l*cp_l)/k_l;
% St    = (cp_l*dt)/h_gl;
% We    = (rho_l*u*u*l)/sig;
% Pe    = Re*Pr;
% gs     = Fr*sqrt(g/gz);
% gs     = 1/(gs^2);
% 
% t = l/u;


%% FC-72 dhir microgravity

% sig    = 9.055e-3;
% sig1    = 8.62e-3;
% 
% g      = 9.8;
% gz     = 9.8*2.5*(10^-7);
% 
% rho_l  = 1636;
% rho_g  = 9.92;
% 
% rho_l1  = 1627;
% rho_g1  = 11.7;
% 
% mu_l   = 4.96e-4;
% mu_g   = 4.96e-4;
% 
% mu_l1   = 4.72e-4;
% mu_g1   = 4.72e-4;
% 
% cp_l   = 1084.0;
% cp_g   = 1084.0;
% 
% cp_l1   = 1092.0;
% cp_g1   = 1092.0;
% 
% k_l    = 5.32e-2;
% k_g    = 5.32e-2;
% 
% k_l1    = 5.26e-2;
% k_g1    = 5.26e-2;
% 
% Tl   = 43.5;
% Tw   = 57;
% Tsat = 44.85;
% Tsat1 = 48.9;
% Rg = 8.314/0.33804;
% A = 1.7e-19;
% 
% dt     = Tw-Tl;
% 
% Ts     = (Tsat-Tl)/dt;
% Ts1    = (Tsat1-Tl)/dt;
% 
% h_gl   = 87096.5;
% h_gl1  = 85687.0;
% 
% l   = sqrt((sig)/(g*(rho_l - rho_g)));
% lam = (2*pi)*sqrt(3*sig/(g*(rho_l - rho_g)));
% u = sqrt(g*l);
% 
% rho_d = rho_l/rho_g;
% mu_d  = mu_l/mu_g;
% cp_d  = cp_l/cp_g;
% 
% k_d     = k_l/k_g;
% alpha_d = k_d/(rho_d*cp_d);
% alpha_l = k_l/(rho_l*cp_l);
% 
% rho_d1 = rho_l1/rho_g1;
% mu_d1  = mu_l1/mu_g1;
% cp_d1  = cp_l1/cp_g1;
% 
% k_d1     = k_l1/k_g1;
% alpha_d1 = k_d1/(rho_d1*cp_d1);
% alpha_l1 = k_l1/(rho_l1*cp_l1);
% 
% R_initial     = 2.5581e-3/l;
% theta_initial = 7.44e-3/l;
% 
% 
% Re    = (rho_l*u*l)/mu_l;
% Fr    = u/sqrt(g*l);
% Pr    = (mu_l*cp_l)/k_l;
% St    = (cp_l*dt)/h_gl;
% We    = (rho_l*u*u*l)/sig;
% Pe    = Re*Pr;
% gs    = Fr*sqrt(g/gz);
% gs    = 1/(gs^2);
% Abar  = A/(mu_l*u*l*l);
% Bbar  = (Tsat*u*u)/(h_gl*dt);
% Cbar  = (k_l*Tsat*sqrt(2*pi*Rg*Tsat))/(2*l*h_gl*h_gl*rho_l);
% 
% Re1   = (rho_l1*u*l)/mu_l1;
% Fr1   = u/sqrt(g*l);
% Pr1   = (mu_l1*cp_l1)/k_l1;
% St1   = (cp_l1*dt)/h_gl1;
% We1   = (rho_l1*u*u*l)/sig1;
% Pe1   = Re1*Pr1;
% gs1   = Fr1*sqrt(g/gz);
% gs1   = 1/(gs1^2);
% 
% t = l/u;
% 
% Domain_x = 0.0190/l;
% Domain_y = 0.0761/l;
% 
% t_d = 0.0088/t;

%% R-11

% sig    = 0.017884;
% 
% g      = 9.8;
% 
% rho_l  = 1479.3;
% rho_g  = 5.8528;
% 
% mu_l = 0.00040362;
% mu_g = 1.0111e-05;
% 
% cp_l   = 0.87954e3;
% cp_g   = 0.53333e3;
% 
% k_l    = 0.087173;
% k_g    = 0.0084020;
% 
% Tb   = 23.0;
% Tsat = 23.7;
% Tw   = 35;
% 
% dt  = Tw-Tb;
% Ts  = (Tsat-Tb)/dt;
% 
% %h_gl   = 144350;
% h_gl   = 181195.4;
% 
% rho_d = rho_l/rho_g;
% mu_d  = mu_l/mu_g;
% cp_d  = cp_l/cp_g;
% CP_d  = (rho_l*cp_l)/(rho_g*cp_g);
% 
% k_d   = k_l/k_g;
% alpha_d = k_d/(rho_d*cp_d);
% alpha_l = k_l/(rho_l*cp_l);
% 
% Pr    = (mu_l*cp_l)/k_l;
% St    = (cp_l*dt)/h_gl;
% 
% % u = (sig/mu_l)*(1/10);
% % l = (10*mu_l)/(rho_l*u);
% 
% % R = 4.0e-4;
% % l = (sqrt(27)/2)*St*rho_d*alpha_l*sqrt((rho_l*R)/sig);
% % t = (9/4)*St*rho_d*alpha_l*((rho_l*R)/sig);
% 
% % l = 1e-3;
% % t = 1e-3;
% % u = l/t;
% 
% l   = sqrt(sig/(g*(rho_l - rho_g)));
% lam = (2*pi)*sqrt(3*sig/(g*(rho_l - rho_g)));
% u = sqrt(g*l);
% 
% % R_initial = 50e-6/l;
% % 
% % theta_initial = 352e-6/l;
% % 
% % bl_thickness = 1.35e-5/l;
% % 
% % Domain_x = 2.5e-3/l;
% % Domain_y = 4.0e-3/l;
% 
% Re    = (rho_l*u*l)/mu_l;
% Fr    = u/sqrt(g*l);
% Bo    = (g*(rho_l-rho_g)*l*l)/sig;
% We    = (rho_l*u*u*l)/sig;
% Pe    = Re*Pr;
% 
% t = l/u;
