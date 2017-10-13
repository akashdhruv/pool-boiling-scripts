function solveMicroLayer_odefun

close all

%%Introduction

% Microlayer solution using MATLAB's inbuilt ode functions

% A. V. Dhruv %

% % Setup
Re    = 1345.1;			% Reynolds number
rho   = 1/1624.4;		% Density ratio rho_v/rho_l
Abar  =  3.0812e-11;		% Dimensionless Hamakar constant
Bbar  =  6.5426e-07;		% Dimensionless number #1
Cbar  =  1.0736e-08;		% Dimensionless number #2

Twall = 1.0;			% Dimensionless wall temperature
Tsat  = 0.0;			% Dimensionless saturation temperature

z0    = zeros(4,1);		% Solution arary

% % Upper boundary conditions
z0(1) = 0.005;
z0(2) = tan(38*pi/180);
z0(3) = Abar/(0.005^3);
z0(4) = 0.0;
r0    = [0.005/tan(38*pi/180),0];

% % Lower boundary conditions
% z0(1) = ((Abar*Bbar)/(Re*(Twall-Tsat)))^(1/3);
% z0(2) = tan(38*pi/180);
% z0(3) = ((Twall-Tsat)*Re)/Bbar;
% z0(4) = 0.0;
% r0    = [0.0,0.005];

% % ODE solution
[r,z] = ode23tb(@microFUN,r0,z0);

% % Post processing
N = length(z(:,1));

q = zeros(N,1); 		% Evaporative heat flux
T = zeros(N,1);			% Interface temperature

for i=1:N

    q(i) = (Twall - Tsat - (Bbar/Re)*z(i,3))/(z(i,1) + (Cbar/rho));
    T(i) = Twall - q(i)*z(i,1);

end

% % Plotting

% Flim thickness
figure
hold on
plot(r,z(:,1),'-k')

% Heat flux
figure
hold on
plot(r,q,'-b')

% Interface temperature
figure
hold on
plot(r,T,'-r')

end


function dzdr = microFUN(r,z)

Re  = 1345.1;			% Reynolds number
Pr  = 1.7386;			% Prandtl number
St  = 0.0116;			% Stefan (Jakob) number
rho = 1/1624.4;			% Density ratio rho_v/rho_l
Pe  = Re*Pr;			% Peclet number
We  = 1.0006;			% Weber number
Fr  = 1.0;			% Froude number

Abar =  3.0812e-11;		% Dimensionless Hamakar number
Bbar =  6.5426e-07;		% Dimensionless number #1
Cbar =  1.0736e-08;		% Dimensionless number #2

Twall = 1.0;			% Wall temperature
Tsat  = 0.0;			% Saturation temperature

dzdr = zeros(4,1);		% Solution vector

dzdr(1) = z(2);
dzdr(2) = ((z(3)-Abar/(z(1)^3))/(Re/We))*((1+z(2)^2)^(3/2));
dzdr(3) = (3*St/Pe)*(z(4)/(z(1)^3));
dzdr(4) = (Tsat - Twall + (Bbar/Re)*z(3))/(z(1) + (Cbar/rho));

end
