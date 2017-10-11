function solveMicroLayer

close all

Re    = 1345.1;
rho   = 1/1624.4;
Abar  =  3.0812e-11;
Bbar  =  6.5426e-07;
Cbar  =  1.0736e-08;

Twall = 1.0;
Tsat  = 0.0;

z0    = zeros(4,1);

z0(1) = 0.005;
z0(2) = tan(38*pi/180);
z0(3) = Abar/(0.005^3);
z0(4) = 0.0;
r0    = [0.005/tan(38*pi/180),0];

% z0(1) = ((Abar*Bbar)/(Re*(Twall-Tsat)))^(1/3);
% z0(2) = tan(38*pi/180);
% z0(3) = ((Twall-Tsat)*Re)/Bbar;
% z0(4) = 0.0;
% r0    = [0.0,0.005];

[r,z] = ode23tb(@microFUN,r0,z0);

N = length(z(:,1));

q = zeros(N,1);
T = zeros(N,1);

for i=1:N

    q(i) = (Twall - Tsat - (Bbar/Re)*z(i,3))/(z(i,1) + (Cbar/rho));
    T(i) = Twall - q(i)*z(i,1);

end

figure
hold on
plot(r,z(:,1),'-k')

figure
hold on
plot(r,q,'-b')

figure
hold on
plot(r,T,'-r')

end


function dzdr = microFUN(r,z)

Re  = 1345.1;
Pr  = 1.7386;
St  = 0.0116;
rho = 1/1624.4;
Pe  = Re*Pr;
We  = 1.0006;
Fr  = 1.0;

Abar =  3.0812e-11;
Bbar =  6.5426e-07;
Cbar =  1.0736e-08;

Twall = 1.0;
Tsat  = 0.0;

dzdr = zeros(4,1);

dzdr(1) = z(2);
dzdr(2) = ((z(3)-Abar/(z(1)^3))/(Re/We))*((1+z(2)^2)^(3/2));
dzdr(3) = (3*St/Pe)*(z(4)/(z(1)^3));
dzdr(4) = (Tsat - Twall + (Bbar/Re)*z(3))/(z(1) + (Cbar/rho));

end