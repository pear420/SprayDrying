% Spray drying steady state problem:

clc
clear all
close all

%% Data:

global MWa R muG D_eff MWw cpG k Dp Dsd P ETA m_dry DHev cpL rhoP g
global AA BB CC

% Milk:
m_milk = 1750/3600; %[kg/s] mass flow rate
TP0 = 303; %[K] inlet temperature
w_fat = 4.76/100; %[-] fat weight fraction
rhoP = 1000; %[kg/m^3] particle
cpL = 4184; %[J/kg] specific heat
DHev = 540*4.184; %[J/kg] water latent enthalpy of evaporation

% Air:
m_dry = 72000/3600; %[kg/s] mass flow rate
TG0 = 403; %[K]inlet temperature
P = 1*1.01325E5; %[Pa] pressure
muG = 2.3E-5; %[kg/m/s] viscosity
cpG = 0.25*4.184; %[J/kg] 
D_eff = 1.8E-5; %[m^2/s] Effective diffusivity
k = 8.10E-6; % Conductivity
MWa = 28.9647*1e-3; %[kg/mol] molar weight 

% Water:
MWw = 28*1e-3; %[kg/mol] molar weight

% Antoine coefficents for water [mmHg],[K]:
AA = 18.3036;
BB = 3816.44;
CC = -46.13;

% Design Specifications: 
Dsd = 5.5; %[m] Diameter of the spray drier
W_target = 0.0005; %[kg_water/kg_fat] target particle moisture

% Universal gas constant:
R = 8.314; %[J/mol/K]

% Gravity acceleration:
g = 9.81; %[m/s^2]

%% 1. Integration:
disp('1.')

% Initial (and constant) particle size:
Dp = 0.2E-3; %[m]
mP0 = rhoP*pi*Dp^3/6; %[kg]

% Number concentration of particles:
ETA = m_milk/mP0; %[-/s]

% Initial particle velocity:
vP0 = 0.3; %[m/s]

% Initial gas velocity:
mG0 = m_dry;  %[kg/s]
rhoG = P*MWa/R/TG0; %[kg/m^3]
vG0 = mG0/rhoG*4/pi/Dsd^2; %[m/s]

% Initial drag velocity:
vS0 = vP0-vG0; %[m/s]

% Initial coordinate along the spray drier:
z0 = 0; %[m]

% Packing:
Y0 = [mP0 mG0 TP0 TG0 vS0 z0];

% Desired final mass of the particle:
mP_fat = w_fat*mP0; %[kg] fat(dry) mass in the particle (constant)
mP_target = mP_fat*(1+W_target); %[kg]

% Options:
options = odeset('Event',@(t,Y)EventFunction(t,Y,mP_target)); 

% Integration:
tspan = [0,20]; %[s]
[t,Y,te,Ye] = ode23s(@ODE,tspan,Y0,options)

% Unpacking:
mP = Y(:,1); %[kg]
mG = Y(:,2); %[kg/s]
TP = Y(:,3); %[K]
TG = Y(:,4); %[K]
vS = Y(:,5); %[m/s]
z = Y(:,6); %[m]

% Other profiles:
for i = 1:length(t)
  [~,vG(i,:),vP(i,:)] = ODE(t,Y(i,:));
end

%% 2. Post-processing:
disp('2.')

% Velocities:
figure
plot(t,[vG, vP, vS],'Linewidth',1)
title('Velocity profiles')
grid on;
xlabel('t[s]')
ylabel('v[m/s]')
legend('v_G','v_P','v_S')

% Temperatures:
figure
plot(t,[TG, TP],'Linewidth',1)
title('Temperature profiles')
grid on;
xlabel('t[s]')
ylabel('T[K]')
legend('T_G','T_P')

% Particle Mass:
figure
p = plot(t,[mP, mP_target*ones(1,length(mP)).'],'Linewidth',1)
p(2).LineStyle = '--';
title('Particle mass profile')
grid on;
xlabel('t[s]')
ylabel('mP[kg]')
legend('mP','target mP')

% Distance:
figure
plot(t,z,'Linewidth',1)
title('Distance')
grid on;
xlabel('t[s]')
ylabel('z[m]')

%% Functions:

function [dY_dt,vG,vP] = ODE(t,Y)
% Governing equation for a plug-flow spray dryer:

global MWa R muG D_eff MWw cpG k Dp Dsd P ETA m_dry DHev cpL rhoP g

% Unpacking:
mP = Y(1); %[kg]
mG = Y(2); %[kg/s]
TP = Y(3); %[K]
TG = Y(4); %[K]
vS = Y(5); %[m/s]
z = Y(6); %[m]

% Section area of the spray drier:
A = pi*Dsd^2/4; %[m^2]

% Gas density:
rhoG = P*MWa/R/TG; %[kg/m^3]

% Gas velocity:
vG = mG/rhoG/A; %[m/s]

% Particle velocity:
vP = vG+vS; %[m/s]

% Transport properties:
Re = rhoG*abs(vS)*Dp/muG;
Pr = muG*cpG/k;
Sc = muG/rhoG/D_eff;
Nu = 2 + 0.4*Re^0.5*Pr^(1/3);
Sh = 2 + 0.4*Re^0.5*Sc^(1/3);
Kp = Sh*D_eff/Dp*MWw/R/TG;
h = Nu*k/Dp;

% Vapor pressure:
Pev = Antoine(TG)*1.01325e5; %[Pa]

% Surface of the particle:
Sp = pi*Dp^2/4;

% Volume of the particle:
Vp = pi*Dp^3/8;

% Water Pressure in the gas phase:
y = (mG-m_dry)/MWw/mG*MWa;
Pw = P*y;

% Friction factor in viscous regime:
f = 24/Re;

% Derivatives:
dY_dt(1) = Kp*Sp*(Pw-Pev);
dY_dt(2) = Kp*Sp*(Pev-Pw)*A*ETA*vP;
dY_dt(3) = (h*Sp*(TG-TP)+dY_dt(1)*DHev)/mP/cpL;
dY_dt(4) = h*Sp*(TP-TG)*A*ETA*vP/mG/cpG;
dY_dt(5) = ((rhoP-rhoG)*Vp*g-0.5*f*rhoG*vS*abs(vS)*pi*Dp^2/4-vS*dY_dt(1))/mP;
dY_dt(6) = vP; 

% Transpose:
dY_dt = dY_dt.';

end

function [value,isterminal,direction] = EventFunction(t,Y,mP_target)
% Event function: stop integration upon reaching the target particle 
% mass (determined from the target moisture):

value = Y(1)-mP_target;
isterminal = 1;
direction = 0;

end

function Pev = Antoine(T)
% Vapour pressure of water:

global AA BB CC

Pev = exp(AA-BB/(T+CC))/760; %[Pa]

end



