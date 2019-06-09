% Dynamic model of a spray dryer:
% 
% ref. C. Palencia, J. Nava, E. Herman, G. C. Rodriguez-Jimenes, 
%      and M. A. Garc�a-Alvarado* - SPRAY DRYING DYNAMIC MODELING 
%      WITH A MECHANISTIC MODEL - DRYING TECHNOLOGY, 20(3), 569�586 (2002)

%% 0.

clc
clear all
close all

%% 1.

global G_beta G_gamma

% Steady-state variables:
% beta - Dry matter in the dispersed phase;
% gamma- Dry matter in the continuous phase;
G_beta = 2.5E-4; %[kg/s] mass velocity
G_gamma = 1.872E-1; %[kg/s] mass velocity
X_beta0 = 4; %[-] Initial water content
X_gamma0 = 1.960E-2; %[-] Initial water content
T_beta0 = 2.7E1+273.15; %[K] Initial temperature
T_gamma0 = 1.23E2+273.15; %[K] Initial temperature

%% 2.

global  Cp_beta Cp_gamma Cp_wv Cp_w D_w_beta D_w_gamma H_wv DHev_w rho_beta rho_gamma epsi R V MW_air MW_w

% Properties used for model solution:
% beta - Dry matter in the dispersed phase;
% gamma- Dry matter in the continuous phase;
Cp_beta = 1657; %[J/kg/K] specific heat 
Cp_gamma = 1000; %[J/kg/K] specific heat
Cp_wv = 1800; %[J/kg/K] specific heat of the water vapour
Cp_w = 4185; %[J/kg/K] specific heat of water
D_w_beta = 6.67E-10; %[m^2/s] mass diffusivity
D_w_gamma =1.8E-5; %[m^2/s] mass diffusivity
H_wv = 2501E3; %[J/kg] enthalpy of the water vapour
DHev_w = 2260.44E3; %[J/kg] enthalpy of evaporation of water
rho_beta = 800; %[kg/m^3] volumetric concentration of dry matter
rho_gamma = 1; %[kg/m^3] volumetric concentration of dry matter
epsi = 0.98; %[-] volume fraction of continuous phase in the chamber
R = 5E-5; %[m] radius of the particles
V = 1.26; %[m^3]
MW_air = 28.84E-3; %[kg/mol]
MW_w = 18E-3; %[kg/mol]

%% 3.

% Integration:
tspan = [0 1E5]; %[s]
Y0 = [X_beta0 X_gamma0 T_beta0 T_gamma0];

% Drying stage:
[ts,Ys] = ode23s(@(t,Y)DryingStage(t,Y,Y0),tspan,Y0)

%% *. Functions:

function dY_dt = DryingStage(t,Y,Y_INLET)
% Spray drying stage:
% beta - Dry matter in the dispersed phase;
% gamma- Dry matter in the continuous phase;

global G_beta G_gamma Cp_beta Cp_gamma Cp_wv DHev_w Cp_w D_w_beta D_w_gamma H_wv rho_beta rho_gamma epsi R V MW_air MW_w

% Unpacking:
X_beta = Y(1); %[-] 
X_gamma = Y(2); %[-]
T_beta = Y(3); %[K]
T_gamma = Y(4); %[K]
X_beta_INLET = Y_INLET(1); %[-] 
X_gamma_INLET = Y_INLET(2); %[-]
T_beta_INLET = Y_INLET(3); %[K]
T_gamma_INLET = Y_INLET(4); %[K]

% Properties:
k_beta = 0.1418+0.00493*X_beta/(1+X_beta); %[W/m/s]
k_gamma = 8.4044E-5*(T_gamma)+4.63E-5; %[W/m/s]
a = 3*(1-epsi)/R; %[1/m]
Sh = 2; %[-]
Nu = 2; %[-]
kc_beta = 3*pi^2*D_w_beta/R; %[W/m/s]
kc_gamma = Sh*D_w_gamma/2/R; %[W/m/s]
h_beta = 3*pi^2*k_beta/R; %[W/m^2/K]
h_gamma = Nu*k_gamma/2/R; %[W/m^2/K]
a_w = 1-exp(-1.71*T_beta^0.3*X_beta^1.013); %[-]
% k1 = exp(12.872);
% k2 = -2.0582;
% k3 = -10.380;
% k4 = 0.07509;
% k5 = -0.0001279;
% a_w = 1-exp(-k1*T_gamma^k2*X_gamma^(k3+k4*T_gamma+k5*T_gamma^2));
A = MW_air/MW_w*X_gamma; %[-]
B = A/(1+A); %[-]
X_gamma_INTERPHASE = MW_w/MW_air*a_w*B/(1-a_w*B); %[-]
X_beta_INTERPHASE = X_beta+kc_gamma*rho_gamma/kc_beta/rho_beta*(X_gamma-X_gamma_INTERPHASE);
T_INTERPHASE = (-rho_beta*kc_beta*(X_beta-X_beta_INTERPHASE)*DHev_w+h_gamma*T_gamma+h_beta*T_beta)/(h_gamma+h_beta);

% Terms:
N_w_beta = rho_beta*kc_beta*(X_beta-X_beta_INTERPHASE);
N_w_gamma = rho_gamma*kc_gamma*(X_gamma_INTERPHASE-X_gamma);
q_beta = h_beta*(T_INTERPHASE-T_beta);
q_gamma = h_gamma*(T_gamma-T_INTERPHASE);
H_beta = (Cp_beta+Cp_w*X_beta)*T_beta; %[J/kg]
H_gamma = Cp_gamma*T_gamma+(H_wv+Cp_wv*T_gamma)*X_gamma; %[J/kg]
H_beta_INLET = (Cp_beta+Cp_w*X_beta_INLET)*T_beta_INLET; %[J/kg]
H_gamma_INLET = Cp_gamma*T_gamma_INLET+(H_wv+Cp_wv*T_gamma_INLET)*X_gamma_INLET; %[J/kg]

% Corrections:
N_w_beta = N_w_gamma;

% Derivatives:
dX_beta_dt = (-N_w_beta*a*V-G_beta*(X_beta-X_beta_INLET))/(rho_beta*(1-epsi)*V);
dX_gamma_dt = (+N_w_gamma*a*V-G_gamma*(X_gamma-X_gamma_INLET))/(rho_gamma*epsi*V);
dT_beta_dt = (+q_beta*a*V-G_beta*(H_beta-H_beta_INLET))/(rho_beta*(1-epsi)*V*(Cp_beta+X_beta*Cp_w));
dT_gamma_dt = (-q_gamma*a*V+N_w_beta*a*V*DHev_w-G_gamma*(H_gamma-H_gamma_INLET))/(rho_gamma*epsi*V*(Cp_gamma+Cp_w*X_gamma));

% Output:
dY_dt = [dX_beta_dt; dX_gamma_dt; dT_beta_dt; dT_gamma_dt]; 

end
