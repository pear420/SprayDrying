% Dynamic model of a spray dryer:
% 
% ref. C. Palencia, J. Nava, E. Herman, G. C. Rodriguez-Jimenes, 
%      and M. A. Garcìa-Alvarado* - SPRAY DRYING DYNAMIC MODELING 
%      WITH A MECHANISTIC MODEL - DRYING TECHNOLOGY, 20(3), 569–586 (2002)

%% 1.

% Steady-state variables:
% beta - Dry matter in the dispersed phase;
% gamma- Dry matter in the continuous phase;
G_beta = 2.5E-4; %[kg/s] mass velocity
G_gamma = 1.872E-1; %[kg/s] mass velocity
X_beta0 = 4; %[-] Initial water content
X_gamma0 = 1.960E-2; %[-] Initial water content
T_beta0 = 2.7E1+273.15; %[K] Initial temperature
T_gamma0 = 1.6E2; %[K] Initial temperature

%% 2.

% Properties used for model solution:
% beta - Dry matter in the dispersed phase;
% gamma- Dry matter in the continuous phase;
Cp_beta = 1657; %[J/kg/K] specific heat 
Cp_gamma = 1000; %[J/kg/K] specific heat
Cp_vw = 1800; %[J/kg/K] specific heat of the water vapour
Cp_w = 4185; %[J/kg/K] specific heat of water
Dw_beta = 6.67E-10; %[m^2/s] mass diffusivity
H_wv = 2501E3; %[J/kg] enthalpy
rho_beta = 800; %[kg/m^3] volumetric concentration of dry matter
epsi = 0.98; %[-] volume fraction of continuous phase in the chamber


%% *. Functions:

function dY_dt = DryingStage(t,Y)
% Spray drying stage:
% beta - Dry matter in the dispersed phase;
% gamma- Dry matter in the continuous phase;

% Unpacking:
X_beta = Y(1); %[-] 
X_gamma = Y(2); %[-]
T_beta = Y(3); %[K]
T_gamma = Y(4); %[K]

% Derivatives:
dX_beta_dt = -kc_beta*a*(X_beta-X_beta_INTERPHASE)/(1-epsi)-G_beta*(X_beta-X_beta_INLET)/rho_beta/(1-epsi)/V;
dX_gamma_dt = kc_gamma*a*(X_gamma_INTERPHASE-X_gamma)/epsi-G_gamma*(X_gamma-X_gamma_INLET)/rho_gamma/epsi/V;
dT_beta_dt = h_beta*a*(T_gamma_INTERPHASE-T_beta)/epsi-Cp_w*T_beta/(Cp_beta+Cp_w*X_beta)*dX_beta_dt-G_beta*((Cp_beta+Cp_w*X_beta)*T_beta-(Cp_beta+Cp_w*X_beta_INLET)*T_INLET)/(rho_beta*(1-epsi)*V*(Cp_beta+Cp_w*X_beta));
dT_gamma_dt = 


end
