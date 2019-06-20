% Dynamic model of a spray dryer:
% 
% ref. C. Palencia, J. Nava, E. Herman, G. C. Rodriguez-Jimenes, 
%      and M. A. Garcìa-Alvarado* - SPRAY DRYING DYNAMIC MODELING 
%      WITH A MECHANISTIC MODEL - DRYING TECHNOLOGY, 20(3), 569–586 (2002)

%% 0.

clc
clear all
close all

%% 1. Data:

global Rp D_w_beta D_w_gamma M_air M_w rho_gamma rho_beta Lambda
global T_amb h_out A_out T_ref H_wv_ref Cp_beta Cp_gamma Cp_w Cp_wv epsi V G_beta G_gamma

% From paper:
Cp_gamma = 1000;        %[J/kg/K]
Cp_wv = 1800;          %[J/kg/K]
Cp_beta = 1657;         %[J/kg/K]
Cp_w = 4185;            %[J/kg/K]
D_w_beta = 6.67E-10;    %[m^2/s]
H_wv_ref = 2501E3;      %[J/kg]
rho_beta = 800;         %[kg/m^3]
epsi = 0.98;            %[-]
Rp = 5E-5;              %[m]
V = 1.26;               %[m^3]
A_out = 5.06;           %[m^2]
h_out = 19.1;           %[W/m^2/K]

% Missing:
Lambda = 2260.44E3;     %[J/kg]
M_air = 28.84E-3;       %[kg/mol]
M_w = 18E-3;            %[kg/mol]
rho_gamma = 1;          %[kg/m^3]
D_w_gamma = 1.8E-5;     %[m^2/s]
T_ref = 273.15;         %[K]
T_amb = 298.15;         %[K]

%% 2. Inlet/Initial values:

% Constant mass flow rates:
G_beta = 2.5E-4;        %[kg/s]
G_gamma = 1.872E-1;     %[kg/s]

% Initial values:
X_beta0 = 4;            %[-]
X_gamma0 = 1.96E-2;     %[-]
T_beta0 = 27+273.15;    %[K]
T_gamma0 = 160+273.15;  %[K]

% % Inlet values:
% X_beta_IN = 4.0;      %[-]
% X_gamma_IN = 1.96E-2;    %[-]
% T_beta_IN = 27+273.15;   %[K]
% T_gamma_IN = 160+273.15; %[K]

%% 3. Integration:

% Time span:
tspan = [0, 3600]; %[s]

% Initial conditions - Packing:
Y0 = [X_beta0 X_gamma0 T_beta0 T_gamma0];
YIN = Y0;

[Nw_beta,Nw_gamma,q_beta,q_gamma] = Fluxes(Y0)
[dY_dt] = ODE(0,Y0,YIN)

% Numerical integration:
options = odeset('Abstol',1E-6,'Reltol',1E-6);
[ts,Ys] = ode23s(@(t,Y)ODE(t,Y,YIN),tspan,Y0,options);

%% 4. Post processing:

% Moisture content:
figure
plot(ts,Ys(:,1:2)); grid on;

% Temperatures:
figure
plot(ts,Ys(:,3:4)-273.15); grid on;

%% *. Functions:

function [dY_dt] = ODE(t,Y,Y_IN)
% System of governing equations for a fully mixed spray dryer:

global T_amb h_out A_out T_ref H_wv_ref Cp_beta Cp_gamma Cp_w Cp_wv ...
    epsi Rp V G_beta G_gamma rho_beta rho_gamma

% Unpacking:
X_beta = Y(1);              %[-]
X_gamma = Y(2);             %[-]
T_beta = Y(3);              %[K]
T_gamma = Y(4);             %[K]

% Unpacking - INLET:
X_beta_IN = Y_IN(1);        %[-]
X_gamma_IN = Y_IN(2);       %[-]
T_beta_IN = Y_IN(3);        %[K]
T_gamma_IN = Y_IN(4);       %[K]

% Material and heat fluxes:
[Nw_beta,Nw_gamma,q_beta,q_gamma] = Fluxes(Y);  
q_out = h_out*(T_gamma-T_amb);                  %[W/m^2/K]
Nw_beta = Nw_gamma;

% Mass-specific enthalpies:
H_beta = (Cp_beta+Cp_w*X_beta)*(T_beta-T_ref); %[J/kg]
H_gamma = (Cp_gamma+X_gamma*Cp_wv)*(T_gamma-T_ref)+H_wv_ref*X_gamma; %[J/kg]
H_gamma = Cp_gamma*(T_gamma-273.15)+(H_wv_ref+Cp_wv*(T_gamma-273.15))*X_gamma; %[J/kg]

H_beta_IN = (Cp_beta+Cp_w*X_beta_IN)*(T_beta_IN-T_ref); %[J/kg]
H_gamma_IN = (Cp_gamma+X_gamma_IN*Cp_wv)*(T_gamma_IN-T_ref)+H_wv_ref*X_gamma_IN; %[J/kg]

% Contact surface for unit of volume:
a = 3*(1-epsi)/Rp; %[1/m]

% Diffusive and convective mass and enthalpic flow rates:
mD_beta =  -Nw_beta*V*a;                    %[kg/s]
mC_beta =  -G_beta*(X_beta-X_beta_IN);      %[kg/s]
mD_gamma = Nw_gamma*V*a;                    %[kg/s]
mC_gamma = -G_gamma*(X_gamma-X_gamma_IN);   %[kg/s]
HD_beta = q_beta*a*V;                       %[W]
HC_beta = -G_beta*(H_beta-H_beta_IN);       %[W]
HD_gamma = -q_gamma*a*V-q_out*A_out;        %[W]
HC_gamma = -G_gamma*(H_gamma-H_gamma_IN);    %[W]

% Derivatives:
dX_beta_dt = (mD_beta+mC_beta)/rho_beta/V/(1-epsi);          %[1/s]
dX_gamma_dt = (mD_gamma+mC_gamma)/rho_gamma/V/epsi;          %[1/s]
dT_beta_dt = (HD_beta+HC_beta)/(Cp_beta+Cp_w*X_beta)/rho_beta/V/(1-epsi);  %[K/s]
dT_gamma_dt = (HD_gamma+HC_gamma)/(Cp_gamma+Cp_wv*X_gamma)/rho_gamma/V/epsi; %[K/s]

% Packing:
dY_dt = [dX_beta_dt dX_gamma_dt dT_beta_dt dT_gamma_dt].';

end

function [Nw_beta,Nw_gamma,q_beta,q_gamma] = Fluxes(Y)
% Material and heat fluxes:

global Rp D_w_beta D_w_gamma M_air M_w rho_gamma rho_beta Lambda

% Unpacking:
X_beta = Y(1);              %[-]
X_gamma = Y(2);             %[-]
T_beta = Y(3);              %[K]
T_gamma = Y(4);             %[K]

% Transport properties - material:
Nu = 2;
Sh = 2;
kc_beta = 3*pi^2*D_w_beta/Rp; %[m/s]
kc_gamma = Sh*D_w_gamma/2/Rp; %[m/s]

% Interphase - material:
a_w = 1-exp(-1.71*T_beta^0.3*X_beta^1.013); %[-]
A = M_air/M_w*X_gamma;                      %[-]
B = A/(1+A);                                %[-]
X_gamma_I = M_w/M_air*a_w*B/(1-a_w*B);      %[-]
K = rho_gamma*kc_gamma/rho_beta/kc_beta;    %[-]
X_beta_I = K*(X_gamma_I-X_gamma)+X_beta;           %[-]

% Fluxes - material:
Nw_beta = rho_beta*kc_beta*(X_beta-X_beta_I);       %[kg/m^2/s]
Nw_gamma = rho_gamma*kc_gamma*(X_gamma-X_gamma_I);  %[kg/m^2/s]

% Transport properties - heat:
k_beta = 0.1418+0.00493*X_beta/(1+X_beta);  %[W/m/s]
k_gamma = 8.4044E-5*T_gamma+4.63E-5;        %[W/m/s]
h_beta = 3*pi^2*k_beta/Rp;                  %[W/m^2/s]
h_gamma = k_gamma*Nu/2/Rp;                  %[W/m^2/s]

% Interphase - heat:
T_I = (h_beta*T_beta+h_gamma*T_gamma-Nw_beta*Lambda)/(h_gamma+h_beta); %[K]

% Fluxes - heat:
q_beta = h_beta*(T_I-T_beta);               %[W/m^2/s]
q_gamma = h_gamma*(T_gamma-T_I);            %[W/m^2/s]

end
