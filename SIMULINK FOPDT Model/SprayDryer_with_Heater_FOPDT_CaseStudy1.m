% Complementary file to SprayDryer_FOPDT.slx

% This case study is based on the article "Tan, L.W., Taip, F.S., & Aziz, 
% N.A. (2009). Simulation and control of spray drying using nozzle atomizer
% spray dryer." with the addition of an electrical heater in series.

%% 0.

clc
clear all
close all

%% 1. Setup:

% Model Name:
modelName = 'SprayDryer_with_Heater_FOPDT';

%% 2. Model variables (base workspace):

% Spray Dryer Transfer function:
K_process = 2;              %[-]    Gain
tau_process = 6.8;          %[min]  Characteristic time
theta_process = 1;          %[min]  Delay

% Electrical heater Transfer function:
K_heater = 5.8;             %[-]    Gain
tau_heater = 1;             %[min]  Characteristic time
theta_heater = 7.3/60;         %[min]  Delay

% Step setpoint change: 
timeSetpointChange = 10; %[min]
initialSetpoint = 0;        %[-]
finalSetpoint = 1;          %[-]

% PID controller:
kP_controller = 1E-2;       %[-]    Proportional Gain
tauI_controller = 0.5E-2;   %[min]  Integration time
tauD_controller = 0;        %[min]  Derivative time        
kI_controller = kP_controller/tauI_controller;    %[-]    Integral Gain
kD_controller = kP_controller*tauD_controller;    %[-]    Differential Gain

%% 3. Simulation:

% Run simulation:
tstop = '360';      %[min]
tstep = '0.01';     %[min]
SolverType = 'Fixed-step'; %[-]
simOut = sim(modelName,'SolverType',SolverType,...
    'FixedStep',tstep,'StopTime',tstop,'SignalLogging',...
    'on','SignalLoggingName','logsout');

%% 4. Post-Processing:

% Recover logged data:
tout = simOut.tout;         %[min]   Time Span

for i=1:simOut.logsout.numElements
    Yout(:,i) = simOut.logsout{i}.Values.Data; %[-]    Data
    leg{i} = simOut.logsout{i}.Values.Name;    %[-]    Name
end

% Plot:
plot(tout,Yout,'LineWidth',0.9); 
xlabel('time [min]'); ylabel('[-]'); legend(leg)
title('Logged Signals'); grid on;


