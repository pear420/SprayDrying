% Complementary file to SprayDryer_FOPDT.slx

% This case study is based on the article "Woun Tan, Lee & Ibrahim, Nordin 
% & Kamil, Raja & Taip, Farah. (2011). Empirical modeling for spray drying 
% process of sticky and non-sticky products. Procedia Food Science. 1. 
% 690-697. 10.1016/j.profoo.2011.09.104." 

%% 0.

clc
clear all
close all

%% 1. Setup:

% Model Name:
modelName = 'SprayDryer_FOPDT';

% Read step experiment data from file:
filenameIn = 'CaseStudy1_StepDataIn.txt';
stepDataIn = dlmread(filenameIn);
filenameOut = 'CaseStudy1_StepDataOut.txt';
stepDataOut = dlmread(filenameOut);

% Find time at which the step change occurs:
[~,i_step] = max(diff(stepDataIn(:,2)));
t_step = stepDataIn(i_step,1);

%% 2. Model variables (base workspace):

% Spray Dryer Transfer function:
K_process = 0.86;              %[-]    Gain
tau_process = 2.44;            %[min]  Characteristic time
theta_process = 0.16;          %[min]  Delay

% Step setpoint change: 
timeSetpointChange = t_step; %[min]
initialSetpoint = 0;         %[-]
finalSetpoint = 1;           %[-]

% PID controller:
kP_controller = 2.68;       %[-]    Proportional Gain
tauI_controller = 5.3084;   %[min]  Integration time
tauD_controller = 0;        %[min]  Derivative time        
kI_controller = kP_controller/tauI_controller;    %[-]    Integral Gain
kD_controller = kP_controller*tauD_controller;    %[-]    Differential Gain

% Temperatures from the sp change experiment:
Tin0 = 165;                 %[�C]   Inlet air temperature before sp change
Tin1 = 175;                 %[�C]   Inlet air temperature after sp change
Tout0 = 86;                 %[�C]   Outlet air temperature before sp change

%% 3. Simulation:

% Run simulation:
tstop = '20';       %[min]
tstep = '0.01';     %[min]
SolverType = 'Fixed-step'; %[-]
simOut = sim(modelName,'SolverType',SolverType,...
    'FixedStep',tstep,'StopTime',tstop,'SignalLogging',...
    'on','SignalLoggingName','logsout');

%% 4. Post-Processing:

% We're interested only in the setpoint change and the open loop response, 
% which is at position 2 of the output dataset:
tout = simOut.tout;         %[min]   Time Span
y_sp_out = simOut.logsout{4}.Values.Data;      %[-]    Data
y_ol_out = simOut.logsout{2}.Values.Data;      %[-]    Data

% Signals:
figure

subplot(2,1,1)
plot(tout,y_sp_out,'LineWidth',0.9); hold on
xlabel('time [min]'); ylabel('[-]');
title('Setpoint signal'); grid on;

subplot(2,1,2)
plot(tout,y_ol_out,'LineWidth',0.9); hold on
xlabel('time [min]'); ylabel('[-]');
title('Open loop signal'); grid on;

% Temperatures:
Tin = y_sp_out*(Tin1-Tin0)+Tin0;               %[�C] Inlet air temperature
Tout = y_ol_out*(Tin1-Tin0)+Tout0;           %[�C] Outlet air temperature

figure

subplot(2,1,1)
plot(tout,Tin,'LineWidth',0.9); hold on
plot(stepDataIn(:,1),stepDataIn(:,2),'LineWidth',0.6,'Marker','+','LineStyle','--'); hold off
xlabel('time [min]'); ylabel('[-]');
title('Inlet temperature'); grid on;
legend('Model','Experimental')

subplot(2,1,2)
plot(tout,Tout,'LineWidth',0.9); hold on
plot(stepDataOut(:,1),stepDataOut(:,2),'LineWidth',0.6,'Marker','+','LineStyle','--'); hold off
xlabel('time [min]'); ylabel('[-]');
title('Outlet temperature'); grid on;
legend('Model','Experimental')