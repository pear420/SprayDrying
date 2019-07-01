% Complementary file to SprayDryer_FOPDT.slx

% This case study is based on the article "Tan, L.W., Taip, F.S., & Aziz, 
% N.A. (2009). Simulation and control of spray drying using nozzle atomizer
% spray dryer."

%% 0.

clc
clear all
close all

%% 1. Setup:

% Model Name:
modelName = 'SprayDryer_DIYcontroller_FOPDT';

%% 2. Model variables (base workspace):

% Spray Dryer Transfer function:
K_process = 2;              %[-]    Gain
tau_process = 6.8;          %[min]  Characteristic time
theta_process = 1;          %[min]  Delay

% Step setpoint change: 
timeSetpointChange = 10; %[min]
initialSetpoint = 0;        %[-]
finalSetpoint = 1;          %[-]

% PID controller:
kPs = [1.7 1.9744 2.6839 3.0600];       %[-]    Proportional Gain
tauIs = [6.8 6.9923 5.3084 3.33];       %[min]  Integration time
tauDs = [0 0 0 0];                      %[min]  Derivative time            %[-]    Differential Gain
ControllerSaturationMax = +Inf;                %[-] 
ControllerSaturationMin = -Inf;                %[-]

%% 3. Simulation:

for j=1:length(kPs)
    
    % Change controller settings:
    kP_controller = kPs(j);       %[-]    Proportional Gain
    tauI_controller = tauIs(j);       %[min]  Integration time
    tauD_controller = tauDs(j);
    kI_controller = kP_controller/tauI_controller;    %[-]    Integral Gain
    kD_controller = kP_controller*tauD_controller;
    
    % Run simulation:
    tstop = '60';       %[min]
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
    figure
    plot(tout,Yout,'LineWidth',0.9);
    xlabel('time [min]'); ylabel('[-]'); legend(leg)
    title('Logged Signals'); grid on;
    
    % Evaluate overshoot:
    OS(j) = max((Yout(:,3)-finalSetpoint)/finalSetpoint);
    
    % Legend
    leg1{j} = sprintf('M%d, OS = %.2f %%',j,OS(j)*100);
    
    % Save closed loop signal:
    y_cl(:,j) = Yout(:,3);
    
end

% Plot:
figure
plot(tout,[Yout(:,2),Yout(:,4), y_cl],'LineWidth',0.9);
xlabel('time [min]'); ylabel('[-]'); 
title('PI-controlled unit - reaction curves'); grid on;
legend([{'Open loop','Setpoint'},leg1])


    


