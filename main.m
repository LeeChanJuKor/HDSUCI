% Initialize workspace
clear all; close all; clc;  
addpath('./src');  

%% Define user-controlled parameters
num_substeps = 5;  % Number of substeps (Choose between 2, 3, 4, or 5)
rho_inf = 0.83734; % Numerical dissipation parameter


%%
% Define system properties
M = 1;   
C = 10;  
K = 100; 

% Define simulation parameters
del_t = 0.1;  
last_t = 10;   
t = 0:del_t:last_t;  

% Compute integration parameters
HDSUCI_params = HDSUCI_time_integration_params(num_substeps, rho_inf);

% Compute numerical solution using HDSUCI time integration
[HDSUCI_u, HDSUCI_v] = HDSUCI_TI(del_t, t, M, C, K, num_substeps, HDSUCI_params);
[HDSUCI_a] = HDSUCI_acc(HDSUCI_u, del_t, t);

% Plot results
figure;
plot(t, HDSUCI_u(1,:), 'r', 'LineWidth', 1.5);  
legend('HDSUCI');  
grid on;

figure;
plot(t, HDSUCI_v(1,:), 'r', 'LineWidth', 1.5);  
legend('HDSUCI');  
grid on;

figure;
plot(t, HDSUCI_a(1,:), 'r', 'LineWidth', 1.5);  
legend('HDSUCI');  
grid on;
