function [u, v] = HDSUCI_TI(del_t, t, M, C, K, num_substeps, HDSUCI_params)
% HDSUCI_TI: s-substep implicit time integration method
%
%
% Inputs:
%   del_t       - Time step size
%   t           - Time vector
%   M           - Mass matrix
%   C           - Damping matrix
%   K           - Stiffness matrix
%   rho_inf     - Numerical damping parameter
%   num_substeps - Number of sub-steps (2, 3, 4, or 5)
%   HDSUCI_params - Struct containing integration parameters
%
% Outputs:
%   u - Displacement vector
%   v - Velocity vector
%   a - Acceleration vector

% Number of time steps
num_time_steps = length(t) - 1;

% Initialize displacement, velocity, and acceleration vectors
u = zeros(1, length(t));
v = zeros(1, length(t));

% Extract the effective matrix coefficients dynamically based on num_substeps
Effective_matrix_coeffs = zeros(1, num_substeps);
for i = 1:num_substeps
    Effective_matrix_coeffs(i) = HDSUCI_params.(['a' num2str(i) num2str(i)]);
end

% Compute Effective Stiffness Matrices and LU decomposition dynamically
Eff_K = cell(1, num_substeps);
for i = 1:num_substeps
    Eff_K{i} = C/(del_t*Effective_matrix_coeffs(i)) + K + M/(del_t^2 * Effective_matrix_coeffs(i)^2);
end

% Time integration loop
for i = 1:num_time_steps
    for j = 1:num_substeps
        % First sub-step
        if j == 1
            g1 = HDSUCI_params.g1;
            a11 = HDSUCI_params.a11;
            F1 = sin(2 * (t(i) + g1 * del_t)); % External force
            Eff_R1 = F1 - (-C*u(:,i)/(a11*del_t) - M*(v(:,i)*a11*del_t + u(:,i))/(a11^2*del_t^2));
            g1_u = Eff_K{j} \ Eff_R1; % Solve for displacement
            g1_v = -(u(:,i) - g1_u)/(a11*del_t); % Compute velocity
            g1_a = -(v(:,i)*a11*del_t + u(:,i) - g1_u)/(a11^2*del_t^2); % Compute acceleration

        elseif j == 2
            % Second sub-step
            g2 = HDSUCI_params.g2;
            a21 = HDSUCI_params.a21;
            a22 = HDSUCI_params.a22;
            b1 = HDSUCI_params.b1;
            b2 = HDSUCI_params.b2;
            F2 = sin(2 * (t(i) + g2 * del_t)); % External force
            Eff_R2 = F2 - (-C*(del_t*a21*g1_v + u(:,i))/(del_t*a22) ...
                - M*(del_t^2*a21*g1_a*a22 + v(:,i)*del_t*a22 + del_t*a21*g1_v + u(:,i))/(del_t^2*a22^2));
            g2_u = Eff_K{j} \ Eff_R2;
            g2_v = -(del_t*a21*g1_v + u(:,i) - g2_u)/(del_t*a22);
            g2_a = -(del_t^2*a21*g1_a*a22 + v(:,i)*del_t*a22 + del_t*a21*g1_v + u(:,i) - g2_u)/(del_t^2*a22^2);
            
            % Update final displacement and velocity if last sub-step
            if j == num_substeps
                u(:,i+1) = u(:,i) + del_t * (b1 * g1_v + b2 * g2_v);
                v(:,i+1) = v(:,i) + del_t * (b1 * g1_a + b2 * g2_a);
                continue;
            end

        elseif j == 3
            % Third sub-step
            g3 = HDSUCI_params.g3;
            a31 = HDSUCI_params.a31;
            a32 = HDSUCI_params.a32;
            a33 = HDSUCI_params.a33;
            b3 = HDSUCI_params.b3;
            F3 = sin(2 * (t(i) + g3 * del_t)); % External force
            Eff_R3 = F3 - (-C*(a31*del_t*g1_v + a32*del_t*g2_v + u(:,i))/(del_t*a33) ...
                - M*(a31*a33*del_t^2*g1_a + a32*a33*del_t^2*g2_a + v(:,i)*del_t*a33 + a31*del_t*g1_v ...
                + a32*del_t*g2_v + u(:,i))/(del_t^2*a33^2));
            g3_u = Eff_K{j} \ Eff_R3;
            g3_v = -(a31*del_t*g1_v + a32*del_t*g2_v + u(:,i) - g3_u)/(del_t*a33);
            g3_a = -(a31*a33*del_t^2*g1_a + a32*a33*del_t^2*g2_a + v(:,i)*del_t*a33 + a31*del_t*g1_v + a32*del_t*g2_v + u(:,i) - g3_u)/(del_t^2*a33^2);
            
            if j == num_substeps
                u(:,i+1) = u(:,i) + del_t * (b1 * g1_v + b2 * g2_v + b3 * g3_v);
                v(:,i+1) = v(:,i) + del_t * (b1 * g1_a + b2 * g2_a + b3 * g3_a);
                continue;
            end

        elseif j == 4
            % Fourth sub-step
            g4 = HDSUCI_params.g4;
            a41 = HDSUCI_params.a41;
            a42 = HDSUCI_params.a42;
            a43 = HDSUCI_params.a43;
            a44 = HDSUCI_params.a44;
            b1 = HDSUCI_params.b1;
            b2 = HDSUCI_params.b2;
            b3 = HDSUCI_params.b3;
            b4 = HDSUCI_params.b4;
            F4 = sin(2 * (t(i) + g4 * del_t));
            Eff_R4 = F4 - (-C*(a41*del_t*g1_v + a42*del_t*g2_v + a43*del_t*g3_v + u(:,i))/(del_t*a44) ...
                - M*(a41*a44*del_t^2*g1_a + a42*a44*del_t^2*g2_a + a43*a44*del_t^2*g3_a ...
                + v(:,i)*del_t*a44 + a41*del_t*g1_v + a42*del_t*g2_v + a43*del_t*g3_v + u(:,i))/(del_t^2*a44^2));
            g4_u = Eff_K{j} \ Eff_R4;
            g4_v = -(a41*del_t*g1_v + a42*del_t*g2_v + a43*del_t*g3_v + u(:,i) - g4_u)/(del_t*a44);
            g4_a = -(a41*a44*del_t^2*g1_a + a42*a44*del_t^2*g2_a + a43*a44*del_t^2*g3_a + v(:,i)*del_t*a44 + a41*del_t*g1_v + a42*del_t*g2_v + a43*del_t*g3_v + u(:,i) - g4_u)/(del_t^2*a44^2);
            if  j == num_substeps
                u(:,i+1) = u(:,i) + del_t * (b1 * g1_v + b2 * g2_v + b3 * g3_v + b4 * g4_v);
                v(:,i+1) = v(:,i) + del_t * (b1 * g1_a + b2 * g2_a + b3 * g3_a + b4 * g4_a);
                continue;
            end

        elseif j == 5
            % Fifth sub-step
            g5 = HDSUCI_params.g5;
            a51 = HDSUCI_params.a51;
            a52 = HDSUCI_params.a52;
            a53 = HDSUCI_params.a53;
            a54 = HDSUCI_params.a54;
            a55 = HDSUCI_params.a55;
            b1 = HDSUCI_params.b1;
            b2 = HDSUCI_params.b2;
            b3 = HDSUCI_params.b3;
            b4 = HDSUCI_params.b4;
            b5 = HDSUCI_params.b5;
            F5 = sin(2 * (t(i+1)));
            Eff_R5 = F5 - (-C*(a51*del_t*g1_v + a52*del_t*g2_v + a53*del_t*g3_v + a54*del_t*g4_v ...
                + u(:,i))/(del_t*a55) - M*(a51*a55*del_t^2*g1_a + a52*a55*del_t^2*g2_a + a53*a55*del_t^2*g3_a ...
                + a54*a55*del_t^2*g4_a + v(:,i)*del_t*a55 + a51*del_t*g1_v + a52*del_t*g2_v ...
                + a53*del_t*g3_v + a54*del_t*g4_v + u(:,i))/(del_t^2*a55^2));
            g5_u = Eff_K{j} \ Eff_R5;
            g5_v = -(a51*del_t*g1_v + a52*del_t*g2_v + a53*del_t*g3_v + a54*del_t*g4_v + u(:,i) - g5_u)/(del_t*a55);
            g5_a = -(a51*a55*del_t^2*g1_a + a52*a55*del_t^2*g2_a + a53*a55*del_t^2*g3_a + a54*a55*del_t^2*g4_a + v(:,i)*del_t*a55 ...
                + a51*del_t*g1_v + a52*del_t*g2_v + a53*del_t*g3_v + a54*del_t*g4_v + u(:,i) - g5_u)/(del_t^2*a55^2);
            if  j == num_substeps
                u(:,i+1) = u(:,i) + del_t * (b1 * g1_v + b2 * g2_v + b3 * g3_v + b4 * g4_v + b5 * g5_v);
                v(:,i+1) = v(:,i) + del_t * (b1 * g1_a + b2 * g2_a + b3 * g3_a + b4 * g4_a + b5 * g5_a);
                continue;
            end
                                           
        end
    end
end

end
