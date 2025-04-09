function a = HDSUCI_acc(u, del_t, t)
% Compute acceleration using finite difference method
% 8th order of accuracy
%
% Inputs:
%   u      - Displacement vector 
%   del_t  - Time step size (scalar)
%
% Outputs:
%   a      - Acceleration vector 

a = zeros(1, length(t));

% Compute acceleration using finite difference method
% 8th order of accuracy
for i = 1:length(t)
    if i <= 5
        % Forward difference approximation
        a(:,i) = (6515/1008*u(:,i) - 4609/140*u(:,i+1) + 5869/70*u(:,i+2) ...
            - 6289/45*u(:,i+3) + 6499/40*u(:,i+4) - 265/2*u(:,i+5) ...
            + 6709/90*u(:,i+6) - 967/35*u(:,i+7) + 3407/560*u(:,i+8) - 761/1260*u(:,i+9)) / del_t^2;
    elseif i >= length(t) - 5
        % Backward difference approximation
        a(:,i) = (6515/1008*u(:,i) - 4609/140*u(:,i-1) + 5869/70*u(:,i-2) ...
            - 6289/45*u(:,i-3) + 6499/40*u(:,i-4) - 265/2*u(:,i-5) ...
            + 6709/90*u(:,i-6) - 967/35*u(:,i-7) + 3407/560*u(:,i-8) - 761/1260*u(:,i-9)) / del_t^2;    
    else
        % Central difference approximation
        a(:,i) = ((-1/560)*u(:,i-4) + (8/315)*u(:,i-3) + (-1/5)*u(:,i-2) + (8/5)*u(:,i-1) ...
            + (-205/72)*u(:,i) + (8/5)*u(:,i+1) + (-1/5)*u(:,i+2) ...
            + (8/315)*u(:,i+3) + (-1/560)*u(:,i+4)) / del_t^2;
    end
end