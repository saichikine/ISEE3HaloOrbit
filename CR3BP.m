%% Sai Chikine - ASEN 5050 Research Project

% CR3BP Equations of motion 

%% 
function [result] = CR3BP(times,r,mu)

% Positions
x = r(1);
y = r(2);
z = r(3);
% Velocities
xdot = r(4);
ydot = r(5);
zdot = r(6);

% Distance from the third body to the larger and smaller body respectively
r1 = sqrt((x+mu)^2 + y^2 + z^2); 
r2 = sqrt((x-1+mu)^2 + y^2 + z^2); 

% Accelerations 
xdotdot = 2*ydot + x -(1-mu)*((x+mu)/(r1^3)) - mu*(x-1+mu)/(r2^3);
ydotdot = -2*xdot + y - (1-mu)*(y/(r1^3)) - mu*(y)/(r2^3);
zdotdot = -(1-mu)*(z)/(r1^3) - mu*(z)/(r2^3); 

% the State transition Matrix for the CRTBP is 
% A(t) = | 0 | I |
%        |g' | 2O|
% where 2O  = [0 1 0;-1 0 ;0 0 0]; and g' is partial deriv matrix

%Pull Jacobian
A2 = Jacobian(x,y,z,xdot,ydot, zdot, mu);

% Multiply the A matrix with the existing STM. 
Phi1 = A2*reshape(r(7:end), 6,[]);
% Reshape to vector. 
Phi2 = reshape(Phi1,1,[]);
% create results vector. 
result = [r(4:6); [xdotdot; ydotdot; zdotdot]; Phi2'];

end