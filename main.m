%% Sai Chikine - ASEN 5050 Research Project
%
% This project aims to replicate the ISEE-3 Halo orbit around Sun-Earth L1.
% 

%%
clear all; close all; clc

%% Basic Constants

mSun = 1.9891e30;% kg
mEarth = 5.9742e24; %kg
mu = mEarth/(mSun + mEarth);
L = 1.495978714e8; %Sun-Earth Distance, km
E_L1 = 1.497610041e6;
gamma1 = E_L1/L;%km, distance between Earth and L1

%% Constants from Koon Lo Marsden Ross

% Sun Earth L1

xL1 = 1 - mu - gamma1;
xL2 = 1.0100752000;

c2 = 4.0610735668;
c3 = 3.0200105081;
c4 = 3.0305378797;

s1 = -8.246608317e-1;
s2 = 1.210985938e-1;
l1 = -1.596560314e1;
l2 = 1.740900800;
deltaconst = 0.29221444425;
kappa = 3.2292680962;

Az = 110000/E_L1; %km
% Ax = sqrt((-l2*(Az/E_L1)^2 - deltaconst)/l1);
Ax = sqrt((-l2*(Az)^2 - deltaconst)/l1);
Ay = kappa*Ax;

wp = 2.086453455;
wv = 2.0152105515;

nu1 = 0;
nu2 = s1*Ax^2 + s2*Az^2;
nu = 1 + nu1 + nu2;
 
T = 2*pi/(wp*nu);

%% Initial Guess (third order)

m = 3; %southern orbit?
deltam = 2-m;

%Constants needed for third order Richardson extrapolation solution
%Used to find approximate initial conditions

%lambda = sqrt((c2 + sqrt(9*c2^2 - 8*c2))/2);
%lambda = -1.46713;
lambda = 2.08646;
d1 = 3*lambda^2/kappa*(kappa*(6*lambda^2 - 1) - 2*lambda);
d2 = 8*lambda^2/kappa*(kappa*(11*lambda^2 - 1) - 2*lambda);

a21 = 3*c3*(kappa^2 - 2)/(4*(1+2*c2));
a22 = 3*c3/(4*(1+2*c2));
a23 = -3*c3*lambda/(4*kappa*d1)*(3*kappa^3*lambda - 6*kappa*(kappa - lambda) + 4);
a24 = -3*c3*lambda/(4*kappa*d1)*(2 + 3*kappa*lambda);

b21 = -3*c3*lambda/(2*d1)*(3*kappa*lambda - 4);
b22 = -3*c3*lambda/d1;

d21 = -c3/(2*lambda^2);

a31 = -9*lambda/(4*d2)*(4*c3*(kappa*a23 - b21) + kappa*c4*(4 + kappa^2)) + ...
    (9*lambda^2 + 1 - c2)/(2*d2)*(3*c3*(2*a23 - kappa*b21) + c4*(2 + 3*kappa^2));
a32 = -9*lambda/(4*d2)*(4*c3*(3*kappa*a24 - b22) + kappa*c4) - ...
    (-3*(9*lambda^2 + 1 - c2))/(2*d2)*(c3*(kappa*b22 + d21 - 2*a24) - c4);

b31 = 3/(8*d2)*8*lambda*(3*c3*(kappa*b21 - 2*a23) - c4*(2 + 3*kappa^2)) + ...
    3/(8*d2)*(9*lambda^2 + 1 + 2*c2)*(4*c3*(kappa*a23 - b21) + kappa*c4*(4 + kappa^2));
b32 = 9*lambda/d2*(c3*(kappa*b22 + d21 - 2*a24) - c4) + ...
    3*(9*lambda^2 + 1 + 2*c2)/(8*d2)*(4*c3*(kappa*a24 - b22) + kappa*c4);

d31 = 3/(64*lambda^2)*(4*c3*a24 + c4);
d32 = 3/(64*lambda^2)*(4*c3*(a23 - d21) + c4*(4 + kappa^2));

s1c = (2*lambda*(lambda*(1+kappa^2) - 2*kappa))^(-1)*(3/2*c3*(2*a21*(kappa^2 - 2) - a23*(kappa^2 + 2) - 2*kappa*b21) - 3/8*c4*(3*kappa^4 - 8*kappa^2 + 8));


%% Initial Guess

xguess = a21*Ax^2 + a22*Az^2 - Ax + (a23*Ax^2 - a24*Az^2) + (a31*Ax^3 - a32*Ax*Az^2);
ydotguess = wp*kappa*Ax + 2*wp*(b21*Ax^2 - b22*Az^2) + 3*wp*(b31*Ax^3 - b32*Ax*Az^2);
zguess = deltam*Az + deltam*d21*Ax*Az*(1-3) + deltam*(d32*Az*Ax^2 - d31*Az^3);
% find rough value of x0, z0 and dy0

% x0 = (1-mu)+(-Ax*dist_normalized/L);
% z0 = Az*dist_normalized/L;
% dy0 = wp*Ay*dist_normalized/L;

x0 = (1-mu) + (-9.697924e-4);
z0 = 7.412574e-4;
dy0 = 0.00909109;
% 
% x0 = (1-mu) + (-9.697919e-4);
% z0 = 7.412574e-4;
% dy0 = 0.00909109;

x0c = (1-mu) + gamma1*(xguess-1);
z0c = zguess*gamma1;
dy0c = ydotguess*gamma1;

%% Plot initial orbit
% Initialize the state transition matrix
STM = reshape(eye(6),36,[]);
X0Guess = [x0c; 0; z0c; 0; dy0c; 0; STM]; %[x, y, z, xdot, ydot, zdot]

% Orbit period computed above
orbitPeriod = T;

% Integrate orbit from initial guess
ode_options = odeset('RelTol',1e-13,'AbsTol',1e-22);
[tGuess, XGuess] = ode45(@CR3BP, [0 orbitPeriod], X0Guess, ode_options, mu);

%% Initial Plot

%3d plot of sun, earth, l1, l2, and guessed orbit
%Do not remove holds, figure breaks if you do
figure(1)
plot3((XGuess(:,1)), XGuess(:,2), XGuess(:,3)); hold on
plot3(0-mu, 0, 0,'ok', 'markerfacecolor', 'y', 'markersize', 10); hold on %Sun
plot3(1-mu, 0, 0,'ok', 'markerfacecolor', 'b', 'markersize',5); hold on %Earth
plot3(xL1, 0, 0,'ok', 'markerfacecolor', 'r', 'markersize', 2); hold on %L1 point
plot3(xL2, 0, 0,'ok', 'markerfacecolor', 'r', 'markersize', 2); hold on %L2 point
title('Initial (Approximate) Halo Orbit');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
grid on;
hold off

% Close up of orbit centered around L1
figure(2)
plot3(xL1*L,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 2); hold on
plot3((XGuess(:,1))*L, XGuess(:,2)*L, XGuess(:,3)*L)
title('Initial (Approximate) Halo Orbit (Close Up)');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
grid on;

%% Differential Correction (single shooting method)

% Stop integrating when y=0 (corresponds to half an orbit)
ode_options = odeset('Events',@Findyzero,'RelTol',1e-13,'AbsTol',1e-13);
X0 = X0Guess;

verbose = 1;

deltaVec = [100; 100];
counter = 0;
while abs(norm(deltaVec)) > 1e-13
    counter = counter+1;
    [tHalfOrbit,XHalfOrbit] = ode45(@CR3BP, [0 Inf], X0, ode_options, mu);

    % State when y=0 (half orbit period)
    dX = XHalfOrbit(end,1:6);
    
    % STM at half orbit period
    dPhi = reshape(XHalfOrbit(end,7:end),6,[]);

    %% dX/dt (X is state)
    xdot = dX(4);
    ydot = dX(5);
    zdot = dX(6);

    r1 = sqrt((dX(1)+mu)^2 + dX(2)^2 + dX(3)^2); %S/C distance to Sun
    r2 = sqrt((dX(1)-1+mu)^2 + dX(2)^2 + dX(3)^2); %S/C distance to Earth

    % Accelerations
    xdotdot = 2*ydot + dX(1) - (1 - mu)*((dX(1) + mu)/(r1^3)) - mu*(dX(1) - 1 + mu)/(r2^3);
    ydotdot = -2*xdot + dX(2) - (1 - mu)*(dX(2)/(r1^3)) - mu*(dX(2))/(r2^3);
    zdotdot = -(1 - mu)*(dX(3))/(r1^3) - mu*(dX(3))/(r2^3); 
    
    % Derivative of new state
    dXdt = [xdot ydot zdot xdotdot ydotdot zdotdot];
    
    %Update matrix to correct deltaVec
    %Only changes x0 and doty0, holds z0 constant
    updateMat = [dPhi(4,1)-dPhi(2,1)*(xdotdot/ydot), dPhi(4,5)-dPhi(2,5)*(xdotdot/ydot);...
        dPhi(6,1)-dPhi(2,1)*(zdotdot/ydot), dPhi(6,5)-dPhi(2,5)*(zdotdot/ydot)];

    deltaVec = inv(updateMat)*[-XHalfOrbit(end,4); -XHalfOrbit(end,6)];
    
    if verbose
        norm(deltaVec)
        fprintf('Iteration counter: %d\n', counter)
    end

    deltaX = [deltaVec(1), 0, 0, 0, deltaVec(2), 0]';

    X0(1:6) = X0(1:6) + deltaX;
end

%% New orbit with periodic initial conditions

%Integrate orbit with solved initial conditions
ode_options = odeset('RelTol',1e-13,'AbsTol',1e-22);
[tFinal, XFinal] = ode45(@CR3BP, [0 tHalfOrbit(end)*2], X0, ode_options, mu); 

% 3d plot of sun, earth, l1, l2, and new periodic orbit
%Do not remove holds, figure breaks if you do
figure(10)
plot3((XFinal(:,1)),XFinal(:,2),XFinal(:,3),'linewidth',1.5); hold on
plot3(0-mu,0,0,'ok', 'markerfacecolor', 'y', 'markersize', 10); hold on %Sun
plot3(1-mu,0,0,'ok', 'markerfacecolor', 'b', 'markersize',5); hold on %Earth
plot3(xL1,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 2); hold on %L1 Point
plot3(xL2,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 2); hold on %L2 Point
title('Periodic Halo Orbit');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
grid on

% Close up view of new periodic orbit centered around L1
figure(20)
plot3(xL1*L,0,0,'ok', 'markerfacecolor', 'r', 'markersize', 2); hold on
plot3((XFinal(:,1))*L, XFinal(:,2)*L, XFinal(:,3)*L)
title('Periodic Halo Orbit (Close Up)');
xlabel('X [km]');
ylabel('Y [km]');
zlabel('Z [km]');
grid on;

%% 
%KLMR Plot
figure(5)
subplot(1,3,1) %y vs x
plot((XFinal(:,1)-1+mu+gamma1)*L/10e4, XFinal(:,2)*L/10e4)
xlim([-7,7]);
xlabel('x');
axis('equal');
ylim([-7,7]);
ylabel('y');

subplot(1,3,2) %z vs x
plot((XFinal(:,1)-1+mu+gamma1)*L/10e4, XFinal(:,3)*L/10e4)
xlim([-7,7]);
xlabel('x');
axis('equal');
ylim([-7,7]);
ylabel('z');

subplot(1,3,3) %z vs y
plot(XFinal(:,2)*L/(10e4), XFinal(:,3)*L/(10e4))
xlim([-7,7]);
xlabel('y');
axis('equal');
ylim([-7,7]);
ylabel('z');

AxMag = (max(XFinal(:,1)*L) - min(XFinal(:,1)*L))/2
AyMag = (max(XFinal(:,2)*L) - min(XFinal(:,2)*L))/2
AzMag = (max(XFinal(:,3)*L) - min(XFinal(:,3)*L))/2
PeriodMag = tFinal(end)