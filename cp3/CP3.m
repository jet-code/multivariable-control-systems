% Alexander Hebert
% ECE 6390
% Computer Project #3 continued from Python

clear; clc;

A = load('A.txt');
B = load('B.txt');
C = load('C.txt');

disp('Original system:')
sys_original = ss(A,B,C,0)

figure(1)
pzmap(sys_original)


L1 = dlmread('Lambda1.txt',' '); % L = capital lambda
Bd1 = dlmread('Bd1.txt',' ');
Cd1 = dlmread('Cd1.txt',' ');
Dstar = dlmread('Dstar.txt',' ');

disp(' ')
disp('Matrix Sign Function (msf) 4th order reduced model:')
sys_msf_2nd = ss(L1,Bd1,Cd1,Dstar)

figure(2)
impulse(sys_original,'b',sys_msf_2nd,'g--')
legend('Original','Matrix Sign Function (msf) 6th order','Location','SouthOutside','Orientation','horizontal')


% L1_eigvals = [...
% -3.491460840764012419e-01+6.344006476918209181e+00j,...
% -3.491460840764012419e-01-6.344006476918209181e+00j,...
% -1.042053905840172989e+00,...
% -2.344831273439368091e-01,...
% -1.000000000000000355e+01,...
% -1.666700000000000070e+00]

% Decide on desired pole/eigenvalue locations (determined by trial and error).
% Complex conjugate pair must have equal real and imaginary parts.
A_eigvals_desired = [-0.5, -1, -2, -40+40j, -40-40j, -45]

% Design controller for reduced order model.
% Compute state-feedback matrix K using Matlab's place command
K = place(L1,Bd1,A_eigvals_desired)

A_control = L1 - Bd1*K;
sys_msf_2nd_control = ss(A_control,Bd1,Cd1,Dstar)

figure(3)
impulse(sys_original,'b',sys_msf_2nd_control,'r--')
legend('Original','Controlled 6th order msf model','Location','SouthOutside','Orientation','horizontal')

% Design steady-state Kalman filter with generic noise profile.
% F = V = W = I

% Steady-state Kalman filter
% Solve ARE for error covariance matrix P
% For Matlab's care() function, R^-1 = W^-1 = I

% For compatibility with class notes:
% A_care = A'
% B_care = C'
% X_care = P
% Q = F*V*F' = I
Q = eye(6);

[P_hat,L,G,report] = care(L1',Cd1',Q);
report

% Steady-state Kalman gain
K_hat = P_hat*Cd1' % W^-1 = I so it is left out of the product

% Overall system with Kalman filter and state-feedback control:
A_co = [L1, -Bd1*K; K_hat*Cd1, (L1 - K_hat*Cd1 - Bd1*K)];
B_co = [Bd1; Bd1];
C_co = [Cd1, Cd1];

sys_co = ss(A_co,B_co,C_co,Dstar)

figure(4)
impulse(sys_original,'b',sys_co,'g--')
legend('Original','6th order model with KF and feedback control','Location','SouthOutside','Orientation','horizontal')


