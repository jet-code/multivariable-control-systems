% Alexander Hebert
% ECE 6390
% CP #1 continued from Python

clc; clear;

A = load('A.txt');
B = load('B.txt')';
C = load('C.txt');

disp('Original system:')
sys_original = ss(A,B,C,0)

A_hat = csvread('A_hat_cp1.csv');
B_hat = csvread('B_hat_cp1.csv');
C_hat = csvread('C_hat_cp1.csv');

disp(' ')
disp('Second Cauer reduced model (2nd order):')
sys_2ndCauer = ss(A_hat,B_hat,C_hat,0)

Ar = csvread('Ar_cp1.csv');
Br = csvread('Br_cp1.csv');
Cr = csvread('Cr_cp1.csv');
Dr = csvread('Dr_cp1.csv');

disp(' ')
disp('Reduced model from residualization:')
sys_residue = ss(Ar,Br,Cr,Dr)

figure(1)
step(sys_original,'b',sys_2ndCauer,'g--')
legend('Original','2nd Cauer','Location','SouthOutside','Orientation','horizontal')

figure(2)
impulse(sys_original,'b',sys_2ndCauer,'g--')
legend('Original','2nd Cauer','Location','SouthOutside','Orientation','horizontal')

% figure(3)
% step(sys_original,'b',sys_residue,'r-.')
% legend('Original','Residualization')
% 
% figure(4)
% impulse(sys_original,'b',sys_residue,'r-.')
% legend('Original','Residualization')

% sinusoidal input
% one positive and one negative sine wave input
[u_sine,t] = gensig('sin',25,50,0.1);  
u_sine_neg = -1*u_sine;             
figure(5)
lsimplot(sys_original,'b',sys_2ndCauer,'r:')
legend('Original system','2nd Cauer 2nd order reduced model','Location','SouthOutside','Orientation','horizontal')

