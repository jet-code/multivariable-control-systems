close all
clear 
clc
fprintf('********************************************************************\n')
disp('-------------------------------------------------------------------')
disp('     (LINEAR MULTIVARIABLE CONTROL SYSTEM Computer Project #3)')
disp('-------------------------------------------------------------------')
fprintf('\n********************************************************************\n\n\n') 
A=[0 1 0 0 0 0 0 0 0 0; 
    0 -.11323 -0.98109 -11.847 -11.847 -63.080 -34.339 -34.339 -27.645 0; 
    324.121 -1.1755 -29.101 0.12722 2.83448 -967.73 -678.14 -678.14 0 -129.29; 
    -127.30 0.46167 11.4294 -1.0379 13.1237 380.079 266.341 266.341 0 1054.85; 
    -186.05 0.67475 16.7045 0.86092 -17.068 555.502 389.268 389.268 0 -874.92; 
    341.917 1.09173 1052.75 756.465 756.465 -29.774 0.16507 3.27626 0 0; 
    -30.748 -0.09817 -94.674 -68.029 -68.029 2.67753 -2.6558 4.88497 0 0; 
    -302.36 -0.96543 -930.96 -668.95 -668.95 26.3292 2.42028 -9.5603 0 0;
    0 0 0 0 0 0 0 0 -1.6667 0; 
    0 0 0 0 0 0 0 0 0 -10.000;]
B=[0 0 0 0 0 0 0 0 1.66667 0; 
    0 0 0 0 0 0 0 0 0 10]'
C=[1 0 0 0 0 0 0 0 0 0; 
    -0.49134 0 -0.63203 0 0 -0.20743 0 0 0 0]
D=zeros(2)

poles=eig(A)
[M,D]=eig(A);

disp('-------------------------------------------------------------------')
disp('(Soln 1.) 2nd Order Reduced Model using sign Approach')
disp('-------------------------------------------------------------------')
% Calulating gamma value
r=trace(D)/10
%Determining Ahat with gamma =-5
Ah=A-r*eye(10);

%Calculating Sign(A_hat)
A0=Ah;
A1=0.5*(A0+inv(A0));
x1=trace(A1)-trace(A0);
A2=0.5*(A1+inv(A1));
x2=trace(A2)-trace(A1);
A3=0.5*(A2+inv(A2));
x3=trace(A3)-trace(A2);
A4=0.5*(A3+inv(A3));
x4=trace(A4)-trace(A3);
A5=0.5*(A4+inv(A4));
x5=trace(A5)-trace(A4);
A6=0.5*(A5+inv(A5));
x6=trace(A6)-trace(A5);
A7=0.5*(A6+inv(A6));
x7=trace(A7)-trace(A6);
A8=0.5*(A7+inv(A7));
x8=trace(A8)-trace(A7);
A9=0.5*(A8+inv(A8));
x9=trace(A9)-trace(A8);
A10=0.5*(A9+inv(A9));
x10=trace(A10)-trace(A9);
A11=0.5*(A10+inv(A10));
x11=trace(A11)-trace(A10);
A12=0.5*(A11+inv(A11));
x12=trace(A12)-trace(A11);
A13=0.5*(A12+inv(A12));
x13=trace(A13)-trace(A12);
A14=0.5*(A13+inv(A13));
x14=trace(A14)-trace(A13);
A15=0.5*(A14+inv(A14));
x15=trace(A15)-trace(A14);
A16=0.5*(A15+inv(A15));
x16=trace(A16)-trace(A15)
sign_Ahat=A16

%Calculating sign+(A_hat) & sign-(A_hat)
sign_pA= 0.5*(eye(10)+sign_Ahat);
tr_sign_pA=trace(sign_pA);
sign_mA=0.5*(eye(10)-sign_Ahat);
tr_sign_mA=trace(sign_mA);

disp('Block Decoupled Model') 
m1=[sign_pA(:,1),sign_pA(:,4),sign_pA(:,5),sign_pA(:,8),sign_pA(:,9),sign_pA(:,10)]
m2=[sign_mA(:,1),sign_mA(:,2),sign_mA(:,5),sign_mA(:,6)]
M1=[m1 m2]
Ad = inv(M1)*A*M1
Bd = inv(M1)*B
Cd = C*M1

disp('Reduced Order Model')
lamda1=Ad(1:6,1:6)
lamda2=Ad(7:10,7:10);
B_d1=Bd(1:6,:)
B_d2=Bd(7:10,:);
C_d1=Cd(:,1:6)
C_d2=Cd(:,7:10);
D_star=C*inv(-A)*B+zeros(2)-C_d1*inv(-lamda1)*B_d1

figure(1)
sys1=ss(A,B,C,0);
step(sys1)
hold on
sys2=ss(lamda1,B_d1,C_d1,D_star);
step(sys2)
legend('Original system','2nd Order Reduced Model')

disp('Poles of Reduced Order Model')
poles2=eig(lamda1)

%Controller Designing
disp('New Poles for controller of Reduced Order Model')
P=[-0.3491+0.3491i -0.3491-0.3491i -0.5445+0.5445i -0.5445-0.5445i -7.2421+7.2421i -7.2421-7.2421i]

K=place(lamda1,B_d1,P)
Ahat=lamda1-B_d1*K

figure(2)
sys2=ss(lamda1,B_d1,C_d1,D_star);
impulse(sys2)
hold on
sys3=ss(Ahat,B_d1,C_d1,D_star);
impulse(sys3)
legend('2nd Order Reduced Model','Controlled system')