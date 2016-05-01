clc; clear;

A = load('A_ex1.txt')
B = load('B_ex1.txt')

n = 4;
m = 2;

A_eigvals_desired = [-0.2   , -0.5   , -5.0566, -8.6659];

[U,S,V] = svd(B)
U_0 = U(:,1:m)
U_1 = U(:,m+1:n)

Z = S*V;
Z = Z(1:m,:)

X = zeros(4);

for j = 1:n
    lambda_j = A_eigvals_desired(j);
    temp = null(U_1'*(A - lambda_j*eye(n)));
    X(:,j) = temp(:,2);
end

K =  place(A,B,A_eigvals_desired);

M = A + B*-1*K;

[M_evec,M_evals] = eig(M);

Bc = [0,0;0,0;1,0;0,1]

Tc1 = Bc'*inv([B,A*B])

Tc = [Tc1;Tc1*A]

Tc_inv = inv(Tc)

Ac = Tc*A*Tc_inv

Tc*B












