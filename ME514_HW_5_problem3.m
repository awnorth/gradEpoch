% ME514 HW Chapter 5
% Andrew North - 3/24/2022
% --------------------------- Problem 3 --------------------------- %
E = 210E9; % Pa
G = 84E9; % Pa
I = 2E-4; % m^4
J = 1E-4; % m^4
A = 1E-2; % m^2

L1 = 3; % m
L2 = 3; % m
L3 = sqrt(9+9); % m

% Element 1
c = 0;
s = 1;
L = L1;
k_G_prime_1 = [(12*E*I)/L^3 0 -(6*E*I)/L^2 ; 0 (G*J)/L 0 ; -(6*E*I)/L^2 0 (4*E*I)/L];
T_G_1 = [1 0 0;0 c s;0 -s c]; %transformation matrix
k_G_1 = T_G_1'*k_G_prime_1*T_G_1;

% Element 3
c = sqrt(2)/2;
s = sqrt(2)/2;
L = L3;
k_G_prime_2 = [(12*E*I)/L^3 0 -(6*E*I)/L^2 ; 0 (G*J)/L 0 ; -(6*E*I)/L^2 0 (4*E*I)/L];
T_G_2 = [1 0 0;0 c s;0 -s c]; %transformation matrix
k_G_2 = T_G_2'*k_G_prime_2*T_G_2;
k_G = (k_G_1 + k_G_2);
F = [-10000;0;0];
d = inv(k_G)*F; % nodal displacements

% ------------ Element Forces 
% Element 1
c = 0;
s = 1;
L = L1;
k_G_prime_1 = [(12*E*I)/L^3 0 (6*E*I)/L^2 -(12*E*I)/L^3 0 (6*E*I)/L^2 ;...
    0 (G*J)/L 0 0 -(G*J)/L 0;...
    (6*E*I)/L^2 0 (4*E*I)/L -(6*E*I)/L^2 0 (2*E*I)/L ;...
    -(12*E*I)/L^3 0 -(6*E*I)/L^2 (12*E*I)/L^3 0 -(6*E*I)/L^2 ;...
    0 -(G*J)/L 0 0 (G*J)/L 0;...
    (6*E*I)/L^2 0 (2*E*I)/L -(6*E*I)/L^2 0 (4*E*I)/L];
a = (1/1E7)*k_G_prime_1;
T_G_1 = [1 0 0 0 0 0;0 c s 0 0 0;0 -s c 0 0 0 ;...
    0 0 0 1 0 0; 0 0 0 0 c s; 0 0 0 0 -s c]; %transformation matrix
d_1 = [0;0;0;d(1,1);d(2,1);d(3,1)];
f_1_1 = T_G_1*d_1;
f_1 = k_G_prime_1*f_1_1;

% Element 2
c = sqrt(2)/2;
s = sqrt(2)/2;
L = L3;
k_G_prime_2 = [(12*E*I)/L^3 0 (6*E*I)/L^2 -(12*E*I)/L^3 0 (6*E*I)/L^2 ;...
    0 (G*J)/L 0 0 -(G*J)/L 0;...
    (6*E*I)/L^2 0 (4*E*I)/L -(6*E*I)/L^2 0 (2*E*I)/L ;...
    -(12*E*I)/L^3 0 -(6*E*I)/L^2 (12*E*I)/L^3 0 -(6*E*I)/L^2 ;...
    0 -(G*J)/L 0 0 (G*J)/L 0;...
    (6*E*I)/L^2 0 (2*E*I)/L -(6*E*I)/L^2 0 (4*E*I)/L];
b = (1/1E7)*k_G_prime_2;
T_G_2 = [1 0 0 0 0 0;0 c s 0 0 0;0 -s c 0 0 0 ;...
    0 0 0 1 0 0; 0 0 0 0 c s; 0 0 0 0 -s c]; %transformation matrix
d_2 = [0;0;0;d(1,1);d(2,1);d(3,1)];
f_2_1 = T_G_2*d_2;
f_2 = k_G_prime_2*f_2_1;


