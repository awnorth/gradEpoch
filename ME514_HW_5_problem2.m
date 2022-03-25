% ME514 HW Chapter 5
% Andrew North - 3/24/2022
% --------------------------- Problem 2 --------------------------- %
E = 200E9; % 200 GPa
A = 1E-2; % m^2
I = 2E-4; % m^4
L = 1; % meter

% Element 1
theta = pi/2;
c = 0; % cos(theta)=0
s = sin(theta); % = 1

k_1 = (E/L)*[
    (A*c^2)+(12*I*s^2)/(L^2) (A-((12*I)/L^2))*c*s (6*I*s)/L ;...
    (A-((12*I)/L^2))*c*s (A*s^2)+(12*I*c^2)/(L^2) -(6*I*c)/L ;...
    (6*I*s)/L -(6*I*c)/L 4*I];

% Element 2
theta = 0;
c = 1; % cos(theta)=0
s = 0; % = 0
k_2 = (E/L)*[(A*c^2)+(12*I*s^2)/(L^2) (A-((12*I)/L^2))*c*s -(6*I*s)/L ;...
    (A-((12*I)/L^2))*c*s (A*s^2)+(12*I*c^2)/(L^2) (6*I*c)/L ;...
    -(6*I*s)/L (6*I*c)/L 4*I];

k = k_1 + k_2;
F = [0;0;40000];
d = inv(k)*F; % ANS - displacements and angle at node 2.
