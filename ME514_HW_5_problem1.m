% ME514 HW Chapter 5
% Andrew North - 3/23/2022
% --------------------------- Problem 1 --------------------------- %
% Frame Problem
% part a)
E = 30E3; % k-psi
A = 10; % in^2
I = 200; % in^4

% Element 1 stiffness matrix
theta1 = atan(3/2); % radians
L1 = sqrt((30*12)^2 + (20*12)^2); % Element 1 length (in)
c = cos(theta1);
s = sin(theta1);
k_1 = [(A*c^2)/L1+(12*I*s^2)/(L1^3) (A/L1-((12*I)/L1^3))*c*s (6*I*s)/L1^2 ;...
    (A/L1-((12*I)/L1^3))*c*s (A*s^2)/L1+(12*I*c^2)/(L1^3) -(6*I*c)/L1^2 ;...
    (6*I*s)/L1^2 -(6*I*c)/L1^2 (4*I)/L1];

% Element 2 stiffness matrix
theta2 = 0; % radians
L2 = 30*12; % Element 2 length (in)
c = cos(theta2);
s = sin(theta2);
k_2 = [(A*c^2)/L2+(12*I*s^2)/(L2^3) (A/L2-((12*I)/L2^3))*c*s -(6*I*s)/L2^2 ;...
    (A/L2-((12*I)/L2^3))*c*s (A*s^2)/L2+(12*I*c^2)/(L2^3) (6*I*c)/L2^2 ;...
    -(6*I*s)/L2^2 (6*I*c)/L2^2 (4*I)/L2];

k = E*(k_1 + k_2); % global stiffness matrix at note 2
d = inv(k)*[0;(-40-10-20);(600.07-1800)]; % node 2 displacements & rotations. Nodes 1 & 3 = 0

% ---------- part b) - find the element forces ----------
% ---------- Element 1 forces ---------- %
theta = theta1;
C = cos(theta);
S = sin(theta);
T = [C S 0 0 0 0;-S C 0 0 0 0;0 0 1 0 0 0;...
     0 0 0 C S 0;0 0 0 -S C 0;0 0 0 0 0 1];
L = L1; 
C1 = (A*E)/L;
C2 = (E*I)/L^3;
k_prime = [C1 0 0 -C1 0 0;0 12*C2 6*C2*L 0 -12*C2 6*C2*L;...
    0 6*C2*L 4*C2*L^2 0 -6*C2*L 2*C2*L^2;-C1 0 0 C1 0 0;...
    0 -12*C2 -6*C2*L 0 12*C2 -6*C2*L;0 6*C2*L 2*C2*L^2 0 -6*C2*L 4*C2*L^2];
d_1 = [0;0;0;d(1,1);d(2,1);d(3,1)];
T_d_1 = T*d_1;
f_prime_1 = k_prime*T_d_1; % effective nodal forces

% Element 1 equivalent local nodal forces
f_0_1 = [-8.32 ; -5.547 ; -600.07 ; -8.32 ; -5.547 ; 600.07];

% Element 1 local nodal forces
f_1 = f_prime_1 - f_0_1;

% ---------- Element 2 forces ---------- %
theta = theta2;
C = cos(theta);
S = sin(theta);
T = [C S 0 0 0 0;-S C 0 0 0 0;0 0 1 0 0 0;...
     0 0 0 C S 0;0 0 0 -S C 0;0 0 0 0 0 1];
L = L2; 
C1 = (A*E)/L;
C2 = (E*I)/L^3;
k_prime = [C1 0 0 -C1 0 0;0 12*C2 6*C2*L 0 -12*C2 6*C2*L;...
    0 6*C2*L 4*C2*L^2 0 -6*C2*L 2*C2*L^2;-C1 0 0 C1 0 0;...
    0 -12*C2 -6*C2*L 0 12*C2 -6*C2*L;0 6*C2*L 2*C2*L^2 0 -6*C2*L 4*C2*L^2];
d_2 = [d(1,1);d(2,1);d(3,1);0;0;0];
T_d_2 = T*d_2;
f_prime_2 = k_prime*T_d_2; % effective nodal forces

% Element 2 equivalent local nodal forces
f_0_2 = [0 ; -20 ; -1800 ; 0 ; -20 ; 1800];

% Element 2 local nodal forces
f_2 = f_prime_2 - f_0_2;

% ---------- part c) - find the reactions ----------
% Element 1 stiffness matrix
theta1 = atan(3/2); % radians
L = sqrt((30*12)^2 + (20*12)^2); % Element 1 length (in)
c = cos(theta1);
s = sin(theta1);
k_1_full = E*[(A*c^2)/L+(12*I*s^2)/(L^3) (A/L-((12*I)/L^3))*c*s -(6*I*s)/L^2 ...
    -((A*c^2)/L+(12*I*s^2)/(L^3)) -(A/L-((12*I)/L^3))*c*s -(6*I*s)/L^2 ;...
    (A/L-((12*I)/L^3))*c*s (A*s^2)/L+(12*I*c^2)/(L^3) (6*I*c)/L^2 ...
    -((A/L-((12*I)/L^3))*c*s) -((A*s^2)/L+(12*I*c^2)/(L^3)) (6*I*c)/L^2 ;...
    -(6*I*s)/L^2 (6*I*c)/L^2 (4*I)/L ...
    (6*I*s)/L^2 -(6*I*c)/L^2 (2*I)/L ;...
    -((A*c^2)/L+(12*I*s^2)/(L^3)) -((A/L-((12*I)/L^3))*c*s) (6*I*s)/L^2 ...
    (A*c^2)/L+(12*I*s^2)/(L^3) (A/L-((12*I)/L^3))*c*s (6*I*s)/L^2 ;...
    -(A/L-((12*I)/L^3))*c*s -((A*s^2)/L+(12*I*c^2)/(L^3)) -(6*I*c)/L^2 ...
    (A/L-((12*I)/L^3))*c*s (A*s^2)/L+(12*I*c^2)/(L^3) -(6*I*c)/L^2 ;...
    -(6*I*s)/L^2 (6*I*c)/L^2 (2*I)/L ...
    (6*I*s)/L^2 -(6*I*c)/L^2 (4*I)/L];

F_1_reactions = k_1_full*d_1 - [0;-10;-600.07;0;-10;600.07];


% Element 2 stiffness matrix
theta2 = 0; % radians
L = 30*12; % Element 2 length (in)
c = cos(theta2);
s = sin(theta2);
k_2_full = E*[(A*c^2)/L+(12*I*s^2)/(L^3) (A/L-((12*I)/L^3))*c*s -(6*I*s)/L^2 ...
    -((A*c^2)/L+(12*I*s^2)/(L^3)) -(A/L-((12*I)/L^3))*c*s -(6*I*s)/L^2 ;...
    (A/L-((12*I)/L^3))*c*s (A*s^2)/L+(12*I*c^2)/(L^3) (6*I*c)/L^2 ...
    -((A/L-((12*I)/L^3))*c*s) -((A*s^2)/L+(12*I*c^2)/(L^3)) (6*I*c)/L^2 ;...
    -(6*I*s)/L^2 (6*I*c)/L^2 (4*I)/L ...
    (6*I*s)/L^2 -(6*I*c)/L^2 (2*I)/L ;...
    -((A*c^2)/L+(12*I*s^2)/(L^3)) -((A/L-((12*I)/L^3))*c*s) (6*I*s)/L^2 ...
    (A*c^2)/L+(12*I*s^2)/(L^3) (A/L-((12*I)/L^3))*c*s (6*I*s)/L^2 ;...
    -(A/L-((12*I)/L^3))*c*s -((A*s^2)/L+(12*I*c^2)/(L^3)) -(6*I*c)/L^2 ...
    (A/L-((12*I)/L^3))*c*s (A*s^2)/L+(12*I*c^2)/(L^3) -(6*I*c)/L^2 ;...
    -(6*I*s)/L^2 (6*I*c)/L^2 (2*I)/L ...
    (6*I*s)/L^2 -(6*I*c)/L^2 (4*I)/L];

F_2_reactions = k_2_full*d_2 - [0;-20;-1800;0;-20;1800];
