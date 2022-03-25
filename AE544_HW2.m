

Ct0 = [0.9385 0.3231 0.1219; -0.3450 0.8618 0.3718; 0.0151 -0.3910 0.9203];
Ctf = [.9768 .0512 -.2079; -.0977 0.9706 -.22; .1905 .2352 .9531];
A = [.9079 -.3702 -.1966; .1951 .7884 -.5835; 0.371 0.4914 0.788];

syms b1 b2 b3 n1 n2 n3
n = [n1;n2;n3]; % inertial frame
b = [b1;b2;b3]; % body frame

b0 = Ct0*n; % body frame at t0
bf = A*b0;
bff = Ctf*n;

syms c(a) s(a) c(b) s(b) c(g) s(g)
R3a = [c(a) s(a) 0; -s(a) c(a) 0; 0 0 1];
R1b = [1 0 0; 0 c(b) s(b); 0 -s(b) c(b)];
R3g = [c(g) s(g) 0; -s(g) c(g) 0; 0 0 1];

Q = R3g*R1b*R3a;
% Q = R3a*R1b*R3g;

apha = atan(A(3,1)/-A(3,2))*180/pi;
beta = acos(A(3,3))*180/pi;
gamma = atan(A(1,3)/A(2,3))*180/pi;
Phi = acos(0.5*(A(1,1) + A(2,2) + A(3,3) - 1))*180/pi; % ____ degrees
e = 1/(2*sin(Phi))*[A(2,3)-A(3,2) ; A(3,1) - A(1,3) ; A(1,2) - A(2,1)];

one = sqrt(e(1,1)^2 + e(2,1)^2 + e(3,1)^2);
