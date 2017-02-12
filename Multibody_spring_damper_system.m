%% ARO406-01 - Advanced Dynamics & Vibrations
% California State Polytechnic University, Pomona
% Ryan Gunawan
clc
clear all

%% Wing Parameters
L = 4;
EI = 2E6;
J = .02; 
F = L^3/EI*[1/24 5/48       % Flexibility Matrix
     5/48 1/3];
 
K = [2 5
    5 16];
                   
Force = [1000
    0];

I = [1 0
    0 1];

%% Case 1: Mass 1 > Mass 2
m1 = 100;
m2 = 1;
M = [m1 0 
    0 m2];


Matrix1 = K*M;

syms T

Matrix1 - T*I;
det(Matrix1 - T*I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigenValue = solve(det(Matrix1 - T*I));
EigenValue(1);
EigenValue(2);

natfreq1 = (48*EI/(EigenValue(1)*L^3))^(1/2);
natfreq2 = (48*EI/(EigenValue(2)*L^3))^(1/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1_z2_1 = -Matrix1(1,2)/((Matrix1(1,1)-EigenValue(1)));
z1_z2_2 = -Matrix1(1,2)/((Matrix1(1,1)-EigenValue(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigenVector1 = [z1_z2_1
    1];
EigenVector2 = [z1_z2_2
    1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correction Factor

alpha1 = (transpose(EigenVector1)*M*EigenVector1)^(-1/2);
alpha2 = (transpose(EigenVector2)*M*EigenVector2)^(-1/2);

OrthoEigenVector1 = alpha1*EigenVector1;
OrthoEigenVector2 = alpha2*EigenVector2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ModalMatrix = [OrthoEigenVector1 OrthoEigenVector2]

syms q1 q2
q = [q1
    q2];

Q = transpose(ModalMatrix)*Force
zfunc = ModalMatrix*q;

w_n1 = natfreq1;
w_n2 = natfreq2;
w_d1 = w_n1*sqrt(1-J^2);
w_d2 = w_n2*sqrt(1-J^2);
syms t
q1 = Q(1)/w_n1^2*(1-exp((-J*w_n1*t)*(cos(w_d1*t)+J/sqrt((1-J^2))*sin(w_d1*t))));
q2 = Q(2)/w_n2^2*(1-exp((-J*w_n2*t)*(cos(w_d2*t)+J/sqrt((1-J^2))*sin(w_d2*t))));


%Plotting
plotfunc = ModalMatrix(1,1)*q1+ModalMatrix(1,2)*q2

figure(1)
ezplot(plotfunc, [0, .4])
title('Displacement of M1 > M2')
xlabel('Time (s)') % x-axis label
ylabel('Displacement (m)') % y-axis label
%% Case 2: Mass 1 == Mass 2
m1 = 100;
m2 = 100;
M = [m1 0 
    0 m2];


Matrix1 = K*M;

syms T

Matrix1 - T*I;
det(Matrix1 - T*I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigenValue = solve(det(Matrix1 - T*I));
EigenValue(1);
EigenValue(2);

natfreq1 = (48*EI/(EigenValue(1)*L^3))^(1/2);
natfreq2 = (48*EI/(EigenValue(2)*L^3))^(1/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1_z2_1 = -Matrix1(1,2)/((Matrix1(1,1)-EigenValue(1)));
z1_z2_2 = -Matrix1(1,2)/((Matrix1(1,1)-EigenValue(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigenVector1 = [z1_z2_1
    1];
EigenVector2 = [z1_z2_2
    1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correction Factor

alpha1 = (transpose(EigenVector1)*M*EigenVector1)^(-1/2);
alpha2 = (transpose(EigenVector2)*M*EigenVector2)^(-1/2);

OrthoEigenVector1 = alpha1*EigenVector1;
OrthoEigenVector2 = alpha2*EigenVector2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ModalMatrix = [OrthoEigenVector1 OrthoEigenVector2]

syms q1 q2
q = [q1
    q2];

Q = transpose(ModalMatrix)*Force
zfunc = ModalMatrix*q;

w_n1 = natfreq1;
w_n2 = natfreq2;
w_d1 = w_n1*sqrt(1-J^2);
w_d2 = w_n2*sqrt(1-J^2);
syms t
q1 = Q(1)/w_n1^2*(1-exp((-J*w_n1*t)*(cos(w_d1*t)+J/sqrt((1-J^2))*sin(w_d1*t))));
q2 = Q(2)/w_n2^2*(1-exp((-J*w_n2*t)*(cos(w_d2*t)+J/sqrt((1-J^2))*sin(w_d2*t))));


%Plotting
plotfunc = ModalMatrix(1,1)*q1+ModalMatrix(1,2)*q2
figure(2)
ezplot(plotfunc, [0, .4])
title('Displacement of M1 = M2')
xlabel('Time (s)') % x-axis label
ylabel('Displacement (m)') % y-axis label
%% Case 3: Mass 1 < Mass 2

m1 = 1;
m2 = 100;
M = [m1 0 
    0 m2];


Matrix1 = K*M;

syms T

Matrix1 - T*I;
det(Matrix1 - T*I);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigenValue = solve(det(Matrix1 - T*I));
EigenValue(1);
EigenValue(2);

natfreq1 = (48*EI/(EigenValue(1)*L^3))^(1/2);
natfreq2 = (48*EI/(EigenValue(2)*L^3))^(1/2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
z1_z2_1 = -Matrix1(1,2)/((Matrix1(1,1)-EigenValue(1)));
z1_z2_2 = -Matrix1(1,2)/((Matrix1(1,1)-EigenValue(2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EigenVector1 = [z1_z2_1
    1];
EigenVector2 = [z1_z2_2
    1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correction Factor

alpha1 = (transpose(EigenVector1)*M*EigenVector1)^(-1/2);
alpha2 = (transpose(EigenVector2)*M*EigenVector2)^(-1/2);

OrthoEigenVector1 = alpha1*EigenVector1;
OrthoEigenVector2 = alpha2*EigenVector2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ModalMatrix = [OrthoEigenVector1 OrthoEigenVector2]     

syms q1 q2
q = [q1
    q2];

Q = transpose(ModalMatrix)*Force
zfunc = ModalMatrix*q;

w_n1 = natfreq1;
w_n2 = natfreq2;
w_d1 = w_n1*sqrt(1-J^2);
w_d2 = w_n2*sqrt(1-J^2);
syms t
q1 = Q(1)/w_n1^2*(1-exp((-J*w_n1*t)*(cos(w_d1*t)+J/sqrt((1-J^2))*sin(w_d1*t))));
q2 = Q(2)/w_n2^2*(1-exp((-J*w_n2*t)*(cos(w_d2*t)+J/sqrt((1-J^2))*sin(w_d2*t))));


%Plotting
plotfunc = ModalMatrix(1,1)*q1+ModalMatrix(1,2)*q2
figure(3)
ezplot(plotfunc, [0, .4])
title('Displacement of M1 < M2')
xlabel('Time (s)') % x-axis label
ylabel('Displacement (m)') % y-axis label

