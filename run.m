%--------------------------------------------------------------------------
% Matlab M-file Project: HyEQ Toolbox @  Hybrid Systems Laboratory (HSL), 
% https://hybrid.soe.ucsc.edu/software
% http://hybridsimulator.wordpress.com/
% Filename: run_ex1_2.m
%--------------------------------------------------------------------------
% Project: Simulation of a hybrid system (bouncing ball)
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%   See also HYEQSOLVER, PLOTARC, PLOTARC3, PLOTFLOWS, PLOTHARC,
%   PLOTHARCCOLOR, PLOTHARCCOLOR3D, PLOTHYBRIDARC, PLOTJUMPS.
%   Copyright @ Hybrid Systems Laboratory (HSL),
%   Revision: 0.0.0.3 Date: 05/20/2015 3:42:00
%close all

global A M L T1 T2 T1_d T2_d q T_t d c

T1 = 0.2;
T2 = 1;
T1_d = 0;
T2_d = 0.2;
A=1;
M=1;
L = 1.1 - exp(-1);
q = 0;
T_t = 3;
d = 0.5;
c = 0.5;

% plant initial conditions
z_0 = 1;

% plant estimate initial conditions for

% Observer with delays
zhat_0 = 0;

% Observer without delays
zhat_nd_0 = 0;

% Auxilary observer variable for analysis
zhat_p_0 = 0;

% State measurment
y_0 = M*z_0;
y_m0 = 0;

% initial conditions for Clocks and Timers

tauP0 = 0;
tauO0 = 0;
taud0 = 2*T2_d+1;
tauN0 = T1+(T2-T1).*rand(1,1);

tP_m0 = 0;
tP_mp0 = 0;

% Protocol Timer
tauT0 = 0;

% Master-Slave Timer
tauM0 = 0;

% Slave-Master Timer
tauS0 = 0;

% Slave Memory Buffer
M_m0 = zeros(6,1);

% Master Memory Buffer
M_s0 = zeros(6,1);

% Logic State Variables
q0 = 0;
qP0 = 0;
p0 = 1;

% 1588 state initial condition
w0 = [tauT0; tauM0; tauS0; qP0; p0; M_m0; M_s0;];

% state initial condition
x0 = [z_0; zhat_0; tauP0; tauO0; tauN0; taud0; q0; y_0; y_m0; tP_m0; zhat_nd_0; zhat_p_0; tP_mp0; w0];

% simulation horizon
TSPAN=[0 10];
JSPAN = [0 1000];

% rule for jumps
% rule = 1 -> priority for jumps
% rule = 2 -> priority for flows
rule = 1;

options = odeset('RelTol',1e-6,'MaxStep',.01);

% simulate
[t,j,x] = HyEQsolver( @f,@g,@C,@D,...
    x0,TSPAN,JSPAN,rule,options,'ode45');

%%

% plot solution
figure(1) % position
clf
subplot(2,1,1), plotHarc(t,j,x(:,1));
grid on
ylabel('$z$','Interpreter','latex','FontSize',20)
subplot(2,1,2), plotHarc(t,j,x(:,2));
grid on
ylabel('$\hat{z}$','Interpreter','latex','FontSize',20)
xlabel('$t$','Interpreter','latex','FontSize',20)

% plot solution
figure(2) % position
clf
subplot(2,1,1), plotHarc(t,j,x(:,1) - x(:,2));
grid on
ylabel('$| \varepsilon |$','Interpreter','latex','FontSize',20)
subplot(2,1,2), plotHarc(t,j,x(:,7));
grid on
ylabel('$q$','Interpreter','latex','FontSize',20)
xlabel('$t$','Interpreter','latex','FontSize',20)

figure(3) % position
clf
subplot(2,1,1), plotHarc(t,j,x(:,3));
grid on
ylabel('$\tau_P$','Interpreter','latex','FontSize',20)
subplot(2,1,2), plotHarc(t,j,x(:,4));
grid on
ylabel('$\tau_O$','Interpreter','latex','FontSize',20)
xlabel('$t$','Interpreter','latex','FontSize',20)

figure(4) % position
clf
subplot(2,1,1), plotHarc(t,j,x(:,5));
grid on
ylabel('$\tau_N$','Interpreter','latex','FontSize',20)
subplot(2,1,2), plotHarc(t,j,x(:,6));
grid on
ylabel('$\tau_{\delta}$','Interpreter','latex','FontSize',20)
xlabel('$t$','Interpreter','latex','FontSize',20)

figure(5) % position
clf
subplot(2,1,1), plotHarc(t,j,x(:,8));
grid on
ylabel('$y$','Interpreter','latex','FontSize',20)
subplot(2,1,2), plotHarc(t,j,x(:,9));
grid on
ylabel('$y_{meas}$','Interpreter','latex','FontSize',20)
xlabel('$t$','Interpreter','latex','FontSize',20)

figure(6) % position
clf
plotHarc(t,j,x(:,1) - x(:,2));
grid on
ylabel('$|\varepsilon|$','Interpreter','latex','FontSize',50)
xlabel('$t,j$','Interpreter','latex','FontSize',50)