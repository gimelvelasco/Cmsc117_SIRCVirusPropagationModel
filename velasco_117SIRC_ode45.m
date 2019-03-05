function velasco_117SIRC_ode45()
%%
%Cmsc 117 FINAL PAPER
%%
%Solving SIRC ODE using ode45%
%%
%%%%%%%%%%%%
%2012-58922
%VELASCO, Gimel David F.
%May 23, 2016
%Cmsc 117
%%%%%%%%%%%%
%%
%%%%%%%%
%This matlab code for solving a system of differential equations uses the
%ode45
%The whole code is based on how ode45 is explained in
%http://www3.nd.edu/~nancy/Math20750/Demos/3dplots/dim3system.html
%%%%%%%%
%%
%CHANGE PARAM
%change the values of the ff variables for the system of ode.
b = 60;             %rate at which external computers are connected to the network
bet = 0.024;         %rate at which a susceptible computer can become infected
mu = 0.25;           %rate at which a computer is removed from the network
gam = 0.07;          %recovery rate of infected computers due to antivirus
sig = 10;           %aveg delay of the alert notification in virus infection
s0 = 80;            %initial number of susceptible computers
i0 = 40;             %initial number of infected computers
r0 = (b*bet)/mu*(mu+gam);   %basic reproduction number
c0 = 120;            %initial number of computers in the network
t_final = 1000;
%%
%Preset Initial Conditions
%20,0.01,0.2,0.4,5
%10,0.02,0.1,0.6,4
%Preset Initial Condition for KA301 Computer Lab
%8 0.01 0.4 0.8 2
%%
%%%%%%%
syms t x
%this block of codes is unused since the function dsolve cant solve the
%system
%eqn1 = diff(S) == b - bet*C*I - mu*S;
%eqn2 = diff(I) == bet*C*I - (mu + gam)*I;
%eqn3 = diff(R) == gam*I - mu*R;
%eqn4 = diff(C) == (1/sig)*(S - C);
%%
%eqn1 = diff(S) == b - bet*x(4)*x(2) - mu*x(1);
%eqn2 = diff(I) == bet*x(4)*x(2) - mu*x(2) - gam*x(2);
%eqn3 = diff(R) == gam*x(2) - mu*x(3);
%eqn4 = diff(C) == (1/sig)*x(1) - (1/sig)*x(4);
%%
%S is x(1)
%I is x(2)
%R is x(3)
%C is x(4)
%%
f = @(t,x)[b - bet*x(4)*x(2) - mu*x(1);bet*x(4)*x(2) - mu*x(2) - gam*x(2);gam*x(2) - mu*x(3);(1/sig)*x(1) - (1/sig)*x(4)];
[t,xa] = ode45(@(t,x) f(t,x),[0 t_final],[s0 i0 r0 c0]);    %S(0) I(0) R(0) C(0)
%%
%Note:
%Conditions for initial values S(0), I(0), R(0) and C(0)
%C <= b/mu
%S+I <= b/mu
%%%%Preset Initial Conditions
%50 0 0.83 100
%20 30 2.86 50
%Preset Initial Condition for KA301 Computer Lab
%15 0 ? 34
%%
%Plot Solution
plot(t,xa)
legend('Susceptible','Infected','Recovered','Computers')
xlabel('Time')
ylabel('SIRC')
r0                      %displays the initial value of R so that the user would know if the desired initial value of r is met
b_div_mu = b/mu        %displays the value of b/mu just to notify the user that the bounded condition of S + I and C is satisfied %just for user ergonomicity
%S Blue
%I Red
%R Yellow
%C Purple
end