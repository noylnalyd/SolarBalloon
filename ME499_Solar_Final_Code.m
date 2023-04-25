%% Final Report
clear all;
clc;
format long g;

% Atmospheric conditions
alt = 10000;
P_atm = 2.65e4;
rho_atm = 4.135e-1;
T_atm = -49.9;
g = 9.776;
M = 2.016;
R = 0.0821;

rho_H2 = (P_atm/101300)*M/(R*(T_atm + 273))

r = 15;
m = rho_atm*((4/3)*pi*r^3) - (rho_H2)*((4/3)*pi*r^3)



%% NASA study
clear all;
clc;
format long g;

% Atmospheric conditions
alt = 30000;
P_atm = 1.197e3;
rho_atm = 1.841e-2;
T_atm = -46.64;
g = 9.776;
M = 4;
R = 0.0821;

rho_H2 = (P_atm/101300)*M/(R*(T_atm + 273))

m = rho_atm*(121762.4403) - (rho_H2)*(121762.4403)


%%
load('Data001.mat','Data001');
x = Data001(:,1);
y = Data001(:,2);
plot(x,y);








