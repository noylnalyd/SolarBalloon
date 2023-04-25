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

% H2 density
rho_H2 = (P_atm/101300)*M/(R*(T_atm + 273))

% Total mass that can be lifted
r = 15;
m_lift = rho_atm*((4/3)*pi*r^3) - (rho_H2)*((4/3)*pi*r^3)

% Cable mass
linear_density = (31.43/1000)*(3280.84)*(1/2.2);
L_cable = 10 %km;
m_cable = linear_density*L_cable

% Balloon mass
rho_balloon = 916; %kg/m3
A_balloon = 4*pi*r^2;
thickness_balloon = 0.002/100 %m;
m_balloon = rho_balloon*A_balloon*thickness_balloon;

% Solar panel mass
m_payload = m_lift - m_cable - m_balloon;
SF = 1.2;
m_panel = m_payload/1.2
A_panel = m_panel/11.66


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

xlabel('Altitude (km)','FontSize', 22);
ylabel('Irradiance (W/(m^2)','FontSize', 22);

