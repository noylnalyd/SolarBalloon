%% Final Report
clear all;
clc;
format long g;

r = [10; 15; 20];
for i=1:length(r)

    % Atmospheric conditions
    alt = [5; 10; 15; 20; 25; 30];
    P_atm = [5.405e4; 2.65e4; 1.211e4; 5.529e3; 2.549e3; 1.197e3];
    rho_atm = [7.364e-1; 4.135e-1; 1.948e-1; 8.891e-2; 4.008e-2; 1.841e-2];
    T_atm = [-17.47; -49.9; -56.5; -56.5; -51.6; -46.64];
    M = 2.016;
    R = 0.0821;
        
    % H2 density
    rho_H2 = (P_atm./101300).*M./(R.*(T_atm + 273))
        
    % Total mass that can be lifted
    
    m_lift = rho_atm.*((4/3).*pi.*r(i)^3) - (rho_H2).*((4/3).*pi.*r(i)^3)
    m_lift = m_lift.*ones(6,1);
        
    % Cable mass
    linear_density = (31.43./1000)*(3280.84)*(1/2.2);
    L_cable = alt %km;
    m_cable = linear_density.*alt
        
    % Balloon mass
    rho_balloon = 916; %kg/m3
    A_balloon = 4*pi*r(i)^2;
    thickness_balloon = 0.002/100 %m;
    m_balloon = rho_balloon*A_balloon*thickness_balloon;
    m_balloon = m_balloon.*ones(6,1);
        
    % Solar panel mass
    m_payload = m_lift - m_cable - m_balloon;
    SF = 1.2;
    m_panel = m_payload./1.2
    A_panel = m_panel./11.66
        
    % Graph
    plot(alt,A_panel)
    hold on;
    xlabel('Altitude (km)','FontSize', 22);
    ylabel('Solar Panel Area (m^2)','FontSize', 22);
    xline(0);
    yline(0);
end
xlim([0,26]);
ylim([0,1800]);

qw{1} = plot(nan, 'Color', "[0 0.4470 0.7410]");
qw{2} = plot(nan, 'Color', "[0.8500 0.3250 0.0980]");
qw{3} = plot(nan, 'Color', "[0.9290 0.6940 0.1250]");
legend([qw{:}], {'20 m diameter','30 m diameter','40 m diameter'}, 'location', 'best')
hold on;

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






