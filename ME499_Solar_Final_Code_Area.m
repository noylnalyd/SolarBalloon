%% Solar Panel Area Calcs
clear all;
clc;
format long g;


% Scaling
load('Data001.mat')

Phoenix_Ground = 364; %W/m^2
Phoenix_Sunny_Ground = 406; %W/m^2
Altitude_Data = Data001(1:end,1);
Ground = Data001(1,2);
Radiation = Data001(1:end,2);
Radiation_Scaling = Radiation/Ground;
Phoenix_Scaled = Phoenix_Sunny_Ground*Radiation_Scaling;

figure(2)
plot(Altitude_Data,Phoenix_Scaled)

yyaxis left
xlabel('Altitude (km)','FontSize', 22);
ylabel('Irradiance (W/(m^2)','FontSize', 22);
xlim([5,20])

yyaxis right
plot(Altitude_Data,Phoenix_Scaled*365*24/10^6)

xlabel('Altitude (km)','FontSize', 22);
ylabel('Yearly Irradiance (MWh/(m^2)','FontSize', 22);
xlim([5,20])
ylim([4.03,4.73])

figure(3)
plot(Altitude_Data,Phoenix_Scaled/Phoenix_Ground*100 - 100)

xlabel('Altitude (km)','FontSize', 22);
ylabel('Increase Over Ground (%)','FontSize', 22);
xlim([5,20])


color = ["r-","b-","m-"]
counter = 1;
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
    plot(alt,A_panel,color(counter),'LineWidth',2)
    hold on;
    counter = counter+1;
    xlabel('Altitude (km)','FontSize', 22);
    ylabel('Solar Panel Area (m^2)','FontSize', 22);
    xline(0);
    yline(0);
    
    xlabel('Altitude (km)','FontSize', 22);
    ylabel('Solar Panel Area (m^2)','FontSize', 22);
    xline(0);
    yline(0);

end
xlim([0,26]);
ylim([0,2000]);

qw{1} = plot(nan, 'r-');
qw{2} = plot(nan, 'b-');
qw{3} = plot(nan, 'm-');
legend([qw{:}], {'20 m diameter','30 m diameter','40 m diameter',}, 'location', 'best')
hold on;
ax = gca;
ax.FontSize = 22;
