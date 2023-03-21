% Economics and Smart Management of Electrical Systems Exam
% Project 2

%% Problem data
% Network import 
parms.Pnom_prosumer= 30; % Pnom prosumer kW

% HVAC data
parms.R = 5; %walls thermal resistance [°C/kW]
parms.C = 0.654; % thermal capacitance of air [kWh/°C]
parms.eta_c = 1.9; % cooling system COP
parms.eta_h = 0.85; % hating system efficiency
parms.Pnom_hvac = 2; % nominal power of the cooling system [kW]

% PV generator data
PV_panel_Area = 1.6275; % Single panel area [m2]
npanels = 60; % number of panels
PV_Area =  npanels*PV_panel_Area; % Total PV area
eta_pv = 0.222; % PV system efficiency
parms.Pnom_PV = PV_Area*eta_pv; % Nominal power at 1kW/m2

% % Battery data
% parms.Pnom_b = 35; % Battery nominal power [kW]
% parms.Eb = 35; % Battery capaciy [kWh]
% parms.eta_ch = 0.95; % Battery charging efficiency [pu]
% parms.eta_dsc = 0.95; % Battery discharging efficiency [pu]

% PEV data
parms.Pnom_pev = 22; % Rated charging power [kW]
parms.Epev = 85; % PEV battery capaciy [kWh]
parms.eta_pev = 0.9; % PEV battery charging efficiency [pu]

% ABP data
parms.n_abp_phases = 4; % number of phases [-]
parms.Eabp = [0.11 0.2 0.07 0.8]'; % phases energy [kWh]
parms.Tabp = [1 1 1 1]'; % phase time execution [time steps]
parms.Dabp = [2 2 2 2]'; % maximal delay among phases [time steps]
parms.Pmax_abp = [0.15 1.6 0.15 1.6]'; % maximal power for each phase [kW]
parms.Pmin_abp = [0 0 0 0]'; % minimal power for each phase [kW]

% Control specifications data
parms.Tsp = 22; % temperature set-point [°C]
parms.Delta  = 2; % temperature regulation tolerance [°C] 
parms.SoCmin = 0.2; % desired minimal battery SoC [pu]
parms.SoCmax = 0.8; % desired maximal battery SoC [pu]
parms.Dts = 1; % control sampling time [h]

% Selection of days
months_days = [31 28 31 30 31 30 31 31 30 31 30 31];
m1 = 6; %starting month
d1 = 10; %staring day
m2 = 6; %end month
d2 = 17; %end day
firstday = 3; % 1=monday etc.
idxs = (sum(months_days(1:m1-1))+(d1-1))*24+1:(sum(months_days(1:m2-1))+d2)*24+1;

fun = 'compute_control_step';
nargin(fun)

% External temperature data 
load t_ex_rome_campus_bio_medico_january_july_2021
% [month,day,hour,forecasted temperature (°C), actual temperature (°C)]
T_ex = T_ex(idxs,:);

% Solar irradiation data
load Ir_rome_campus_bio_medico_january_june_2021
% [month,day,hour,forecasted Ir (°C), actual Ir (°C)]
Ir = Ir(idxs,:);

load residential_load_january_october_30 % nominal load 30kW
% [hour,forecasted load (kWh), actual load (kWh)]
Pi_tot = Pi_tot(idxs,:);

% Prices data 
load  pun_january_november_2021
% [date,hour,price [€/MWh]]
pun = pun(idxs,3)/1000; % exported power price[€/kWh]

parms.tdSHPV = 119/1000; % PV shared-consumption incentive (tr+tMiSE) [€/kWh]

% Imported energy cost (ARERA price)
load energy_prices_january_july_2021
priceF1 = prices(m1,1); % F1 price [eur/kWh]
priceF2 = prices(m1,2); % F2 price [eur/kWh]
priceF3 = prices(m1,3); % F3 price [eur/kWh]

load Pul_162_days.mat

% assuming working days
cday_lv = [ones(6,1)*priceF3;
        priceF2;
        ones(11,1)*priceF1;
        ones(4,1)*priceF2;
        ones(2,1)*priceF3];
cday_s = [ones(6,1)*priceF3;
          ones(17,1)*priceF2;
          priceF3];
cday_d = ones(24,1)*priceF3;
cweek = [repmat(cday_lv,5,1);cday_s;cday_d];

cweek = [cweek(24*(firstday-1)+1:end);cweek(1:24*(firstday-1))];

c = repmat(cweek,2,1);

%% Control parameters
parms.T = 24; % control time horizon
%% Simulation parameters
Tf = 24*7; % Simulation Time
T0 = parms.Tsp-1; % Initial temperature
%SoC0 = 0.5; % Initial battery state of charge

%% Offline computation
parms.alpha = exp(-parms.Dts/(parms.R*parms.C));
parms.beta = 1 - parms.alpha;

%% Initialization
Pc = zeros(Tf,1);        % shared PV power [kW] 
Pe = zeros(Tf,1);        % exported PV power [kW] 

Ppv = zeros(Tf,1);       % PV generation [kW]
Ppev = zeros(Tf,1);      % Pev Charging Power [kW]

Pabp = zeros(Tf,1);      % Pabp Power Consumption [kW]

Pi_0 = zeros(Tf,1);      % total prosumer load [kW]
Pi = zeros(Tf,1);        % total consumer load [kW]
P_ul = zeros(Tf,1);      % Uncontrollable load Power [kW]

Pcool = zeros(Tf,1);     % air cooling power [kW] 
Ph = zeros(Tf,1);        % air heating power [kW]

T = zeros(Tf+1,1);       % Internal air temperature [°C]
T(1) = T0;               % Initial Internal air temperature [°C]

% SoC = zeros(Tf+1,1);     % Battery state of charge [°C]     
% SoC(1) = SoC0;           % Initial Battery state of charge [°C]

SoCpev = zeros(Tf+1,1);  % PEV Battery state of charge [°C]
SoCpev(1) = 1;           % Initial PEV Battery state of charge [°C]

UR_hvac = zeros(Tf,1);   % User Requirements for HVAC
UR_pev = zeros(Tf,1);    % User Requirements for PEV
UR_abp = zeros(Tf,1);    % User Requirements for ABP 

d1 = 0;                  % day counter for PEV (to define user requirements)
d2 = 0;                  % day counter for ABP (to define user requirements)

rng(18)

% intialize variables to manage ABP
abp_varsk.tabp_donek =  zeros(parms.n_abp_phases,1); 
abp_varsk.s_abpk = zeros(parms.n_abp_phases,1); 
abp_varsk.d_abpk = zeros(parms.n_abp_phases,1);
abp_varsk.Eabp_donek =  parms.Eabp; 
abp_varsk.Tabp_donek =  parms.Tabp; 

k_start_pev = randi(5)+18; % first pev recharging starting time
k_start_abp = randi(5)+18; % abp pev recharging starting time

for k=1:Tf

    % get measurement
    Tk = T(k);
    SoCpevk = SoCpev(k);
    
    % get forecasts

    ck = c(rem(k-1,24*7)+1:rem(k-1,24*7)+1+parms.T-1);
    pek = pun(k:k+parms.T-1);
    T_ex_forecast_k = [T_ex(k,5);
                       T_ex(k+1:k+parms.T-1,4)];
    Pi_tot_forecast_k = Pi_tot(k:k+parms.T-1,2);
    P_PV_forecast_k = parms.Pnom_PV*Ir(k:k+parms.T-1,4);
    P_ul_forecast_k= P_u_loads(k:k+parms.T-1,2);
    
    % gest user requirements
    UR_hvac_k = ones(parms.T,1);
    UR_pev_k = zeros(parms.T,1);
    SoCpev_objk = 1;

    if k >= k_start_pev+d1*24
        UR_pev_k(1:k_start_pev+d1*24+12-k) = ones(k_start_pev+12+d1*24-k,1); % complete the cycle in 12 hours
        if k == k_start_pev+d1*24
            SoCpev(k) = rand;
            SoCpevk = SoCpev(k);
        end
        if k_start_pev+12+d1*24-k == 0
            d1 = d1+1;
            k_start_pev = randi(5)+18;
        end
    end

    UR_abp_k = zeros(parms.T,1);
    if k >= k_start_abp+d2*24
        UR_abp_k(1:k_start_abp+d2*24+12-k) = ones(k_start_abp+12+d2*24-k,1); % complete the cycle in 12 hours
        if k==k_start_abp+d2*24
            abp_varsk.Eabp_donek = zeros(parms.n_abp_phases,1);
            abp_varsk.Tabp_donek = zeros(parms.n_abp_phases,1); 
        end
        if k_start_abp+12+d2*24-k == 0
            d2 = d2+1;
            k_start_abp = randi(5)+18;
        end
    end
    if abp_varsk.Tabp_donek(end)==parms.Tabp(end) % abp cycle termined 
        abp_varsk.tabp_donek = zeros(parms.n_abp_phases,1); 
        abp_varsk.s_abpk = zeros(parms.n_abp_phases,1);
        abp_varsk.d_abpk = zeros(parms.n_abp_phases,1); 
    end
    
   
    
    % compute control 
    [Pek,hatPck,Pi_0k,Ppev(k),Pcool(k),Ph(k),Pabp(k),abp_varsk] = compute_control_step(parms,pek,ck,Pi_tot_forecast_k,T_ex_forecast_k,P_PV_forecast_k,P_ul_forecast_k,Tk,SoCk,SoCpevk,UR_hvac_k,UR_pev_k,SoCpev_objk,UR_abp_k,abp_varsk);
    
    Ppv(k) = parms.Pnom_PV*Ir(k,5);

    % simulate real system
    Pi(k) = Pi_tot(k,3);
    P_ul(k)= P_u_loads(k,3);
    Pexk = Ppv(k)- Pcool(k) - Ph(k)-Ppev(k)- P_ul(k)-Pabp(k);

    Pi_0(k)= -min(0,Pexk);
    Pe(k)= max(0,Pexk);
    Pc(k)=min(Pe(k),Pi(k));

    T(k+1) = parms.alpha*T(k)-parms.beta*parms.R*(parms.eta_c*Pcool(k)-parms.eta_h*Ph(k))+parms.beta*T_ex(k,5);
    SoCpev(k+1) = SoCpev(k) + parms.Dts/parms.Epev*parms.eta_pev*Ppev(k);
    
    % update ABP vars
    abp_varsk.Tabp_donek =  abp_varsk.Tabp_donek +  abp_varsk.d_abpk;
    abp_varsk.tabp_donek =  abp_varsk.tabp_donek +  abp_varsk.t_abpk;
    abp_varsk.Eabp_donek =  abp_varsk.Eabp_donek + parms.Dts*Pabp(k)*abp_varsk.d_abpk;
    
    % save current user requirements
    UR_hvac(k) = UR_hvac_k(1);
    UR_pev(k) = UR_pev_k(1);
    UR_abp(k) = UR_abp_k(1);
end



%% Plot results
close all

% Temperatures
figure(1)
subplot(3,1,1)
plot(0:Tf,ones(Tf+1,1)*parms.Tsp,':','LineWidth',2)
hold on
plot(0:Tf,ones(Tf+1,1)*(parms.Tsp+parms.Delta),'LineWidth',1.5)
plot(0:Tf,ones(Tf+1,1)*(parms.Tsp-parms.Delta),'LineWidth',1.5)
plot(0:Tf,T,'k','LineWidth',1.5)
plot(0:Tf-1,T_ex(1:Tf,5),'--','LineWidth',1.5)
plot(0:Tf-1,T_ex(1:Tf,4),':','LineWidth',1.5)
xlabel('Time [h]')
ylabel('Temperature [°C]')
xlim([0 Tf])
grid on
legend('Set-point','Max Temp','Min Temp','Internal Temp','External Temp','Forecasted External Temp')

% Power HVAC
subplot(3,1,2)
stairs(0:Tf-1,Pcool,'LineWidth',2)
hold on
stairs(0:Tf-1,Ph,':','LineWidth',1.5)
plot(0:Tf-1,parms.Pnom_hvac*ones(Tf,1),'LineWidth',1.5)
xlabel('Time [h]')
ylabel('Power [kW]')
xlim([0 Tf-1])
grid on
legend('Cooling Power','Heating Power','P^{nom}')

% UR HVAC
subplot(3,1,3)
stairs(0:Tf-1,UR_hvac,'LineWidth',1.5)
xlabel('Time [h]')
xlim([0 Tf-1])
grid on
legend('HVAC User Requirements')


% PEV SOC
figure(2)
subplot(3,1,1)
stairs(0:Tf,SoCpev.*[0;UR_pev],'LineWidth',1.5)
xlabel('Time [h]')
ylabel('SoC_{pev} [pu]')
xlim([0 Tf])
ylim([-0.1 1.1])
grid on

% PEV Power
subplot(3,1,2)
stairs(0:Tf-1,Ppev,'LineWidth',1.5)
hold on
plot(0:Tf-1,parms.Pnom_pev*ones(Tf,1),'LineWidth',1.5)
xlabel('Time [h]')
ylabel('Power [kW]')
xlim([0 Tf-1])
grid on
legend('PEV Recharging Power','Charging limit')

% UR PEV
subplot(3,1,3)
stairs(0:Tf-1,UR_pev,'LineWidth',1.5)
xlabel('Time [h]')
xlim([0 Tf-1])
grid on
legend('PEV User Requirements')


% Powers
figure(3)
hold on
plot(0:Tf-1,-(Pcool+Ph),'LineWidth',1.5)
plot(0:Tf-1,Ppv,'LineWidth',2.5)
plot(0:Tf-1,Pc,'LineWidth',1.5)
plot(0:Tf-1,-Ppev,'LineWidth',1.5)
plot(0:Tf-1,-Pabp,'LineWidth',1.5)
plot(0:Tf-1,-Pi_0,'k','LineWidth',2,'Color','g')
plot(0:Tf-1,Pe,'k','LineWidth',2,'Color','b')
grid on
box on
legend('HVAC Power Consumption','PV Generation','Power shared','PEV charging power','ABP Power Consumption','Imported power','Exported power')
xlim([0 Tf-1])
xlabel('Time [h]')
ylabel('Power [kW]')

% ABP Power
figure(4)

subplot(2,1,1)
stairs(0:Tf-1,Pabp,'LineWidth',1.5)
xlabel('Time [h]')
ylabel('Power [kW]')
xlim([0 Tf-1])
grid on
legend('ABP Power Consumption')

% UR ABP
subplot(2,1,2)
stairs(0:Tf-1,UR_abp,'LineWidth',1.5)
xlabel('Time [h]')
xlim([0 Tf-1])
grid on
legend('ABP User Requirements')

% % Costs
figure(5)

subplot(3,1,1)
stairs(0:Tf-1,-Pi.*c(1:Tf)*parms.Dts,'LineWidth',2.5)
hold on
stairs(0:Tf-1,Pc*parms.tdSHPV*parms.Dts,'LineWidth',2.5)
hold on
stairs(0:Tf-1,Pe.*pun(1:Tf)*parms.Dts,'LineWidth',2.5)
hold on
stairs(0:Tf-1,-Pi_0.*c(1:Tf)*parms.Dts+Pc*parms.tdSHPV*parms.Dts+Pe.*pun(1:Tf)*parms.Dts,'k','LineWidth',1);
xlim([0 Tf-1])
grid on
ylabel('Hourly Energy Cost [€]')
legend('Imported energy cost','Self-consumed CER Power','Exported energy revenue', 'Total savings')


subplot(3,1,2)
stairs(0:Tf-1,cumsum(-Pi.*c(1:Tf)*parms.Dts+Pc*parms.tdSHPV*parms.Dts+Pe.*pun(1:Tf)*parms.Dts),'LineWidth',1.5)
xlim([0 Tf-1])
grid on
ylabel('Total Energy Cost [€]')

subplot(3,1,3)
stairs(0:Tf-1,c(1:Tf),'LineWidth',1.5)
hold on
stairs(0:Tf-1,ones(Tf,1)*parms.tdSHPV,'LineWidth',1.5)
xlim([0 Tf-1])
grid on
xlabel('Time [h]')
ylabel('[€/kWh]')
legend('Energy Price','Self-consumed PV tariff discount')


%% MPC step
function [Pek,Pck,Pi_0k,Ppevk,Pcoolk,Phk,Pabpk,abp_varsk] = compute_control_step(parms,pek,ck,Pi_tot_forecast_k,T_ex_forecast_k,P_PV_forecast_k,P_ul_forecast_k,T_k,SoCpev_k,UR_hvac_k,UR_pev_k,SoCpev_objk,UR_abp_k,abp_varsk1)
%% AMPL SETUP
setupAMPL %MATLAB path setup
ampl = AMPL('./AMPL'); %open AMPL session
ampl.setOption('solver', 'cplex'); %set the solver
ampl.setOption('log_file', 'logfile.log');%make log file
ampl.read('Progetto2nob.mod') %read the ampl file (.mod)

%Communication of parameters to AMPL model

% time index
j = (1:parms.T)';
j = num2cell(j);
J = ampl.getSet('J');
J.setValues(j);


% abp phases index
i = (1:parms.n_abp_phases)';
i = num2cell(i);
I = ampl.getSet('I');
I.setValues(i);


% constant parameters
%Pnom_i = ampl.getParameter('Pnom_i');
%Pnom_i.setValues(parms.Pnom_i);
C = ampl.getParameter('C');
C.setValues(parms.C);

Pnom_prosumer= ampl.getParameter('Pnom_prosumer');
Pnom_prosumer.setValues(parms.Pnom_prosumer);

Pnom_hvac = ampl.getParameter('Pnom_hvac');
Pnom_hvac.setValues(parms.Pnom_hvac);

alpha = ampl.getParameter('alpha');
alpha.setValues(parms.alpha);

beta = ampl.getParameter('beta');
beta.setValues(parms.beta);

eta_c = ampl.getParameter('eta_c');
eta_c.setValues(parms.eta_c);

eta_h = ampl.getParameter('eta_h');
eta_h.setValues(parms.eta_h);

R = ampl.getParameter('R');
R.setValues(parms.R);


Tsp = ampl.getParameter('Tsp');
Tsp.setValues(parms.Tsp);

Delta = ampl.getParameter('Delta');
Delta.setValues(parms.Delta);

SoCmax = ampl.getParameter('SoCmax');
SoCmax.setValues(parms.SoCmax);

SoCmin = ampl.getParameter('SoCmin');
SoCmin.setValues(parms.SoCmin);

Pnom_pev = ampl.getParameter('Pnom_pev');
Pnom_pev.setValues(parms.Pnom_pev);

Epev = ampl.getParameter('Epev');
Epev.setValues(parms.Epev);

SoCpev_obj = ampl.getParameter('SoCpev_obj');
SoCpev_obj.setValues(SoCpev_objk);

eta_pev = ampl.getParameter('eta_pev');
eta_pev.setValues(parms.eta_pev);

Eabp = ampl.getParameter('Eabp');
Eabp.setValues(parms.Eabp);

Tabp = ampl.getParameter('Tabp');
Tabp.setValues(parms.Tabp);

Dabp = ampl.getParameter('Dabp');
Dabp.setValues(parms.Dabp);

Pmax_abp = ampl.getParameter('Pmax_abp');
Pmax_abp.setValues(parms.Pmax_abp);

Pmin_abp = ampl.getParameter('Pmin_abp');
Pmin_abp.setValues(parms.Pmin_abp);

tdSHPV = ampl.getParameter('tdSHPV');
tdSHPV.setValues(parms.tdSHPV);

Dts = ampl.getParameter('Dts');
Dts.setValues(parms.Dts);

% time varying parameters
pe = ampl.getParameter('pe');
pe.setValues(pek);

c = ampl.getParameter('c');
c.setValues(ck);

P_ul_forecast= ampl.getParameter('P_ul_forecast');
P_ul_forecast.setValues(P_ul_forecast_k);

Pi_tot_forecast = ampl.getParameter('Pi_tot_forecast');
Pi_tot_forecast.setValues(Pi_tot_forecast_k);

Tex_forecast = ampl.getParameter('Tex_forecast');
Tex_forecast.setValues(T_ex_forecast_k);

P_PV_forecast = ampl.getParameter('P_PV_forecast');
P_PV_forecast.setValues(P_PV_forecast_k);

Tk  = ampl.getParameter('Tk');
Tk.setValues(T_k);


SoCpevk = ampl.getParameter('SoCpevk');
SoCpevk.setValues(SoCpev_k);
 
UR_hvac = ampl.getParameter('UR_hvac');
UR_hvac.setValues(UR_hvac_k);

UR_pev = ampl.getParameter('UR_pev');
UR_pev.setValues(UR_pev_k);

UR_abp  = ampl.getParameter('UR_abp');
UR_abp.setValues(UR_abp_k);

s_abp_k1  = ampl.getParameter('s_abp_k1');
s_abp_k1.setValues(abp_varsk1.s_abpk);

d_abp_k1  = ampl.getParameter('d_abp_k1');
d_abp_k1.setValues(abp_varsk1.d_abpk);

Tabp_done  = ampl.getParameter('Tabp_done');
Tabp_done.setValues(abp_varsk1.Tabp_donek);

t_abp_done  = ampl.getParameter('t_abp_done');
t_abp_done.setValues(abp_varsk1.tabp_donek);

Eabp_done  = ampl.getParameter('Eabp_done');
Eabp_done.setValues(abp_varsk1.Eabp_donek);


% AMPL SOLUTION 
ampl.solve;

% Get control trajectory
Pc=ampl.getVariable('Pc'); 
Pc=Pc.getValues;
Pc=Pc.getColumnAsDoubles('Pc.val');

Pe=ampl.getVariable('Pe'); 
Pe=Pe.getValues;
Pe=Pe.getColumnAsDoubles('Pe.val');

Pcool=ampl.getVariable('Pcool'); 
Pcool=Pcool.getValues;
Pcool=Pcool.getColumnAsDoubles('Pcool.val');

Ph=ampl.getVariable('Ph'); 
Ph=Ph.getValues;
Ph=Ph.getColumnAsDoubles('Ph.val');


Pi_0=ampl.getVariable('Pi_0'); 
Pi_0=Pi_0.getValues;
Pi_0=Pi_0.getColumnAsDoubles('Pi_0.val');

Ppev=ampl.getVariable('Ppev'); 
Ppev=Ppev.getValues;
Ppev=Ppev.getColumnAsDoubles('Ppev.val');

Pabp=ampl.getVariable('Pabp_tot'); 
Pabp=Pabp.getValues;
Pabp=Pabp.getColumnAsDoubles('Pabp_tot.val');

s_abp=ampl.getVariable('s_abp'); 
s_abp=s_abp.getValues;
s_abp=s_abp.getColumnAsDoubles('s_abp.val');
s_abp=reshape(s_abp,parms.n_abp_phases,parms.T); % Consider that the order is inverted as the usual assumed for matrices

d_abp=ampl.getVariable('d_abp'); 
d_abp=d_abp.getValues;
d_abp=d_abp.getColumnAsDoubles('d_abp.val');
d_abp=reshape(d_abp,parms.n_abp_phases,parms.T); 

t_abp=ampl.getVariable('t_abp'); 
t_abp=t_abp.getValues;
t_abp=t_abp.getColumnAsDoubles('t_abp.val');
t_abp=reshape(t_abp,parms.n_abp_phases,parms.T); 


% Apply receding horizon principle

Pcoolk = Pcool(1); 
Phk = Ph(1); 
Pi_0k = Pi_0(1); 
Ppevk = Ppev(1); 
Pck = Pc(1); 
Pek = Pe(1);
Pabpk = Pabp(1);
abp_varsk = abp_varsk1;
abp_varsk.s_abpk = s_abp(:,1); 
abp_varsk.d_abpk = d_abp(:,1);  
abp_varsk.t_abpk = t_abp(:,1);


ampl.close(); % close the AMPL Session

end





