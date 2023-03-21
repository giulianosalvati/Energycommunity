# Sets
set J ordered;
set I ordered;

#Parameters
# Network import (da rivedere)

param Pnom_prosumer;

# HVAC data
param C;
param R;
param eta_c;
param eta_h;
param Pnom_hvac;

# PV generator data 

param Pnom_PV;

# Battery data
param Pnom_b;
param Eb;
param eta_ch;
param eta_dsc;

# PEV data

param Pnom_pev;
param Epev;
param eta_pev;

# ABP 

param Eabp{I};
param Tabp{I};
param Dabp{I};
param Pmax_abp{I};
param Pmin_abp{I};

# Control specifications data
param Tsp;
param Delta;
param SoCmin;
param SoCmax;
param Dts;

# PV shared consumption 

param tdSHPV;

param alpha;
param beta;

param P_PV_forecast{J};
param Pi_tot_forecast{J};
param P_ul_forecast{J};

param pe{J};
param c{J};
	
param Tk;
param SoCpevk;
param SoCk;

param t_abp_done{I};
param s_abp_k1{I};
param d_abp_k1{I};
param Eabp_done{I};
param Tabp_done{I};

param Tex_forecast{J};
param UR_hvac{J};
param UR_pev{J};
param UR_abp{J};
param SoCpev_obj;


# Variables

var Pcool{J} >= 0;
var Ph{J} >= 0;
var Pch{J} >= 0;
var Pdsc{J} >= 0;
var Ppev{J} >= 0;
var Pe{J} >= 0;
var Pc{J} >= 0;
var Pi_0{J} >= 0;
var Pabp{J,I} >= 0;
var Pabp_tot{J} >= 0;
var T{J};
var SoC{J};
var d_c{J} binary;
var d_h{J} binary;
var d_b{J} binary;
var d_pev{J} binary;
var d_ie{J} binary;
var d_abp{J,I} binary;
var s_abp{J,I} binary;
var t_abp{J,I} binary;



# Objective function

minimize total_cost : sum {j in J} (-tdSHPV*Dts*Pc[j]-pe[j]*Dts*Pe[j]+ c[j]*Dts*Pi_0[j]);

#Constraints for HVAC

subject to con_hvac_1 {j in J: ord(j)=1}:       T[j] == alpha*Tk-beta*R*(eta_c*Pcool[j]-eta_h*Ph[j])+beta*Tex_forecast[j];
subject to con_hvac_2 {j in J: ord(j)>1}:       T[j] == alpha*T[j-1]-beta*R*(eta_c*Pcool[j]-eta_h*Ph[j])+beta*Tex_forecast[j];
subject to con_hvac_3 {j in J}:                 T[j] <= (Tsp+Delta)*UR_hvac[j];
subject to con_hvac_4 {j in J}:                 T[j] >= (Tsp-Delta)*UR_hvac[j];
subject to con_hvac_5 {j in J}:                 Pcool[j]<= d_c[j]*Pnom_hvac;
subject to con_hvac_6 {j in J}:                 Ph[j]<= d_h[j]*Pnom_hvac;
subject to con_hvac_7 {j in J}:                 d_h[j] + d_c[j]<= UR_hvac[j];

#Constraints for battery

subject to con_battery_1  {j in J: ord(j)=1}:   SoC[j] == SoCk + Dts/Eb*(eta_ch*Pch[j]-1/eta_dsc*Pdsc[j]);
subject to con_battery_2  {j in J: ord(j)>1}:   SoC[j] == SoC[j-1] + Dts/Eb*(eta_ch*Pch[j]-1/eta_dsc*Pdsc[j]);
subject to con_battery_3  {j in J}:             SoC[j] <= SoCmax;
subject to con_battery_4  {j in J}:             SoC[j] >= SoCmin;
subject to con_battery_5  {j in J}:             Pch[j]<= d_b[j]*Pnom_b;
subject to con_battery_6  {j in J}:             Pdsc[j]<= (1-d_b[j])*Pnom_b;
subject to con_battery_7  {j in J}:				Pch[j] <= P_PV_forecast[j];   # Carica batteria deve essere rinnovabile


#Contraints for PEV
subject to con_pev_1:                           sum {j in J} (Dts*eta_pev*Ppev[j]) == Epev*(SoCpev_obj-SoCpevk);
subject to con_pev_2 {j in J}:                  Ppev[j]<= d_pev[j]*Pnom_pev; 
subject to con_pev_3 {j in J}:                  d_pev[j] <= UR_pev[j];

#Constraints for ABP
subject to con_abp_1 {i in I}:                            sum {j in J} (Dts*Pabp[j,i]) = Eabp[i]-Eabp_done[i]; 
subject to con_abp_2 {i in I}:                            sum {j in J} (d_abp[j,i]) = Tabp[i]-Tabp_done[i];
subject to con_abp_3 {j in J, i in I}:                    Pabp[j,i]<= d_abp[j,i]*Pmax_abp[i];
subject to con_abp_4 {j in J, i in I}:                    Pabp[j,i]>= d_abp[j,i]*Pmin_abp[i];
subject to con_abp_5 {j in J, i in I}:                    d_abp[j,i]+s_abp[j,i]<=1;
subject to con_abp_6 {j in J, i in I: ord(j)=1}:          d_abp_k1[i]-d_abp[j,i]-s_abp[j,i]<=0;
subject to con_abp_7 {j in J, i in I: ord(j)>1}:          d_abp[j-1,i]-d_abp[j,i]-s_abp[j,i]<=0;
subject to con_abp_8 {j in J, i in I: ord(j)=1}:          s_abp_k1[i]-s_abp[j,i]<=0;
subject to con_abp_9 {j in J, i in I: ord(j)>1}:          s_abp[j-1,i]-s_abp[j,i]<=0;
subject to con_abp_10 {j in J, i in I: ord(i)>1}:          d_abp[j,i]-s_abp[j,i-1]<=0;
subject to con_abp_11 {j in J, i in I: ord(i)>1}:         t_abp[j,i] == s_abp[j,i-1]-d_abp[j,i]-s_abp[j,i];
subject to con_abp_12 {i in I: ord(i)>1}:                 sum {j in J} (t_abp[j,i]) <= Dabp[i]-t_abp_done[i]; 
subject to con_abp_13 {j in J, i in I}:                   d_abp[j,i] <= UR_abp[j];
subject to con_abp_14 {j in J}:                           Pabp_tot[j] == sum {i in I} (Pabp[j,i]); 



#Constraints for system

subject to con_cer_1 {j in J}:                Pi_0[j] - Pe[j] == Pcool[j]+ Ph[j]+Ppev[j]+ P_ul_forecast[j]+ Pabp_tot[j] +Pch[j] - Pdsc[j] - P_PV_forecast[j];
subject to con_cer_2 {j in J}:                Pc[j] <= Pe[j];
subject to con_cer_3 {j in J}:                Pc[j] <= Pi_tot_forecast[j];
subject to con_cer_4 {j in J}:             	  Pi_0[j] <= d_ie[j]* Pnom_prosumer;
subject to con_cer_5 {j in J}:				  Pe[j] <= (1- d_ie[j])* Pnom_prosumer;
