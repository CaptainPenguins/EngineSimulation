function [T_cc, Tad_Count] = ChamberTemperature(P_cc, n, T_ad_guess, Tad_Count)
%%%%% Created by J. Osborne, Oct 2013

%%% THIS PROGRAM WILL CALCULATE THE COMBUSTION CHAMBER TEMPERATURE OF A
%%% PARAFFIN/ALUMINUM HYBRID ROCKET ENGINE AS A FUNCTION OF
%%% COMBUSTION CHAMBER PRESSURE P_cc (Pa) AND PRODUCTS MOLAR FRACTIONS (n). 
%%%% THERMOPHYSICAL PROPERTIES OF PARAFFIN (APPROXIMATED AS C25H52) ARE 
%%%% TAKEN AT STANDARD ATMOSPHERIC CONDITIONS (PRESSURE = 1ATM) AND SCALING FACTORS
%%% ARE APPLIED TO THE TEMPERATURE VALUE TO OBTAIN IT AT THE COMBUSTION
%%% CHAMBER PRESSURES. CALCULATIONS ARE BASED ON THE CONSTANT-PRESSURE
%%% ADIABATIC FLAME TEMPERATURE

P_calc = 1*101325; %1atm in Pa for initial calculations

tolerance = 1e-5; %Convergence criteria

%%% Thermophysical Properties %%%

h_comb_paraffin = -41400 * (25*12.01 + 52*1.01); %Enthalpy of combustion of paraffin, kJ/kmol
h_comb_CO2 = -394000; %Enthalpy of combustion of CO2
h_comb_H2O = -286000; %Enthalpy of combustion of H2O

h_f_CO2 = -393546; %CO2 at 298 in kJ/kmol
h_f_H2O = -241845; %H2O at 298 in kJ/kmol
h_f_Al = 0; %Aluminum, elementary means 0 enthalpy of formation
h_f_N2O = 82050; %N2O at 298 in kJ/kmol
h_f_N2 = 0; %Elementary
h_f_Al2O3 = -1669800; %Aluminum oxide at 298 in kJ/mol

Cp_paraffin = 2.52*(25*12.01 + 52*1.01); %Heat capacity of Paraffin in kJ/kmol-K, NOTE: We need to have this value at multiple temperatures, but the NIST database is down at the time of writing this code
Cp_Al2O3 = 79.04; %Heat capacity of Aluminum Oxide in kJ/kmol-K, NOTE: We need to have this value at multiple temperatures, but the NIST database is down at the time of writing this code
Cp_Al = 880*26.98/1000; %Heat capacity of aluminum in kJ/kmol-k. NOTE: We need to have this value at multiple temperatures, but the NIST database is down at the time of writing this code
Cp_N2O = 892/1000*44.04; %Heat capacity of nitrous oxide in kJ/kmol-K. NOTE: We need to have this value at multiple temperatures, but the NIST database is down at the time of writing this code

BP = importdata('ParaffinBoilingPoint.txt'); 
Pcc_ref = BP.data(:,1); T_ref = BP.data(:,2); 

SH = importdata('SpecificHeat.txt'); 
Cp_T = SH.data(:,1); Cp_CO2_ref = SH.data(:,2); Cp_H2O_ref = SH.data(:,3); Cp_N2_ref = SH.data(:,4); 

T_i = interp1(Pcc_ref*101325, T_ref, P_calc) + 273;  %Evaporation temperature of paraffin at P_calc

% Calculate Enthalpy of Formation of Paraffin %

h_f_paraffin = (25 * h_comb_CO2 + 52/2*h_comb_H2O - h_comb_paraffin); %Enthalpy of formation of paraffin, kJ/kmol

H_reactants = 1*(h_f_paraffin + Cp_paraffin*(T_i - 298)) + n(1) * (h_f_Al + Cp_Al * (T_i - 298)) + n(2) * (h_f_N2O + Cp_N2O * (T_i - 298)); %Total enthalpy of reactants at evaporation temperature of paraffin

Q = n(3) * h_f_CO2 + n(4) * h_f_H2O + n(5) * h_f_N2 + n(6) * h_f_paraffin + n(7) * h_f_Al + n(8) * h_f_Al2O3; 
R = n(3) * interp1(Cp_T, Cp_CO2_ref, T_ad_guess) + n(4) * interp1(Cp_T, Cp_H2O_ref, T_ad_guess) + n(5) * interp1(Cp_T, Cp_N2_ref, T_ad_guess) + n(6) * Cp_paraffin + n(7) * Cp_Al + n(8) * Cp_Al2O3;

T_ad = (H_reactants - Q)/R + 298;

while abs((T_ad - T_ad_guess)/(T_ad_guess)) > tolerance
    
    Tad_Count = Tad_Count + 1;
    
    T_ad_guess = T_ad;
    
    Q = n(3) * h_f_CO2 + n(4) * h_f_H2O + n(5) * h_f_N2 + n(6) * h_f_paraffin + n(7) * h_f_Al + n(8) * h_f_Al2O3;
    
    R = n(3) * interp1(Cp_T, Cp_CO2_ref, T_ad_guess) + n(4) * interp1(Cp_T, Cp_H2O_ref, T_ad_guess) + n(5) * interp1(Cp_T, Cp_N2_ref, T_ad_guess) + n(6) * Cp_paraffin + n(7) * Cp_Al + n(8) * Cp_Al2O3;
    
    T_ad = (H_reactants - Q)/R + 298;
    
end

T_cc = (0.0373*log(P_cc/101325)+1)*T_ad; %Fudge factor to account for combustion chamber pressure, as T_ad was calculated based on standard pressure (1atm). Equation was taken from the relation between stoichiometric adiabatic flame temperature of C12H22 at different pressures



