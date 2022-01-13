%D Partridge, University of Exeter, Sept 2021

clc
clear all 
close all 

Pars = []

Pars(1,5) = 0.35   %Updraft velocity (ms-1)
Pars(1,9) = 285    %Temperature (K)

Pars(1,1) = 1e+008 %Number mode 1 (m^-3)
Pars(1,6) = 5e+007 %Number mode 2 (m^-3)

Pars(1,2) = 0.0200 %Mean radius mode 1 (microns)
Pars(1,7) = 0.0700 %Mean radius mode 2 (microns)

Pars(1,3) = 1.4000 %GSD mode 1 
Pars(1,8) = 1.4000 %GSD mode 2 

Pars(1,4) = 0.9000 %Soluble volume fraction (soluble component = (nh4)hso4, insoluble = BC)

%Compute effective kappa
density_sol   = 1770;       %Kg/m3 (nh4)2so4
density_insol = 2000;       %Kg/m3 BC
molar_mass_sol = 0.13216;   %(nh4)2so4
ionic_dissociation_sol = 3; %(nh4)2so4
rho_water = 1000;
Mv_water  = 0.018;

volfrac_sol = Pars(1,7);
volfrac_insol = 1- Pars(1,7);
ns_sol = (ionic_dissociation_sol * 1 / (molar_mass_sol/density_sol));
ns = ns_sol;
kappa_sol = ns * Mv_water / rho_water;
kappa = kappa_sol * volfrac_sol;

% Write the parameters to a structure (all in SI units)
dt = .25;                      %Timestep
input = struct(...
    'n_bins', 250, ...           %Number of bins
    'T',  Pars(1,9), ...         %Initial Temperature (K)
    'RH', .999, ...              %Initial Relative Humidity (1)
    'p',  1.00e+05, ...          %Initial Pressure (Pa)
    'w',  Pars(1,5), ...         %Updraft velocity (m/s)
    'kappa', kappa, ... 
    'n_tot', [Pars(1,1), Pars(1,6)], ...  %N mode 1,2
    'meanr', [Pars(1,2) * 1e-6, Pars(1,7) * 1e-6], ...  %Mean radius mode 1,2
    'gstdv', [Pars(1,3), Pars(1,8)], ...  %GSD mode 1,2
    'dt', dt, ... 
    'nt', ceil(1000 / Pars(1,5) / dt), ... % max output height (but actually stops at S_max)
    'sigma', 0.072 ... % surface tension coefficient for water
);

% Run the forward model 
output = drops_py(input);
   
% Define the return vector for comparison with the measurements
N_Parcel = [output.N_act, output.S_max]; 
N_Parcel = N_Parcel(:);

Pars
input
output
