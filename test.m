% D Partridge, University of Exeter, Sept 2021
% S Arabas, University of Illinois & Jagiellonian University

clc
clear all 
close all 


density_sol   = 1770;       %Kg/m3 (nh4)2so4
molar_mass_sol = 0.13216;   %(nh4)2so4
ionic_dissociation_sol = 3; %(nh4)2so4
rho_water = 1000;
Mv_water  = 0.018;
ns = (ionic_dissociation_sol * 1 / (molar_mass_sol / density_sol));
kappa_sol = ns * Mv_water / rho_water;
dt = .25;                      %Timestep

%------------------------------------------------------------------------
%Define input variables for bi-modal aerosol having internal mixture.

Pars = zeros(1, 13);

%Mode 1
Pars(1,1) = 1e8; % Number concentration of mode 1 (m^-3)
Pars(1,2) = .02; % Mean radius of mode 1 (microns)
Pars(1,3) = 1.4; % Geometric standard deviation of mode 1 

%Mode 2
Pars(1,4) = 5e7; % Number concentration of mode 2 (m^-3)
Pars(1,5) = .07; % Mean radius of mode 2 (microns)
Pars(1,6) = 1.4; % Geometric standard deviation of mode 2

%Chemistry: Volume fractions. These (in order) in ICPM are set as: h2so4, (nh4)hso4, (nh4)2so4, OC, BC, dust and seasalt
%We represent chemistry as an internal mixture of 1 soluble compound, and 1 insoluble (BC).

Pars(1,7) = 0.07;  % Volume fraction of soluble compound 1: (nh4)2so4
Pars(1,8) = 1;     % Mass accomodation coefficient
Pars(1,9) = 72.8;  % Surface tension

%Meteorology
Pars(1,10) = 0.35; % Vertical velocity (ms^-1)
Pars(1,11) = 285;  % Cloud base temperature (K)
Pars(1,12) = 1e5;  % Cloud base pressure (Pa)
Pars(1,13) = 99.9; % Initial parcel RH (%)

%------------------------------------------------------------------------

input = struct(...
    'n_bins', 250, ...           %Number of bins
    'T',  Pars(1,11), ...        %Initial Temperature (K)
    'RH', Pars(1,13) / 100, ...  %Initial Relative Humidity (1)
    'p',  Pars(1,12), ...        %Initial Pressure (Pa)
    'w',  Pars(1,10), ...        %Updraft velocity (m/s)
    'kappa', {{kappa_sol * Pars(1,7)}}, ... 
    'n_tot', {{[Pars(1,1), Pars(1,4)]}}, ...  %N mode 1,2 [m^{-3}] at initial density
    'meanr', {{[Pars(1,2) * 1e-6, Pars(1,5) * 1e-6]}}, ...  %Mean radius mode 1,2 [m]
    'gstdv', {{[Pars(1,3), Pars(1,6)]}}, ...  %GSD mode 1,2
    'dt', dt, ... 
    'nt', ceil(1000 / Pars(1,10) / dt), ... % max output height (but actually stops at S_max)
    'sigma', Pars(1,9) / 1000, ... % surface tension coefficient for water [J/m2]
    'MAC', Pars(1,8) ...
)
output = drops_py(input)
   
%------------------------------------------------------------------------
%Define input variables for 4 external modes (each mode represented by internal aerosol mixture).
%Prior range for each parameter to come from UKESM1 GCM. 
%Each of these variables will be defined prior to the line %Initialise PySDM"

Pars = zeros(1, 22);

%Soluble Mode 1
Pars(1,1) = 1.1e8;   % Number concentration of mode 1 (m^-3)
Pars(1,2) = .4;    % Mean radius of mode 1 (microns)
Pars(1,3) = 1.1;   % Geometric standard deviation of mode 1 

%Soluble Mode 2
Pars(1,4) = 1.2e8;   % Number concentration of mode 2 (m^-3)
Pars(1,5) = .3;    % Mean radius of mode 2 (microns)
Pars(1,6) = 1.2;   % Geometric standard deviation of mode 2

%Soluble Mode 3
Pars(1,7) = 1.3e8;   % Number concentration of mode 3 (m^-3)
Pars(1,8) = .2;    % Mean radius of mode 3 (microns)
Pars(1,9) = 1.3;   % Geometric standard deviation of mode 3

%Soluble Mode 4
Pars(1,10) = 1.4e8;  % Number concentration of mode 4 (m^-3)
Pars(1,11) = .1;   % Mean radius of mode 4 (microns)
Pars(1,12) = 1.4;  % Geometric standard deviation of mode 4

%Chemistry: Volume fractions. These (in order) in ICPM are set as: h2so4, (nh4)hso4, (nh4)2so4, OC, BC, dust and seasalt
%We represent chemistry as an internal mixture of 1 soluble compound, and 1 insoluble (BC) for each mode.
%I may keep this simple and use the same fraction for each modes, but it is
%easier to decide that later and still have the parameters set up as
%follows now.

Pars(1,13) = .1;     % Volume fraction of soluble compound 1 in mode 1: (nh4)2so4
Pars(1,14) = .2;     % Volume fraction of soluble compound 1 in mode 2: (nh4)2so4
Pars(1,15) = .3;     % Volume fraction of soluble compound 1 in mode 3: (nh4)2so4
Pars(1,16) = .4;     % Volume fraction of soluble compound 1 in mode 4: (nh4)2so4

Pars(1,17) = 1;      % Mass accomodation coefficient
Pars(1,18) = 72.8;   % Surface tension

%Meteorology
Pars(1,19) = 1.;     % Vertical velocity (ms^-1)
Pars(1,20) = 290.;   % Cloud base temperature (K)
Pars(1,21) = 90000.; % Cloud base pressure
Pars(1,22) = 99.;    % Initial parcel RH (%)

%------------------------------------------------------------------------

input = struct(...
    'n_bins', 250, ...           %Number of bins
    'T',  Pars(1,20), ...        %Initial Temperature (K)
    'RH', Pars(1,22) / 100, ...  %Initial Relative Humidity (1)
    'p',  Pars(1,21), ...        %Initial Pressure (Pa)
    'w',  Pars(1,19), ...        %Updraft velocity (m/s)
    'kappa', {{ ...
        kappa_sol * Pars(1,13), ...
        kappa_sol * Pars(1,14), ...
        kappa_sol * Pars(1,15), ...
        kappa_sol * Pars(1,16), ...
     }}, ... 
    'n_tot', {{ ...
        [Pars(1,1)], ...
        [Pars(1,4)], ...
        [Pars(1,7)], ...
        [Pars(1,10)], ...
    }}, ...  %N mode [m^{-3}] at initial density
    'meanr', {{ ...
        [Pars(1,2) * 1e-6] ...
        [Pars(1,5) * 1e-6] ...
        [Pars(1,8) * 1e-6] ...
        [Pars(1,11) * 1e-6] ...
    }}, ...  %Mean radius mode [m]
    'gstdv', {{ ...
        [Pars(1,3)] ...
        [Pars(1,6)] ...
        [Pars(1,9)] ...
        [Pars(1,12)] ...
    }}, ...  %GSD
    'dt', dt, ... 
    'nt', ceil(1000 / Pars(1,19) / dt), ... % max output height (but actually stops at S_max)
    'sigma', Pars(1,18) / 1000, ... % surface tension coefficient for water [J/m2]
    'MAC', Pars(1,17) ...
)
output = drops_py(input)

% Define the return vector for comparison with the measurements
N_Parcel = [output.N_act, output.S_max]; 
N_Parcel = N_Parcel(:);
