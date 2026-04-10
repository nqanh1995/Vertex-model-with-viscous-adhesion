%% run avm
addpath("lees_edwards_functions/")
addpath("matlab_functions/")

nsteps = 1000; %%% number of simulation steps
interval = 10; %%% data saving frequency
seed=1;
Npts = 100; %%% number of cells in the system
p0 = 3.7; %%% p0 value for all cells
KA = 1; %%% area elasticity
dt = 0.05;
tau_T1 = 10; %%% cooldown time between T1 transition on an edge
L_thresh = 0.05; %%% edge length threshold for T1 transition
v0 = 0.05;
dr = 0.5;
B = 1; %%% cell-cell viscous coefficient \xi
avm_viscous(nsteps,interval,seed,Npts,p0,KA,dt,tau_T1,L_thresh,v0,dr,B)

%% run simple shear avm
addpath("lees_edwards_functions/")
addpath("matlab_functions/")

nsteps = 1000; %%% number of simulation steps
interval = 10; %%% data saving frequency
seed=1;
Npts = 100; %%% number of cells in the system
p0 = 3.7; %%% p0 value for all cells
KA = 1; %%% area elasticity
dt = 0.05;
tau_T1 = 10; %%% cooldown time between T1 transition on an edge
L_thresh = 0.05; %%% edge length threshold for T1 transition
v0 = 0.05;
dr = 0.5;
B = 1; %%% cell-cell viscous coefficient \xi
strain_amp = 0.005; %%% harmonic strain amplitude
omega = 1;  %%% shearing frequency
stress_only= 1; %%% option to save only the stress data to save storage, 1 to save only stress values, 0 to save configuration and stress

harmonic_oscillatory_shear(strain_amp,nsteps, interval, seed, Npts, p0, KA, dt, L_thresh, v0, dr, omega, B,stress_only)

