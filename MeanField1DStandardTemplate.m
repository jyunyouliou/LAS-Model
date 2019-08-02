%% This is a standard template by a 1-dimensional rate model. 
% 
% This script calls function 'gen_param.m' 
% 
% Version: 2017/04/22, Jyun-you Liou

%% Network size
n = [500,1]; % Number of neurons 

%% Set parameters for each of the equations
% gen_param is a parameter generation function gen_param(mu,sigma,n,flag)
% mu: mean of the distribution
% sigma: standard deviation of the distribution
% n: network partition 
% flag: what type of distribution it is
% ------------- Equation 1 ---- Membrane potential ------------------------
V = gen_param(-58,0,n,0); % Initial membrane potential
C=gen_param(100,0,n,1); % Cell capacitance, NeuroElectro median 100 pF, highly variable
% Leaky conductance
g_L=gen_param(4,0,n,1); % NeuroElectro: mean input resistence = 300 M
E_L=gen_param(-58,0,n,0); % NeuroElectro: Resting membrane potential -62 mV
% Spike generator - threshold noise model
T_refractory = gen_param(5,0,n,1);% ms, the duration of AP + absolutely refractory period 
                                  % in spiking model, during T_refractory, phi will be set to +inf
f_max = 1./T_refractory;% kHz, highest possible firing rate
f0 = gen_param(0.001,0,n,1); % kHz (because time unit is ms)
beta = gen_param(2.5,0,n,1); % mV, high beta = low threshold noise
f = @(u) f_max./(1+exp(-u./beta)); % unit: kHz

% Synaptic input filter and reversal potential
tau_syn.E = gen_param(15,0,n,0); % Excitatory synaptic time constant 
tau_syn.I = gen_param(15,0,n,0); % Inhibitory synaptic time constant
E_Esyn= gen_param(0,0,n,0); % Reversal potential of excitatory synapses 

% Supplementary dynamics - short-term plasticity

% Short-term plasticity variables, reference: Science 2008, Mongillo et. al.
% Larry & Misha's models are actually equivalent if you set u constant in Misha's model
% It's just Larry's output is x, but Misha's output is u*x, for detailed explanation, 
% please read the word file 
flag_STP = false; % Whether STP is allowed or not
tau_D = gen_param(0,0,n,1); % ms, put it 0 to disable depression
tau_F = gen_param(0,0,n,1); % ms, put it 0 to disable faciliation
U = gen_param(0.2,0,n,1); % 1-U is actually f_D in Larry's model if u is constant
    % U stands for portion of calcium influx, set tau_F = 0 to terminate
    % facilitation process
    % The model actually requires you to think about u-(left limit) & u+(right limit) 
    % u stands for calcium concentration, and it has upper limit 1
    % the amount of vesicle release depends on u+ (limit from the right)
    % but the equation, u' = -u + U(1-u)dirac, actually describes u- (limit from the left)
    % so in their paper, the equation u is actually u+, which u+ = u-+U(1-u-) 

% ------------- Equation 2 ---- Dynamical spiking threshold ---------------
% General form of dynamic threshold is dphi = f(phi)
phi = gen_param(-45,0,n,0); % mV, initial condition of threshold
phi_0 = gen_param(-45,0,n,0); % mV, NeuroElectro, median approximately -45
phi_inf = @(u) phi_0+u.*60; % u: ratio of rate/f_max
tau_phi = gen_param(100,0,n,1); % ms

% ------------- Equation 3 ---- Chloride dynamic parameters ---------------
% 1) So assume diameter = 15 micro, 40% available space of cytosol, then 50% free water, then
% 2) Use John Rinzel's parameter, diameter = 10 micro, round
% 3) border 20 micro, isotetrahedrum
Vd_Cl = 0.25 .* sqrt(2) / 12 * 20^3 /1000 .* ones(n); % Unit: pL (10^-12L, which is 10^-15m^3, which is 10^3 microm^3, so if you use micro, remember to /1000)
                                                  % Reference: Thus, the amount of restricted water is thought to be of the order of 50% of all water in the cytoplasm (Luby-Phelps, 2000; Fullerton & Cameron, 2007).
Cl_ex = 110; % Reference: J Neurophysiol. 1988b;60:195–124. (Berglund et al. 2006; Glykys et al. 2014). 6~14
Cl_in_eq = gen_param(6,0,n,1); % The equilibruim intracellular chloride concentration, 
tau_Cl = gen_param(5000,0,n,1); % Chloride clearance time constant
Cl_in = Cl_in_eq;

% ------------- Equation 4 ---- slow afterhyperpolarization ---------------
tau_K = gen_param(5000,0,n,1); % Can also be function
g_K_max = gen_param(40,0,n,1); % maximal g_K when firing rate is f_max
E_K = gen_param(-90,0,n,0);
g_K = zeros(n);