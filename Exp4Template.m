SpikingModel1DStandardTemplate;

% Parameter adjustment
beta = gen_param(1.5,0,n,1); 
f = @(u) f0+fs.*exp(u./beta); % unit: kHz

E_L = ones(n)*-60;
V = E_L;