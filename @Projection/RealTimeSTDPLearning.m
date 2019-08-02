function RealTimeSTDPLearning(P,dt)
% RealTimeSTDPLearning(Projection,dt)
%
% Modify projection strength real-time
%
% Jyun-you Liou 2017/04/17
%
% Final update: remove bugs & allow BatchSTDPLearning to call
% RealTimeSTDPLearning while rate model is used 


%% Check the projection settings do allow the synapses to learn
if isempty(P.STDP)
    disp('You need to set up STDP learning rule by setting Projection.STDP');return;
end

if ~strcmpi(P.Method,'multiplication')
    error('Projection.Method needs to be multiplication.');
end

%% Catch the NeuralNetwork objects & their Recorders
Os = P.Source; % Source network
Ot = P.Target; % Target network


%% Update the filtered S variable for STDP
% Pre-synaptic part
if isa(Os,'SpikingNetwork')
    S_pre = double(Os.S.S); % How many spikes are there (consider it as Dirac function)
    P.S.LTP = P.S.LTP.*exp(-dt./P.STDP.tau_LTP) + S_pre./P.STDP.tau_LTP;
elseif isa(Os,'RateNetwork')
    S_pre = dt.*Os.R.f;
    P.S.LTP = P.S.LTP.*exp(-dt./P.STDP.tau_LTP) + S_pre./P.STDP.tau_LTP;       
end

% Post-synaptic part
if isa(Ot,'SpikingNetwork')
    S_post = double(Ot.S.S);
    P.S.LTD = P.S.LTD.*exp(-dt./P.STDP.tau_LTD) + S_post./P.STDP.tau_LTD;    
elseif isa(Ot,'RateNetwork')
    S_post = dt.*Ot.R.f;
    P.S.LTD = P.S.LTD.*exp(-dt./P.STDP.tau_LTD) + S_post./P.STDP.tau_LTD;       
end


%% Calculate how much change spiking produced given soft/hard bound (Wmin & Wmax) is considered 
if strcmpi(P.STDP.BoundType,'soft') % I did not write it in exponential integration form so that it is easier to understand
    dW = (P.STDP.dLTP .* (P.STDP.Wmax - P.STDP.W) .* (S_post * P.S.LTP(:)')  - ... % LTP-part
          P.STDP.dLTD .* (P.STDP.W - P.STDP.Wmin) .* (P.S.LTD(:) * S_pre')) * dt; % LTD-part
elseif strcmpi(P.STDP.BoundType,'hard')
    dW = (P.STDP.dLTP .* (P.STDP.Wmax > P.STDP.W) .* (S_post * P.S.LTP(:)')  - ... % LTP-part
          P.STDP.dLTD .* (P.STDP.W > P.STDP.Wmin) .* (P.S.LTD(:) * S_pre')) * dt; % LTD-part                
end

P.STDP.W = P.STDP.W + dW;

if P.STDP.Homeostasis
    P.STDP.W = bsxfun(@rdivide, P.STDP.W, sum(P.STDP.W)); 
end

end