function [T_TC, T_PreTerm, T_End, N_spk] = MeasureSeizureActivity(rate, varargin)
% [T_TC, T_PreTerm, T_End, N_spk] = MeasureSeizureActivity(rate, varargin)
%
% Input: rate - an n_neuron by n_time matrix, unit: normalized firing rate
%
% Output: T_TC - Tonic to clonic transition timing
%         T_PreTerm - Pretermination timing 
%         T_End - Termination time
% 
% All the time unit here is 'time index', which can only assume integers
%
% Measure the phasic transition timing 
%
% Jyun-you Liou 2017/04/08
p = inputParser;
addParameter(p,'center',63);
addParameter(p,'threshold',0.1);
addParameter(p,'termination_threshold',0.1);
addParameter(p,'T_stimulation',5000); % Until what time the external stimulation is still on
addParameter(p,'T_start',6000); % Since what time it is considered successful seizure initiation
parse(p,varargin{:});


% Decide if seizure initiation is successful
if all(rate(p.Results.center,p.Results.T_start:end)) < p.Results.threshold
    disp('Seizure initiation failed.');
    T_TC = nan;
    T_PreTerm = nan;
    T_End = nan;
    return;
end

% Remove the external stimulation part
rate = rate(:,p.Results.T_stimulation:end);

% Find tonic to clonic transition
[Pk,PreTerm_Idx] = findpeaks(rate(p.Results.center,:),'MinPeakProminence',p.Results.threshold);
if numel(PreTerm_Idx) <= 1
    disp('Tonic only seizure');
    T_TC = inf;
else
    T_TC = PreTerm_Idx(2); % Because the first peak must be transient
end

% Find termination
for i = size(rate,2):-1:1
    if any((rate(:,i)) > p.Results.termination_threshold)
        T_End = i;
        break;
    end 
end

% Find annihilation of wavefront
if T_TC ~= inf && ~isnan(T_TC) % Only do this if TC transition happens
    % Check if pre-term exists first    
    r_avg = mean(rate);
    [Pk,Tp,w,prom] = findpeaks(-r_avg);
    [~,Max_Idx] = max(r_avg);    
    [~,PreTerm_Idx] = max(prom);
    T_PreTerm = Tp(PreTerm_Idx);  
    
    [Pk,Tp] = findpeaks(rate(p.Results.center,T_PreTerm:end),'MinPeakProminence',p.Results.threshold/10);
    N_spk = numel(Pk);
    % This should the characteristics of pre-termination phase
    if isempty(Pk) || ...
       any(diff(Pk) > 0) || ...
       any(diff(diff(Tp)) < 0)
        disp('The pretermination phase may be wrong, suggest to review by eyes.');
    end
else
    T_PreTerm = nan;
end
end