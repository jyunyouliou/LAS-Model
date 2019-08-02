function [ P_E, P_I1, P_I2 ] = StandardRecurrentConnection( O, varargin )
% [ P_E, P_I1, P_I2 ] = StandardRecurrentConnection( O, varargin )
%
% This is for building standard recurrent connection.  
%
% Input: O: needs to be a NeuralNetwork 
%        varargin:  'Adjust' - default: false
%                              If the simulation space has aperiodic
%                              boundary, cells at the boundary will have
%                              different maximal amount of synaptic inputs
%                              than other cells.  If you turn this on, it
%                              will address this problem by amplifying 
%                              their synaptic inputs.
%
% P_E: Spatially localized recurrent excitation 
% P_I1: Spatially localized recurrent inhibition
% P_I2: Non-localized recurrent inhibition
%
% PS: I don't want to put this function into NeuralNetwork folder because
% this function is just commonly called for this epilepsy study.
% 
% Jyun-you Liou, 2017/04/30
p = inputParser;
addParameter(p,'Adjust',false);
parse(p,varargin{:});
Adj = p.Results.Adjust;

% Build recurrent excitation & configure it
P_E = Projection(O,O,'Type','E','Topology','linear');
Sigma_E = diag(O.n) * 0.02; % percentage of the field 
Kernelize(P_E, @(x) mvnpdf(x,[0 0],Sigma_E.^2), 'KerSize', ceil(2.5*diag(Sigma_E)));
if Adj;AdjustWeight(P_E);end % Adjust strength at space border; 
P_E.WPost = P_E.WPost * 100; % Projection strength 

% Build recurrent inhibition & configure it
P_I1 = Projection(O,O,'Type','I','Topology','linear');
Sigma_I = diag(O.n) * 0.03; % percentage of the field  
Kernelize(P_I1, @(x) mvnpdf(x,[0 0],Sigma_I.^2), 'KerSize', ceil(2.5*diag(Sigma_I)));
if Adj;AdjustWeight(P_I1);end % Adjust strength at space border; 
P_I1.WPost = P_I1.WPost * 250; % Projection strength

% Build the global recurrent inhibition & configure it
P_I2 = Projection(O,O,'Type','I','Method','function');
P_I2.W = @(x) sum(x(:))/prod(O.n); % uniform distribution
if Adj;AdjustWeight(P_I2);end % Adjust strength at space border; 
P_I2.WPost = P_I2.WPost * 50; % Projection strength

end

