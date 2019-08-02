classdef SpikingNetwork < NeuralNetwork
    % SpikingNetwork is a sub-class of NeuralNetwork.  It inherits all 
    % methods and properties from NeuralNetwork.  In addition, it defines a
    % spiking-associated property: .S ,which is a structure that records
    % spike-related dynamic variables, such as time from the last spike 
    % (which is required to implement refractory period), short-term 
    % plasticity variables ... etc.  
    
    properties
        % Spiking related dynamical variables
        S % S needs to be a structure, whose fieldnames should be defined 
          % by the specific type of Models you use 
          % Commonly used for spiking models: S, dT ... etc.         
    end
    
    methods
    end
    
end

