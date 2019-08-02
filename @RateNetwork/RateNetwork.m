classdef RateNetwork < NeuralNetwork
    % RateNetwork is a sub-class of NeuralNetwork.  It inherits all methods
    % and properties from NeuralNetwork.  In addition, it defines a
    % rate-associated property: .R ,which is a structure that records
    % rate-related dynamic variables.  
    
    properties
        % Firing rate and related dynamical variables
        R % R needs to be a structure, whose fieldnames should be defined 
          % by the specific type of Models you use    
    end
    
    methods
    end
    
end

