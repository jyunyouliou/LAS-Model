function Project(P)
% Project(P)
%
% Project the information saved in P.Value (effective 
% pre-synaptic # of spikes), calculate its projection
% (Method, Topology, WPre, W), and save the projection results 
% to P.Target.Input.Type according to to the Projection.Type
% 
% This function is called by @NeuralNetwork.Project
%
% Jyun-you Liou 2017/05/04

%% Projection
for i = 1:numel(P)
    
    % Strength of each pre-synaptic neuron
    if ~isempty(P.WPre) && any(P.WPre(:) ~= 1)                 
        P(i).Value = P(i).WPre .* P(i).Value; 
    end
    
    % Projection, according it how synaptic projection is
    % distributed, use the correspondent method
    switch lower(P(i).Method)
        case 'convolution'
            P(i).Value = conv2_field(P(i).Value, ...
                                     P(i).W, ...
                                     P(i).Target.n, ...
                                     P(i).Topology);
        case 'multiplication'
            P(i).Value = P(i).W * P(i).Value(:); 
            P(i).Value = reshape(P(i).Value,P(i).Target.n);
        case 'function'
            P(i).Value = P(i).W(P(i).Value);
    end
    
    % Strength of each post-synaptic receptors                
    if ~isempty(P.WPost) && any(P.WPost(:) ~= 1) 
        P(i).Value = P(i).Value .* P.WPost;
    end
end

end