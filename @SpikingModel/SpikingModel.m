classdef SpikingModel < SpikingNetwork
    % Class definition for a spiking model 
    % 
    % Model = SpikingModel(template)
    %
    % Input: template - a filepath that leads to a .m script file which
    %                   specifies parameters of the model
    %
    % Default setting for properties inherited from NeuralNetwork
    % Input: 'I' & 'E'
    % S: 'S', 'dT', 'x', 'u'
    % 
    % Author: Jyun-you Liou
    % Final update: 2016/12/22
    
    properties
        % The model inherit from NeuralNetwork with the following
        % properties: n,t,Input,Proj,Ext,Recorder,UserData
        
        % The model's intrinsic parameters
        param     
        
        % The 4 dynamical variables - if you want to change the name or #
        % of dynamical variables, please remember to update the Constant 
        % property 'VarName' 
        V % averaged membrane potential
        phi % spiking threshold
        Cl_in % intracellular chloride concentration 
        g_K % slow AHP conductance  
    end
    
    properties (Constant)
        VarName = {'V','phi','Cl_in','g_K'}
    end
    
    methods % Simulation functions
        function IndividualModelUpdate(O,dt) 
            % IndividualModelUpdate(Model,dt)
            %
            % Perform the first section of dt-step simulation. This 
            % internal function is called by Update(Model,dt).
            %
            % Jyun-you Liou 2016/12/30
     
            % Collect input from 'ExternalInput'
            I_ext = 0;
            for i = numel(O.Ext):-1:1 
                % You need to evaluate from the back because ExternalInput
                % may be destroyed after expiration time. (.Tmax)
                I_ext = Evaluate(O.Ext(i),O.t,dt) + I_ext;
            end
            
            % Collect input from 'Projection'
            for i = 1:numel(O.Proj.In) % So the maximum of O.Input is W0 * f_max
                O.Input.(O.Proj.In(i).Type) = O.Input.(O.Proj.In(i).Type) + ...
                                              O.Proj.In(i).Value ./ O.param.tau_syn.(O.Proj.In(i).Type);
            end
            
            % ----- Update equation 1 (membrane potential) -----
            g_sum = O.param.g_L + ...
                    O.Input.E ./ O.param.f_max + ...
                    O.Input.I ./ O.param.f_max + ...
                    O.g_K ./ O.param.f_max; % from outside

            % V_inf - this equation still needs 'I'
            E_Cl = 26.7*log(O.Cl_in./O.param.Cl_ex); % Calculate E_Cl (inhibition effectiveness) 
            V_inf = (O.param.g_L .* O.param.E_L  + ...
                     O.Input.E ./ O.param.f_max .* O.param.E_Esyn + ...
                     O.Input.I ./ O.param.f_max .* E_Cl + ...  
                     O.g_K ./ O.param.f_max .* O.param.E_K + ...
                     I_ext) ./ g_sum;  
                 
            tau_V_eff = O.param.C./g_sum;

            O.V = V_inf + (O.V - V_inf).*(exp(-dt./tau_V_eff));   
            
            % ----- Update equation 2 (threshold dynamics) -----
            O.phi = O.param.phi_0 + (O.phi - O.param.phi_0).*(exp(-dt./O.param.tau_phi)) + ... % Exponential decay part
                    O.param.delta_phi(O.phi-O.param.phi_0) .* O.S.S; % Dirac pulse part

            % ----- Update equation 3 (chloride dynamics) ------ 
            T_AP = O.param.dt_AP - O.S.dT; % How much time spent for action potential
            T_AP = T_AP .* (T_AP>0);
            T_non_AP = dt - T_AP; % How much time not spent for action potential
            Veff = T_AP.*O.param.V_AP./dt + ...
                   T_non_AP.*O.V./dt; % effective membrane potential 
            Faraday = 96500;               
            Cl_in_inf = (O.param.tau_Cl./O.param.Vd_Cl./Faraday.*O.Input.I.*(Veff-E_Cl) + O.param.Cl_in_eq);
            O.Cl_in = Cl_in_inf + (O.Cl_in - Cl_in_inf).*(exp(-dt./O.param.tau_Cl));   

            % ----- Update equation 4 (sAHP dynamics) ----------
            O.g_K = O.g_K .* exp(-dt ./ O.param.tau_K) + O.param.g_K_max.*O.S.S./O.param.tau_K;    
           
            % ######### Generate spikes #########
            % Notice in numeric simulation, dirac function should take a value so that delta * dt = 1
            Spike = O.param.f(O.V - O.phi)*dt > rand(O.n); % Spike generation
            Spike(O.S.dT < O.param.T_refractory) = 0; % Remove the spike if still refractory
            O.V(Spike) = O.param.V_reset(O.V(Spike)); % Reset membrane potential
            O.S.S = Spike; % Save it
            
            % ######### Calculate all S-derived variables ######### 
            % The reason why update of S-derived variables needs to happen 
            % in the same cycle is because Dirac function cause step from 
            % the right limit, not the left limit. 
            % Refractory period: .dT - when the neurons spike last time
            O.S.dT = O.S.dT + dt;
            O.S.dT(O.S.S) = 0;
  
            % Short-term plasticity variables: .x & .u can affect effective
            % output strength
            %
            % Save the effective firing rate into its projections 
            if O.param.flag_STP && all(O.param.tau_D(:)) && all(O.param.U(:))
                % Within the parenthesis now becomes the RATIO compare the strength to a fresh spike                
                [O.Proj.Out.Value] = deal(((O.S.x.*O.S.u)./O.param.U) .* O.S.S); % divided by dt to preserve .Proj.Value 's unit as rate
                if all(O.param.tau_F(:)) % Short-term facilitation
                    O.S.x = 1 - (1-O.S.x).*exp(-dt./O.param.tau_D) - O.S.u.* O.S.x .* O.S.S; % x:percentage of synaptic vesicle reserve
                    O.S.u = O.param.U - (O.param.U-O.S.u).*exp(-dt./O.param.tau_F) + O.param.U.*(1-O.S.u).*O.S.S; % u:percentage of synaptic vesicles to be used 
                else
                    O.S.x = 1 - (1-O.S.x).*exp(-dt./O.param.tau_D) - O.param.U.*O.S.x.*O.S.S;
                end
            else
                [O.Proj.Out.Value] = deal(O.S.S); % Send to .Proj.Value 'rate'
            end
            
            % Filter the synaptic input (required for spiking network)          
            O.Input.E = O.Input.E.*exp(-dt./O.param.tau_syn.E);
            O.Input.I = O.Input.I.*exp(-dt./O.param.tau_syn.I);
            
            % Update time
            O.t = O.t + dt;
        end
    end
   
    %% Construct the object
    methods (Static)
        function O = SpikingModel(template)
            % Model = SpikingModel(template)
            %
            % Final update: 2016/12/12
            % Load parameters & initial conditions    
            eval(template);
            VarList = who;
            for i = 1:numel(VarList)
                if any(strcmp(VarList{i},O.VarName)) % Initial conditions of dynamical variables
                    O.(VarList{i}) = eval(VarList{i});
                elseif strcmp(VarList{i},'n') % Size of the network
                    O.n = eval(VarList{i});                    
                else % Parameters
                    O.param.(VarList{i}) = eval(VarList{i});
                end              
            end
            
            %% Set firing-associated dynamical variables 
            O.S.S = false(O.n); % Logical value, decide whether there is a spike or not
            O.S.dT = zeros(O.n); % The previous spike            
            O.S.x = ones(O.n); % short-term plasticity variables
            O.S.u = O.param.U; % short-term plasticity variables
            
            %% Set acceptable input types 
            O.Input.E = zeros(O.n);
            O.Input.I = zeros(O.n);
        end
    end
    
    
end

