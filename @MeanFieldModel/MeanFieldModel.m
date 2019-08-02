classdef MeanFieldModel < RateNetwork
    % Class definition for a mean-field model (This is a rate approach of 
    % SpikingModel)
    % 
    % Model = MeanFieldModel(template)
    %
    % Input: template - a filepath that leads to a .m script file which
    %                   specifies parameters of the model
    %
    % Default setting for properties inherited from NeuralNetwork
    % Input: 'I' & 'E'
    % S: 'S', 'dT', 'x', 'u'
    % 
    % Author: Jyun-you Liou
    % Final update: 2017/01/18
    
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
                if all(O.param.tau_syn.(O.Proj.In(i).Type) == 0) % Instant synapse, no synaptic filtering
                    O.Input.(O.Proj.In(i).Type) = O.Input.(O.Proj.In(i).Type) + ...
                                                  O.Proj.In(i).Value;
                else % Apply synaptic filtering
                    O.Input.(O.Proj.In(i).Type) = O.Input.(O.Proj.In(i).Type) + ...
                                                  O.Proj.In(i).Value ./ O.param.tau_syn.(O.Proj.In(i).Type);
                end
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
            O.phi = O.phi + (O.param.phi_inf(O.R.f./O.param.f_max)-O.phi).*(1-exp(-dt./O.param.tau_phi));
                % In rate model, use rate instead of Dirac pulse method
                % (spiking) because it causes more numeric error
            
            % ----- Update equation 3 (chloride dynamics) ------ 
            Faraday = 96500;               
            Cl_in_inf = (O.param.tau_Cl./O.param.Vd_Cl./Faraday.*O.Input.I.*(O.V-E_Cl) + O.param.Cl_in_eq);
            O.Cl_in = Cl_in_inf + (O.Cl_in - Cl_in_inf).*(exp(-dt./O.param.tau_Cl));   

            % ----- Update equation 4 (sAHP dynamics) ----------
            O.g_K = O.g_K .* exp(-dt ./ O.param.tau_K) + O.param.g_K_max.*O.R.f.*dt./O.param.tau_K;    

            % ######### Generate rate #########
            O.R.f = O.param.f(O.V-O.phi); % Unit: kHz            
            
            % ######### Calculate all R-derived variables ######### 
  
            % ######### Short-term plasticity variables ########
            % Calculate effective strength (# of spikes) by considering 
            % .x & .u (See Mongillo et. al. 2008 Science)
            %
            % Save the effective # of spikes to its projection 
            if O.param.flag_STP && all(O.param.tau_D(:)) && all(O.param.U(:))
                % Within the parenthesis now becomes the RATIO compare the strength to a fresh spike                
                [O.Proj.Out.Value] = deal(((O.R.x.*O.R.u)./O.param.U) .* O.R.f); % divided by dt to preserve .Proj.Value 's unit as rate
                if all(O.param.tau_F(:)) % Short-term facilitation
                    O.R.x = 1 - (1-O.R.x).*exp(-dt./O.param.tau_D) - O.R.u.* O.R.x .* O.R.f; % x:percentage of synaptic vesicle reserve
                    O.R.u = O.param.U - (O.param.U-O.R.u).*exp(-dt./O.param.tau_F) + O.param.U.*(1-O.R.u).*O.R.f; % u:percentage of synaptic vesicles to be used 
                else
                    O.R.x = 1 - (1-O.R.x).*exp(-dt./O.param.tau_D) - O.param.U.*O.R.x.*O.R.f;
                end
            else
                [O.Proj.Out.Value] = deal(O.R.f.*dt); % Send to .Proj.Value # of spikes
            end
            
            % Filter the synaptic input
            if all(O.param.tau_syn.(O.Proj.In(i).Type) == 0) 
                % Instant synapse, no synaptic filtering, just clear the value
                O.Input.E = 0;
                O.Input.I = 0;
            else
                % Apply synaptic filtering
                O.Input.E = O.Input.E.*exp(-dt./O.param.tau_syn.E);
                O.Input.I = O.Input.I.*exp(-dt./O.param.tau_syn.I);
            end
            
            % Update time
            O.t = O.t + dt;
        end
    end
  
    %% Construct the object
    methods 
        function O = MeanFieldModel(template)
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
            O.R.f = zeros(O.n); % firing rate         
            O.R.x = ones(O.n); % short-term plasticity variables
            O.R.u = O.param.U; % short-term plasticity variables
            
            %% Set acceptable input types 
            O.Input.E = zeros(O.n);
            O.Input.I = zeros(O.n);
        end
    end
    
    
end

