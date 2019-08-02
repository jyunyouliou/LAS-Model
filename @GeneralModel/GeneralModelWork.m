classdef GeneralModel < RateNetwork
    % Class definition for a general model (neural field)
    % 
    % Model = GeneralModel(template)
    %
    % Input: template - a filepath that leads to a .m script file which
    %                   specifies parameters of the model
    %
    % Author: Jyun-you Liou
    % Final update: 2016/12/28
    
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
        z % effectiveness of local inhibition 
        U % resting membrane potential 
    end
    
    properties (Constant)
        VarName = {'V','phi','z','U'} % For quick inquiry about what dynamic variables this type of model has
    end
    
    methods % Simulation functions        
        function IndividualModelUpdate(O,dt) 
            % IndividualModelUpdate(Model,dt)
            %
            % Perform the first section of dt-step simulation. This 
            % internal function is called by Update(Model,dt).
            %
            % Jyun-you Liou 2016/12/28
            
            % Collect input from 'ExternalInput'
            I_ext = 0;
            for i = numel(O.Ext):-1:1 
                % You need to evaluate from the back because ExternalInput
                % may be destroyed after expiration time. (.Tmax)
                I_ext = Evaluate(O.Ext(i),O.t,dt) + I_ext;
            end
            
            % Collect input from 'Projection'
            for i = 1:numel(O.Proj.In)
                O.Input.(O.Proj.In(i).Type) = O.Input.(O.Proj.In(i).Type) + O.Proj.In(i).Value;
            end
            
            % Update equation 1 (membrane potential)
            V_inf = O.U + O.Input.E - O.z.*O.Input.I + I_ext;
            O.V = V_inf + (O.V - V_inf).*exp(-dt./O.param.tau_V);      
            
            % Update equation 2 (threshold dynamics)
            phi_inf = O.param.phi_0 + O.param.Delta_phi.*O.S.f;
            O.phi = phi_inf + (O.phi - phi_inf).*exp(-dt./O.param.tau_phi);      
            
            % Update equation 3 (local inhibition effectiveness)
            O.z = O.param.z_inf(O.Input.I) + (O.z - O.param.z_inf(O.Input.I)).*exp(-dt./O.param.tau_z);
            
            % Update equation 4 (resting membrane potential dynamics)
            U_inf = O.param.Delta_U.*O.S.f + O.param.U0;
            O.U = U_inf + (O.U - U_inf).*exp(-dt./O.param.tau_U);    
            
            % ######### Generate firing #########
            O.S.f = O.param.sigmoid(O.V-O.phi); % (rate)

            % Save effective firing strength into .Proj (inherited from NeuralNetwork)
            [O.Proj.Out.Value] = deal(O.S.f); % You need 'deal' because .Proj can be multiple
                                              % Here we directly projected
                                              % 'normalized firing rate'
            
            % Erase all fields of Input structure 
            O.Input = structfun(@(x) 0,O.Input,'Un',false);
            
            % Update time
            O.t = O.t + dt;
        end
    end
    
    %% Construct the object
    methods (Static)
        function O = GeneralModel(template)
            % Model = GeneralModel(template)
            %
            % Final update: 2016/12/22
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
            O.S.f = zeros(O.n);
            
            %% Set acceptable input types
            O.Input.E = zeros(O.n);
            O.Input.I = zeros(O.n);
        end
    end
    
    
end

