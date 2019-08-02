classdef ExternalInput < matlab.mixin.Copyable % ExternalInput inherit shallow copy method
    % Ext = ExternalInput
    %
    % Class ExternalInput provides additional input to NeuralNetwork.
    % ExternalInput is a copyable handle object.  You can call its method
    % 'copy' to make a shallow copy of an instance of this class.
    %
    % Properties:
    %
    % .Target - an instance of NeuralNetwork where the ExternalInput is
    %           delivered, only 1 NeuralNetwork can be targeted by one 
    %           ExternalInput instance.  Copy the original instance if you 
    %           want to build several of them 
    % .Random - a strcture that encodes information about the random part 
    %           of ExternalInput.  Please see doc ExternalInput for more 
    %           information 
    % .Deterministic - the deterministic part of ExternalInput.  
    % .Tmax - the expiration time of this instance.  Once simulation of the
    %         target NeuralNetwork goes beyond this time, the instance will
    %         be automatically destroyed 
    % .UserData - For user-defined functions
    %                  
    % Jyun-you Liou, 2017/04/01
    
    properties
        % .Target is an instance of NeuralNetwork which the ExternalInput 
        % exerts its effect on  
        Target
        
        % .Random is a structure that has the following fields: .tau_x, 
        % .tau_t, sigma, and Iz.  This allows you to generate noise by OU
        % process
        %
        % .tau_x - spatial constant for creating pink noise in space   
        %          by Ornstein–Uhlenbeck process (scalar)
        % .tau_t - time constant for creating pink noise in time 
        %          by Ornstein–Uhlenbeck process (scalar)
        % .sigma - amplitude of the noise (scalar)
        % .Iz - the initial condition of OU process.  It can be set to a 
        %       standard normal distribution  saved for enforcing so that
        %       you don't need to wait several time constants until the 
        %       process achieves its steady state 
        %
        % If you want just white noise, set tau_x & tau_t both to empty.
        Random = struct('tau_x',[], ...
                        'tau_t',[], ...
                        'sigma',0, ...
                        'Iz',0);

        % .Deterministic is a function handle that gives the deterministic 
        % part of ExternalInput.  The function needs to have two inputs. 
        % The first input is position and the second input is time.
        % The function requires to work in a vectorized manner.  What does
        % this mean? Here is the illustration
        % 
        % @(x,t) (t>0)*(x(:,1)>0)
        %
        % x is a n_neuron by 2 matrix, with its first row indicates its
        % position along the first dimesion (unit: # of neurons, not 
        % normalized length).  The second row is its position along the
        % second dimension.
        Deterministic 
        
        % .Tmax indicates until what time this ExternalInput is valid. 
        % Over this time, the ExternalInput instance will be destroyed
        % Not absolutely required.  If left blank, the instance of 
        % ExternalInput will never be destoryed.  However, to faciliate
        % simulation speed, I suggest to set Tmax so that it won't need to
        % be evaluated every simulation cycle.
        Tmax 
        
        % .UserData is for any additional data that needs to be stored.
        UserData 
    end
    
    methods % Evaluate an ExternalInput
        function I = Evaluate(Ext,t,dt)
            % I = Evaluate(Ext,t,dt)
            % 
            % Exert the effect of an external input device at time t with
            % time step dt.
            %
            % Inputs: Ext - an instance of 'ExternalInput'
            %         t - evaluate at time t
            %         dt - the time step in the simulation
            %
            % Get its projection's information
            O = Ext.Target;
            n = O.n;
            [x2,x1] = meshgrid(1:n(2),1:n(1));
            I = 0;
            % Check these ExternalInput have not expired
            CheckValidity(Ext);
            if ~isvalid(Ext);return;end
            
            % Deterministic Part
            if ~isempty(Ext.Deterministic)
                I_d = Ext.Deterministic([x1(:),x2(:)],t);
                I_d = reshape(I_d,n);
                I = I + I_d;
            end

            % Random part - Ornstein–Uhlenbeck process
            % Space
            z = randn(n);
            if ~isempty(Ext.Random.tau_x) 
                z = SpaceFilter(z,Ext.Random.tau_x);
            end
            function z = SpaceFilter(z,tau_x) 
                dx = 1;
                for i = 1:numel(tau_x)
                    if tau_x(i)>0
                        z = permute(z,[i sort(setdiff(1:numel(size(z)),i))]);
                        nz = size(z);
                        z = reshape(z,size(z,1),[]);
                        z = filter(sqrt(2*dx/tau_x(i)),[1, (dx-tau_x(i))/tau_x(i)],z, ...
                                   (1-sqrt(2*dx/tau_x(i)))*z(1,:));                        
                                   % Initial condition needs to be set to avoid 'warm up' phenomenon
                        z = reshape(z, nz);
                        z = ipermute(z,[i sort(setdiff(1:numel(size(z)),i))]);
                    end
                end          
            end
            
            % Time - Ornstein–Uhlenbeck process
            if ~isempty(Ext.Random.tau_t) && Ext.Random.tau_t > 0
                % I may need to be initiated if it is empty to prevent
                % warm-up phenomenon
                if all(Ext.Random.Iz(:)) == 0
                    z0 = randn(n);
                    if ~isempty(Ext.Random.tau_x)
                        z0 = SpaceFilter(z0,Ext.Random.tau_x);
                    end
                    Ext.Random.Iz = z0;
                end
                Ext.Random.Iz = Ext.Random.Iz .* exp(-dt ./ Ext.Random.tau_t) + ... % The autoregression part
                               sqrt(2*dt/Ext.Random.tau_t) .* z; % Innovation part 
            else
                Ext.Random.Iz = z;
            end
                
            % Calculate result
            I = I + Ext.Random.sigma*Ext.Random.Iz;
        end
    end

    methods % Property set functions
        function Ex = set.Target(Ex,O)
            if ~isa(O,'NeuralNetwork')
                error('The target of ExternalInput needs to be an instance of NeuralNetwork or its subclasses');
            end
            % Ensure the instance of ExternalInput is properly registered            
            Ex.Target = O;
        end
    end
    
    methods % Housekeeping functions
        function CheckValidity(E)
            % CheckValidity(ExternalInput)
            %
            % Check whether this ExternalInput instance is still within its
            % effective time range.  If not, delete it to save computation
            % power.
            %
            % Jyun-you Liou 2017/01/09
            if numel(E)>1
                for i = 1:numel(E)
                    CheckValidity(E(i));
                end
                return;
            end
            O = E.Target;
            if ~isempty(E.Tmax) && O.t > E.Tmax
                delete(E)
                disp('An ExternalInput instance has expired & deleted.')
            end
        end
    end
    
    methods % Destructor functions
        function delete(E)
            % Make sure the deleted ExternalInput will be removed from the
            % registration list of its .Target NeuralNetwork 
            %
            % Jyun-you Liou 2017/01/09
            if isvalid(E.Target)
                O = E.Target;
                Sel = O.Ext == E;
                O.Ext(Sel) = [];
            end
        end
    end

    methods (Static = true) % Static methods for create some stimulation patterns
        Func = RingPattern
        Func = FocalPattern
    end
end

