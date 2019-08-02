classdef Projection < handle
    % Projection(NeuralNetwork1, NeuralNetwork2, 'parameter', parameter)
    %
    % Connection is a Class that is used to link two models together.  
    % You can also create recurrent connection by putting NeuralNetwork1 &
    % NeuralNetwork2 as the same object.
    %
    % Inputs: NeuralNetwork1, NeuralNetwork2 - Needs to be instances of 
    %         Class NeuralNetwork
    %
    % 'parameter'
    %
    % Method - 'multiplication', 'convolution', or 'function' 
    %          This parameter determines how the projection from Oi
    %          to Oj is computed.
    %          Default: 'convolution'
    % WPre - Projection strength, indexed by Oi neurons.
    % WPost - Post-synaptic strength - indexed by Oj neurons 
    % W - Connectivity information, if Method is 'convolution', 
    %     W will be treated as convolution kernel
    %     If Method is 'multiplication', then W is weight matrix
    % Topology - Only effective in 'convolution' condition. 
    %            'linear' or 'circular'
    % Type - A character string describing the type of this projection.
    %        Commonly used: 'E' or 'I'.  NeuralNetwork will identify this 
    %        tag and load the Projection's .Value accordingly to compute
    % STDP - Control if STDP learning is on or off and how STDP is going to
    %        be implemented.
    %        It's a structure with the following fields: Enabled, W, 
    %        tau_LTP, tau_LTD, dLTP, dLTD, Wmin, Wmax, BoundType, 
    %        Homeostasis.  
    %        Please see the classdef file itself to undestand how these
    %        parameters work.
    % Value - The value of this projection. 
    %         While you update a presynaptic NeuralNetwork, the
    %         NeuralNetwork will put its effective strength into this
    %         property.  
    %         While you Project this Projection, the Value will be computed
    %         so that post-synaptic NeuralNetwork can load as its input.
    % S - Filtered pre-synaptic spike train/rate for STDP learning
    %     .LTP - for the potentiation part of STDP
    %     .LTD - for the depression part of STDP
    % Jyun-you Liou, 2017/04/29 (Change W0 to WPre)
    
    properties
        Source % Source (Pre-synaptic) NeuralNetwork
        Target % Target (Post-synaptic) NeuralNetwork
        Method = 'convolution'
               % 'multiplication' or 'convolution' or 'function'
               %  This parameter determines how the projection from Source 
               %  to Target is computed.
        WPre = 1 % Strength of the projection, indexed by pre-synaptic neurons
        WPost = 1% Strength of the projection, indexed by post-synaptic neurons
        W = 1
          % Connectivity information, if Method is 'convolution', W 
          % will be treated as convolution kernel.  If Method is 
          % 'multiplication', then W is weight matrix.  If Method is 
          % 'function', then this is a function handle that while Project 
          % method is implemented, Proj.Value = Proj.W(Proj.WPre.*Proj.Value)
        Topology = 'linear'
                 % - Only effective when method is 'convolution'. 
                 %   Can assume value as 'linear' or 'circular'
        Type = 'E'
             % Determine what type of projection this is, commonly used 
             % 'E' or 'I'.  This information is required to allow .Target 
             % NeuralNetwork to properly load .Value as their synaptic 
             % input.
        STDP = struct('Enabled',false, ...   % Whether the real time STDP learning to going to be used or not
                      'W',[], ...            % The STDP weight, effective weight = Proj.W .* Proj.STDP.W
                      'tau_LTP',15, ...      % The LTP side time constant
                      'tau_LTD',15, ...      % The LTD side time constant
                      'dLTP',10e-3, ...      % The LTP side strength
                      'dLTD',10e-3, ...      % The LTD side strength
                      'Wmin',0.5, ...        % The minimal weight
                      'Wmax',2, ...          % The maximal weight
                      'BoundType','hard', ...% 'hard' or 'soft', decide how Wmin & Wmax limit synaptic strength change 
                      'Homeostasis',false);  % Whether whole projection needs to be normalized or not
        Value = 0 % When you invoke the method 'Project', 
                  % .Value will be updated from pre-synaptic to post-synaptic format  
        S = struct('LTP',[], ...
                   'LTD',[])
               % Filtered pre-synaptic spike train/rate for STDP learning
               % .LTP - for the potentiation part of STDP
               % .LTD - for the depression part of STDP
        UserData % Stored userdata
                 % .Wfun - the function form of W 
    end
    
    properties (Dependent)
    end
    
    methods % Functions for building kernel/weight matrix
        function K = Kernelize( P, func, varargin )
            % K = Kernelize( Proj, func, varargin )
            % 
            % Set Projection into a convolution projection.
            %
            % Inputs: Proj - An instance of Projection
            %         func - The kernel function.  It should only have one 
            %                input (relative position)
            %
            %        Options
            %        'N' - Stochasticity.  If N == 0, the weight will be 
            %              determined by its expected value
            %              If N is a positive integer, N repetitive samples
            %              will be drawn from the distribution and
            %              normalized.  default: 0
            %        'KerSize' - The Projection range of the Kernel. 
            %                    It will generate a 2*KerSize - 1 matrix as
            %                    the kernel for numeric simulation. 
            %                    Default: the size of post-synaptic network
            %                    KerSize needs to be at least 1
            %       
            % Outputs: K - Sampled kernel function (distribution)
            %
            % Jyun-you Liou 2016/12/29
            p = inputParser;
            p.CaseSensitive = false;
            addParameter(p,'N',0,@(x) x>=0 & mod(x,1)==0); % How many samples you want to draw from the distribution
            addParameter(p,'KerSize',P.Target.n); % only effective in Kernel case, and this will ignore all other things
            parse(p,varargin{:});
            % Check whether func has the right format
            if nargin(func)~=1;error('The kernel function can only have one input (relative position).');end
            % Sample from the function 'func'
            [x2,x1] = meshgrid(-p.Results.KerSize(2)+1:p.Results.KerSize(2)-1, ...
                               -p.Results.KerSize(1)+1:p.Results.KerSize(1)-1);
            K = func([x1(:),x2(:)]); % Sample it
            K = K/sum(K(:)); % Normalization
            if p.Results.N > 0  % When stochastic resampling is requested
                K = ResamplePMF( K, p.Results.N );
            end
            K = reshape(K,size(x1)); % Change it back to its shape
            % Set the corresponding parameters in P
            P.Method = 'convolution';
            P.W = K;
            P.STDP.Enabled = false; % STDP learning rule is not compatible with convolution
            P.UserData.Wfunc = func;% Save the information            
        end
        
        function W = KernelToMultiplication(P)
            % W = KernelToMultiplication(Proj)
            % 
            % Change the computation method of P from convolution to 
            % multiplication.  It will allow the projection to learn by
            % STDP rule.  Notice this operation is irreversible. (You can
            % only do Kernel -> Multiplication not vice versa)
            %
            % Inputs: Proj - An instance of Projection
            %
            % Outputs: W - the convolution matrix of the kernel
            %          Also, since Projection is a handle class, Proj.W 
            %          will be changed and Proj.Method will be properly set 
            %
            % Jyun-you Liou 2016/12/09
            P.UserData.Kernel = P.W;% Save the old data
            N = P.Target.n; % n & m are parameters to crop out the no-existing projections 
            M = size(P.W); 
            m = (M-1)/2;            
            W = convmtx2(P.W,N); % This is the raw 2D convolution matrix
            if strcmpi(P.Topology,'linear') % Crop the connectivity matrix 
                selector = false(N(1)+M(1)-1,N(2)+M(2)-1); % This is the size of raw convolution result
                selector(m(1)+1:end-m(1),m(2)+1:end-m(2)) = true; % Only select the middle part
                W = W(selector,:);
            elseif strcmpi(P.Topology,'circular') % Wrap the connectivity matrix
                Idx = 1:size(W,1); % Each row of W is operating independently to produce the vectorized result 
                [Idx1,Idx2] = ind2sub([N(1)+M(1)-1,N(2)+M(2)-1],Idx); % This is the corresponding index of conv(X,K)
                Idx1 = Idx1 - m(1);
                Idx2 = Idx2 - m(2);
                Idx1 = mod(Idx1-1,N(1))+1;
                Idx2 = mod(Idx2-1,N(2))+1;            
                Idx = sub2ind(N,Idx1(:),Idx2(:));
                w = zeros(prod(N)); % Start to wrap around
                for i = 1:prod(N)
                    w(i,:) = sum(W(Idx==i,:),1);
                end
                W = w;
            end
            % Set the corresponding parameters in P
            P.Method = 'multiplication';
            P.W = W;
        end
        
        function W = DiscretizeWeightFunction( P, func, varargin )
            % W = DiscretizeWeightFunction( P, func, varargin )
            % 
            % Sample from an function 'func' to create weight information.
            % Depending on the P.Method informtion, func would be treated 
            % as a kernel function ('convolution') or a weight function 
            % ('multiplication')
            %
            % Notice 'func' here is a function handle that takes two input
            % arguments (pre-synaptic network position and post-synaptic 
            % network position)
            %
            % Inputs: P - An instance of Projection
            %         func - a function handle, indexed by
            %                its position at pre & post-synaptic network
            %        Options
            %        'N' - Stochasticity.  If k == 0, the weight will be 
            %              determined by its expected value
            %              If k > 0 & is an integer, it will be considered 
            %              # of samples getting from pre-synaptic 
            %              projection's  distribution.
            p = inputParser;
            p.CaseSensitive = false;
            addParameter(p,'N',0,@(x) x>=0 & mod(x,1)==0); % For stochasticity
            parse(p,varargin{:});            
            % Check whether func has the right format
            if nargin(func)~=2;error('The weight function need to have two inputs specifying pre & post-synaptic neuron positions.');end
            [I1,I2] = ind2sub(P.Target.n,1:prod(P.Target.n));
            [J1,J2] = ind2sub(P.Source.n,1:prod(P.Source.n));
            I = [I1(:),I2(:)];
            J = [J1(:),J2(:)];
            for i = size(I,1):-1:1
                for j = size(J,1):-1:1
                    W(i,j) = func(I(i,:),J(i,:));
                end
            end
            % Normalization - make W column-wise sum to 1 and put the
            % strength to WPre
            P.WPre = sum(W);
            W = bsxfun(@rdivide,W,sum(W));
            if p.Results.N > 0  % When stochastic resampling is requested
                for i = 1:size(W,2)
                    W(:,i) = ResamplePMF( W(:,i), p.Results.N );
                end
            end
            % Change it to sparse matrix if > 90% of the entries are zero
            if sum(W(:)==0) > 0.9*numel(W(:))
                W = sparse(W);
            end
            % Set the corresponding parameters in P
            P.Method = 'multiplication';
            P.W = W;
        end
    end

    methods % Spike-timing dependent plasticity learnings
        RealTimeSTDPLearning(P,dt) % Learn real time by filtering method
        dW = BatchSTDPLearning(P,varargin) % Learn at once from Recorder data 
    end

    methods % Quick change of weight distributions
        P = AdjustWeight( P, varargin )        
    end
    
    methods % Execute projection
        Project(P)
    end
    
    methods % Set functions, mostly just for checking the user's input makes sense
        function set.Method(P,arg)
            if ~any(strcmpi(arg,{'convolution','multiplication','function'}))
                error('Now I only support convolution, multiplication, or function as a valid Method argument');
            end
            P.Method = arg;
        end
        
        function set.Topology(P,arg)
            if ~any(strcmpi(arg,{'linear','circular'}))
                error('Now I only support linear or circular as a valid Topology argument');
            end
            P.Topology = arg;
        end        
    end
    
    methods % Connectivity analyzing functions
        [W,STAT] = AnalyzeConnectivityByDistance(P) 
    end
    
    methods (Static)
        function P = Projection(Oi,Oj,varargin)
            if ~isa(Oi,'NeuralNetwork') || ~isa(Oj,'NeuralNetwork')
                error('Inputs need to belong to class or subclasses of NeuralNetwork.');
            end
            P.Source = Oi;
            P.Target = Oj;
            Oi.Proj.Out = [Oi.Proj.Out;P];
            Oj.Proj.In = [Oj.Proj.In;P];            
            n = P.Source.n;
            P.WPre = ones(n);
            m = P.Target.n;
            P.STDP.W = ones(prod(m),prod(n));
            P.S.LTP = zeros(n);
            P.S.LTD = zeros(m);
            % Save varargin into properties
            n_var = numel(varargin);
            for i = 1:(n_var/2)
                PropertyName = varargin{(i-1)*2 + 1};
                PropertyValue = varargin{(i-1)*2 + 2};
                P.(PropertyName) = PropertyValue;
            end
        end
    end
    
end

