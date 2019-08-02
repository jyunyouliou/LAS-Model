classdef Recorder < handle
    % Recorder is a class that records the simulation results of
    % NeuralNetwork.  To create a Recorder that is associated with an 
    % instance of NeuralNetwork, use CreateRecorder(NeuralNetwork,Capacity)
    % 
    % To write to a Recorder, use WriteToRecorder(NeuralNetwork). 
    % 
    % Jyun-you Liou, 2016/12/25
    
    properties (AbortSet = true, SetObservable=true)
        Idx = 1 % the index for the next saving
        Capacity % Capacity of this recorder        
    end
    
    properties
        % .Parent keeps track where the data are recorded from
        % There are two possible classes of a Recorder's Parent
        %    - NeuralNetwork
        %    - Recorder
        % If .Parent is a NeuralNetwork, it indicates the data is
        % directly recorded from a simulating NeuralNetwork.  
        % If .Parent is a Recorder, it indicates the data is from
        % compressing a finer scale Recorder        
        Parent 
        
        % .Var is a structure used for recording dynamic variables
        % Default: it will record all dynamic variables from the 
        %          NeuralNetworks it records from.  
        % You can add more variables to record by AddVar(Recorder,'VarName')
        Var 
        
        T % Time record, an 1 by Capacity vector
        SBuffer  % For incoming spiking data, an 1 by nt cell array 
        
        % .S is for recording spikes.  
        % First you need to transfer data from .SBuffer to .S by method 
        % 'TransferSpikeRecord'.  
        % .S is a 2 by n_spike numeric matrix. The first row is time, and 
        % the second row is the neuron's index.
        S 
        
        % .Graph is a structure that stores the associated real-time graph.
        % It has two fields: .Figure (for .Var) and .Raster (for .S)
        Graph = struct('Figure',[], ...
                       'Raster',[]);
                   
        Listener % The listeners associated with this Recorder
        
        % .Children stores the compressed version of this recorder.
        % The original NeuralNetwork can be recursively trace .Parent, or 
        % simply use .Ancestor to find it 
        Children 
        
        UserData % For user-defined notes
    end
    
    properties (Dependent)
        % .VarName perform a quary looking for names of variables the 
        % recorder is recording.  Default: the dynamic variables of the 
        % .Ancestor NeuralNetwork           
        VarName 
        
        % .Ancestor will invoke a recursive function to go through .Parent
        % all the way to the original NeuralNetwork        
        Ancestor 
    end
    
    methods % Record
        function Record(R) % Record from R.Parent
            % Record(Recorder) - Record from its Parent NeuralNetwork
            O = R.Parent;
            R.T(R.Idx) = O.t;
            for N = R.VarName
                if any(strcmp(N{:},fieldnames(O)))
                    % Recoding the dynamic variables
                    R.Var.(N{:})(:,R.Idx) = O.(N{:})(:);
                else
                    % Recording the user-defined variables
                    RecordExtraVar(R,N{:});
                end
            end
            if isa(R.Parent,'SpikingNetwork') % Check if it is a SpikingNetwork
                R.SBuffer{1,R.Idx} = find(O.S.S);
            end
            R.Idx = R.Idx + 1;
        end
    end
    
    methods % Trasnfer spike record from temporary .SBuffer to permanent .S
        function TransferSpikeRecord(R)
            if numel(R)>1
                for i = 1:numel(R)
                    TransferSpikeRecord(R(i))
                end
                return;
            end
            if isempty(R.SBuffer);disp('No spike has been recorded');return;end
            s = R.SBuffer(:)'; % Make sure the concatenation can work as a row
            s = cellfun(@(x) x(:)',s,'un',false);
            n = cellfun(@numel,s);
            I = find(n);
            if isempty(I);return;end; % No spikes detected at all.
            for i = numel(I):-1:1
                tp{1,i} = repmat(R.T(I(i)),[1 n(I(i))]);
            end
            tp = cell2mat(tp);
            s = cell2mat(s);
            R.S = [R.S,[tp;s]];
            R.SBuffer = cell(1,R.Capacity);
        end
    end
    
    methods % Data compression - see 'Compress' & 'Refresh' in the class folder
    end

    methods % For AddVar (Additional variables you want to record)
        function AddVar(R,varargin)
            % AddVar(Recorder,'ExtraVarName1','ExtraVarName2', ...)
            % 
            % AddExtraVar is a method for Recorder class.  It allows the 
            % user to freely add any variable required to be recorded if 
            % the corresponding methods are defined.  
            % 
            % Now all additional Var needs to be time series data 
            % (numeric matrix).  The next expansion will allow the users to
            % record point process data.
            %
            % Jyun-you Liou, 2017/03/27
            for i = 1:numel(varargin)
                if any(strcmp(varargin{i},R.VarName))
                    warning(['The variable ' varargin{i} ' has been defined.  No need to redefine it.']);
                    continue;
                end
                R.Var.(varargin{i}) = nan(prod(R.Parent.n),R.Capacity);
            end
        end
    end
    
    methods % Plotting methods  - see 'plot' & 'vplot' in separated files
    end
    
    methods % For dependent properties
        function VarName = get.VarName(R)
            VarName = fieldnames(R.Var);
            VarName = VarName(:)';
        end
        
        function S = get.Ancestor(R)
            S = R.Parent;
            while ~isa(S,'NeuralNetwork')
                S = S.Parent;
            end
        end
    end
    
    methods % Listener functions        
        function CapacityChange(R) % Triggered while .Capacity is changed 
            C = R.Capacity; % Target Capacity (Remember it is PostSet trigger)
            c = numel(R.T); % Original Capacity
            if c > C
                warning('The new capacity is smaller than original capcity, this operation may cause loss of data.');
                ANS = input('Are you sure to proceed? If you do, ensure you have safely transfer the original data first. (1/0)');
                if ANS
                    R.T(C+1:c) = [];
                    R.SBuffer(:,C+1:c) = [];
                    for m = R.VarName
                        R.Var.(m{:})(:,C+1:c) = [];
                    end
                    if R.Idx > C;warning('Writing Idx has been changed.');R.Idx = C+1;end 
                else
                    R.Capacity = c;
                end
            elseif c < C
                R.T(c+1:C) = nan; 
                R.SBuffer{end,C} = [];
                for m = R.VarName
                    R.Var.(m{:})(:,c+1:C) = nan;
                end
            end
        end
        
        function AutoCompressCheck(R)
            if R.Idx > R.Capacity
                disp('Capacity is full ... Compressing data to Recorder.Children');
                Compress(R);
                Refresh(R);
            end
        end
    end
    
    methods % Load from the recorder
        function Load(R,Idx)
            % Load(Recorder, Idx)
            % 
            % Load from the Recorder to its .Ancestor NeuralNetwork from
            % the data recorded at index 'Idx'.  If Idx is not given, it
            % will automatically load from the Idx = 1.
            %
            % Jyun-you Liou 2017/01/07
            if nargin < 2 || isempty(Idx)
                Idx = 1;
            end
            if numel(R)>1
                for i = 1:numel(R)
                    Load(R,Idx);
                end
                return;
            end
            % Start of the loading function
            O = R.Parent;
            for N = R.VarName
                O.(N{:}) = R.Var.(N{:})(:,Idx);
                O.(N{:}) = reshape(O.(N{:}),O.n);
            end
            O.t = R.T(Idx);
            R.Idx = Idx+1;
        end
    end
    
    methods (Static) % Recorder creating function
        function R = Recorder(O,N)
            % Record = Recorder(NeuralNetwork,Capacity)
            %
            % Input: NeuralNetwork - instances of class NeuralNetwork
            %                        ,which includes GeneralModel and 
            %                        SpikingModel
            %        Capacity - number of savings you can do. 
            %
            % Final update: 2016/12/12
            if nargin == 0;return;end
            if numel(O) > 1;error('To create multiple recorders for a bunch of NeuralNetwork, please use CreateRecorder(NeuralNetwork,Capacity)');end
            if ~isa(O,'NeuralNetwork');error('Recorder can only be created and associated with a NeuralNetwork.');end
            % Register them together
            O.Recorder = R; 
            R.Parent = O;
            % Capacity & time vector
            if nargin < 2;N=ceil(10e6/prod(O.n));end
            R.Capacity = N;
            R.T = nan(1,N);  
            % Look at the type of NeuralNetwork and assign the
            % corresponding way to record it.
            % Dynamical variables
            for i = O.VarName % Dynamical variables
                R.Var.(i{:}) = nan([prod(O.n),N]);                    
            end
            % Firing variables - Recorder only needs to remember .SBuffer (spikes)
            %                    because: .f is a function of V & phi
            %                             .dT, .u, .x, are just filtered
            %                             version of .SBuffer
            R.SBuffer = cell(1,N); % Each cell record what neurons fire.  The 
                             % index of the cell array indicates its index 
                             % in .T vector 
            % Default listener: if .Idx is larger than capacity, it will
            % automatically trigger method 'Compress' 
            R.Listener = addlistener(R,'Idx','PostSet',@(~,evd) AutoCompressCheck(evd.AffectedObject));
            R.Listener(2) = addlistener(R,'Capacity','PostSet',@(~,evd) CapacityChange(evd.AffectedObject));         
        end 
    end
end