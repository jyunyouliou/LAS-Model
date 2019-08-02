classdef NeuralNetwork < handle & matlab.mixin.Heterogeneous
    % NeuralNetwork is a Superclass of fundamental building unit of this 
    % simulation package.  To build a functional NeuralNetwork, you need to
    % call any of the subtypes of NeuralNetwork to that the computation
    % methods, ... etc. will be loaded.  Directly call NeuralNetwork won't
    % give a functional NeuralNetwork.
    %
    % Properties that are shared by all types of neural networks are stored
    % and described here, including
    %
    % ##### (Basic information) #####
    % .n - size of the network
    % .t - current time of the network
    %
    % ##### (Association with instances of other classes) #####
    % .Proj - Record the information about the connection topology between 
    %         different NeuralNetworks.  .Proj is a structure with two 
    %         fields - .In and .Out - which include instances of Class 
    %         Projection 
    % .Ext - Additional input, must be an instance of Class ExternalInput
    % .Recorder - an instance of Class Recorder that can be written
    %
    % ##### Input/Output variables ##### (Changing during simulation)
    % .Input - a structure with its field names indicate the type of 
    %          post-synaptic channel/conductance/method it is using. 
    %          Most commonly Input.I & Input.E code for inhibitory and
    %          excitatory inputs.  For each update, Input will be consumed
    %          by the NeuralNetwork, and the results will be projected to
    %          other networks and save in .Input again.    
    % 
    % ##### Supplementary properties #####
    % .Graph - the associated realtime plot figure
    % .UserData
    %
    % ##### Dependent properties ##### 
    % .dim - dimensionality of the neural network, most commonly 1 or 2
    %
    % Both GeneralModel & SpikingModel inherit from this class.
    % 
    % Jyun-you Liou, 2017/11/14 
    %
    % Update note: dim function was incorrect previously
    
    properties (SetObservable=true)
        % Basic description of the network
        n % network size
        t = 0 % current time
        
        % Inter-model connections (not necessary for this paper)
        Input % A structure that records additional E and I projections 
              % the neural fields project to it.  This field will be
              % assimilated into the field and cleared every update cycle
              % Input needs to be a structure, and its fieldnames should 
              % match the 'Type' of Projection it receives
        Proj = struct('In',[], ...
                      'Out',[]);
              % Instances of 'Projection' class.  In & Out mark the
              % network's use those projections as input or output.
              % Each entry corresponds to one group of projections.
              % Please see 'Projection' to see its properties & methods.
               
        % External input
        Ext % An instance of Class ExternalInput 
        
        % Recorder - an instance of Class Recorder, recording the 
        %            simulation results of a NeuralNetwork
        %            Make it by Recorder(NeuralNetwork, Capacity)
        Recorder % Please see Recorder to see its properties and methods      
        
        % Associated figure handle for real time plot
        Graph
        
        % UserData - for User to save note ... etc. 
        UserData        
    end
    
    properties (Dependent)
        dim
    end
    
    methods % For dependent variables
        function d = get.dim(O)
            s = size(O.n);
            if numel(s) > 2
                d = numel(s);
            else
                if prod(O.n) == max(O.n)
                    d = 1;
                else
                    d = 2;
                end
            end
        end
    end

    methods % Simulation functions
        function Update(O,dt) % Vectorized function for UpdateSingle
            % Update(NeuralNetwork,dt)
            %
            % Simulate NeuralNetwork models with a time step dt.
            % NeuralNetwork can be a vector (containing several instances 
            % of NeuralNetwork).  Depends on the specific subtype of the
            % NeuralNetwork, the function calls its corresponding method to
            % execute its update one by one. Namely, every SubClass of 
            % NeuralNetwork should have a method 'IndividualModelUpdate'
            % and 'Project' to be called.
            %
            % Jyun-you Liou 2016/12/22
            if dt <= 0;warning('The model is running backward.');end            
            % Step 1: Pass through the field box
            for i = 1:numel(O)
                IndividualModelUpdate(O(i),dt);
            end
            % Step 2: Calculate its output and save it for additional input
            for i = 1:numel(O) % Each NeuralNetwork
                for j = 1:numel(O.Proj.Out) % Each Projection of a NeuralNetwork
                    Project(O(i).Proj.Out(j));
                end
            end
            % Step 3: Unsupervised learning 
            for i = 1:numel(O) % Each NeuralNetwork
                for j = 1:numel(O.Proj.Out) % Each Projection of a NeuralNetwork
                    if O(i).Proj.Out(j).STDP.Enabled % If the switch is true, then learn.
                        RealTimeSTDPLearning(O(i).Proj.Out(j),dt);
                    end
                end
            end            
        end        
    end

    methods % Plot functions
        function f = plot(O)
            % plot(NeuralNetwork)
            %
            % Plot current state of NeuralNetwork. Depends on the specific 
            % subtype of the NeuralNetwork, the function calls its 
            % corresponding 'IndividualModelPlot' method to execute one by 
            % one. 
            %
            % Jyun-you Liou 2016/12/22
            for i = numel(O):-1:1
                f(i) = IndividualModelPlot(O);
            end
        end
        
        function f = IndividualModelPlot(O) 
            % Class NeuralNetwork's default plotting function
            % figure_handle = IndividualModelPlot(NeuralNetwork)
            %
            % NeuralNetwork can only be a scalar instance.
            %
            % IndividualModelPlot will produce a figure with n axes, with
            % each axes corresponds to each dynamic variable.  When the
            % figure already exists, it will find the axes automatically
            % and update the plot.
            %
            % For 1 dimensional model, it will produce a variable versus
            % space plot.  AttachHotKey will allow you to manipulate the
            % variables. 
            %
            % For 2 dimensional model, it will produce a heat map over the
            % two dimensional space.  Up to now, the function to manipulate
            % the dynamic variables are not built, yet.
            %
            % For future developer - if you want to define your own plot
            % function, please define it in the subclass classdef 
            % file/folder.
            %
            % Jyun-you Liou, 2017/03/13
            
            % Determine if the figure has been created
            if isempty(O.Graph) || ~isvalid(O.Graph)
                f = figure('Unit','normalized','Position',[0.25 0.05 0.5 0.85]);
                O.Graph = f;
                for j = numel(O.VarName):-1:1 % In this way, Ax(1) will be the most top axes
                    Ax(j) = subplot(numel(O.VarName),1,j);
                    Ax(j).NextPlot = 'replacechildren';                              
                    Ax(j).Tag = O.VarName{j}; % Put a tag here so that you can findobj to adjust your plot later
                    title(O.VarName{j});                    
                end
            else
                f = O.Graph;
                % Find the right axes back according to their .Tag
                for j = numel(O.VarName):-1:1
                    Ax(j) = findobj(f,'Tag',O.VarName{j},'Type','Axes');
                end
            end
            
            % Determine the dimension of the model then choose a
            % corresponding method to plot.         
            if any(O.n==1) % one-dimensional case
                for j = 1:numel(O.VarName)
                    if isempty(Ax(j).Children)
                        plot(Ax(j),O.(O.VarName{j}),'Tag',['Data_' O.VarName{j}]);
                    else
                        Obj = findobj(Ax(j).Children,'Tag',['Data_' O.VarName{j}]);
                        Obj.YData = O.(O.VarName{j});
                    end
                end               
            elseif numel(O.n) == 2 && ~any(O.n==1) % two-dimensional case
                for j = 1:numel(O.VarName)
                    if isempty(Ax(j).Children)
                        imagesc(Ax(j),O.(O.VarName{j}),'Tag',['Data_' O.VarName{j}]);
                    else
                        Obj = findobj(Ax(j).Children,'Tag',['Data_' O.VarName{j}]);                        
                        Obj.CData = O.(O.VarName{j});
                    end
                end  
            end
        end        
    end
    
    methods % Recorder function
        function R = CreateRecorder(O,N)
            % Record = CreateRecorder(NeuralNetwork,Capacity)
            %
            % Input: NeuralNetwork - instances of class NeuralNetwork
            %                        ,which includes GeneralModel and 
            %                        SpikingModel
            %        Capacity - number of savings you can do.
            %
            % Final update: 2016/12/26
            for i = 1:numel(O)
                if nargin < 2 % Use 5% of RAM to record 
                    [~,sV] = memory;
                    N = 0.05*sV.PhysicalMemory.Available/8; % Available slots
                    N=N/ prod(O(i).n);
                    N=1000*floor(N/1000); % Make sure the Capacity is counted in 'thousands'
                end
                O(i).Recorder = Recorder(O(i),N); %#ok<CPROPLC>
            end
            R = vertcat(O.Recorder);
        end

        function WriteToRecorder(O)
            % WriteToRecord(NeuralNetwork)
            %
            % Write current data to .Recorder ,this method calls Recorder.
            % Record
            %
            % Final update: 2016/12/23
            for i = 1:numel(O)
                % If .Recorder has not been set up, create it                
                if isempty(O(i).Recorder)
                    CreateRecorder(O(i));
                end
                for j = 1:numel(O(i).Recorder)
                    Record(O(i).Recorder(j));
                end
            end
        end        
    end

    methods % Build connections
        function P = Link(Oi,Oj,varargin)
            % P = Link(Oi,Oj,'parameter',parameter)            
            % 
            % Build a projection from Oi to Oj.  
            %
            % This function is exactly the same as Projection(Oi,Oj)
            %
            % Optional inputs: Please see help Projection
            P = Projection(Oi,Oj,varargin{:});
        end
    end

    methods % House-keeping functions such as rounding-off to prevent small numbers
        function RoundOff(O,p)
            % RoundOff(NeuralNetwork,precision)
            %
            % Rounding off the numerics of dynamic variables 
            % 
            % O - instances of NeuralNetwork
            % p - precision. # of numbers behind decimals you want to preserve 
            %     default: 6 digits
            if nargin < 2
                p = 6;
            end            
            for i = 1:numel(O)
                % Round off dynamic variables
                Var = O(i).VarName;
                for j = 1:numel(Var)
                    O.(Var{j}) = round(O.(Var{j}) * 10^p) / 10^p;
                end
                % Round off rate/spiking-derived variables
                if isfield(O,'S')
                    FieldNames = fieldnames(O(i).S);
                    for j = FieldNames'
                        O(i).S.(j{:}) =  round(O(i).S.(j{:}) * 10^p) / 10^p;
                    end
                elseif isfield(O,'R')
                    FieldNames = fieldnames(O(i).R);
                    for j = FieldNames'
                        O(i).R.(j{:}) =  round(O(i).R.(j{:}) * 10^p) / 10^p;
                    end                    
                end
            end
        end
    end    

    methods % Graphic interface functions, including HotKeys etc.
        F = AttachHotKey( O, varargin )
    end
    
end