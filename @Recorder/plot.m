function H = plot(R, dN) 
    % H = plot(Recorder, dN)
    %
    % Plot all information about the recorder.  
    %
    % Input: Recorder - an instance of Class Recorder
    %        dN - the number of steps you want to plot from current
    %             time 
    %
    % The first figure will display all dynamical variables,
    % recorded in .Var.  Horizontal axis is time and vertical axis 
    % is cell index.  The figure handle will be saved in
    % R.Graph.Figure
    %
    % The second figure is only displayed if the NeuralNetwork of 
    % this Recorder is a spiking model. The data needs to be
    % transfered from .SBuffer to .S first by TransferSpikeRecord 
    % first.  Again, the horizontal axis is time, and the vertical
    % axis is cell index.  The vertical axis will be saved in 
    % R.Graph.Raster.
    %
    % Jyun-you Liou, 2017/05/07 ( Allow raster-only condition )

    % The figure for displaying dynamical variables
    if ~isempty(R.VarName)
        if isempty(R.Graph.Figure) || ~isvalid(R.Graph.Figure) % If none of the figure has been established
            R.Graph.Figure = figure('Unit','Normalized','Position',[0 0.3 0.8 0.4]);
            N_ax = numel(fieldnames(R.Var));
            for i = N_ax:-1:1
                ax(i) = subplot(1,N_ax,i);
                ax(i).YDir = 'normal';title(R.VarName{i});ax(i).NextPlot = 'replacechildren';
                if i ~= 1;ax(i).YTick = [];end
            end
        else % If the figure has been established, just grab them back
            figure(R.Graph.Figure); % move to the main figure for plotting
            ax = R.Graph.Figure.Children;                
        end

        % Decide what range of data needs to be pulled out
        if nargin == 1 || isempty(dN) % If dN is not given, plt the whole thing until current step
            I = 1:R.Idx-1;
        else % if dN is given
            I = (R.Idx-dN):(R.Idx-1);
            I(I < 1) = [];
        end

        % Pull out data and imagesc them to each axes
        for i = 1:numel(ax) 
            CData = R.Var.(R.VarName{i})(:,I);
            Tvec = R.T(I);
            % Check whether the image has been established
            if isempty(ax(i).Children) 
                imagesc(ax(i),Tvec,[1 size(CData,1)],CData);
            else
                set(ax(i).Children,'XData', Tvec, 'CData', CData);
            end
            ylim(ax(i),[1,size(CData,1)]);
            if ~all(isnan(Tvec)) % To prevent error while ploting a completely nascent Recorder
                xlim(ax(i),[min(Tvec) max(Tvec)]);
            end
        end
    end
    
    % The figure for displaying raster plot if it is a spiking model
    if ~isempty(R.S) 
        if isempty(R.Graph.Raster) || ~isvalid(R.Graph.Raster)
            R.Graph.Raster = figure('Unit','Normalized','Position',[0.8 0.3 0.2 0.4]);
            ax = axes;xlabel('time');title('Raster plot');
            ax.NextPlot = 'replacechildren';
        end
        ax = R.Graph.Raster.Children; % Move plotting to the raster figure
        X = R.S(1,:);
        Y = R.S(2,:);
        if nargin > 1 && ~isempty(dN)
            Selector = (X>min(Tvec)) & (X<max(Tvec));
            X = X(Selector);
            Y = Y(Selector);
        end
        Y = repmat(Y,[3,1]);
        X = repmat(X,[3,1]);
        Y = bsxfun(@plus,Y,[-0.5;0.5;0]);
        Y(3,:) = nan;        
        if isempty(ax.Children) % Check whether the raster has existed
            plot(ax,X(:),Y(:));
        else
            set(ax.Children,'XData',X(:),'YData',Y(:));
        end
    end
    if nargout > 0
        H = R.Graph;
    end
end
