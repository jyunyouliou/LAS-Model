function H = vplot(R,idx,varargin)
    % H = vplot(R,idx,varargin)
    %
    % Plot a single cell's membrane potential against time.
    %
    % R - an instance of Recorder
    % idx - index of the cell that you want to plot
    %       If you want to plot several channels at the same time,
    %       please concatenate the channel indices along the first
    %       dimension
    %       If you want to specify the cell's index by pointing at
    %       its coordinate, you need to put it in a row. The 
    %       program will try to convert the coordinate to an
    %       index based on information in R.Parent.n 
    %       For example: [2 3] - plot the cell at coordinate [2,3]
    %                    [2;3] - plot the cells of index 2 and 3
    %     
    % varargin - 'VarName' - the name of the variable that encodes
    %                        voltage, default: 'V'
    %            'VAP' - the action potential peak (unit: voltage)
    %                    default: +40
    %
    % Output: H - the graphic handle (Line) of the results.
    % 
    % Jyun-you Liou 2017/03/27
    p = inputParser;
    addParameter(p,'VarName','V');
    addParameter(p,'VAP',40);
    parse(p,varargin{:});
    if size(idx,2) > 1 % User want to specify its location by coordinates
        idx = mat2cell(idx, size(idx,1), 1+0*idx(1,:));            
        idx = sub2ind(R.Ancestor.n,idx{:});
    end 
    figure;axes;hold on;
    for i = numel(idx):-1:1
        Tvec = R.T;
        V = R.Var.(p.Results.VarName)(idx(i),:);
        % Adding spiking data
        if ~isempty(R.S)
            Selection = R.S(2,:) == idx(i); 
            tvec = R.S(1,Selection);
            v = tvec*0 + p.Results.VAP;
            Tvec = [tvec(:);Tvec(:)];
            V = [v(:);V(:)];
        end
        [Tvec,SortIdx] = sort(Tvec);
        H(i) = plot(Tvec,V(SortIdx));
    end
end
