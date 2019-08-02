function Rcom = Compress(R,r)
    % Rcom = Compress(R,r)
    %
    % Compress the recorder data from R to R.Children
    % If R.Children is empty, a new Recorder will be created to be
    % put into R.Children
    %
    % Input: R - an instance of Recorder
    %        r - compression ratio (optional)
    %            If this argument is left blank. The program will
    %            automatically try to find a divisor that is
    %            closest to 10 as compression ratio.
    % 
    % Jyun-you Liou, final update: 2016/12/24

    % Call recursive functions if there are several instances in R
    if numel(R) > 1
        for i = numel(R):-1:1
            if nargin == 2
                Rcom(i) = Compress(R,r);
            else
                Rcom(i) = Compress(R);
            end
        end
        return;
    end            

    % Check whether R.Children has been created
    if isempty(R.Children) || ~isvalid(R.Children)
        Original_Recorder = R.Ancestor.Recorder; % We don't want to change the hierarchy
        R.Children = Recorder(R.Ancestor,R.Capacity);
        ExtraVarName = setdiff(R.VarName,R.Children.VarName);
        if ~isempty(ExtraVarName)
            AddVar(R.Children,ExtraVarName{:});
        end
        R.Children.Parent = R;                
        R.Ancestor.Recorder = Original_Recorder;
    end
    Rcom = R.Children;

    % Check whether the compression ratio 'r' can avoid
    % intepolation.              
    N = R.Capacity;            
    if nargin < 2 % Automatic suggesting compression ratio (I avoid divisors so that users do not need symbolic toolbox)
        pf = factor(N);         % Prime factors of n
        upf = unique(pf);       % Unique
        d = 1;
        for f = upf
            d = d(:)*(f.^(0:1:sum(pf == f)));
        end
        r_goal = 10;
        [~,idx_choice] = min(abs(d(:)-r_goal));
        r = d(idx_choice);
        disp(['Compression ratio is set to be ' num2str(r)]);
    end
    if mod(N,r);warning('Compress ratio needs to divisible to Capacity to prevent error or loss of data.');end

    % Compress time vector
    Tvec = R.T(r:r:end);            
    n = numel(Tvec);
    Rcom.T(Rcom.Idx:Rcom.Idx+n-1) = Tvec;

    % Compress each dynamical variables & extrafields
    % I avoid to use the function 'smooth' so that the user does
    % not need to have Curve Fitting Toolbox
    width = r + mod(r,2) - 1; % moving average window size, force it to be odd            
    for Name = R.VarName(:)'
        x = R.Var.(Name{:})'; % Original data
        x = filter(ones(width,1)/width,1,x); % moving average, improved from 'smooth' function
        xbegin = cumsum(x(1:width-2,:)); % Shrink the window at the data beginning & ends
        xbegin = bsxfun(@rdivide, xbegin(1:2:end,:), (1:2:(width-2))');
        xend = cumsum(x(R.Capacity:-1:R.Capacity-width+3,:));
        xend = bsxfun(@rdivide, xend(end:-2:1,:),(width-2:-2:1)');
        x = [xbegin;x(width:end,:);xend];
        x = x(r:r:end,:)'; % downsample
        Rcom.Var.(Name{:})(:, Rcom.Idx:Rcom.Idx+n-1) = x; % reshape and vertcat
    end
    Rcom.Idx = Rcom.Idx + n;


    % Deal with spiking data if there is spiking records
    TransferSpikeRecord(R);
    Rcom.S = [Rcom.S,R.S];
end