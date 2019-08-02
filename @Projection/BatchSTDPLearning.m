function dW = BatchSTDPLearning(P,varargin)
% dW = BatchSTDPLearning(Projection, 'parameter', parameter)
%
% Perform learning at once instead of real-time learning.
%
% Learning strategy: for spiking model - we do point-process cross
%                                        correlation learning
%                    for rate model - we do realtime learning 
%
% Input: Projection - an instance of Class Projection
%        Notice, Projection.Source & Projection.Target both need to have 
%        .Recorder that contains their activity history for learning.
%
% 'parameter', parameter pairs
%
% PreIdx: Presynaptic index used for learning
% PostIdx: Postsynaptic index used for leanring
% dIdx: maximal index difference for learning
%            
% Jyun-you Liou, 2017/04/17
%
% Final update notes: allow PreIdx, PostIdx, and dIdx selection so that
% un-necessary learning can be skipped.

p = inputParser;
addParameter(p,'PreIdx',[]);
addParameter(p,'PostIdx',[]);
addParameter(p,'dIdx',[]);
parse(p,varargin{:});

%% Check the projection settings do allow the synapses to learn            
if isempty(P.STDP)
    disp('You need to set up STDP learning rule by setting Projection.STDP');return;
end

if ~strcmpi(P.Method,'multiplication')
    error('Projection.Method needs to be multiplication.');
end

if ~isa(P.Source.Recorder,'Recorder') || ~isa(P.Target.Recorder,'Recorder')
    error('The .Source and .Target of this projection does not have history record.');
end

%% Catch the NeuralNetwork objects & their Recorders

Os = P.Source; % Source network
Rs = Os.Recorder; % Source recorder
Ot = P.Target; % Target network
Rt = Ot.Recorder; % Target recorder

%% Set up what neurons are going to learn
if p.Results.PreIdx
    PreIdx = p.Results.PreIdx;
else
    PreIdx = 1:prod(Os.n);
end

if p.Results.PostIdx
    PostIdx = p.Results.PostIdx;
else
    PostIdx = 1:prod(Ot.n);
end


%% ------ Prepare rate / spiking data -------
% Source network first
if isa(Os,'RateNetwork')
    % You need to have an evenly sampled Recorder
    Rs_Tgrid = Rs.T(~isnan(Rs.T)); % Remember to remove nan (unwritten data)
    dt = unique(diff(Rs_Tgrid));
    if abs(std(dt)/mean(dt))>10e-6
        error('Not an evenly sampled Recorder data');
    else
        dt = mean(dt);
    end                
    pre = Os.param.f(Rs.Var.V(:,1:Rs.Idx-1) - Rs.Var.phi(:,1:Rs.Idx-1));
elseif isa(Os,'SpikingNetwork')
    for i = prod(Os.n):-1:1
        if isempty(Rs.S);error('No spike has been formally recorded into source network .S');end
        pre{i,1} = Rs.S(1,Rs.S(2,:)==i);
    end
end

% Target network
if Ot == Os % If they are the same, then great! We don't need to recalculate
    post = pre;
else % If they are not the same
    if isa(Ot,'RateNetwork')
        % The two Recorder need to have exactly the same time record                  
        if any(Rs.T(1:Rs.Idx-1) ~= Rt.T(1:Rt.Idx-1)) % Use Idx here so that unwritten part (nan) will be ignored
            error('.Source & .Target.Recorder do not have the same time grid.');
        end
        post = Ot.param.f(Rt.Var.V(:,1:Rt.Idx-1) - Rt.Var.phi(:,1:Rt.Idx-1));
    elseif isa(Ot,'SpikingNetwork')
        for i = prod(Ot.n):-1:1
            if isempty(Rt.S);error('No spike has been formally recorded into target network .S');end                        
            post{i,1} = Rt.S(1,Rt.S(2,:)==i);
        end
    end    
end

%% Learning 
if isa(Os,'RateNetwork') && isa(Ot,'RateNetwork')
    % If both Ot & Os are RateNetwork - then do filter learning, which is
    % usually faster    
    pre_filter = 0;
    post_filter = 0;
    W_sub = P.STDP.W(PostIdx,PreIdx);
    for i = 1:numel(Rs_Tgrid)
        pre_filter = pre_filter .* exp(-dt./P.STDP.tau_LTP) + dt.*pre(PreIdx,i)./P.STDP.tau_LTP;
        post_filter = post_filter .* exp(-dt./P.STDP.tau_LTD) + dt.*post(PostIdx,i)./P.STDP.tau_LTD;
        if strcmpi(P.STDP.BoundType,'soft') % I did not write it in exponential integration form so that it is easier to understand
            dW = (P.STDP.dLTP .* (P.STDP.Wmax - W_sub) .* (dt * post(PostIdx,i) * pre_filter')  - ... % LTP-part
                  P.STDP.dLTD .* (W_sub - P.STDP.Wmin) .* (post_filter * pre(PreIdx,i)' * dt)) .* dt; % LTD-part
        elseif strcmpi(P.STDP.BoundType,'hard')
            dW = (P.STDP.dLTP .* (P.STDP.Wmax > W_sub) .* (dt * post(PostIdx,i) * pre_filter')  - ... % LTP-part
                  P.STDP.dLTD .* (W_sub > P.STDP.Wmin) .* (post_filter * pre(PreIdx,i)' * dt)) .* dt; % LTD-part                
        end
        W_sub = W_sub + dW;        
        if mod(100*i,numel(Rs_Tgrid))<100
            disp(['Learning ... ' num2str(100*(i/numel(Rs_Tgrid))) ' %'])
        end
    end
    P.STDP.W(PostIdx,PreIdx) = W_sub;
    
else % Either one or both are SpikingNetwork
    % Warning information about learning from spiking network
    disp('To learn for spiking network, all the Recorder.S needs to be integers.');
    disp('Please adjust the time and spatial unit accordingly.')
    disp('');
    disp('If there is one rate network involved,');
    disp('change the rate network''s recorder time unit to make every entry of Recorder.T integer, too.')
    dt = 1;
 
    % Construct the STDP learning kernel
    % We need to separate the potentiation & depression parts 
    % because we want to be able to implement soft bound         
    % We want to make the integral of K_LTP and K_LTD's integral = 1 so
    % that K_LTP still preserve rate unit
    K_LTP = @(t) (t>0).*exp(-t./P.STDP.tau_LTP)./P.STDP.tau_LTP;  
    K_LTD = @(t) (t<0).*exp(+t./P.STDP.tau_LTD)./P.STDP.tau_LTD;
    n_LTP = ceil(2.5*P.STDP.tau_LTP/dt);
    n_LTD = ceil(2.5*P.STDP.tau_LTD/dt);
    
    % Do cross-correaltion
    X_LTP = zeros(prod(Ot.n),prod(Os.n));
    X_LTD = zeros(prod(Ot.n),prod(Os.n));
    XC_LTP = cell(prod(Ot.n),prod(Os.n));
    XC_LTD = cell(prod(Ot.n),prod(Os.n));
    
    for i = PostIdx
        for j = PreIdx
            % In some situation, you don't need to calculate
            % cross-correaltion to save time.
            if p.Results.dIdx % If some 'skipping rules' are specified, skip them
                if abs(j - i)>p.Results.dIdx
                    XC_LTP{i,j} = 0;
                    XC_LTD{i,j} = 0;
                    continue;
                end
            end
            if Os == Ot % If it is recurrent connection, only calculate half of the matrix
                if j < i
                    continue;
                end
            end
            
            % LTP part
            if iscell(pre) && ~iscell(post)
                XC_LTP{i,j} = pcxcorr(pre{j},post(i,:),n_LTP);
            elseif ~iscell(pre) && iscell(post)
                XC_LTP{i,j} = fliplr(pcxcorr(post{i},pre(j,:),n_LTP));
            else % Both are SpikingNetwork
                XC_LTP{i,j} = pxcorr(post{i},pre{j},n_LTP);               
            end

            % LTD part            
            if n_LTP == n_LTD
                XC_LTD{i,j} = XC_LTP{i,j};
            else
                if iscell(pre) && ~iscell(post)
                    XC_LTD{i,j} = pcxcorr(pre{j},post(i,:),n_LTD);
                elseif ~iscell(pre) && iscell(post)
                    XC_LTD{i,j} = fliplr(pcxcorr(post{i},pre(j,:),n_LTD));
                else % Both are SpikingNetwork
                    XC_LTD{i,j} = pxcorr(post{i},pre{j},n_LTD);
                end
            end
        end
        disp(['Learning ... ' num2str(100*(1-i/size(post,1))) ' %'])
    end

    if Os == Ot
        for i = size(post,1):-1:1
            for j = size(pre,1):-1:1
                if j < i
                    XC_LTP{i,j} = fliplr(XC_LTP{j,i});
                    XC_LTD{i,j} = fliplr(XC_LTD{j,i});
                end
            end
        end
    end    
    
    % Calculate cross-correlation .* kernel
    for i = size(post,1):-1:1
        for j = size(pre,1):-1:1
            X_LTP(i,j) = sum(K_LTP((-n_LTP:n_LTP)*dt).* XC_LTP{i,j}*dt);
            X_LTD(i,j) = sum(K_LTD((-n_LTD:n_LTD)*dt).* XC_LTD{i,j}*dt);  
        end
    end

    %%
    if strcmpi(P.STDP.BoundType,'soft') % I did not write it in exponential integration form so that it is easier to understand
        dW = (P.STDP.dLTP .* (P.STDP.Wmax - P.STDP.W) .* X_LTP  - ... % LTP-part
              P.STDP.dLTD .* (P.STDP.W - P.STDP.Wmin) .* X_LTD) * dt; % LTD-part
    elseif strcmpi(P.STDP.BoundType,'hard')
        dW = (P.STDP.dLTP .* (P.STDP.Wmax > P.STDP.W) .* X_LTP  - ... % LTP-part
              P.STDP.dLTD .* (P.STDP.W > P.STDP.Wmin) .* X_LTD) * dt; % LTD-part                
    end
    P.STDP.W = P.STDP.W + dW;
    
    % Warning - if bound is reached
    if any(P.STDP.W(:)>P.STDP.Wmax)
        warning('Upper bound is reached.  Suggest either to decrease learning rate or decrease batch size.');
        warning('Please reset .W to fullfill its upper/lower bounds manually if necessary.');
    elseif any(P.STDP.W(:)<P.STDP.Wmin)
        warning('Lower bound is reached.  Suggest either to decrease learning rate or decrease batch size.');
        warning('Please reset .W to fullfill its upper/lower bounds manually if necessary.');        
    end    
end

%% Homeostasis by pre-synaptic projections
if P.STDP.Homeostasis
    disp('Homeostasic is applied.');
    P.STDP.W = bsxfun(@rdivide, P.STDP.W, sum(P.STDP.W)); 
end            

end 



%% pxcorr for point process cross-correlation
% Code source: https://github.com/VincentToups/matlab-utils/blob/1db2242f4813b07fe763faa150d67385a132ef8f/chronux/spikesort/utility/datatools/private/CORE_pxcorr.m
% http://procyonic.org/
% Original author: Jonathan Vincent Toups
function [Z,lag_idx] = pxcorr(x,y,maxlag)
    %   Z = PXCORR(X,Y,MAXLAG), where X and Y are point process data 
    %   returns the length (2*MAXLAG+1) vector
    %   Z such that Z(lag(i)) = #(X(j)-Y(k) == i), |i| <= MAXLAG
    %
    %   CONDITIONS
    %   ----------
    %   X and Y must contain sorted integer values. (Ascending)
    %   X and Y must be row vectors of type DOUBLE.
    x = sort(x);
    y = sort(y);

    lag_idx = -maxlag:maxlag;
    Z = zeros(1,2*maxlag+1);

    % If either of the data is empty, just return empty 
    if isempty(x) || isempty(y);return;end
    % We do this efficiently by stepping along X and keeping track of those
    % indices in Y that are near by (i.e., within a distance of MAXLAG).
    limit = length(x);
    a = 1;  c = 1;   

    for b = 1:length(y)
        while((y(b)-x(a)) > maxlag),        % move left bracket until less than MAXLAG before x(b)
                a=a+1;   if (a > limit), return; end;
        end
        if (c < a), c = a; end;             % catch up the right bracket if it was passed by
        if (c <= limit)
            while((x(c)-y(b)) <= maxlag),   % move right bracket until more than MAXLAG after x(b)
                c=c+1;   if (c > limit), break; end;
            end
        end

        offset = -y(b)+maxlag+1;            % add 'em up
        for bb = a:(c-1)
            ind = x(bb)+offset;
            Z(ind) = Z(ind) + 1;
        end
    end

    return;

end

function Z = pcxcorr(y,x,d_index) 
    % This is for one input is point process and another one is a time series
    % y - event data index
    % x - time series data
    Z = zeros(d_index*2+1,1);
    for m = 1:numel(y)
        Z = Z + x(y(m)-d_index:y(m)+d_index);
    end
    Z = Z/numel(y);  
end