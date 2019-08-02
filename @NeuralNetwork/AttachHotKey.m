function F = AttachHotKey( O, varargin )
% F = AttachHotKey(NeuralNetwork, 'parameter', parameter)
%
% This function will attach several Callback functions to the 
% NeuralNetwork's Graph (figure).  Several hot keys will be created
% to allow users to control the experimental flow.
%
% Hot Keys
% 'space' - pause & resume
% 'q' - quit simulation
% 's' - set into stimulation mode
% 'z' - set into zap mode, only effective in 1D models
% 'c' - set into clamp mode, only effective in 1D models 
% 'd' - set into delta mode, only effective in 1D models
% 'n' - nullify all stimulation/clamp interventions
% 'r' - restart from the beginning
% 'p' - set parameter changes
% 'f' - flag.  This allows you to define any special intervention you want 
%       to do during experiments.  You need to define its action in
%       Experiment files.  Flag will be saved in Figure.UserData.Flag
% 'uparrow' and 'downarrow' - tune the parameters up or down
% 
% Parameters (not necessary to set since all of them have default values) 
% 'Active' - 1/0, control whether the simulation is going on or not
%            in 'Experiment.m' scripts
% 'Pause' - 1/0, whether the experiment is paused or not
% 'Mode' - 'stimulation' or 'clamp', how your mouse click will
%          affect the simulation system 
% 'StimulationRule' - The artificial stimulation waveforms (this 
%                     will create a new transient ExternalInput 
%                     instance to your NeuralNetwork)
% 'ZapRule' - This will produce a oscillatory square waves that persist
%             
% 'ClampRule' - This will clamp the dynamic parameter to to the value where
%               you click your mouse
% 'DeltaRule' - The perturbation waveforms for the dynamical variables
%               Effective when your 'Mode' is 'delta'
% 'ParameterRule' - The rule applied to change NeuralNetwork.param
%                   when you press 'uparrow' or 'downarrow'
% Jyun-you Liou, 2017/04/01

% For multiple inputs
if numel(O)>1
    for i = numel(O):-1:1
        F(i) = AttachHotKey(O(i), varargin);
    end
    return;
end

% Start of HotKey definitions
p = inputParser;
addParameter(p,'Active',true); % Whether the simulation is active or not
addParameter(p,'Pause',false); 
addParameter(p,'ParameterRule',[]);
addParameter(p,'Flag',false);
addParameter(p,'Mode','stimulation'); % What intervention you are giving to the network
                                      % Options: 'stimulation', 'zap', 'delta', 'clamp'
% Set Stimulation parameters
amp_stimulus = 1000; % amplitude of stimulation
t_stimulus = [0 10 20]; % Biphasic stimulation
sigma_stimulus = 10.^2; % physical range of stimulation
addParameter(p,'StimulationRule',@(x,t) amp_stimulus*((t_stimulus(2)>=t&t>t_stimulus(1))-(t_stimulus(3)>=t&t>t_stimulus(2))) ...
                                        * mvnpdf(x,0,sigma_stimulus*eye(2))/mvnpdf([0 0],0,sigma_stimulus*eye(2)));
% Set up Zap rules                                    
amp_zap = 1000; % amplitude of stimulation
t_cycle = 5; % Biphasic, so frequency = 1/t_cycle/2;
sigma_zap = 10.^2; % physical range of stimulation
addParameter(p,'ZapRule',@(x,t) amp_zap * ...
                                (-sign(mod(floor(t/t_cycle),2)-0.5)) * ...
                                mvnpdf(x,0,sigma_zap*eye(2))/mvnpdf([0 0],0,sigma_zap*eye(2)));
                                    
% Set up Delta rules
sigma_delta = 10.^2;
Amplify_factor = 0.5;
addParameter(p,'DeltaRule',@(x) Amplify_factor * mvnpdf(x,0,sigma_delta*eye(2))/mvnpdf([0 0],0,sigma_delta*eye(2)));

% Set up Clamp rules
sigma_clamp = 10.^2;
addParameter(p,'ClampRule',@(x) mvnpdf(x,0,sigma_clamp*eye(2))/mvnpdf([0 0],0,sigma_clamp*eye(2)));
parse(p,varargin{:});

% Find the .Graph
if isempty(O.Graph) || ~isvalid(O.Graph);plot(O);end
F = O.Graph;
F.UserData = p.Results; % Put the hotkey parameters into Figure.UserData

% Hot Key settings
F.KeyPressFcn = {@HotKeyCallback,O};
    function HotKeyCallback(src,evd,O)
        % src is the figure handle
        switch evd.Key
            case 'space' % Pause & Continue
                if src.UserData.Pause;uiresume(src);else;uiwait(src);end
            case 'q' % Quit
                src.UserData.Active = false;
            case 's' % Set the interaction mode into stimulation mode
                disp('Change to artificial stimulation mode. Click the axes to give stimulus.');
                src.UserData.Mode = 'stimulation';
            case 'z' % Set the interaction mode into zap mode 
                disp('Change to artificial stimulation mode. Click the axes to give stimulus.');
                src.UserData.Mode = 'zap';
            case 'd' % Set the interaction mode into delta (perturbation) mode
                disp('Change to dynamic variable delta mode.  Click the axes to perturb dynamic variables');
                src.UserData.Mode = 'delta';
            case 'c' % Set the interaction mode into clamp mode
                if O.dim > 1
                    disp('Clamping is now only available for 1D model.');
                else
                    src.UserData.Mode = 'clamp';
                end
            case 'n' % Nullify all interventions
                disp('Stop interfering the system.')
                src.UserData.Mode = [];
            case 'r' % Restart from the beginning
                if isempty(O.Recorder);warning('There is no Recorder can be found at all.');end
                Load(O.Recorder);disp('Initial data has been loaded.');
            case 'p' % Change parameters                       
                ParamName = fieldnames(O.param);
                Param = ParamName(structfun(@(x) isnumeric(x),O.param)); % You can only change numeric parameters
                for m = 1:numel(Param) % Magnitude of change
                    Param{m,2} = mean(O.param.(Param{m,1})(:))*0.05;
                end
                [Param{:,3}] = deal(false); % Whether it is going to change
                g = figure;
                T = uitable(g,'Data',Param, ...
                          'ColumnName',{'Parameter','Magnitude (Default 5%)','Enabled'}, ...
                          'ColumnEditable',[false true true], ...
                          'Position',[0.1*g.Position(3) 0.2*g.Position(4) 0.8*g.Position(3) 0.7*g.Position(4)]);
                uicontrol('String','OK', ...
                          'Callback',@(s,~) uiresume(src), ...
                          'Position',[0.1*g.Position(3) 0.05*g.Position(4) 0.8*g.Position(3) 0.1*g.Position(4)]);
                uiwait(src);
                src.UserData.ParameterChange = T.Data(cell2mat(T.Data(:,3)),1:2);
                close(g);
            case 'f'
                src.UserData.Flag = true;
            case 'uparrow' % Adjust parameters up
                C = src.UserData.ParameterChange;
                if isempty(C)
                    disp('Please set parameter change rule first.  Press HotKey p')
                else
                    for m = 1:size(C,1)
                        O.param.(C{m,1}) = O.param.(C{m,1}) + C{m,2};
                        disp(['Parameter ' C{m,1} ' has been increased. Now average = ' num2str(mean(O.param.(C{m,1})(:)))]);                                
                    end
                end
            case 'downarrow' % Adjust parameters up
                C = src.UserData.ParameterChange;
                if isempty(C)
                    disp('Please set parameter change rule first.  Press HotKey p')
                else
                    for m = 1:size(C,1)
                        O.param.(C{m,1}) = O.param.(C{m,1}) - C{m,2};
                        disp(['Parameter ' C{m,1} ' has been decreased. Now average = ' num2str(mean(O.param.(C{m,1})(:)))]);                                    
                    end
                end
        end
    end

% Mouse Callback functions (Give stimulation, delta, or clamp) at 'axes'
set(F.Children,'ButtonDownFcn',{@MouseClickCallback,O});                
    function MouseClickCallback(src,evd,O)
        P = evd.IntersectionPoint; % Intervention position (Px, Py, Pz)
        P = P(1:2);
        P2_original = P(2); % This will only be used for 'clamp' mode 
        P = round(P); % Make the stimulation/delta/clamp right on some cell/population
        if O.dim == 1;P(2)=1;end % Make sure there is no loss of simulation/delta/clamp power in 1D model                
        Fig = src.Parent; % The parent figure
        Polarity = any(strcmp(Fig.SelectionType,{'normal','open'})) - strcmp(Fig.SelectionType,'alt'); % Right or left mouse click
        % Anything that you want to control from mouse click
        switch lower(Fig.UserData.Mode) % The mode is saved in 'Figure'
            case 'stimulation' % This will add an ExternalInput object (handle)
                disp('Give stimulation as an external input.');
                % Generate an ExternalInput instance
                Iext = ExternalInput;
                Iext.Target = O;
                O.Ext = [O.Ext;Iext];
                % Construct the ExternalInput's waveform
                t_offset = O.t; % You can't directly use O.t because O is a handle!
                Iext.Deterministic = @(x,t) Polarity*(Fig.UserData.StimulationRule(bsxfun(@minus,x,P),t-t_offset));
                Iext.Tmax = t_offset+100;
            case 'zap'
                disp('Giving zap as external input.  Click again to stop it.');
                % Generate an ExternalInput instance
                Iext = ExternalInput;
                Iext.Target = O;
                O.Ext = [O.Ext;Iext];
                % Construct the ExternalInput's waveform
                t_offset = O.t; % You can't directly use O.t because O is a handle!
                Iext.Deterministic = @(x,t) Polarity*(Fig.UserData.ZapRule(bsxfun(@minus,x,P),t-t_offset));
                Iext.Tmax = inf;
                set(Fig.Children,'ButtonDownFcn',{@MouseClickCallbackFinishZap,Iext});    
            case 'clamp' % This will directly clamp the dynamical variables to the where you clicked
                disp(['Clamping parameters: ' src.Title.String]);
                [x2,x1] = meshgrid(1:O.n(2),1:O.n(1));                
                ClampFactor = Fig.UserData.ClampRule(bsxfun(@minus,[x1(:) x2(:)],P));
                ClampFactor = reshape(ClampFactor,O.n);
                O.(src.Title.String) = ClampFactor .* P2_original + (1 - ClampFactor) .* O.(src.Title.String);                
            case 'delta' % This will directly perturb the dynamical variables
                disp(['Perturbing parameters: ' src.Title.String]);
                [x2,x1] = meshgrid(1:O.n(2),1:O.n(1));
                PerturbFactor = Fig.UserData.DeltaRule(bsxfun(@minus,[x1(:) x2(:)],P));
                PerturbFactor = reshape(PerturbFactor,O.n);
                O.(src.Title.String) = (1 + Polarity * PerturbFactor) .* O.(src.Title.String);
        end
    end

    function MouseClickCallbackFinishZap(src,~,Iext)
        Fig = src.Parent; % This parent figure
        Iext.Tmax = -inf;
        disp('Finish zapping.');
        set(Fig.Children,'ButtonDownFcn',{@MouseClickCallback,Iext.Target});
    end
end