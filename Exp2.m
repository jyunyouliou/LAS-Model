% Figure 2 is a one dimensional model showing all phasic transitions
% 
% 'PanelFlag' is a inherited variable from Fig2 to determine what 
% parameters to change & how to end the experiments
%
% If 'PanelFlag' is not given, it will just run as Panel A 
%
% Final update: Jyun-you Liou, 2017/04/22
if ~exist('PanelFlag','var')
    PanelFlag = 'A';
end

% Lay out the field
O = MeanFieldModel('Exp2Template');

%% Build recurrent projections
[ P_E, P_I1, P_I2 ] = StandardRecurrentConnection( O );

%% External input
Ic = 200;
stim_t = [2 5]; % second
stim_x = [0.1 0.15]; % normalized spatial unit
O.Ext = ExternalInput;
O.Ext.Target = O;
O.Ext.Deterministic = @(x,t) ( (O.n(1)*stim_x(2))>x(:,1) & x(:,1)>(O.n(1)*stim_x(1)) ) .* ...
                             ( (stim_t(2)*1000) > t & t > (stim_t(1)*1000)) .* ...
                             Ic; % x: position (neuron index), t: ms, unit: pA 

%% Simulation settings 
dt = 1; % ms
R = CreateRecorder(O,100000); % The 2nd argument is Recorder.Capacity 
T_end = R.Capacity - 1; % simulation end time.  

%% Realtime plot setting
flag_realtime_plot = 1; % whether you want to see simulation result real time or not
T_plot_cycle = 1000; % How often updating the figures
if flag_realtime_plot
    AttachHotKey(O);
    f = plot(O);drawnow;        
    ylim(f.Children(1),[-80 -20]); % V
    ylim(f.Children(2),[-80 -20]); % phi
    ylim(f.Children(3),[0 40]); % Chloride concentration
    ylim(f.Children(4),[0 mean(O.param.g_K_max .*O.param.f_max)]); % g_K
    ylabel(f.Children(4),['X ' num2str(1/mean(O.param.f_max(:)))]);
    set(f.Children,'YLimMode','manual','XLim',[0 max(O.n)]);         
end

%% Additional adjustment based on 'PanelFlag'
% Adjust parameters according to PanelFlag (for early termination of the experiment)
T_no_activity = 0; % Use to automatically turn off the experiment
rate_baseline = mean(O.param.f(O.V-O.phi));
switch PanelFlag
    case 'A'
        T_min = 100000; % Minimum simulation time 
    case 'D'
        T_min = 10000; % Minimum simulation time      
        O.param.tau_Cl = ones(O.n)*Condition(n_trial); % Panel D: vary chloride clearance speed
    case 'E'
        T_min = 10000; % Minimum simulation time 
        O.param.g_K_max = ones(O.n)*Condition(n_trial); % Panel E: vary potassium conductance strength
end

%% Simulation
while 1 % You need to press 'q' to escape this while loop     
    % Termination criteria 
    if ~O.Graph.UserData.Active || O.t >= T_end
        break;
    elseif O.t > T_min % Turn off the experiment if there is no activity
        if all(O.param.f(O.V-O.phi)/O.param.f_max(1)<0.01)
            T_no_activity = T_no_activity + 1;
        else
            T_no_activity = 0;
        end
        if T_no_activity > 5000
            break
        end
    end
    
    % Record & update 
    WriteToRecorder(O); 
    Update(O,dt);
    
    % Real time plotting
    % use mod < dt instead of == 0 can deal with floating number errors
    if mod(O.t,T_plot_cycle) < dt && flag_realtime_plot 
        f=plot(O);drawnow
    end   
end
