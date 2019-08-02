%% Spiking model - onset LFP pattern depends on inhibitory projection distribution
% This is a spiking model that qualitatively replicate the finding of mean
% field model.  While changing local inhibition kernel ratio, this can
% produce forward traveling waves.
%
% Final update: Jyun-you Liou, 2017/04/29
if ~exist('PanelFlag','var')
    PanelFlag = 'A';
end

%% Lay out the field
O = SpikingModel('Exp4Template');

%% Build recurrent projections
[ P_E, P_I1, P_I2 ] = StandardRecurrentConnection( O );

if PanelFlag == 'B'
    P_I1.WPost = 150; % Local recurrent inhibition strength    
    P_I2.WPost = 150; % Non-local recurrent inhibition strength     
end

%% External input
Ic = 80;
stim_x = [0.05 0.1];
stim_t = [2 5];
O.Ext = ExternalInput;
O.Ext.Target = O;
O.Ext.Deterministic = @(x,t) ((stim_x(2)*O.n(1))>x(:,1) & x(:,1)>(stim_x(1)*O.n(1))) .* ...
                              ((stim_t(2)*1000)> t & t > (stim_t(1)*1000)) .* ...
                              Ic; % x: position, t: ms, unit: mV

%% Simulation settings 
dt = 1; % ms
R = CreateRecorder(O,50000); % The 2nd argument is Recorder.Capacity 
T_end = R.Capacity - 1; % simulation end time.  
AddVar(R,'EPSC'); % Record PSCs so that you can calculate LFP proxy later
AddVar(R,'IPSC');

%% Realtime plot setting
flag_realtime_plot = 1; % whether you want to see simulation result real time or not
T_plot_cycle = 1000; % How often updating the figures
if flag_realtime_plot 
    AttachHotKey(O);
    f = plot(O);drawnow;        
    ylim(f.Children(1),[-80 0]); % V
    ylim(f.Children(2),[-80 0]); % phi
    ylim(f.Children(3),[0 40]); % Cl_in
    ylim(f.Children(4),[0 mean(O.param.g_K_max .*O.param.f_max)]); % g_K
    ylabel(f.Children(4),['X ' num2str(1/mean(O.param.f_max(:)))]);
    set(f.Children,'YLimMode','manual','XLim',[0 max(O.n)]);        
end


%% Simulation
while 1  % You need to press 'q' to escape this while loop
    % Termination criteria 
    if ~O.Graph.UserData.Active || O.t > T_end
        break;
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