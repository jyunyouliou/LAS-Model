% Figure 1 is a 2D mean field model, showing the qualitative simularity 
% between the model and real data (You need to zoom in to show a subset of 
% the recordings so that the two looks alike)

% Lay out the field
O = MeanFieldModel('Exp1Template');

% Build the round mask (round boundary condition)
[xx,yy] = meshgrid( -(O.n(2)/2)+0.5:(O.n(2)/2)-0.5, -(O.n(1)/2)+0.5:(O.n(1)/2)-0.5);
xx = xx/O.n(2); % Normalized spatial unit
yy = yy/O.n(1); % Normalized spatial unit
rr = sqrt(xx.^2 + yy.^2); % Normalized spatial unit
mask = rr < 0.5; % The easiest way to apply mask is to redefine its activation function    
f_original = O.param.f;
O.param.f = @(v) mask.*f_original(v); % Now neurons outside the boundary can not fire

%% Build recurrent projections
[ P_E, P_I1, P_I2 ] = StandardRecurrentConnection( O );

%% External input
Ic = 200;
stim_t = [2 5]; % second
stim_x = [0.5 0.1]; % normalized spatial unit
stim_r = 0.05; % normalized spatial unit
O.Ext = ExternalInput;
O.Ext.Target = O;
O.Ext.Deterministic = @(x,t) ( (sqrt(((x(:,1)-1)/(O.n(2)-1)-stim_x(1)).^2 + ((x(:,2)-1)/(O.n(1)-1)-stim_x(2)).^2)) < stim_r) .* ...
                             ( (stim_t(2)*1000) > t & t > (stim_t(1)*1000)) .* ...
                             Ic; % x: position (neuron index), t: ms, unit: pA 

%% Simulation settings 
dt = 1; % ms
write_cycle = 10; % Every 10 ms save the result once to avoid asking for too much RAM 
R = CreateRecorder(O,10000); % The 2nd argument is Recorder.Capacity 
T_end = write_cycle*R.Capacity - 1*write_cycle; % simulation end time.                           

%% Realtime plot setting
flag_realtime_plot = 1; % whether you want to see simulation result real time or not
T_plot_cycle = 1000; % How often updating the figures
if flag_realtime_plot
    AttachHotKey(O); % So that you can press 'q' to stop the simulation any time
    f = plot(O);drawnow;
    Ax = findobj(f,'Type','Axes');
    xlim(Ax,[0,size(O.V,2)]+0.5);
    ylim(Ax,[0,size(O.V,1)]+0.5);
    caxis(Ax(1),[-80 -20]); % V
    caxis(Ax(2),[-80 -20]); % phi
    caxis(Ax(3),[0 50]); % Chloride concentration
    caxis(Ax(4),[0 mean(O.param.g_K_max(:).*O.param.f_max(:))]); % g_K      
    for j = 4:-1:1
        cbar(j) = colorbar(Ax(j),'EastOutside');
    end
    ylabel(cbar(4),['X ' num2str(1/mean(O.param.f_max(:)))]);
end

%% Simulation
while 1 % You can press 'q' to escape this while loop     
    % End mechanism
    if ~O.Graph.UserData.Active || O.t >= T_end
        break;
    end
    
    % Record & update 
    if mod(O.t,write_cycle) < dt
        WriteToRecorder(O); 
    end
    Update(O,dt);

    % Real time plotting
    % use mod < dt instead of == 0 can deal with floating number errors
    if mod(O.t,T_plot_cycle) < dt && flag_realtime_plot 
        f=plot(O);drawnow
    end    
end
