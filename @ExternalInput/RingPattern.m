function [E, Func] = RingPattern
% [E, Func] = RingPattern
%
% A GUI-based method to create a ring pattern function.  It is used to
% create a ring like stimulation pattern.  Output E is an instance of 
% ExternalInput with its deterministic part Func.  Func function can also
% be used to clamp some dynamic parameters.
%
% Reminder: when you incorporate E into a NeuralNetwork, remember to set 
% its E.Tmax to NeuralNetwork.t + E.Tmax.
%
% Jyun-you Liou, 2017/04/01
disp('Select three points to determine a circle.')
[X,Y] = ginput(3);
set(gca,'NextPlot','add');
Pt = plot(X,Y,'.');hold on;
Pt.MarkerSize = 20;
% It is actually a simple math exercise you can do to think why the circle
% can be determine this way.  If you can not figure it out, go to see this
% http://www.ambrsoft.com/TrigoCalc/Circle3D.htm
M =   [X(1).^2+Y(1).^2, X(1), Y(1), 1; ...
       X(2).^2+Y(2).^2, X(2), Y(2), 1; ... 
       X(3).^2+Y(3).^2, X(3), Y(3), 1];
A = det(M(:,2:4));
B = -det(M(:,[1 3 4]));
C = det(M(:,[1 2 4]));
D = -det(M(:,1:3));
x0 = -B/A/2; % Center of the circle
y0 = -C/A/2;
r = sqrt((B^2 + C^2 - 4*A*D)/ (4*A^2)); % Radius of the circle

% Plot the stimulation ring for the user.
theta = linspace(0,2*pi,360);
L = plot(x0 + r*cos(theta),y0+r*sin(theta),'g--');
L.LineWidth = 4;
% Stimulation details
Q = {'Stimulation width (unit: # of neurons)', ...
     'Stimulation strength (pA)', ...
     'Stimulation duration (ms)', ...
     'Stimulation sinusoidal frequency (Hz) - if empty or 0, it will be just a square wave'};
A = {'2', ...
     '300', ...
     '20', ...
     ''};
A=inputdlg(Q,'Stimulation setting',1,A);
A = str2double(A);
dr = A(1);
strength = A(2);
Tmax = A(3);
if isnan(A(4))
    freq = 0;
    delta_theta = pi/2;
else
    freq = A(4);
    delta_theta = 0;
end
% Construct ExternalInput instance
Func = @(x,t) (((r+dr/2) > sqrt((x(:,1)-y0).^2 + (x(:,2)-x0).^2)) & ...
               (sqrt((x(:,1)-y0).^2 + (x(:,2)-x0).^2) > (r-dr/2))) * ...
               strength * sin(2*pi*freq*t + delta_theta);
E = ExternalInput;
E.Deterministic = Func;
E.Tmax = Tmax;
delete(Pt);
delete(L);
end