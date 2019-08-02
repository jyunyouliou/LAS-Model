function [E, Func] = FocalPattern
% [E, Func] = FocalPattern
%
% A GUI-based method to create a focal pattern function.  It is used to
% create a ring like stimulation pattern.  Output E is an instance of 
% ExternalInput with its deterministic part Func.  Func function can also
% be used to clamp some dynamic parameters.
%
% Reminder: when you incorporate E into a NeuralNetwork, remember to set 
% its E.Tmax to NeuralNetwork.t + E.Tmax.
%
% Jyun-you Liou, 2017/04/01
disp('Please determine stimulation center.');
[x0,y0] = ginput(1);
set(gca,'NextPlot','add');
Pt = scatter(x0,y0,'filled');Pt.SizeData = 50;Pt.CData = [0 1 0];
Q = {'Stimulation radius (unit: # of neurons)', ...
     'Stimulation strength (pA)', ...
     'Stimulation duration (ms)', ...
     'Stimulation sinusoidal frequency (kHz) - if empty or 0, it will be just a square wave', ...
     'Shape (sharp or Gaussian)'};
A = {'5', ...
     '300', ...
     '20', ...
     '', ...
     'sharp'};
A = inputdlg(Q,'Stimulation setting',1,A);
Type = A{5};
A = str2double(A);
r = A(1);
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
switch Type
    case 'Gaussian'
        Func = @(x,t) normpdf(sqrt((x(:,1)-y0).^2 + (x(:,2)-x0).^2),0,r)/normpdf(0) * ...
                      strength * sin(2*pi*freq*t + delta_theta);
    case 'sharp' 
        Func = @(x,t) (sqrt((x(:,1)-y0).^2 + (x(:,2)-x0).^2) < r) * ...
                       strength * sin(2*pi*freq*t + delta_theta);
end
E = ExternalInput;
E.Deterministic = Func;
E.Tmax = Tmax;
delete(Pt);
end