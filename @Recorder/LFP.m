function L = LFP( R, varargin )
% L = LFP(Recorder, 'parameter', parameter)
%
% This is a method function of Recorder.  L is the calculated local field 
% potential by Mazzoni's method (PLoS Computational Biology 2015)
%
% Recorder needs to have EPSC and IPSC recorded already
%
% 'Parameters' 
% 
% alpha: dominance of IPSC for LFP generation, default: 1.65
% dT: the temporal delay of EPSC signals, default is 6 (suppose dt = 1 ms)
%     If dt is not 1 ms, please modify it accordingly.
% sigma: default: 0 = no contribution of neighboring cells
%
% Jyun-you Liou, 2017/06/25

p = inputParser;
addParameter(p,'alpha',1.65);
addParameter(p,'dT',6);
addParameter(p,'sigma',0.025);
parse(p,varargin{:});


if ~all(ismember({'EPSC','IPSC'},R.VarName))
    error('This Recorder does not have EPSC and IPSC recorded.  Use AddVar(Record,''EPSC'',''IPSC'') to record them and start experiment again.');
end


L = R.Var.EPSC;
L = circshift(L,p.Results.dT);
L = L - p.Results.alpha * R.Var.IPSC;
L(:,1:p.Results.dT) = nan;

O = R.Ancestor;
if O.dim > 1
    warning('I have only written the convolution part of LFP simulation for 1-dimensional model.');
    return;
end

N = max(O.n);
kernel = exp(0:-1:-N / N / p.Results.sigma);
kernel = [kernel(end:-1:2),kernel];
kernel = kernel(:);
NAN = isnan(L);
L(NAN) = 0;
L = conv2(L,kernel,'same');
end

