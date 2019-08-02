function [W,STAT] = AnalyzeConnectivityByDistance(P)
% [W,D,STAT] = AnaylzeConnectivityByDistance(P)
%
% Analyze the connectivity matrix by distance.  P is an instance of class 
% 'Projection'.  P.method should be 'multiplication'.  
%
% Outputs: W - The connectivity matrix, circularly shifted so that the
%              first row is auto-projection
%          STAT - a structure showing some basic statistics according to 
%                 distance such as 'mean' and 'variance'.
% 
% Jyun-you Liou, final update: 2017/06/25
W = P.W;
[~,n] = size(W);
for i = 1:n
    W(:,i) = circshift( W(:,i), [1-i, 0] );
end
STAT.mean = mean(W,2);
STAT.var = var(W,0,2);
end