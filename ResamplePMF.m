function [ dist, x ] = ResamplePMF( Dist, N )
% Resample from a discrete distribution
%
% Dist: the distribution (probability mass function)
% N: number of samples drawn from the distribution
%
% dist: the sample distribution
% x: sample category sequence (only computed if inquired)
%
% Jyun-you Liou, 2016/12/09

S = size(Dist);
Dist = Dist(:) / sum(Dist(:)); % Normalization
D_cumsum = cumsum([0;Dist]);
dist = histcounts(rand(N,1),D_cumsum);
if nargout == 2
    x = nan(N,1);
    idx_event = find(dist);
    pointer = 0;
    for i = idx_event(:)'
        x(pointer + (1:dist(i))) = i;
        pointer = pointer + dist(i);
    end
    x = x(randperm(N));
end
dist = dist / sum(dist(:));
dist = reshape(dist,S);
end

