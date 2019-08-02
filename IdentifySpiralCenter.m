function [R,P] = IdentifySpiralCenter(theta,varargin)
% [R,P] = IdentifySpiralCenter(theta,varargin)
%
% Use phase singularity method to identify where spiral wave center is.
%
% Input: theta: a numeric matrix, with each entry as a phase (-pi to pi) of
%               at that point 
%        varargin:
%        'Range1' - where to search along the first dimension
%        'Range2' - where to search along the second dimension
%        'mask' - a logical matrix with size as theta, where to search.
%
% Output: R: a logical matrix with the same size as theta, indicating where
%            a path integral is 2pi or -2pi
%         P: an n_center by 2 numetric matrix.  Each row indicates where
%            ths center is 
% 
% Reference: An Experimentalist's Approach to Accurate Localization of Phase Singularities during Reentry
%            PMID: 11219507
%
% Final update: Jyun-you Liou, 2017/04/19

p = inputParser;
addParameter(p,'Range1',2:size(theta,1)-1);
addParameter(p,'Range2',2:size(theta,2)-1);
addParameter(p,'mask',true(size(theta)));
parse(p,varargin{:});

% Construct the 3 by 3 box that allows circular integration
seq = [1 2 3;
       8 inf 4; 
       7 6 5];
[~, seq] = sort(seq(:));
seq = seq(1:8);

% Search range
M1 = sort(p.Results.Range1);
M2 = sort(p.Results.Range2);

% Remove border
M1 = setdiff(M1,[1 size(theta,1)]);
M2 = setdiff(M2,[1 size(theta,2)]);

R = zeros(size(theta));

for m1 = M1
    for m2 = M2
        theta_local = theta(m1+(-1:1),m2+(-1:1));
        theta_local = theta_local(seq);
        theta_local(end+1) = theta_local(1);
        dtheta_local = diff(theta_local);
        dtheta_local(dtheta_local > pi) = dtheta_local(dtheta_local > pi)- 2*pi;
        dtheta_local(dtheta_local < -pi) = dtheta_local(dtheta_local < -pi)+ 2*pi;
        R(m1,m2) = sum(dtheta_local);
    end
end
R(~p.Results.mask) = 0;
[I1,I2] = find(abs(R)>pi);

n_center = numel(I1)/4;

if n_center < 1
    disp('Can not find any center');
    P = [];return;
end

for i = 1:n_center
    P1 = min(I1);
    P2 = min(I2(I1 == P1));
    P(i,:) = [P1 + 0.5, P2 + 0.5];
    sel = (ismember(I1,[P1 P1+1]) & ismember(I2,[P2 P2+1]));
    I1(sel) = [];
    I2(sel) = [];
end

end