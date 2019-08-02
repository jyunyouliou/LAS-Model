function Output = conv2_field(R,ker,nT,topology)
% 2D convultion with circular topology
% Input: R - presynaptic source layer activity, matrix
%        ker - kernel, needs to have size of [2m(1)+1,2m(1)+1]
%        Nt - size of target layer, if not given = size(R)
%        topology - linear / circular (deafult: linear)
%
% Output: R_output - Circular convolution of R and ker
%
% This function is called by NeuralField methods
% Author: Jyun-you Liou
% Final update: 2016/05/23 
if nargin < 4
    topology = 'linear';
end

if nargin < 3
    nT = size(R);
end

% check kernel size
if ~mod(size(ker),2)
    error('Kernel size needs to be odd')
end
n = size(R);

% Check if n - nt can be evenly divided, if not, zero-pad the input layer
odd = mod((n-nT),2);
if any(odd)
    R(end+odd(1),end+odd(2)) = 0;
    R = (R + circshift(R,odd))/2;
end
n = size(R);

% Convolution
O = conv2(R,ker);

% Check size of O matrix & change O into a correct size
nO = size(O);
nRound = ceil(nO ./ nT);
if ~all(nO == nT)
    O( nRound(1) * nT(1)  ,  nRound(2) * nT(2) ) = 0;
end

dshift = (nO-nT)/2;
O = circshift(O,-dshift);

if strcmpi(topology,'circular') 
    Output = zeros(nT);
    for iter1 = 1:nRound(1)
        for iter2 = 1:nRound(2)
            Output = Output + O((iter1-1)*nT(1) + (1:nT(1)),(iter2-1)*nT(2) + (1:nT(2)));
        end
    end
elseif strcmpi(topology,'linear')
    Output = O(1:nT(1),1:nT(2));    
end


end

% Test code
% R = eye(10);Ker = [0 1 0;1 0 1;0 1 0];nT = [10 10];
% Output = conv2_circular(R,Ker,nT);imagesc(Output)