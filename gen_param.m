function [ output ] = gen_param(mu,sigma,m,flag)
% Distribution generator
% [ output ] = gen_param(mu,sigma,m,flag)
% m: size of the tensor
% flag = 0: normal (default), flag = 1: gamma, flag = 2: lognormal
%
% This function is called by all neural field template script
% Author: Jyun-you Liou
% Final update: 2016/05/18
if nargin < 4
    flag = 0;
end

if sigma == 0
    output = mu * ones(m);
elseif sigma > 0
    switch flag 
        case {'normal',0}
            output = mu + sigma .* randn(m);
        case {'gamma',1}
            output = gamrnd((mu.^2)./(sigma.^2),(sigma.^2)./(mu),m);
        case {'lognormal',2}
            MU = @(MEAN,DELTA) log(MEAN) - log(1+DELTA.^2) / 2;
            STD = @(MEAN,DELTA) sqrt(log(1+DELTA.^2));
            output = exp(MU(mu,sigma./mu) + STD(mu,sigma./mu).*randn(m));
    end
else
    error('Standard deviation needs to be positive')
end


end

