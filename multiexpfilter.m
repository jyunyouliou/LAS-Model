function x = multiexpfilter(x,input,tau,dt)
    % This function does multi-stage exponential filtering
    % x: n1 by n2 by nstage tensor
    % input: n1 by n2 matrix, unit: kHz
    % tau: n1 by n2 by n2 positive tensor, if any tau = 0, it will produce instant
    %      response without any filter effect
    % dt: dt
    %
    % This function is called by Class NeuralField - method - update
    % Author: Jyun-you Liou
    % Final update: 2016/05/18
    if ~all(tau(:))
        x = repmat(input,[1,1,size(x,3)]);
    else
        for m = 1:size(x,3)                 
            if m == 1
                x(:,:,m) = x(:,:,m).*exp(-dt./tau(:,:,m)) + input.*dt./tau(:,:,m);
            else
                x(:,:,m) = x(:,:,m).*exp(-dt./tau(:,:,m)) + x(:,:,m-1).*dt./tau(:,:,m);
            end           
        end                  
    end    
end