function X = ChasePeak(T,X,varargin)
% Chase where traveling waves meet 
p = inputParser;
addParameter(p,'dTmax',10); % maximal temporal lag allowed
addParameter(p,'dXmax',1); % search size
parse(p,varargin{:});
dXmax = p.Results.dXmax;
dTmax = p.Results.dTmax;

dT = nan(size(T)+2*dXmax);
dT(dXmax+1:end-dXmax,dXmax+1:end-dXmax) = T;

X = X + dXmax;

while true
    dTlocal = dT(X(1)+(-dXmax:dXmax), X(2)+(-dXmax:dXmax)) - dT(X(1),X(2));
    dTlocal(dTlocal>dTmax) = nan;
    [~,local_idx] = max(dTlocal(:));
    [I,J] = ind2sub(size(dTlocal),local_idx);
    I = I - dXmax - 1;
    J = J - dXmax - 1;
    if I==0 && J==0;break;end
    X = X + [I,J];
end

end