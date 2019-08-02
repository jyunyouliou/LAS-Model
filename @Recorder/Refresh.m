function Refresh(R)
    % Refresh(R)
    %
    % Recorder Class specific function 
    %
    % Empty the current Recorder, all dynamic variable record fields, spike
    % recording (S) will be erased.  Writing index will be set back to 1.
    %
    % Jyun-you Liou, 2017/03/22
    if numel(R)>1
        for i =1:numel(R)
            Refresh(R(i));
        end
        return;
    end
    R.T = nan*R.T;
    for Name = R.VarName
        R.Var.(Name{:}) = nan*R.Var.(Name{:});
    end
    R.SBuffer = cell(1,R.Capacity);
    R.S = [];
    R.Idx = 1;
end
