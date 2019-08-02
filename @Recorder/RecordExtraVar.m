function RecordExtraVar( R, VarName )
% RecordExtraVar(Recorder) is a method of handle class object 'Recorder'.
% In addition to the standard recording variables (dynamical variables), 
% you can freely add more variables to record by adding new fields in the 
% structure .Var. 
%
% Now the options include: 'EPSC' & 'IPSC'
% 
% How to use it?  First, you need to AddVar(R,'Parameter'), by adding
% the available extra variables you want to record.  Now only EPSC and IPSC
% are available.  Then, just keep Record(R), the fields will just function
% as a regular variables.  It will also compress and plot as a regular
% variable.
%
% Now only time series data can be added.  In the future, I may expand this
% to include point process data.
%
% Jyun-you Liou 2017/03/27


O = R.Ancestor;

switch VarName
    case 'EPSC'
        % Due to poor program flow design, unfortunately we need to
        % recalculate conductance
        % Step 1A: Get g (condutance)
        g = O.Input.E;
        % Step 1B: Add the recurrent projection
        for i = 1:numel(O.Proj.In) 
            if O.Proj.In(i).Type == 'E'
            g = g + O.Proj.In(i).Value ./ O.param.tau_syn.E;
            end
        end
        % Step 2: Calculate EPSCs
        R.Var.EPSC(:,R.Idx) = g .* (O.param.E_Esyn-O.V) ./ O.param.f_max;
    case 'IPSC'
        % Due to poor program flow design, unfortunately we need to
        % recalculate conductance
        % Step 1A: Get g (condutance)
        g = O.Input.I;
        % Step 1B: Add the recurrent projection
        for i = 1:numel(O.Proj.In) 
            if O.Proj.In(i).Type == 'I'                
            g = g + O.Proj.In(i).Value ./ O.param.tau_syn.I;
            end
        end
        % Step 2: Calculate PSCs
        E_Cl = 26.7*log(O.Cl_in./O.param.Cl_ex);                        
        R.Var.IPSC(:,R.Idx) = g .* (E_Cl-O.V) ./ O.param.f_max;
end

end

