function P = AdjustWeight( P, varargin )
% P = AdjustWeight( P, varargin )
%
% In the case with space boundary, if you want to keep your projection
% toward the bounder having the same maximal value as other spot, you can
% use AdjustWeight to adjust the result.
%
% This will change .WPost.  So now if all of your input strength = 1, you
% will have output all = 1.  Please further adjust .WPost to achieve the
% strength you want to have.
%
% Final update: Jyun-you Liou, 2017/04/30
p = inputParser;
parse(p,varargin{:});

Value_safe_copy = P.Value;
Os = P.Source;
P.Value = ones(Os.n);
Project(P);
w = 1./P.Value;
P.WPost = P.WPost .* w;

P.Value = Value_safe_copy;

end

