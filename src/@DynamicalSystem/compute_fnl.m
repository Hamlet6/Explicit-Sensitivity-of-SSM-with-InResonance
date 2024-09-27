function fnl = compute_fnl(obj,x,xd)
% COMPUTE_FNL We compute the nonlinear internal force in a second-order
% mechanical system. Currently, we do not treat velocity dependent
% nonlinearity.

assert(obj.order == 2, ' fnl can only be computed for second-order systems')

fnl = zeros(obj.n,1);
switch obj.Options.Intrusion
    case 'full'
        for j = 1:length(obj.fnl)
            if size(obj.fnl(j).ind,2) == obj.N % check if the nonlinearity is velocity dependent as well
                z = [x;xd];
            else
                z = x;
            end
            fnl = fnl + expand_multiindex(obj.fnl(j),z);
        end

    case 'semi'
        error('implementation for semi-intrusive computation is not available')

    case 'none'
        fnl = obj.fnl_non(x);

    otherwise
        error('options for compute fnl are none, semi or full')

end