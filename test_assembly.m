function [Res,Tang] = test_assembly(mat_size, iter )
%returns zero res. and Identity tangent.
% for testing the main script
switch iter 
    case 0 
        Res = ones(mat_size,1);
        Tang = diag(ones(mat_size,1));
    case 1 
        Res = zeros(mat_size,1);
        Tang = diag(ones(mat_size,1));
    otherwise 
        error('iteration number is >=2')
end
end

