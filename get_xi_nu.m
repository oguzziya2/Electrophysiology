function [xi, nu] = get_xi_nu(quad_iterator,n_quad,quad_rule)

if ~exist('quad_rule','var')
    % default quad rule is gauss
    error("quad rule does not exist")
end

if (quad_rule=="lumped") 
    %quadrature points are the nodes.
    switch n_quad
        case 4
            points= [-1 -1; 1 -1; 1 1 ; -1 1 ] ;
            point = points(quad_iterator, :);
            xi= point (1); nu= point(2);
            
        case 9
            points= [-1 -1; 0 -1; 1 -1; ...
                     -1  0; 0  0; 1  0; ...
                     -1  1; 0  1; 1  1] ;
            point = points(quad_iterator, :);
            xi= point (1); nu= point(2);
        otherwise
            error('this quadrature rule is not implemented')
    end
    
elseif (quad_rule=="higher-order")
    %quadrature points are between that of gauss and nodes.
    % details are in p446 of Hughes FEM book 
    switch n_quad
        case 4
            points= (1+1/sqrt(3))/2*[-1 -1; 1 -1; 1 1 ; -1 1 ] ;
            point = points(quad_iterator, :);
            xi= point (1); nu= point(2);
            
        case 9
            error('this quadrature rule is not implemented')
            
        otherwise
            error('this quadrature rule is not implemented')
    end
elseif (quad_rule=="gauss")
    %quadrature points are the gauss points.
    switch n_quad
        case 4
            points= 1/sqrt(3)*[-1 -1; 1 -1; 1 1 ; -1 1 ] ;
            point = points(quad_iterator, :);
            xi= point (1); nu= point(2);
            
        case 9
            points= sqrt(3/5)*[-1 -1; 0 -1; 1 -1; ...
                -1  0; 0  0; 1  0; ...
                -1  1; 0  1; 1  1] ;
            point = points(quad_iterator, :);
            xi= point (1); nu= point(2);
            
        otherwise
            error('this quadrature rule is not implemented')
    end
else
    error('this quadrature rule is not implemented')
end

end

