function [dN_dxi] = get_shape_fnc_derv(quad_iterator,n_quad,...
    element_node_number, quad_rule)

if ~exist('quad_rule','var')
    % default quad rule is gauss
    error("quad rule does not exist")
end

[xi, nu] = get_xi_nu(quad_iterator,n_quad, quad_rule);

switch element_node_number
    case 9
    %derivative of shape function wrt xi
        dN_dxi(1,1)=1/4*(2*xi-1)*(nu-1)*nu;
        dN_dxi(2,1)=1/4*(2*xi+1)*(nu-1)*nu;
        dN_dxi(3,1)=1/4*(2*xi+1)*(nu+1)*nu;
        dN_dxi(4,1)=1/4*(2*xi-1)*(nu+1)*nu;
        dN_dxi(5,1)=1/2*(-2*xi)*(nu-1)*nu;
        dN_dxi(6,1)=1/2*(2*xi+1)*(1-nu.^2);
        dN_dxi(7,1)=1/2*(-2*xi)*(nu+1)*nu;
        dN_dxi(8,1)=1/2*(2*xi-1)*(1-nu.^2);
        dN_dxi(9,1)=(1-nu.^2)*(-2*xi);

    %derivative of shape function wrt eta
        dN_dxi(1,2)=1/4*(2*nu-1)*(xi-1)*xi;
        dN_dxi(2,2)=1/4*(2*nu-1)*(xi+1)*xi;
        dN_dxi(3,2)=1/4*(2*nu+1)*(xi+1)*xi;
        dN_dxi(4,2)=1/4*(2*nu+1)*(xi-1)*xi;
        dN_dxi(5,2)=1/2 *(2*nu-1)*(1-xi.^2);
        dN_dxi(6,2)=1/2 *(-2*nu)*(xi+1)*xi;
        dN_dxi(7,2)=1/2 *(2*nu+1)*(1-xi.^2);
        dN_dxi(8,2)=1/2 *(-2*nu)*(xi-1)*xi;
        dN_dxi(9,2)=(1-xi.^2)*(-2*nu);

    case 8
%         dN_dxi=zeros(element_node_number,2); 
%         N_xi= zeros(element_node_number,1);
%         N_nu= zeros(element_node_number,1);

        dN_dxi(5,1)= 1/2 * (-2*xi) * (1-nu);
        dN_dxi(6,1)= 1/2 * (1-nu.^2) * (1);
        dN_dxi(7,1)= 1/2 * (-2*xi) * (1+nu);
        dN_dxi(8,1)= 1/2 * (1-nu.^2) * (-1);
        dN_dxi(1,1)= 1/4 * (-1) * (1-nu) - 1/2 * (dN_dxi(5,1)+dN_dxi(8,1));
        dN_dxi(2,1)= 1/4 * (1) * (1-nu)  - 1/2 * (dN_dxi(5,1)+dN_dxi(6,1));
        dN_dxi(3,1)= 1/4 * (1) * (1+nu)  - 1/2 * (dN_dxi(6,1)+dN_dxi(7,1));
        dN_dxi(4,1)= 1/4 * (-1) * (1+nu) - 1/2 * (dN_dxi(7,1)+dN_dxi(8,1));
        
        dN_dxi(5,1)= 1/2 * (1-xi.^2) * (-1);
        dN_dxi(6,1)= 1/2 * (-2*nu) * (1+xi);
        dN_dxi(7,1)= 1/2 * (1-xi.^2) * (1);
        dN_dxi(8,1)= 1/2 * (-2*nu) * (1-xi);
        dN_dxi(1,1)= 1/4 * (1-xi) * (-1) - 1/2 * (dN_dxi(5,1)+dN_dxi(8,1));
        dN_dxi(2,1)= 1/4 * (1+xi) * (-1) - 1/2 * (dN_dxi(5,1)+dN_dxi(6,1));
        dN_dxi(3,1)= 1/4 * (1+xi) * (1)  - 1/2 * (dN_dxi(6,1)+dN_dxi(7,1));
        dN_dxi(4,1)= 1/4 * (1-xi) * (1)  - 1/2 * (dN_dxi(7,1)+dN_dxi(8,1));
        
        %   dN_dxi=[N_xi  N_nu];
        
    case  4
        %derivative of shape function wrt xi
        dN_dxi(1,1)= 1/4 * (-1) * (1-nu);
        dN_dxi(2,1)= 1/4 * (1) * (1-nu);
        dN_dxi(3,1)= 1/4 * (1) * (1+nu);
        dN_dxi(4,1)= 1/4 * (-1) * (1+nu);

        %derivative of shape function wrt eta
        dN_dxi(1,2)= 1/4 * (1-xi) * (-1);
        dN_dxi(2,2)= 1/4 * (1+xi) * (-1);
        dN_dxi(3,2)= 1/4 * (1+xi) * (1);
        dN_dxi(4,2)= 1/4 * (1-xi) * (1);
         
    otherwise
        error('shape function derivatives for this node numbers is not implemented')
end
        

end

