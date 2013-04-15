function [J] = jacob3D(p,b,f)
%Jacobian of the 2D->3D "reprojection" funciton
    u_l = p(1);
    v_l = p(2);
    u_r = p(3);
    v_r = p(4);
    J = [-b*u_r/(u_l - u_r)^2, 0, b*u_l/(u_l - u_r)^2, 0;
        -b/2*(v_l + v_r)/(u_l - u_r)^2, b/2*1/(u_l - u_r), (b/2)*(v_l+v_r)/(u_l - u_r)^2, (b/2)*1/(u_l-u_r); 
        -b*f/(u_l-u_r)^2, 0,  b*f/(u_l - u_r)^2, 0;
    ];
end

