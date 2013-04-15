function [R_hat, T_hat] = motion_est_matrix(P, C, Sigma_P, Sigma_C)
%MOTION_EST_MATRIX 
%Estimate the motion based on Matrix weights
%See AER1514 Lecture 15 for details

%P and C are 3xN matrices where every ith column represents the 3D location
%of a landmark in the previous and current frames

%Sigma_P and Sigma_C are 3x3xN 3D matrices that represent the co-variance 
%matrices of the ith landmark in the current and previous frames


mcross = @(u) [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0]; 


%No of landmarks
no_lm = length(P(1,:));

%Initialize other variables
error = 10;


[R, T] = motion_est_initial(P, C, Sigma_P, Sigma_C);


i = 1;
MAX_ITER = 200;
M = [T; 0; 0; 0];
lambda = 0;

while error > 10^-4 && i < MAX_ITER

    LS = zeros(6);
    RS = zeros(6,1);
    e = zeros(3, no_lm);
    E = zeros(3,6,no_lm);
    W = zeros(3,3,no_lm);
    
    
    for j=1:no_lm
        
        e(:,j) = C(1:3,j) - R*(P(1:3,j) - T);
        E(:,:,j) = [R -mcross(R*(P(1:3,j) - T))]; 
        w_inner = R*Sigma_P(:,:,j)*R' + Sigma_C(:,:,j);
        W(:,:,j) = inv(w_inner);
        
        %Build up Ej'*Wj*Ej and -Ej'*W'ej; the left and right sides
        %Use LM optimization with lambda = 5
        LS = LS + E(:,:,j)'*W(:,:,j)*E(:,:,j) + lambda*diag(diag(E(:,:,j)'*W(:,:,j)*E(:,:,j)));
        RS = RS + E(:,:,j)'*W(:,:,j)*e(:,j);   
    end

    
    %Calculate objective function (using old m)
    obj_fn = 0;
    for j=1:no_lm
        t = C(1:3,j) - R*(P(1:3,j) - T);
        obj_fn = obj_fn + 0.5*(t'*W(:,:,j)*t);
    end
    
   %Update M
    M = LS\(-RS);
    
    %Check to see if objective function has decreased and update lambda
    %Based on http://cronos.rutgers.edu/~meer/TEACH/ADD/spareLM.pdf
%     obj_fn_trial = 0;
%     phi = M(4:6);
%     eps = M(1:3);
%     R_trial = (eye(3) - mcross(phi))*R;
%     T_trial = T + eps;
%     for j=1:no_lm
%         t =C(1:3,j) - R_trial*(P(1:3,j) - T_trial);
%         obj_fn_trial = obj_fn_trial + 0.5*(t'*W(:,:,j)*t);
%     end
%     if obj_fn_trial < obj_fn
%        lambda = lambda/2;
%     else
%        lambda = lambda*2;
%     end
   % Perform line search 
    for alpha = 1:-0.01:0.01
            phi = alpha*M(4:6);
            eps = alpha*M(1:3);
            R_trial = (eye(3) - mcross(phi))*R;
            T_trial = T + eps;
        obj_fn_trial = 0;
        for j=1:no_lm
            t =C(1:3,j) - R_trial*(P(1:3,j) - T_trial);
            obj_fn_trial = obj_fn_trial + 0.5*(t'*W(:,:,j)*t);
        end
        if obj_fn_trial < obj_fn
            M = alpha*M;
            break;
        end
    end
    
    phi = M(4:6);
    eps = M(1:3);
    
    R = (eye(3) - mcross(phi))*R;
    T = T + eps;
    error = norm(M);
    i=i+1;
end
%disp(i);
%disp(alpha);
R_hat = R;
T_hat = T;
end

