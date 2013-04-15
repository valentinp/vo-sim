% Least squares estimation for motion
function [R_hat, T_hat] = motion_est_initial(P, C, Sigma_P, Sigma_C)

weights = zeros(1, length(P(1,:)));
for i=1:length(P(1,:))
    weights(i) = 1/(det(Sigma_P(:,:,i)) + det(Sigma_C(:,:,i)));
end

w = sum(weights);

Qc = zeros(3,1);
Qp = zeros(3,1);
A = zeros(3,3);

for i=1:length(P(1,:))
    Qc = Qc + weights(i)*C(1:3,i)/w;
    Qp = Qp + weights(i)*P(1:3,i)/w;
end
for i=1:length(P(1,:))
    A = A + weights(i)*(C(1:3, i) - Qc)*(P(1:3,i) -Qp)'/w;
end
W = A/w;

[V,S,U] = svd(W);

R_hat = V*[1 0 0; 0 1 0; 0 0 det(U)*det(V)]*U';
T_hat = -R_hat'*Qc + Qp;
end


