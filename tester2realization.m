function [rho, M] = tester2realization(T,d)
% d = [din dout]
% T_i = T(:,:,i) in I O

din  = d(1);
dout = d(2);
No   = size(T,3);
M = zeros(din*dout,din*dout,No);


sigma = PartialTrace(sum(T,3),2,[din dout])/dout;

Phi = din*MaximallyEntangled(din);

% rho in I AUX
rho = kron(eye(din),sqrtm(sigma))*Phi*kron(eye(din),sqrtm(sigma))';

% M_i = M(:,:,i) in AUX O
for i=1:No
    M(:,:,i) = kron(inv(sqrtm(sigma)),eye(dout))*T(:,:,i)*kron(inv(sqrtm(sigma)),eye(dout))';
end