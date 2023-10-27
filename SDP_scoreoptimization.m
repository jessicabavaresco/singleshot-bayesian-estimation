function [T,score] = SDP_scoreoptimization(Xi,d,minmax)
% set minmax =  1 for maximization
% set minmax = -1 for minimization

d_in  = d(1);
d_out = d(2);
No    = size(Xi,3);

yalmip('clear');

T = sdpvar(d_in*d_out,d_in*d_out,No,'hermitian','complex');

F = [trace(sum(T,3))==d_out,
    sum(T,3)==kron(PartialTrace(sum(T,3),2,[d_in d_out]),eye(d_out)/d_out)
    ];

score = 0;
for i=1:No
    F = F + [T(:,:,i)>=0];
    score = score + real(trace(T(:,:,i)*Xi(:,:,i)));
end


solution = solvesdp(F,-minmax*score,sdpsettings('solver','sedumi','verbose',0,'cachesolvers',1));


T = double(T);
score = double(score);

