function [T,score] = script_phaseestimation(method)
% set method = 1, 2, or 3, to choose desired method

d = 3;

No_initial = 2;
No_final   = 10;

Sz = eye(d);

for i=1:d
    Sz(i,i) = (i-1);
end

score = zeros(No_final-No_initial+1,1);
T     = cell(No_final-No_initial+1,1);

for No=No_initial:No_final
    
    No = No
    
    %%%%%%%%% M1 %%%%%%%%%
    if method==1
        Nh = No;  
        
    %%%%%%%%% M2 and M3 %%%%%%%%%   
    else
        Nh = 1000; 
    end
    
    estimators = 0:2*pi/(No):2*pi;      % uniform distribution of the estimator 
    discretization = 0:2*pi/(Nh):2*pi;  % uniform discretization of theta 
    
    theta_i = estimators(1:No);         % will use only first No values, not including 2pi
    theta_k = discretization(1:Nh);     % will use only first Nh values, not including 2pi

    p(1:Nh,1) = 1/Nh;       % uniform prior
    
%       p = normpdf(theta_k,pi,1);  % gaussian prior
%       p = p./sum(p);         % gaussian prior normalized
    
    Ck = zeros(d^2,d^2,Nh); 
    for k=1:Nh
        Ck(:,:,k) = kraus2choi(expm(-1i*theta_k(k)*Sz));    % channels
    end
    
    Xi = zeros(d^2,d^2,No);
    r = zeros(No,Nh);
    for i=1:No
        for k=1:Nh
            r(i,k) = cos((theta_k(k)-theta_i(i))/2)^2;          % reward function
            Xi(:,:,i) = Xi(:,:,i) + p(k)*r(i,k)*Ck(:,:,k);
        end
    end

    [T{No-No_initial+1,1},score(No-No_initial+1,1),~] = score_optimization(Xi,[d d],1);
    
    
    %%%%%%%%% M3 (seesaw) %%%%%%%%% 
    if method==3    
        
        T_temp = T{No-1,1};
        
        gap = 1;
        precision = 10^(-6);
        rounds = 0;
        old_score = score(No-1,1);
        flag_value = 0;
        
        while gap>precision
            
            rounds = rounds + 1;
                       
            %%%%% step 1 %%%%%
            
            [estimators] = estimator_optimization(p,T_temp,Ck,theta_k);
            
            theta_i = estimators;
            
            Xi = zeros(d^2,d^2,No);
            r  = zeros(No,Nh);
            for i=1:No
                for k=1:Nh
                    r(i,k) = cos((theta_k(k)-theta_i(i))/2)^2;
                    Xi(:,:,i) = Xi(:,:,i) + p(k)*r(i,k)*Ck(:,:,k);
                end
            end
            
            %%%%% step 2 %%%%%
            
            [T_temp,score_temp,flag] = score_optimization(Xi,[d d],1);
            if flag.problem~=0
                sdp_problem = flag.problem
                info        = flag.info
                flag_value = 1;
                %pause
            end
            
            gap       = abs(old_score - score_temp);
            old_score = score_temp;
            %pause
        end
        
       T{No-No_initial+1,1} = T_temp;
       score(No-No_initial+1,1) = score_temp;            
    end
       
end

end

function [T,score,solution] = score_optimization(Xi,d,minmax)
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

end
 
function [estimators] = estimator_optimization(p,T,Ck,theta_k)
%%%% estimator optimization for phase estimation %%%%

Nh = max(size(p));
No = size(T,3);

post = zeros(Nh,No);
avg_sin = zeros(No,1);
avg_cos = zeros(No,1);
estimators = zeros(No,1);

for k = 1:Nh
    for i = 1:No
        post(k,i) = p(k)*real(trace(T(:,:,i)*Ck(:,:,k)));
    end
end
post = post./sum(post,1);

for i = 1:No
    for k = 1:Nh
        avg_sin(i,1) = avg_sin(i,1) + post(k,i)*sin(theta_k(k));
        avg_cos(i,1) = avg_cos(i,1) + post(k,i)*cos(theta_k(k));
    end
    if avg_cos(i,1)>=0 
       estimators(i,1) = atan(avg_sin(i,1)/avg_cos(i,1));
    else 
       estimators(i,1) = atan(avg_sin(i,1)/avg_cos(i,1)) + pi;
    end
    if estimators(i,1)<0
        estimators(i,1) = estimators(i,1) + 2*pi;
    end
end


end 



