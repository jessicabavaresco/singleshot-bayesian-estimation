function [T,score] = script_thermometry(method,PPT)
% set method = 1, 2, or 3, to choose desired method
% set PPT = 0 for standard problem and PPT = 1 for no-entanglement
% approximation

d = 2;

Tmin = 0.1;     % minimum temperature
Tmax = 2;       % maximum temperature

Nt = 100;       % number of time steps

No_initial = 2;
No_final   = 2;

time = 0:(1/(Nt-1)):1;

score = zeros(No_final-No_initial+1,Nt);
T     = cell(No_final-No_initial+1,Nt);

for No=No_initial:No_final
   
    No = No
    
    %%%%%%%%% M1 %%%%%%%%%
    if method==1
        Nh = No;  
        
    %%%%%%%%% M2 and M3 %%%%%%%%%   
    else
        Nh = 1000; 
    end
   
    theta_i(1:No) = Tmin:(Tmax-Tmin)/(No-1):Tmax;
    theta_k(1:Nh) = Tmin:(Tmax-Tmin)/(Nh-1):Tmax;
    
     p(1:Nh,1) = 1/Nh;      % uniform prior
     
     r = zeros(No,Nh);
     for i=1:No
         for k=1:Nh
             r(i,k) = (theta_i(i)-theta_k(k))^2;    % cost function
         end
     end
        
    for t = 1:Nt
        
        Ck = zeros(d^2,d^2,Nh);
        for k=1:Nh
            Ck(:,:,k) = ChoiOperatorThermo(theta_k(k),0.1,2,time(t));
        end
        
        Xi = zeros(d^2,d^2,No);
        for i=1:No
            for k=1:Nh
                Xi(:,:,i) = Xi(:,:,i) + p(k)*r(i,k)*Ck(:,:,k);
            end
        end
        
        [T{No-No_initial+1,t},score(No-No_initial+1,t),~] = score_optimization(Xi,[d d],-1,PPT);
        
        
        %%%%%%%%% M3 (seesaw) %%%%%%%%%
        if method==3
            
            T_temp = T{No-1,t};
            
            gap = 1;
            precision = 10^(-6);
            rounds(t,1) = 0;
            old_score = score(No-1,t);
            flag_value = 0;
            
            while gap>precision
                
                rounds(t,1) = rounds(t,1) + 1;
                
                %%%%% step 1 %%%%%
                
                [estimators] = estimator_optimization(p,T_temp,Ck,theta_k);
                
                theta_i = estimators;
                
                Xi = zeros(d*d,d*d,No);
                r = zeros(No,Nh);
                for i=1:No
                    for k=1:Nh
                        r(i,k) = (theta_k(k)-theta_i(i))^2;
                        Xi(:,:,i) = Xi(:,:,i) + p(k)*r(i,k)*Ck(:,:,k);
                    end
                end
                
                %%%%% step 2 %%%%%
                
                [T_temp,score_temp,flag] = score_optimization(Xi,[d d],-1,PPT);
                
                if flag.problem~=0
                    sdp_problem = flag.problem
                    info        = flag.info
                    flag_value = 1;
                    %pause
                end
                
                gap         = abs(old_score - score_temp);
                old_score   = score_temp;
                %pause
            end
            
            T{No-No_initial+1,t} = T_temp;
            score(No-No_initial+1,t) = score_temp;
            
        end

    end   
    
end

end

function [T,score,solution] = score_optimization(Xi,d,minmax,PPT)
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
    if PPT==1
        F = F + [PartialTranspose(T(:,:,i),2,[d_in d_out])>=0];
    end
    score = score + real(trace(T(:,:,i)*Xi(:,:,i)));
end


solution = solvesdp(F,-minmax*score,sdpsettings('solver','sedumi','verbose',0,'cachesolvers',1));


T = double(T);
score = double(score);

end

function [estimators] = estimator_optimization(p,T,Ck,theta_k)
%%%% estimator optimization for thermometry %%%%

Nh = max(size(p));
No = size(T,3);

post = zeros(Nh,No);
estimators = zeros(No,1);

for k = 1:Nh
    for i = 1:No
        post(k,i) = p(k)*real(trace(T(:,:,i)*Ck(:,:,k)));
    end
end
post = post./sum(post,1);

for i = 1:No
    for k = 1:Nh
        estimators(i,1) = estimators(i,1) + post(k,i)*theta_k(k);
    end
end

end

function J = ChoiOperatorThermo(T,eps,g,t)

% explicit form of Choi state
N = 1/(exp(eps/T) - 1);
gamma = g*(2*N+1);
a = (N+N*exp(-gamma*t)+1)/(2*(2*N+1));
b = ((1-exp(-gamma*t))*(N+1))/(2*(2*N+1));
c = (N-N*exp(-gamma*t))/(2*(2*N+1));
d = (exp(-gamma*t)*(N+1)+N)/(2*(2*N+1));

diag_terms = [a b c d];
off_diag_term = exp(-(gamma*t)/2)/2;

J = diag(diag_terms);
J(1,4) = off_diag_term;
J(4,1) = off_diag_term';

% change choi states to choi matrices + swap parties for consistency with the notes
J = 2 * J;
J = Swap(J);

end

