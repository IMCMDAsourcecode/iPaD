function [beta, mm, m, nin]= LassoSolver(x, y, lam, x_var, beta_old, mm, m, nin)
%% Solve the lasso problem using the coordinate descent algorithm

%%
[n, p] = size(x);
beta = beta_old;  % warm start with beta_old 

yhat = x*beta_old;%zeros(n,1);

tol = 1e-02;
iter_outer = 0;
iter_inner = 0;


%%
while 1
    %loops on all variable
    dlx = 0.0;
    for j = 1:p
        beta(j) = shrinkage(x(:,j)'*(y-yhat) + x_var(j)*beta_old(j), lam)/x_var(j);   
        diff_beta = beta_old(j) - beta(j);
        diffabs = abs(diff_beta);
        if diffabs > 0
            yhat = yhat - x(:,j)*diff_beta;  
            beta_old(j) = beta(j);            
            % record j if j is not in active set
            if (mm(j) == 0)
               nin = nin + 1;
               mm(j) = nin;
               m(nin) = j;
            end            
            dlx = max([dlx,diffabs]); % max diff in coeff
        end       
    end
    iter_outer = iter_outer+1;
    if dlx<tol
        break;
    end
    
    %loops on the active set
    while (1)
        dlx = 0.0;
        iter_inner = iter_inner+1;
        for k = 1:nin
            j = m(k);
            beta(j) = shrinkage(x(:,j)'*(y-yhat) + x_var(j)*beta_old(j), lam)/x_var(j);   
            diff_beta = beta_old(j) - beta(j);
            diffabs = abs(diff_beta);
            if (diffabs > 0)
                yhat = yhat - x(:,j)*diff_beta;  
                beta_old(j) = beta(j);
                
                dlx = max([dlx,diffabs]); % max diff in coeff
            end
            
        end
        if (dlx<tol)
            break; % converged on the active set
        end
    end
    
end

%fprintf('Lasso solved (outer loops: %d; inner loops: %d ).\n',iter_outer, iter_inner);
end

function y = shrinkage(a, kappa)
    y = max(0, a-kappa) - max(0, -a-kappa);
end
