function [X, funcVal] = update_X_miss(Y1, Y2, B1, B2, X, opts)

% min ||Y-XB||^2_F
% st. ||X^2_j||^2 <= 1

Y = [Y1 Y2];
B = [B1 B2];
miss_idx = isnan(Y);

[n, q] = size(Y);
p = size(B,1);

mu = norm(B)^2;
BB = B*B';

funcVal_out = [];


for iter_out = 1:opts.max_iter
    Y_tilde = X*B;
    Y(miss_idx) = Y_tilde(miss_idx);
    funcVal = [];
    XX = X;
    YB = Y*B';
    t=1;
    for iter = 1:opts.max_iter
        X_old = X;
        
        % compute gradient
        gradient = X*BB - YB;
        
        XX = XX - gradient/mu;
        
        X = XX;
        
        %% projection
        for j = 1:p
            normXj = norm(X(:,j),'fro');
            if normXj>1
                X(:,j) = X(:,j)/normXj;
            end
        end
        
        r = Y-X*B;
        
        funcVal_1 = sum(r(:).^2);
        
        funcVal = cat(1,funcVal,funcVal_1);
        % check convergence
        if iter > 1
            if abs(funcVal(iter)-funcVal(iter-1))/abs(funcVal(iter-1)) <= opts.tol
                break
            elseif funcVal(iter) > funcVal(iter-1)
                t = 1; % restart
            end
        end
        
        % for next iteration
        t_old = t;
        t = (1+sqrt(1+4*t^2))/2;%nesterov method
        XX = X + (t_old-1)/t*(X-X_old);
    end
    r = Y - X*B;
    funcVal_out_1 = sum(r(~miss_idx).^2);
    funcVal_out = cat(1,funcVal_out,funcVal_out_1);
    if iter_out > 1
        if abs((funcVal_out(iter_out) - funcVal_out(iter_out-1))/funcVal_out(iter_out)) <= opts.tol
            break;
        end
    end
    
    if sum(sum(miss_idx)) == 0
        break;
    end
end
funcVal = funcVal_out(end);
end
