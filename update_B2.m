function [B2, mm, m, nin] = update_B2(Y2, X, lambda, B2, mm, m, nin, L2_prior, opts)

%% Update the drug-pathway association coefficient matrix

%% Input values

% Y2:       The drug sensitivity profile matrix

% X:        The pathway activity level matrix

% lambda:   The penalty parameter

% B2:       The initial value for the drug-pathway association coefficient
% matrix

% mm, m, nin:   Several helping variables for the lasso algorithm
% mm:   record whether the variable in the active set
% m:    record active variables
% nin:  initialize the number of active variables

% L2_prior:   The indicator matrix for drug-pathway associations. Note that this
% matrix will be used as prior knowledge

% opts:     A set of option parameters

%% Output values

% B2:       The updated drug-pathway association coefficient matrix
% mm, m, nin:   See 'Input values' for details

%%

[~, d] = size(Y2);
miss_idx = isnan(Y2);

for i = 1:d
    alpha = B2(~L2_prior(:,i), i);
    alpha0 = B2(L2_prior(:,i), i);
    X_p = X(~miss_idx(:, i), ~L2_prior(:,i));
    X_np = X(~miss_idx(:, i), L2_prior(:,i));
    XtX = X_np'*X_np;
    x_var = sum(X_p.^2, 1);
    
    Y = Y2(~miss_idx(:, i), i);
    Y_tilde = Y - X_np*alpha0;
    for j = 1:opts.max_iter
        [new_alpha, mm(:, i), m(:, i), nin(i)] = LassoSolver(X_p, Y_tilde, lambda, x_var, alpha, mm(:, i), m(:, i), nin(i));
        new_alpha0 = (XtX+lambda*eye(size(XtX, 1)))\(X_np'*(Y-X_p*new_alpha));
        
        diff_alpha0 = abs((new_alpha0 - alpha0)./new_alpha0);
        if(isempty(diff_alpha0))
            break;
        end
        if(max(diff_alpha0) < opts.tol)
            break;
        end
        Y_tilde = Y - X_np*new_alpha0;
        alpha = new_alpha;
        alpha0 = new_alpha0;
    end

    B2(~L2_prior(:,i), i) = new_alpha;
    B2(L2_prior(:,i), i) = new_alpha0;
end

end
