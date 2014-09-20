function [X] = initialize_X(Y1, L1, tol)

%% Initialize pathway activity level matrix

%% Input values

% Y1:   The transcription data matrix
% L1:   The indicator matrix for gene-pathway relationships
% tol:  The precision tolerance for convergence

%% Output values

% X:    The initialized X matrix

%% Normalizing Data

g = size(Y1, 2);
for i = 1:g
    Y1(:,i) = (Y1(:,i)-nanmean(Y1(:,i)))./nanstd(Y1(:,i));
end


%%
max_iter = 10000;

Y = Y1;
B = zeros(size(L1));
B(L1) = 1;

[n, q] = size(Y);
p = size(B,1);
BB = B*B';

X = randn(n, p);
X = X./repmat(sqrt(sum(X.^2)), n, 1);
XX = X;

mu = norm(B)^2;

YB = Y*B';
t=1;

funcVal = [];

for iter = 1:max_iter
    X_old = X;
    
    % compute gradient
    gradient = X*BB - YB;
   
    XX = XX - gradient/mu;
    
    X = XX;
    
    % projection
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
        if abs(funcVal(iter)-funcVal(iter-1))/abs(funcVal(iter-1)) <= 1e-6
			if max(abs((X_old(:) - X(:))./X(:))) < tol
            	break
			end
        elseif funcVal(iter) > funcVal(iter-1)
            t = 1; % restart
        end
    end
    
    % for next iteration
    t_old = t;
    t = (1+sqrt(1+4*t^2))/2;%nesterov method
    XX = X + (t_old-1)/t*(X-X_old);
end

end
