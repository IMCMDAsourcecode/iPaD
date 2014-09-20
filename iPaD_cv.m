function [cv_results] = iPaD_cv(Y1, Y2, L1, L2_prior, init_x, lamseq, nfold, opts)

%% Perform cross-validation to select the penalty parameter - lambda

%% Input values:

% Y1:   The transcription data matrix

% Y2:   The drug sensitivity data matrix

% L1:   The indicator matrix for pathway-gene relationships

% L2_prior:   The indicator matrix for drug-pathway associations. Note that this
% matrix will be used as prior knowledge

% init_x:   The initial X matrix. 

% lamseq:   A sequence lambdas. 

% nfold:    The number of folds for cross-validation. The default is
% 10-fold.

% opts:     A set of option parameters.

%% Output values:

% cv_results.lamseq:    A sequence of lambdas used in the cross-validation
% cv_results.cv_err:    The residual some of squares (RSS) for each lambda
% cv_results.cv_lam:    The lambda value with the smallest RSS.

%% Options
if nargin < 8
    opts = [];
    opts.tol = 1e-4;
    opts.max_iter = 500;
    opts.nfold = 10;
end

if ~isfield(opts, 'tol')
    opts.tol = 1e-4;
end
if ~isfield(opts, 'max_iter')
    opts.max_iter = 500;
end
if ~isfield(opts, 'nfold')
    opts.nfold = 10;
end
%%

[n, g] = size(Y1);
[~, d] = size(Y2);
[p, ~] = size(L1);

%% Normalizing Data
for i = 1:g
    Y1(:,i) = (Y1(:,i)-nanmean(Y1(:,i)))./nanstd(Y1(:,i));
end

for i = 1:d
    Y2(:,i) = (Y2(:,i)-nanmean(Y2(:,i)))./nanstd(Y2(:,i));
end


%% Split data into folds
fold_size = floor(sum(sum(~isnan(Y2)))/nfold);
fold_ind = zeros(n, d);
fold = 1;
count = 0;
for i = 1:n
    for j = 1:d
        if isnan(Y2(mod(j+i-2, n)+1, j))
            continue;
        end
        fold_ind(mod(j+i-2, n)+1, j) = fold;
        count = count + 1;
        if count == fold_size
            count = 0;
            fold = fold + 1;
        end
    end
end

%% Cross-Validation
cv_err = Inf(nfold, length(lamseq));
opts.verbose = 0;
for ifold = 1:nfold
	fprintf('Fold %d...\n', ifold);
    Y2_poked = Y2;
    Y2_poked(fold_ind==ifold) = NaN;
    fold_result = iPaD(Y1, Y2_poked, L1, L2_prior, init_x, lamseq, opts);
    for i = 1:length(lamseq)
        Y2_tilde = fold_result.xseq{i}*fold_result.b2seq{i};
        cv_err(ifold, i) = nansum((Y2_tilde(fold_ind==ifold) - Y2(fold_ind==ifold)).^2);
    end
end
cv_err = sum(cv_err);
cv_lam = lamseq(find(cv_err == min(cv_err), 1));
%%

cv_results.cv_err = cv_err;
cv_results.lamseq = lamseq;
cv_results.cv_lam = cv_lam;

end
