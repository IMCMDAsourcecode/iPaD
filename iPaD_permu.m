function [permu_results] = iPaD_permu(Y1, Y2, L1, L2_prior, init_x, lam, npermu, opts)

%% Perform permutation test to calculate the p-values for the drug-pathway associations

%% Input values:

% Y1:       The transcription profile matrix

% Y2:       The drug sensitivity profile matrix

% L1:       The indicator matrix for gene-pathway relationship matrix

% L2_prior:      The indicator matrix for drug-pathway associations. Note that this
% matrix will be used as prior knowledge

% init_x:   The initial X matrix

% lam:      The lambda value used in the permutation test

% npermu:   The number of permutations

% opts:     A set of option parameters


%% Output values:

% permu_results.b2:     The B2 matrix estimated based on the original data

% permu_results.pval_exact: The p-values calcuated using the empirical null
% distribution

% permu_results.pval_approx: The p-values calculated using the approximated
% null distribution

% permu_results.b2_permu:   The B2 matrices calculated from the permuted
% data 


%% Options
if nargin < 8
    opts = [];
    opts.tol = 1e-4;
    opts.max_iter = 500;
end

if ~isfield(opts, 'tol')
    opts.tol = 1e-4;
end
if ~isfield(opts, 'max_iter')
    opts.max_iter = 500;
end

opts.verbose = 0;
%%
[n, g] = size(Y1);
[~, d] = size(Y2);
[p, ~] = size(L1);

%% Run with original data
results = iPaD(Y1, Y2, L1, L2_prior, init_x, lam, opts);
b2 = results.b2seq{1};

b2_permu = single(zeros(npermu, p, d));

fprintf('Number of permutations done:  0');
nback = 1;
for ipermu = 1:npermu
    Y2_permu = Y2(randsample(1:n, n),:);
	Y2_permu = Y2_permu - repmat(nanmean(Y2_permu), n, 1);
    results = iPaD(Y1, Y2_permu, L1, L2_prior, init_x, lam, opts);
    b2_permu(ipermu,:,:) = results.b2seq{1};
	if floor(log10(ipermu)+1) ~= nback
		fprintf(' ');
		nback = floor(log10(ipermu)+1);
	end
	fprintf(repmat('\b', 1, nback));
	fprintf('%d', ipermu);
end
fprintf('\n');

pval_exact = zeros(p, d);
pval_approx = zeros(p, d);
for i = 1:p
    for j = 1:d
        pval_exact(i, j) = mean(abs(b2_permu(:,i,j)) >= abs(b2(i,j)));
		nz_frac =  mean(b2_permu(:,i,j)~=0);
		mu = mean(b2_permu(b2_permu(:,i,j)~=0,i,j));
		sigma = std(b2_permu(b2_permu(:,i,j)~=0,i,j));
		pval_approx(i, j) = nz_frac*(normcdf(-abs(b2(i,j)), mu, sigma) + normcdf(-abs(b2(i,j)), -mu, sigma)) + (1-nz_frac)*(b2(i,j)==0);
		if sum(b2_permu(:,i,j)~=0) < 10
			pval_approx(i, j) = pval_exact(i, j);
		end
    end
end

permu_results.b2 = b2;
permu_results.pval_exact = pval_exact;
permu_results.pval_approx = pval_approx;
permu_results.b2_permu = b2_permu;


end
