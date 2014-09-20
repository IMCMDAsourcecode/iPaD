%% Set the random number seed
seed=1;
rng(seed);

%% Load the pooled NCI-60 data
Y1 = load('NCI60_Y1.txt');
Y2 = load('NCI60_Y2.txt');
L1 = load('NCI60_L1.txt') == 1;
% Note that this is the drug-pathway associations that are used as prior knowledge. Set the matrix to all zeros if no prior knowledge are used.
L2_prior = load('NCI60_L2.txt') == 1; 



%% Do some quality control of the data
[Y1, L1, Y2, L2_prior, gene_idx, drug_idx, pathway_idx] = quality_control(Y1, L1, Y2, L2_prior);



%% Set the options for running the program
opts.tol = 1e-4; % precision tolerance for convergence 
opts.init_tol = 1e-14; % precision tolerance for convergence while initializing X
opts.max_iter = 500; % the maximum number of iterations allowed
opts.verbose = 1; % Do you want to print a lot of information while running the program
opts.b2_max_nz_frac = 0.5; % The maximum fraction of zeros in B2 matrix that are allowed to be non-zeros while solving a sequence of lambdas.

%% Initialize X
init_x = initialize_X(Y1, L1, opts.init_tol);

%% Run iPaD for a deceasing sequence of lambdas
tic;
first_run_results = iPaD(Y1, Y2, L1, L2_prior, init_x, [], opts);
toc;

lamseq = [first_run_results.lamseq(1:5)/(0.9^5) first_run_results.lamseq];

%% Ten-fold cross-validation
tic;
cv_results = iPaD_cv(Y1, Y2, L1, L2_prior, init_x, lamseq, 10, opts);
toc;

cv_lam = cv_results.cv_lam;

%% Run the permutation test (10,000 permutations)
tic;
permu_results = iPaD_permu(Y1, Y2, L1, L2_prior, init_x, cv_lam, 10, opts);
toc;

%% Save the results
save(strcat('NCI60_result_', num2str(seed), '.mat'), 'first_run_results', 'cv_results', 'permu_results');
