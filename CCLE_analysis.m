%% Set the random number seed
seed=1;
rng(seed);

%% Load the CCLE data
Y1 = load('CCLE_Y1.txt');
Y2 = load('CCLE_Y2.txt');
L1 = load('CCLE_L1.txt') == 1;
% Note that no known drug-pathway associations will not be used as prior
% knowledge as the CCLE_L2_prior.txt are all zeros.
L2_prior = load('CCLE_L2_prior.txt') == 1;


%% Take the log-transformation
Y1 = log(Y1);
Y2(Y2 == 0) = min(Y2(Y2 ~= 0 & ~isnan(Y2)))*0.9;
Y2 = log(Y2);



%% Do some quality control of the data
[Y1, L1, Y2, L2_prior, gene_idx, drug_idx, pathway_idx] = quality_control(Y1, L1, Y2, L2_prior);



%% Set the options for running the program
opts.tol = 1e-4; % precision tolerance for convergence 
opts.init_tol = 1e-10; % precision tolerance for convergence while initializing X
opts.max_iter = 500; % the maximum number of iterations allowed
opts.verbose = 1; % Do you want to print a lot of information while running the program
opts.b2_max_nz_frac = 0.2; % The maximum fraction of zeros in B2 matrix that are allowed to be non-zeros while solving a sequence of lambdas.


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

%% Run the permutation test (2,000 permutations)
tic;
permu_results = iPaD_permu(Y1, Y2, L1, L2_prior, init_x, cv_lam, 2000, opts);
toc;

%% Save the results
save(strcat('CCLE_result_', num2str(seed), '.mat'), 'first_run_results', 'cv_results', 'permu_results');
