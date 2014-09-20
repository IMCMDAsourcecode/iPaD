function [iPaD_results] = iPaD(Y1, Y2, L1, L2_prior, init_x, lamseq, opts)

%% The main function for iPaD

%% Input values:

% Y1:   The transcription data matrix

% Y2:   The drug sensitivity data matrix

% L1:   The indicator matrix for pathway-gene relationships

% L2_prior:   The indicator matrix for drug-pathway associations. Note that this
% matrix will be used as prior knowledge

% init_x:   The initial X matrix. iPaD will initialize X automatically if empty

% lamseq:   A sequence lambdas. iPaD will choose a sequence of lambdas
% automatically if empty

% opts:     A set of option parameters.


%% Output values:

% iPaD_results.lamseq: A sequence of lambda values used in this run

% iPaD_results.b2seq:  A sequence of B2 matrices corresponding to the 
% lambdas 

% iPaD_results.b1seq:  A sequence of B1 matrices corresponding to the
% lambdas

% iPaD_results.xseq:   A sequence of X matrices corresponding to the
% lambdas

% iPaD_results.L2_order: The ranking of the entries in B2 matrix by which
% they become non-zeros
%%

%% Options
if nargin < 5
    opts = [];
    opts.tol = 1e-4;
    opts.init_tol = 1e-6;
    opts.max_iter = 1000;
end

if ~isfield(opts, 'tol')
    opts.tol = 1e-4;
end
if ~isfield(opts, 'init_tol')
    opts.init_tol = 1e-6;
end
if ~isfield(opts, 'max_iter')
    opts.max_iter = 1000;
end
if ~isfield(opts, 'verbose')
    opts.verbose = 0;
end
if ~isfield(opts, 'b2_max_nz_frac')
    opts.b2_max_nz_frac = 0.7;
end
%%


%%
L1 = L1 == 1;
L2_prior = L2_prior == 1;
[n, g] = size(Y1);
[~, d] = size(Y2);
[p, ~] = size(L1);


mm = zeros(p,d); 
m = zeros(p,d); 
nin = zeros(d); 

ordered_idx = [];

%% Normalizing Data
for i = 1:g
    Y1(:,i) = (Y1(:,i)-nanmean(Y1(:,i)))./nanstd(Y1(:,i));
end

for i = 1:d
    Y2(:,i) = (Y2(:,i)-nanmean(Y2(:,i)))./nanstd(Y2(:,i));
end

%% Initialization
lams = [];
b2s = {};
b1s = {};
xs = {};

b2 = zeros(p, d);
if isempty(init_x)
    init_x = initialize_X(Y1, L1, opts.init_tol);
end

iPaD_results.init_x = init_x;

%% Alternating optimization algorithm

if isempty(lamseq)
    %% If a sequence lambdas were not pre-specified
    lambda = max(max(abs(init_x'*Y2)))*2;
    for i = 1:opts.max_iter
        x = init_x;
        obj = Inf;
        for j = 1:opts.max_iter
            [b2, mm, m, nin] = update_B2(Y2, x, lambda, b2, mm, m, nin, L2_prior, opts);
            b1 = update_B1(Y1, x, L1);
            [x, funcVal] = update_X(Y1, Y2, b1, b2, x, opts);
            new_obj = funcVal + lambda*sum(b2(L2_prior).^2) + lambda*sum(abs(b2(~L2_prior)));
            if (abs(new_obj-obj)/new_obj < opts.tol)
                break;
            end
            obj = new_obj;
        end
        
        if sum(sum(b2(~L2_prior)~=0))~=0
            lams = [lams lambda];
            b2s = [b2s b2];
            b1s = [b1s b1];
            xs = [xs x];
        end
        
        to_add = setdiff(find(b2~=0), [ordered_idx; find(L2_prior)]);
        if (~isempty(to_add))
            ordered_idx = [ordered_idx; to_add];
        end
        if (mean(b2(~L2_prior) ~= 0) > opts.b2_max_nz_frac)
            break;
        end
        if opts.verbose == 1
            fprintf('lambda = %f\t # non-zero coefficients in b2 = %d\n', lambda, sum(b2(~L2_prior) ~= 0));
        end
        if sum(b2(~L2_prior) ~= 0) == 0
            lambda = lambda*0.9;
        else
            lambda = lambda*0.9;
        end
    end
else
    %% If the a sequence of lambdas were pre-specified
    for i = 1:length(lamseq)
        obj = Inf;
        x = init_x;
        lambda = lamseq(i);
        for j = 1:opts.max_iter
            [b2, mm, m, nin] = update_B2(Y2, x, lambda, b2, mm, m, nin, L2_prior, opts);
            b1 = update_B1(Y1, x, L1);
            [x, funcVal] = update_X(Y1, Y2, b1, b2, x, opts);
            new_obj = funcVal + lambda*(sum(b2(L2_prior).^2)) + lambda*sum(abs(b2(~L2_prior)));
            if (abs(new_obj-obj)/new_obj < opts.tol)
                break;
            end
            obj = new_obj;
        end
        lams = [lams lambda];
        b2s = [b2s b2];
        b1s = [b1s b1];
        xs = [xs x];
        
        to_add = setdiff(find(b2~=0), [ordered_idx; find(L2_prior)]);
        if (~isempty(to_add))
            ordered_idx = [ordered_idx; to_add];
        end
        if opts.verbose == 1
            fprintf('lambda = %f\t # non-zero coefficients in b2 = %d\n', lambda, sum(b2(~L2_prior) ~= 0));
        end
    end
end



L2_order = Inf(p, d);
L2_order(ordered_idx) = 1:length(ordered_idx);

%% Output the results
iPaD_results.lamseq = lams;
iPaD_results.b2seq = b2s;
iPaD_results.b1seq = b1s;
iPaD_results.xseq = xs;
iPaD_results.L2_order = L2_order;

end
