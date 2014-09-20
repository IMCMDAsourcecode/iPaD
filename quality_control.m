function [Y1, L1, Y2, L2_prior, gene_idx, drug_idx, pathway_idx] = quality_control(Y1, L1, Y2, L2_prior)
%% Do the following two quality control steps:
% 1. Remove genes or drugs that have too few unique values because in such
% cases these genes and drugs are not informative
% 2. Merge pathways that have identical member genes

%% 

gene_idx = 1:size(Y1, 2);
drug_idx = 1:size(Y2, 2);
pathway_idx = 1:size(L1, 1);

rm_idx = [];
for i = 1:size(Y1, 2)
    if length(unique(Y1(:,i))) <= 2 
        rm_idx = [rm_idx i];
    end
end
Y1(:, rm_idx) = [];
L1(:, rm_idx) = [];
gene_idx(rm_idx) = [];


rm_idx = [];
for i = 1:size(Y2, 2)
    if length(unique(Y2(:,i))) <= 2
        rm_idx = [rm_idx i];
    end
end
Y2(:, rm_idx) = [];
L2_prior(:, rm_idx) = [];
drug_idx(rm_idx) = [];



%% 

[unique_L1, ia, ~] = unique(L1, 'rows', 'stable');
L1 = unique_L1;
pathway_idx = pathway_idx(ia);

end
