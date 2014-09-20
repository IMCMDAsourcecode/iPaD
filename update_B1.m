function [B1] = update_B1(Y1, X, L1)

%% Update the pathway-gene relationship coefficient matrix

%% Input values:
% Y1:   The transcription data matrix
% X:    The pathway activity level matrix
% L1:   The indicator matrix for pathway-gene relationships

%% Output values:

% B1:   The updated pathway-gene relationship coefficient matrix


%%
[p, g] = size(L1);
B1 = zeros(p, g);

miss_idx = isnan(Y1);
for i = 1:g
    B1(L1(:,i),i) = pinv(X(~miss_idx(:,i), L1(:,i))) * Y1(~miss_idx(:,i), i);
end

end
