function [ cost ] = cost_fun( p, x, x1, x2 )
%COST_FUN written by Liang
% INPUTS:
    % 1. p --> params (for all ims)
    % 2. x --> vector of 1s with the current stim set to 0
    % 3. x1 --> Vector of zeros length(n_stim) with all the fam stim as 1s
    % 4. x2 --> same as x1 but with unfam stim marked as 1s
% 
% OUTPUTS:
    % 1. cost 
    % made up of cost1, cost2, cost3
    % cost1 - mse of the per dim variance bw fam and unfam images
    % cost2 - difference of mean variance bw fam and unfam images
    % cost3 - set to either 1/(min_p*0.001) or 0 depending on whether the min p-val is less than 0.05 or not
    % p-val: vector of pvals consisting of 
    %     single p-val for whether distributions of pairwise distances between fam and unfam params are different
    %     n dim vector (n = nuumber of dimensions) of pvals testing if dsitribution of vals in each dimension for fam and unfam are different
% Description written by vwadia Nov2024

p1 = p(x&x1,:); % fam params
p2 = p(x&x2,:); % all but 1 unfam params

% variance - per dim
var1 = var(p1);
var2 = var(p2);

cost1 = mean((var1-var2).^2); % mse
cost2 = abs(mean(var1) - mean(var2)); % total variance

% pairwise distance distribution
d1 = pdist(p1);
d2 = pdist(p2);
[~, p_pd] = kstest2(d1,d2); % p-value of whether the d1 and d2 come from the same distribution
% cost3 = eval_distribution_difference( d1, d2, 20);


% distribution each dimension
n = size(p1,2); % fam images
ps = zeros(n,1);
for i = 1:n
    [~, ps(i)] = kstest2(p1(:,i),p2(:,i));
end

% compute cost3
p_all = [p_pd; ps];
p_min = min(p_all);

if p_min<0.05
    cost3 = 1/(p_min+0.001);
else
    cost3 = 0;
end

% total cost
cost = cost1 + cost2 + cost3;

end

