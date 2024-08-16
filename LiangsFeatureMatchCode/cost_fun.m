function [ cost ] = cost_fun( p, x, x1, x2 )
%COST_FUN 

p1 = p(x&x1,:);
p2 = p(x&x2,:);

% variance
var1 = var(p1);
var2 = var(p2);

cost1 = mean((var1-var2).^2); % mse
cost2 = abs(mean(var1) - mean(var2)); % total variance

% pairwise distance distribution
d1 = pdist(p1);
d2 = pdist(p2);
[~, p_pd] = kstest2(d1,d2);
% cost3 = eval_distribution_difference( d1, d2, 20);


% distribution each dimension
n = size(p1,2);
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

