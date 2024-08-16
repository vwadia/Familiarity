% 

clear

%parpool(32)

% dim_ind = [1:5 21:25];
dim_ind = [1:10 21:30]; % his fam SA features had only 10 shape and 10 app dims

%% monkey faces
load ./params_1k_120d.mat
params = bsxfun(@rdivide, params, std(params)); % normalize features

unfam_ind = 1:1000;
allfam_ind = 1001:1036;

%%

params = params(:, dim_ind);
n_stim = size(params,1);

x_fam = false(n_stim,1); x_fam(allfam_ind) = true;
x_unfam = false(n_stim,1); x_unfam(unfam_ind) = true;

%%

for k = [33 34 35] %[23 24 25 26]
    
    n_fam_min = k;

x = true(n_stim,1);
n_fam = length(allfam_ind);
n_unfam = length(unfam_ind);


%%

n = 2000;

grad = zeros(n_stim,1);

for i = 1:n
      
    % estimate grad
    parfor j=1:n_stim
        xd = x;
        xd(j) = ~xd(j);
        grad(j) = cost_fun( params, xd, x_fam, x_unfam ) - cost;
    end
    
    if n_fam<=n_fam_min
        grad(x&x_fam) = inf; % finds the indices where both x and x_fam = 1
    end
    [grad_min, ind] = min(grad);
    
    if grad_min<0
        x(ind) = ~x(ind);
        cost = cost + grad(ind);
        
        n_fam = sum(x&x_fam);
        n_unfam = sum(x&x_unfam);
        fprintf('[%d] cost %.5f, Nfam%d, Nunfam%d,\n', i, cost, n_fam, n_unfam);
    
    elseif n_unfam<=n_fam
        break;
        
    else
        break;
    end
    
end

%%
ind_unfam_sub = find(x&x_unfam);
ind_fam_sub = find(x&x_fam);
output_file_name = sprintf('match_by_grad_Nfam%d_Nunfam%d.mat', n_fam, n_unfam);
save(output_file_name, 'ind_unfam_sub', 'ind_fam_sub');
 
end