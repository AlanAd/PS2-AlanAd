close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
cons = k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 

ret = cons .^ (1 - sigma) / (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret(cons < 0) = -Inf;

%%%% Iteration
dis = 1; tol = 1e-06; % tolerance for stopping 
v_guess = zeros(1, num_k);
while dis > tol
    % compute the utility value for all possible combinations of k and k':
    value_mat = ret + beta * repmat(v_guess, [num_k 1]);
    
    % find the optimal k' for every k:
    [vfn, pol_indx] = max(value_mat, [], 2);
    vfn = vfn';
    
    % what is the distance between current guess and value function
    dis = max(abs(vfn - v_guess));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess = vfn;
end

g = k(pol_indx); % policy function

plot(k,vfn)
figure
plot(k,g)


