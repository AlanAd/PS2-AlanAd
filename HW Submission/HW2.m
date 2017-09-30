close all
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;

%%% Set up Technology Shocks
a=1.1;
b=0.678;
A= [a b];

pi = [0.977 0.023 ; 0.074 0.926]; 

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Good State Return Function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
consh = A(1,1)*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 

reth = consh .^ (1 - sigma) / (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
reth(consh < 0) = -Inf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bad State Return Function %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
consl = A(1,2)*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 

retl = consl .^ (1 - sigma) / (1 - sigma); % return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
retl(consl < 0) = -Inf;

%%%% Iteration
disH = 1; disL=1; tol = 1e-06; % tolerance for stopping 
v_guessH = zeros(1, num_k);
v_guessL = zeros(1, num_k);

while (disH > tol) && (disL > tol)
    % compute the utility value for all possible combinations of k and k':
    
    expvalH = pi(1,1)*v_guessH + pi(1,2)*v_guessL;
    expvalL = pi(2,1)*v_guessH + pi(2,2)*v_guessL;
    
    value_math = reth + beta * repmat(expvalH, [num_k 1]);
    value_matl = retl + beta * repmat(expvalL, [num_k 1]);
     
    % find the optimal k' for every k:
    [vfn, pol_indxH] = max(value_math, [], 2);
    vfnH = vfn';
    
     % find the optimal k' for every k:
    [vfn, pol_indxL] = max(value_matl, [], 2);
    vfnL = vfn';
    
    % what is the distance between current guess and value function
    disH = max(abs(vfnH - v_guessH));
    disL = max(abs(vfnL - v_guessL));
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guessH = vfnH;
    v_guessL = vfnL;
end

gH = k(pol_indxH); % policy function
gL = k(pol_indxL); % policy function

plot(k,vfnH,k, vfnL)
title('Value Functions');
legend('High A', 'Low A', 'Location','southeast');

figure
plot(k,gH, k, gL)
title('Policy Functions');
legend('High A', 'Low A', 'Location','southeast');

%%% Simulate Output

%%% This calibrations gets standard deviation of output =1.8% per quarter.
a=1.0;
b=0.908;
A= [a b];

Av=zeros(1000,1);

Av(1,1)=A(1,1);

%%% Set Seed
rng(9989);
p = rand(1,1000);
K=zeros(1000,1);
K(1,1)=10;

for i= 2:1000 
    [c,index]=min(abs(K(i-1,1)-k));
    
    if Av(i-1,1)==A(1,1)
     if p(i-1)>= pi(1,1)
      Av(i,1)=A(1,1);
      K(i,1)=gH(index);
    
     else
      Av(i,1)= A(1,2);
      K(i,1)= gL(index);
     end
    else
        if p >= pi(2,2)
        Av(i,1)=A(1,2);
       K(i,1)=gL(index);
     else
        Av(i,1)= A(1,1);
        K(i,1)=gH(index);
        end
    end
end
    
Y=zeros(1000,1);
c=linspace(1,1000,1000);

for i = 1:1000
    Y(i,1)=Av(i,1)*(K(i,1)^alpha);

    end


plot(c,Y)
title('Output (quarterly) with \sigma = 1.8%');


