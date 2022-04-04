%% GDA algorithm convergence

clear all
clc
addpath(genpath('./utils/'));
warning off

global c d var

xhat_1 = rand;
xhat_0 = rand;

beta = 0.5;
alpha = 0.5;

conditions = [];

c = 1;
d = 1;
var=1;


a1 = 0.1;

a2 = 0.01;

memory = [];

theta = [alpha;beta]

flag = 0;

k=1

while (flag==0) && (k<=5e3)
    
    %% Gradient descent
    xhat = [xhat_0; xhat_1];
    
    
    g = grad_GD(theta(1),theta(2),xhat_0,xhat_1);
    
    xhat_new = xhat - a2*g;
    

    
    %% Projected Gradient Ascent
    v = theta + a1*grad_PGA(theta(1),theta(2),xhat_0,xhat_1);
    theta_new=max(0,min(v,1));
    
    
    xhat = xhat_new;
    xhat_0 = xhat(1);
    xhat_1 = xhat(2);
    
    theta = theta_new;
    alpha = theta(1);
    beta = theta(2);
    
    
    %% Stopping criteria
    [Delta1,Delta2] = FirstNashEquilibriumChecker(theta(1),theta(2),xhat_0,xhat_1);
    
    %    conditions = [conditions; max(Delta1,Delta2)];
    
    conditions = [conditions; Delta1 Delta2];
    
    flag=( max(Delta1,Delta2)<=10^-6) ; % Check whether this is a epsilon first nash equilibrium
    
    k = k+1
end

%% Plot
figure

[theta' xhat']

k

semilogy(conditions)




