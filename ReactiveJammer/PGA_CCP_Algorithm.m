%% PGA_CCP convergence

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

c = 1; % communication cost
d = 1; % jamming cost
var=1; % variance


a = 0.1; % step size

theta = [alpha;beta];

flag = 0;

k=1;

while (flag==0) && (k<=1e4)
    
    
    %% projected gradient ascent PGA
    %   v = theta + (t/sqrt(k))*proj_sub_grad(theta(1),theta(2),xhat_0,xhat_1);
    v = theta + a*grad_PGA(theta(1),theta(2),xhat_0,xhat_1);
    
    theta_new=max(0,min(v,1)); % Projection
    
    theta = theta_new;
    
    alpha = theta(1);
    
    beta = theta(2);
    %% convex-concave procedure CCP
    
    A = [2*(1-alpha) 0; 0 2*(beta+alpha)];
    
    g = grad_CCP(alpha,beta,xhat_0,xhat_1);
    
    xhat_new = pinv(A)*g;
    
    xhat = xhat_new;
    
    xhat_0 = xhat(1);
    
    xhat_1 = xhat(2);
    
    %% Stopping criteria
    [Delta1,Delta2] = FirstNashEquilibriumChecker(alpha,beta,xhat_0,xhat_1);
    
    %conditions = [conditions; max(Delta1,Delta2)];
    
    conditions = [conditions; Delta1 Delta2];
    
    
    flag=(max(Delta1,Delta2)<=10^-6); % Check whether this is a epsilon first nash equilibrium
    
    k=k+1;
end

%% Plot
figure

[theta' xhat']  %pairs of epsilon-FNE

k

semilogy(conditions)


