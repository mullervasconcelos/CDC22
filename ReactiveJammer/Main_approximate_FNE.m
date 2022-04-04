%% PGA-CCP algorithm

clear all
clc
addpath(genpath('./utils/'));

global c d var

%initialization
xhat_1 = rand;
xhat_0 = rand;

alpha = 0.5;
beta = 0.5;


c = 1; % communication cost
d = 1; % jamming cost
var=1; % variance

t = 0.1; % step size

theta = [alpha;beta];

flag = 0;

while flag==0
    %% projected gradient ascent PGA
    %   v = theta + (t/sqrt(k))*grad_PGA(theta(1),theta(2),xhat_0,xhat_1);
    v = theta + t*grad_PGA(theta(1),theta(2),xhat_0,xhat_1);
    
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
    
    flag=((Delta1<=10^-5)&&(Delta2<=10^-5)); % Check whether this is an approximate first nash equilibrium
    
    
end

%% Saddle point
[theta' xhat']




