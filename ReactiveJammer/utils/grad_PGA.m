function [g] = grad_PGA(alpha,beta,xhat_0,xhat_1)

global c d

%xhat_1 = 0.2;
%xhat_0 = -0.2;

%beta = 0.1;
%alpha = 0.4;

%c = 0.5;
%d = 0.3;


f1 = @(x) ((x-xhat_1).^2-(x-xhat_0).^2-d).*(beta*(x-xhat_1).^2+c-d*beta > alpha*(x-xhat_1).^2 + (1-alpha)*(x-xhat_0).^2-d*alpha).*exp(-x.^2/2)/sqrt(2*pi);

f2 = @(x) ((x-xhat_1).^2-d).*(beta*(x-xhat_1).^2+c-d*beta <= alpha*(x-xhat_1).^2 + (1-alpha)*(x-xhat_0).^2-d*alpha).*exp(-x.^2/2)/sqrt(2*pi);

g1 = integral(f1,-Inf,Inf);

g2 = integral(f2,-Inf,Inf);

g = [g1; g2];