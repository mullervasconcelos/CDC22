function [J] = payoff_jammer(theta)

global xhat0 xhat1 c d

alpha = theta(1);
beta = theta(2);

f = @(x) min((1-alpha)*(x-xhat0).^2 + alpha*(x-xhat1).^2 - d*alpha,beta*(x-xhat1).^2+c-d*beta).*exp(-x.^2/2)/sqrt(2*pi);

J = -integral(f,-Inf,Inf);

end