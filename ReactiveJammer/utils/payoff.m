function [J] = payoff(xhat0, xhat1, alpha, beta, c, d)

f = @(x) min((1-alpha)*(x-xhat0).^2 + alpha*(x-xhat1).^2 - d*alpha,beta*(x-xhat1).^2+c-d*beta).*exp(-x.^2/2)/sqrt(2*pi);

J = integral(f,-Inf,Inf);

end