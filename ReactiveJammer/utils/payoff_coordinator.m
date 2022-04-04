function [J] = payoff_coordinator(theta)

global alpha beta c d

xhat0 = theta(1);
xhat1 = theta(2);

f = @(x) min((1-alpha)*(x-xhat0).^2 + alpha*(x-xhat1).^2 - d*alpha,beta*(x-xhat1).^2+c-d*beta).*exp(-x.^2/2)/sqrt(2*pi);

J = integral(f,-Inf,Inf);

end