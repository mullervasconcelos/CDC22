function [g] = grad_GD(alpha,beta,xhat_0,xhat_1)

global c d var


f1 = @(x) (-2*(1-alpha)*(x-xhat_0)).*(beta*(x-xhat_1).^2+c-d*beta > alpha*(x-xhat_1).^2 + (1-alpha)*(x-xhat_0).^2-d*alpha).*exp(-x.^2/(2*var))./sqrt(2*pi*var);

f2 = @(x) (-2*beta*(x-xhat_1)).*(beta*(x-xhat_1).^2+c-d*beta <= alpha*(x-xhat_1).^2 + (1-alpha)*(x-xhat_0).^2-d*alpha).*exp(-x.^2/(2*var))./sqrt(2*pi*var);

f3 = @(x) (-2*alpha*(x-xhat_1)).*(beta*(x-xhat_1).^2+c-d*beta > alpha*(x-xhat_1).^2 + (1-alpha)*(x-xhat_0).^2-d*alpha).*exp(-x.^2/(2*var))./sqrt(2*pi*var);

g1 = integral(f1,-Inf,Inf);

g2 = integral(f2,-Inf,Inf);

g3 = integral(f3,-Inf,Inf);

g = [g1; g2+g3];