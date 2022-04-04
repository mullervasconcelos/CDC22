function [Delta1,Delta2] = FirstNashEquilibriumChecker(alpha,beta,xhat_0,xhat_1)

g=grad_PGA(alpha,beta,xhat_0,xhat_1);

f=@(a,b) g'*([a;b]-[alpha;beta]);

Delta1=max([f(0,0),f(0,1),f(1,0),f(1,1)]);

Delta2=norm([2*(1-alpha)*xhat_0;2*(alpha+beta)*xhat_1]- grad_CCP(alpha,beta,xhat_0,xhat_1));
