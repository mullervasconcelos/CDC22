%% Optimal jamming probabilities $\beta^\star$ as a function of $\sigma^2$

clc
clear
c=1;
d=1;

n=500;
step=0.01;
beta=zeros(n,1);
for i=1:n
    
    var=i*step;
    
    g=@(x) x.^2.*exp(-x.^2/(2*var))./sqrt(2*pi*var);
    flag=(2*integral(g,sqrt(c),Inf)>=d);
    
    tol=1e-8; % Tolerance
    if flag==0
        beta(i)=0;
    else
        f=@(b) 2*integral(g,sqrt(c/(1-b)),Inf)-d;
        [beta(i),f_opt,stepNum]=goldenOpt(f,0,1,tol);
    end
end
%% Plot
figure
plot([step:step:n*step],beta)
xlabel('$\sigma^2$','interpreter','latex')
ylabel('$\beta^\star$','interpreter','latex')
grid on
set(gca,'xtick',[0.5:0.5:5],'ytick',[0.1:0.1:1])

function [x_opt,f_opt,stepNum] = goldenOpt(f,a,b,Theta_error)

r=(sqrt(5)-1)/2;
a1=b-r*(b-a);
a2=a+r*(b-a);
stepNum=0;
while abs(b-a)>Theta_error
    stepNum=stepNum+1;
    f1=feval(f,a1);
    f2=feval(f,a2);
    if f1>0
        a=a1;
    end
    if f2<0
        b=a2;
    end
    a1=b-r*(b-a);
    a2=a+r*(b-a);
end
x_opt=(a+b)/2;
f_opt=feval(f,x_opt);
end
