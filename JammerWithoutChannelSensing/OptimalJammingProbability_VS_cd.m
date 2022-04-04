%% Optimal jamming probability for the  jammer without channel sensing $\beta^\star$ as a function of $c$ and $d$

clc
clear
c=1;
var=1;

n=100;
step=0.01;
beta=zeros(n,n);
for i=1:n
    i
    d=i*step;
    for j=1:n
        c=j*step;
        
        g=@(x) x.^2.*exp(-x.^2/(2*var))./sqrt(2*pi*var);
        flag=(2*integral(g,sqrt(c),Inf)>=d);
        
        tol=1e-6; % Tolerance
        if flag==0
            beta(i,j)=0;
        else
            f=@(b) 2*integral(g,sqrt(c/(1-b)),Inf)-d;
            [beta(i,j),f_opt,stepNum]=goldenOpt(f,0,1,tol);
        end
    end
end
%% Plot
figure
% imagesc([step:step:n*step],[step:step:n*step],beta')
surf([step:step:n*step],[step:step:n*step],beta')
% colormap(gray)
ax=gca;
ax.YDir = 'normal'
xlabel('$d$','interpreter','latex')
ylabel('$c$','interpreter','latex')
zlabel('$\varphi^\star$','interpreter','latex')
% grid on
set(gca,'xtick',[0:0.2:1],'ytick',[0:0.2:1],'ztick',[0:0.2:1])

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
