N=input('\ Enter N :\');
M=input('\ Enter M : \');
TF=input('\Enter upper bound interval :\');
m=M;
n=N;
syms x;
for i=1:n
    for j=0:m
        subplot(n,1,i);
        func=legendreP(j,(((2*n)/TF)*x)-((2*i)-1));
        xlim([0,1]);
        f=fplot(func,[((i-1)/n)*TF,(i/n)*TF]);
        set(f,'LineWidth',1.5);
        set(f,'color','b');
        hold on;
        %%
        r=randi([0,1],1,3);
        col=r;
        f1=fplot((0*x),[(i/n)*TF,TF]);
        set(f1,'LineWidth',1.5);
        set(f1,'color','b');
        grid on;
        hold on;
        %%
        f2=fplot((0*x),[0,((i-1)/n)*TF]);
        set(f2,'LineWidth',1.5);
        set(f2,'color','b');
        grid on;
        hold on;
    end
end