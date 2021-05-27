%===============Solve Optimal Control by block pulse Function===============
%% Minimize J=1/2*integral(0^1)(x^2+u^2)dt subject to
%%  x'(t)=-x(t)+u(t)
%%  x(0)=1
%% ---------------Example from Orthogonal Functions in systems and Control by Datta K.B & Mohan B.M
%% Page 233 exe 5 ---------------
%% =o=o
m=input('\Enter m :\');
l=0;
u=1;
h=u/m;
%% ===============Creat E vector===============
E=ones(m,1);
%% ===============Creat D Matrix===============
D=h*eye(m);
%% ===============Creat P Matrix===============
P=diag(h/2*ones(1,m));
for i=m+1:m:m*m
    for j=i:m+1:m*m
    P(j)=h;
    end
end
%% ===============Solve the System Equations===============
%% ----------Creat Coefficient Matrix----------
A=((3/2)*P*D*P')+D+(D*P')+(P*D);
%% ----------Creat RHS Vector----------
b=(-1)*((2*P*D*E)+(D*E));
a=A\b;
% a(end)=[];
%% ----------Creat x(t)----------
x=(a'*P)+E';
%% ===============Plot Real solution===============
subplot(2,1,1);
%%
syms t;
f=cosh(sqrt(2)*t)+(-0.98)*(sinh(sqrt(2)*t));
ez=ezplot(f,[l,u]);
grid on;
grid minor;
set(ez,'LineWidth',2,'color',[0.5,0,0.77]);
title('(cosh \surd 2 t )-(0.98 sinh \surd 2 t)');
xlabel('t');
ylabel('x(t)');
%% ===============Plot Real ana Approximation solutions===============
subplot(2,1,2);
syms t;
f=cosh(sqrt(2)*t)+(-0.98)*(sinh(sqrt(2)*t));
ez=ezplot(f,[l,u]);
grid on;
grid minor;
set(ez,'LineWidth',2,'color',[0.5,0,0.77]);
hold on;
for i=1:m
    x1=[(i-1)*h i*h];
    y1=[x(i),x(i)];
    plot(x1,y1,'LineWidth',2.5,'color','b');
    %%
    x2=[(i-1)*h (i-1)*h];
    y2=[0 x(i)];
    plot(x2,y2,'--','LineWidth',1,'color','b');
    %%
    x3=[i*h i*h];
    y3=[0 x(i)];
    plot(x3,y3,'--','LineWidth',1,'color','b');
    %%
    xlabel('t');
    ylabel('x(t)');
    title('Real and Approximation Solutions');
    hold on;
    drawnow;
end
legend('Analytic Solution','Approximation Solution');