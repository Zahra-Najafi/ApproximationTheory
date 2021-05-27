%===============Solve Optimal Control by Legendre Polynomials===============
%% Minimize J=1/2*integral(0^1)(x^2+u^2)dt subject to
%%  x'(t)=-x(t)+u(t)
%%  x(0)=1
%% ---------------Example from Orthogonal Functions in systems and Control by Datta K.B & Mohan B.M
%% Page 233 exe 5 ---------------
%% =o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=
m=input('\Enter m :\');
l=0;
u=1;
h=u/m;
%% ===============Creat e vector===============
E=eye(m);
e=E(:,1);
%% ===============Creat D & P Matrix===============
%%==========Creat P Matrix==========
%%==========P is three diagonal matrix==========
p=diag([1,zeros(1,m-1)]);
%disp(p);
d1=zeros(1,m);
for i=1:length(d1)
    d1(i)=(-1)*(1/((2*i)-1));
end
Di=d1;
d2=(-1)*d1;
d1(1)=[];
%disp(d1);
d2(end)=[];
p1=diag(d1,-1);
%disp(p1);
p2=diag(d2,1);
%disp(p2);
p=1/2*(p1+p+p2);
%disp(p);
%%==========Creat D Matrix is diagonal Matrix==========
D=diag((-1)*Di);
%% ===============Solve the System Equations===============
%% ----------Creat Coefficient Matrix----------
A=(2*p*D*p')+D+(D*p')+(p*D);
%% ----------Creat RHS Vector----------
b=(-1)*((2*p*D*e)+(D*e));
a=A\b;
x=((a'*p)+e');
%%
%%==========Legendre polynomials==========
syms t;
px=cell(1,m);
for i=1:m
    px{i}=legendreP((i-1),((2*t)-1));
end
%% ===============Plot Real solution===============
subplot(2,1,1);
syms t;
f=cosh(sqrt(2)*t)+(-0.98)*(sinh(sqrt(2)*t));
ez=ezplot(f,[l,u]);
grid on;
grid minor;
set(ez,'LineWidth',2,'color',[0.5,0,0.77]);
title('(cosh \surd 2 t )-(0.98 sinh \surd 2 t)');
xlabel('t');
ylabel('x(t)');
%% ===============Plot Real and Approximation solutions===============
subplot(2,1,2);
pl=0;
for i=1:m
pl=x(i)*px{i}+pl;
end
%disp(pl);
f=cosh(sqrt(2)*t)+(-0.98)*(sinh(sqrt(2)*t));
ez=ezplot(f,[l,u]);
grid on;
grid minor;
hold on;
set(ez,'LineWidth',3,'color',[0.5,0,0.77]);
X=l:0.001:u;
pll=subs(pl,X);
plot(X,pll,'--','LineWidth',2,'color','k');
grid on;
hold on;
title(['Approximation Solution for m =',num2str(m)]);
xlabel('t');
ylabel('x(t)');
legend('Analytic Solution','Approximation Solution');