%%------------Optimal Control Problem By Legendre polynomials-------------
%%----------integral 0-1 (x'^2+tx')dt----------
%%----------By Approximation Legendre polynomials----------
m=input('\n Enter m:\n');
a1=input('\n Enter LowerBound:\n');
b1=input('\n Enter UpperBound:\n');
%%==========Plot Real Solutions==========
X=a1:0.00001:b1;
Y=((-1/4)*(X.^2))+((1/2)*X);
plot(X,Y,'LineWidth',3,'color','k');
grid on;
hold on;
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
%disp(D);
I=eye(m);
%%==========Creat d Vector==========
d=(1/2)*[1 ;1;zeros(m-2,1)];
%disp(d);
%%==========Solve System Equations==========
%%==========Creat Main Matrix==========
%%====This matrix is Partitioned Matrix 2*2====
A=[2*D I(:,1);I(1,:) 0];
%disp(A);
%%=====Creat Right Hand Side Vector(RHS=b)=====
b=[(-1)*D*d;1/4];
a=A\b;
%disp(a);
a(end)=[];
%disp(a);
a=a'*p;
%disp(a);
%%==========Legendre polynomials==========
syms x;
px=cell(1,m);
for i=1:m
    px{i}=legendreP((i-1),((2*x)-1));
end
%%==========Plot Real solution==========
subplot(3,1,1);
X=a1:0.00001:b1;
Y=((-1/4)*(X.^2))+((1/2)*X);
plot(X,Y,'LineWidth',2,'color','k');
grid on;
hold on;
title('Real Solution');
%%==========Plot Real & Approximation Solutions==========
subplot(3,1,2)
X=a1:0.00001:b1;
Y=((-1/4)*(X.^2))+((1/2)*X);
plot(X,Y,'LineWidth',2.5,'color','k');
grid on;
grid minor;
hold on;
pl=0;
for i=1:m
pl=a(i)*px{i}+pl;
end
%disp(pl);
X=a1:0.001:b1;
pl=subs(pl,X);
plot(X,pl,'--','LineWidth',2,'color',[0 0.2 0.8]);
title('Real and Approximation Solutions');
legend('RealSolution','Approximation solution');
%%==========Plot Approximation solution with Legendre Polynomials==========
subplot(3,1,3);
pl=0;
for i=1:m
pl=a(i)*px{i}+pl;
end
%disp(pl);
X=a1:0.001:b1;
pl=subs(pl,X);
plot(X,pl,'--','LineWidth',2,'color',[0 0.2 0.8]);
grid on;
hold on;
title(['Approximation Solution for m =',num2str(m)]);