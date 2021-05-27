%%------------Optimal Control Problem By Legendre polynomials-------------
%%----------integral 0-1 (x'^2+tx'+x^2)dt----------
%%----------By Approximation Legendre polynomials----------
m=input('\n Enter m:\n');
a=input('\n Enter LowerBound:\n');
b=input('\n Enter UpperBound:\n');
%%----------Plot Real Solutions----------
X=a:0.00001:b;
Y=(-1/4*(X.^2))+(1/2*X);
plot(X,Y,'LineWidth',2,'color','k');
grid on;
hold on;
%%---------- Creat P Matrix----------
%%----------P is three diagonal matrix----------
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
%%----------Creat D Matrix is diagonal Matrix----------
D=diag((-1)*Di);
%disp(D);
I=eye(m);
%%----------Creat d Vector----------
d=(1/2)*[1 ;1;zeros(m-2,1)];
%disp(d);
%%----------Solve System Equations---------
%%----------Creat Main Matrix---------
%%-----This matrix is Partitioned Matrix 2*2-----
A=[2*D I(:,1);I(1,:) 0];
%disp(A);
%%------Creat Right Hand Side Vector(RHS=b)-----
b=[(-1)*D*d;1/4];
a=A\b;
%disp(a);
a(end)=[];
%disp(a);
a=a'*p;
%disp(a);
%%---------- Legendre polynomials----------
syms t;
Leg=cell(1,m);
L0=(0*t)+1;
L1=t;
Leg{1}=L0;
Leg{2}=L1;
k=3;
while (k<=m)
    L=((((2*(k-2))+1)/((k-2)+1))*t*L1)-((k-2)/((k-2)+1)*L0);
    Leg{k}=L;
    L0=L1;
    L1=L;
    k=k+1;
end
%%----------Plot Legendre polynomials----------
f=0;
j=1;
while (j<=m)
    f=f+((a(j)*Leg{j}));
    j=j+1;
end
% disp(f);
ez=ezplot(f,[0,1]);
set(ez,'color','r');
title(['for m =',num2str(m)]);
hold on;
legend('RealSolution','Approximation solution')