%%--------Integral Operational Matrix for Block Pulse Functions by Expansion of the function--------
m=input('\n Enter m:\n');
a=input('\n Enter LowerBound:\n');
b=input('\n Enter UpperBound:\n');
%%==========Plot solution of differential equations==========
%%==========y"(t)-y(t)=t+1 & y(0)=0 & y'(0)=1 & 0=<t<1==========
X=a:0.00001:b;
Y=X+1;
plot(X,Y,'LineWidth',2,'color','k');
grid on;
hold on;
%%==========Creat Operational Matrix==========
h=b/m;
p=diag(h/2*ones(1,m));
for i=m+1:m:m*m
    for j=i:m+1:m*m
    p(j)=h;
    end
end
%%==========Creat E Vector==========
e=ones(m,1);
%%==========Creat d Vector==========
d=zeros(m,1);
for k=1:m
 d(k)=((2*k)-1)*(h/2);   
end
b=d-(p'*e);
A=(eye(m)+(p*p))'\b;
y=A'*p;
%%==========Plot Block Pulse==========
for r=1:m
    x1=[(r-1)*h r*h];
    y1=[y(r),y(r)];
    plot(x1,y1,'LineWidth',2.5,'color','m');
    %%
    x2=[(r-1)*h (r-1)*h];
    y2=[0 y(r)];
    plot(x2,y2,'--','LineWidth',1,'color','m');
    %%
    x3=[r*h r*h];
    y3=[0 y(r)];
    plot(x3,y3,'--','LineWidth',1,'color','m');
    %%
     xlabel('t');
    ylabel('y(t)');
    title(['Block Pulse Function for m =',num2str(m)]);
    hold on;
    drawnow;
end