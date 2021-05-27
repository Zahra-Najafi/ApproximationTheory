%%------------Optimal Control Problem By Block Pulse Function-------------
%%----------integral 0-1 (x'^2+tx'+x^2)dt----------
m=input('\n Enter m:\n');
a=input('\n Enter LowerBound:\n');
b=input('\n Enter UpperBound:\n');
%%----------Plot Real Solutions----------
X=a:0.00001:b;
Y=(-1/4*(X.^2))+(1/2*X);
plot(X,Y,'LineWidth',2,'color','k');
grid on;
hold on;
%%----------Solve a System of equations----------
%%-----P Matrix-----
h=b/m;
p=diag(h/2*ones(1,m));
for i=m+1:m:m*m
    for j=i:m+1:m*m
    p(j)=h;
    end
end
%%----- E Vector-----
e=ones(m,1);
%%-----d Vector-----
d=zeros(m,1);
for k=1:m
 d(k)=((2*k)-1)*(h/2);   
end
%%------Creat Right Hand Side Vector(RHS=b)-----
b=[(-1)*(h*eye(m)*d);1/4];
%disp(b);
%%-----Creat Main Matrix-----
%%-----This matrix is Partitioned Matrix 2*2-----
A=[(2*h*eye(m))+2*p*h*eye(m)*p' h*e;h*e' 0];
%disp(A);
%%-----solve-----
a=A\b;
n=length(a);
a(n,:)=[];
%disp(a);
a=a'*p;
%disp(a);
%%-----Plot by Block Pulse Function-----
for r=1:m
    x1=[(r-1)*h r*h];
    y1=[a(r),a(r)];
    plot(x1,y1,'LineWidth',2.5,'color','b');
    %%
    x2=[(r-1)*h (r-1)*h];
    y2=[0 a(r)];
    plot(x2,y2,'--','LineWidth',1,'color','b');
    %%
    x3=[r*h r*h];
    y3=[0 a(r)];
    plot(x3,y3,'--','LineWidth',1,'color','b');
    %%
    xlabel('t');
    ylabel('y(t)');
    title(['Block Pulse Function for m =',num2str(m)]);
    hold on;
    drawnow;
end