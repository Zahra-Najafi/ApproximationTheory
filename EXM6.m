%%----------Integral Operational Matrix for convolution integral----------
Maxm=input('\n Enter m:\n');
a=input('\n Enter LowerBound:\n');
b=input('\n Enter UpperBound:\n');
%%--------------------------Computing f_i----------------------------------
s=a;
e=b;
for m=1:1:Maxm 
%%==========Plot f(t)==========
X=s:0.001:e;
Y=sin(2.*X);
plot(X,Y,'LineWidth',2,'color','k');
hold on;
h=b/m;
%%==========Main Loop===========
for i=0:h:b-h
    lb=i;
    ub=lb+h;
    fi=(1/h)*integral(@(x)(sin(2.*x)),lb,ub);
    col=[1 0.2 0.5];
    x1=[ub ub];
    y1=[0 fi];
    plot(x1,y1,'--','LineWidth',1,'color',col);
    x2=[lb lb];
    y2=[0 fi];
    plot(x2,y2,'--','LineWidth',1,'color',col);
    x3=[lb ub];
    y3=[fi fi];
    xlabel('t');
    ylabel('f(t)');
    title(['Block Pulse Function for m =',num2str(m)]);
    plot(x3,y3,'LineWidth',2.5,'color',col);
    grid on;
    hold on;
end
hold off;
drawnow;
end