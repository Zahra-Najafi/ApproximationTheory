n=input('\n Enter n: \n');
syms x;
L0=(0*x)+1;
axis([-2 2 -2 2]);
ez1=ezplot(L0,[-1,1]);
set(ez1,'color','b');
grid on;
hold on;
L1=x;
ez2=ezplot(L1,[-1,1]);
set(ez2,'color','g');
hold on;
for i=1:n
    L=((((2*n)+1)/(n+1))*x*L1)-(n/(n+1)*L0);
    ez=ezplot(L,[-1,1]);
    L0=L1;
    L1=L;
    hold on;
    N=num2str(n);
    title(['Legendre polynomials for n=1 to',N])
    drawnow;
end
