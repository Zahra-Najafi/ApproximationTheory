%---------------Plot Chebyshev Polynomials of the First Kind---------------
n=input('\Enter n : \');
syms x;
axis([-1.5 1.5 -2 2]);
T0=(0*x)+1;
T1=x;
ez=ezplot(T0,[-1,1]);
set(ez,'LineWidth',1.2);
set(ez,'color','b');
grid on;
hold on;
ez1=ezplot(T1,[-1,1]);
set(ez1,'LineWidth',1.2);
set(ez1,'color','y');
for i=2:n
    Tn=(2*x*T1)-T0;
    T0=T1;
    T1=Tn;
    ez2=ezplot(Tn,[-1,1]);
    set(ez2,'LineWidth',1.2);
    hold on;
    %%
    xlabel('x')
    ylabel('T_n(x)')
    title('Chebyshev polynomials of the first kind');
end
