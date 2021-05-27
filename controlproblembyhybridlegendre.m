%===============Solve Optimal Control by Hybrid Legendre Polynomials===============
%% Minimize J=1/2*integral(0^1)(x^2+u^2)dt subject to
%%  x'(t)=-x(t)+u(t)
%%  x(0)=1
%% ---------------Example from Orthogonal Functions in systems and Control by Datta K.B & Mohan B.M
%% Page 233 exe 5 ---------------
%% =o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=o=
%%==========GET INPUTS==========
N=input('\Enter order of Block Pulse Function: \');
M=input('\Enter order of Legendre Polynomials: \');
%% ==========COMPUTING OPRATIONAL MATRIX==========
l=0;
u=1;
TF=1;
h=TF/N;
P0=diag(1/2*ones(1,N));
for i=N+1:N:N*N
    for j=i:N+1:N*N
   P0(j)=1;
    end
end
P1=[P0 zeros(N,N*M);zeros(N*M,N) zeros(N*M,N*M)];
L=(N*M)+N-2;
Cell=cell(1,M);
for i=0:M-1
p2=1/(2*((2*i)+1))*eye(N);
Cell{i+1}=p2;
end
%disp(Cell);
P2=[zeros(N*M,N) blkdiag(Cell{1:end});
   zeros(N,N) zeros(N,N*M)];
%disp(P2);
%%
Cell=cell(1,M);
for i=1:M
p3=(-1)*(1/(2*((2*i)+1)))*eye(N);
Cell{i}=p3;
end
%disp(Cell);
P3=[zeros(N,N*M) zeros(N,N);
   blkdiag(Cell{1:end}) zeros(N*M,N)];
%disp(P3);
%%
P=P1+P2+P3;
%disp(P);
%%
P=(TF/N)*P;
%% ==========COMPUTING D MATRIX==========
d=cell(1,M+1);
for i=0:M
    d{i+1}=1/((2*i)+1)*eye(N);
end
D=blkdiag(d{1:end});
%disp(D);
D=(TF/N)*D;
%% ==========COMPUTING e VECTOR==========
I=eye(N*(M+1));
e=I(:,1);
%% ===============Solve the System Equations===============
%% ----------Creat Coefficient Matrix----------
A=(2*P*D*P')+D+(D*P')+(P*D);
%% ----------Creat RHS Vector----------
b=(-1)*((2*P*D*e)+(D*e));
c=A\b;
x=((c'*P)+e');
%%==========CREAT P(X)==========
px=cell(N*(M+1),1);
%%
 syms t;
for i=0:M
    for j=1:N
        px{(N*i)+j}=legendreP(i,((2*N)/TF)*t-(2*j)+1);
    end
end
%% ===============Plot Real solution===============
subplot(3,1,1);
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
subplot(3,1,2);
syms t;
f=cosh(sqrt(2)*t)+(-0.98)*(sinh(sqrt(2)*t));
ez=ezplot(f,[l,u]);
grid on;
grid minor;
hold on;
set(ez,'LineWidth',3,'color',[0.5,0,0.77]);
pl=0;
for i=1:N*(M+1)
pl=(x(i)*px{i})+pl;
end
X=l:0.001:u;
pll=subs(pl,X);
plot(X,pll,'--','LineWidth',2,'color','k');
hold on;
xlabel('t');
ylabel('x(t)');
title('Real and Approximation Solutions');
legend('Analytic Solution','Approximation Solution');
%% ==========PLOT Approximation Solution==========
subplot(3,1,3);
%% Approximation Solutions
X=l:0.001:u;
pl=0;
for i=1:N*(M+1)
pl=(x(i)*px{i})+pl;
end
pll=subs(pl,X);
plot(X,pll,'--','LineWidth',2,'color','k');
grid on;
grid minor;
hold on;
xlabel('t');
ylabel('x(t)');
title(['Approximation Solution for m=',num2str(M),'and n=',num2str(N)]);
disp('Approximation solutions is :')
disp(pl);