%%------------Optimal Control Problem By Hybrid-Legendre Function-------------
%-----Minimize J(y)=integral(0^1)(1/2*y'-y*g(x))dx-----
%-----With y(0)=0 ,y'(0)=0 ,y'(1)=0-----
%-----g(x)=-1 if 0<=0<1/4 & g(x)=3 if 1/4<=x<1/2 & g(x)=-1 if 1/2<=x<=1----
%%==========GET INPUTS==========
N=input('\Enter order of Block Pulse Function: \');
M=input('\Enter order of Legendre Polynomials: \');
%%==========COMPUTING OPRATIONAL MATRIX==========
TF=1;
h=TF/N;
P0=diag(1/2*ones(1,N));
for i=N+1:N:N*N
    for j=i:N+1:N*N
   P0(j)=1;
    end
end
%disp(P0);
%%
P1=[P0 zeros(N,N*M);zeros(N*M,N) zeros(N*M,N*M)];
%disp(P1);
%%
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
%%==========COOMPUTING D MATRIX==========
if M==0
    D=eye(4);
else
d=cell(1,M+1);
for i=0:M
    d{i+1}=1/((2*i)+1)*eye(N);
end
D=blkdiag(d{1:end});
end
%disp(D);
D=(TF/N)*D;
%%==========COMPUTING V VECTOR==========
V1=zeros(1,(M+1)*N);
V1(1:N)=1/N;
%disp(V1);
%==========V2=V(1/4)==========
V2=zeros(1,(M+1)*N);
V2(1)=1/N;
%disp(V2);
%==========V3=V(2/4)=V(1/2)==========
V3=zeros(1,(M+1)*N);
V3(1:2)=1/N;
%disp(V3);
%==========CREAT B VECTOR==========
if M==0
    B0=[1 0 0 0];
else
B0=zeros(1,(M+1)*N);
B0(1:N:end)=1;
B0(N+1)=-1;
end
%disp(B0);
B1=zeros(1,(M+1)*N);
B1(N:N:end)=1;
%disp(B1);
%==========SOLVE THE SYSTEMS EQUATIONS==========
%----------CREAT THE RIGHT HAND SIDE VECTOR----------
b=[(-1)*((4*P*V2')-(4*P*V3')+(P*V1')) ; zeros(2,1)];
%disp(b);
%----------CREAT THE COEFFICIENT MATRIX----------
A=[D B0' B1';B0 0 0;B1 0 0];
%disp(A);
c=A\b;
%disp(c);
c(end-1:end)=[];
%disp(c);
%%==========COMPUTING Y(X)==========
y=c'*P;
%disp(y);
%%==========CREAT P(X)==========
px=cell(N*(M+1),1);
%%
 syms x;
for i=0:M
    for j=1:N
        px{(N*i)+j}=legendreP(i,((2*N)/TF)*x-(2*j)+1);
    end
end
%disp(px);
%%==========PLOT REAL SOLUTION==========
subplot(3,1,1)
%xlim=([-0.5,1.5]);
%ylim=([-1,1]);
X1=0:0.0001:0.25;
Y1=(1.*X1.*X1)/2;
plot(X1,Y1,'LineWidth',2,'color','g');
grid on;
grid minor;
hold on;
X2=1/4:0.0001:0.5;
Y2=((-3.*X2.*X2)/2)+X2-(1/8);
plot(X2,Y2,'LineWidth',2,'color',[0,0.75,0.3]);
hold on;
X3=1/2:0.0001:1;
Y3=((1.*X3.*X3)/2)-X3+(3/8);
plot(X3,Y3,'LineWidth',2,'color',[0,0.53,0.4]);
title('Real Solution');
%%==========PLOT REAL AND APPROXIMATION SOLUTIONS==========
subplot(3,1,2)
xlim=([-0.5,1.5]);
ylim=([-1,1]);
X1=0:0.0001:0.25;
Y1=(1.*X1.*X1)/2;
plot(X1,Y1,'LineWidth',2,'color','g');
grid on;
grid minor;
hold on;
X2=1/4:0.0001:0.5;
Y2=((-3.*X2.*X2)/2)+X2-(1/8);
plot(X2,Y2,'LineWidth',2,'color',[0,0.75,0.3]);
hold on;
X3=1/2:0.0001:1;
Y3=((1.*X3.*X3)/2)-X3+(3/8);
plot(X3,Y3,'LineWidth',2,'color',[0,0.53,0.4]);
%%Approximation Solutions
YP=cell(N,1);
for i=1:N-1
Y=0;
for j=i:N:N*(M+1)
   Y=(y(j)*px{j})+Y;
end
YP{i}=Y;
end
YP{end}=YP{N-1};
%disp(YP);
for i=1:N
    XP=((i-1)/N)*TF:0.0001:(i/N)*TF;
    PL=plot(XP,subs(YP{i},XP),'--','LineWidth',2.5);
    set(PL,'color','k');
    hold on;
end
title('Real and Approximation Solutions');
%%==========PLOT Approximation Solution==========
subplot(3,1,3);
YP=cell(N,1);
for i=1:N-1
Y=0;
for j=i:N:N*(M+1)
   Y=(y(j)*px{j})+Y;
end
YP{i}=Y;
end
YP{end}=YP{N-1};
%disp(YP);
for i=1:N
    XP=((i-1)/N)*TF:0.0001:(i/N)*TF;
    grid on;
    grid minor;
    PL=plot(XP,subs(YP{i},XP),'--','LineWidth',2.5);
    set(PL,'color',[0,0.66,0.3]);
    hold on;
end
title(['Approximation Solution for m=',num2str(M),'and n=',num2str(N)]);