%%------------Optimal Control Problem By Hybrid-Legendre Function-------------
%-----Minimize J(y)=integral(0^1)(1/2*y'-y*g(x))dx-----
%-----With y(0)=0 ,y'(0)=0 ,y'(1)=0-----
%-----g(x)=-1 if 0<=0<1/4 & g(x)=3 if 1/4<=x<1/2 & g(x)=-1 if 1/2<=x<=1----
%%==========PLOT REAL SOLUTIONS==========
N=input('\Enter order of Block Pulse Function: \');
M=input('\Enter order of Legendre Polynomials: \');
xlim=([-0.5,1.5]);
ylim=([-1,1]);
X1=0:0.0001:0.25;
Y1=1/2.*X1.*X1;
plot(X1,Y1,'LineWidth',2,'color','g');
grid on;
grid minor;
hold on;
X2=1/4:0.0001:0.5;
Y2=(-3/2.*X2.*X2)+X2-(1/8);
plot(X2,Y2,'LineWidth',2,'color',[0,0.75,0.3]);
hold on;
X3=1/2:0.0001:1;
Y3=(1/2.*X3.*X3)-X3+(3/8);
plot(X3,Y3,'LineWidth',2,'color',[0,0.53,0.4]);
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
d=cell(1,M+1);
for i=0:M
    d{i+1}=1/((2*i)+1)*eye(N);
end
D=blkdiag(d{1:end});
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
B0=zeros(1,(M+1)*N);
B0(1:N:end)=1;
B0(N+1)=-1;
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