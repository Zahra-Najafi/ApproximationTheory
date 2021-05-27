i=input('\n Enter i :\n');
axis([0 1.5 -1.5  1.5]);
grid on;
    for j=1:1:2^i
        r=(2^i)+j-1;
        %%----------Rationalised Haar----------
        J=zeros(1,3);
        k=1;
        for u=1:-1/2:0
        Ju=(j-u)/(2^i);
        J(k)=Ju;
        k=k+1;
        end
        drawnow;
        x1=[0 1];
        y1=[0 0];
        plot(x1,y1,'LineWidth',2,'color','k');
        hold on;
        grid on;
        %%
        x1=[J(1) J(2)];
        y1=[1 1];
        plot(x1,y1,'LineWidth',2,'color','b');
        %%
        x2=[J(2),J(2)];
        y2=[1,-1];
        plot(x2,y2,'--','LineWidth',2,'color','b');
        %%
        x3=[J(1),J(1)];
        y3=[0,1];
        plot(x3,y3,'--','LineWidth',2,'color','b');
        %%
        x4=[J(3),J(3)];
        y4=[0,-1];
        plot(x4,y4,'--','LineWidth',2,'color','b');
        %%
        x5=[J(2) J(3)];
        y5=[-1 -1];
        plot(x5,y5,'LineWidth',2,'color','b');
        hold on;
        drawnow;
        hold off;
        title(['RH for r=',num2str(r)])
    end