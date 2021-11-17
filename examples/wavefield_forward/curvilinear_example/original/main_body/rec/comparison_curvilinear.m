clear all;
close all;
dt=10^-3;
nt=3*10^3;

i=1;
tt=load('./rec_1.mat');
v1_con=tt.data;
tt=load('./curvilinear/rec_1.mat');
v1_cur=tt.data;

shift=1;
v1_con(shift:end)=v1_con(1:end-shift+1);

tt=load('./rec_2.mat');
v2_con=tt.data;
tt=load('./curvilinear/rec_2.mat');
v2_cur=tt.data;

shift=1;
v2_con(shift:end)=v2_con(1:end-shift+1);

tt=load('./rec_1.mat');
v3_con=tt.data;
tt=load('./curvilinear/rec_3.mat');
v3_cur=tt.data;

shift=1;
v3_con(shift:end)=v3_con(1:end-shift+1);

figure;
subplot(1,3,1)
ax=plot(v1_con,dt:dt:dt*nt,'color','red','linewidth',3);
set(gca,'Ydir','reverse');
xlabel('v1 [m/s]');
ylabel('t [s]');
ylim([.8,3])
hold on;
ax2=plot(v1_cur,dt:dt:dt*nt,'color','blue','linewidth',3);
subplot(1,3,2)
ax=plot(v2_con,dt:dt:dt*nt,'color','red','linewidth',3);
set(gca,'Ydir','reverse');
xlabel('v2 [m/s]');
ylabel('t [s]');
ylim([.8,3])
hold on;
ax2=plot(v2_cur,dt:dt:dt*nt,'color','blue','linewidth',3);
subplot(1,3,3)
ax=plot(v3_con,dt:dt:dt*nt,'color','red','linewidth',3);
set(gca,'Ydir','reverse');
xlabel('v3 [m/s]');
ylabel('t [s]');
ylim([.8,3])
hold on;
ax2=plot(v3_cur,dt:dt:dt*nt,'color','blue','linewidth',3);
legend([ax,ax2],'Conventional','Curvilinear','Location','southeast');