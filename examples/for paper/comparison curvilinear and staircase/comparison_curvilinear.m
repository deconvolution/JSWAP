clear all;
close all;
dt=10^-3;
nt=10^3;

i=1;
tt=load('../staircase/main_body/rec/rec_1.mat');
v1_con=tt.data(1:nt);
tt=load('../curvilinear/main_body/rec/rec_1.mat');
v1_cur=tt.data(1:nt);

shift=1;
v1_con(shift:end)=v1_con(1:end-shift+1);

tt=load('../staircase/main_body/rec/rec_2.mat');
v2_con=tt.data(1:nt);
tt=load('../curvilinear/main_body/rec/rec_2.mat');
v2_cur=tt.data(1:nt);

shift=1;
v2_con(shift:end)=v2_con(1:end-shift+1);

tt=load('../staircase/main_body/rec/rec_1.mat');
v3_con=tt.data(1:nt);
tt=load('../curvilinear/main_body/rec/rec_3.mat');
v3_cur=tt.data(1:nt);

shift=1;
v3_con(shift:end)=v3_con(1:end-shift+1);

figure;
subplot(1,3,1)
ax=plot(v1_con,dt:dt:dt*nt,'color','red','linewidth',3);


set(gca,'Ydir','reverse');
xlabel('v1 [m/s]');
ylabel('t [s]');
%ylim([.8,3])
hold on;
ax2=plot(v1_cur,dt:dt:dt*nt,'color','blue','linewidth',3);


subplot(1,3,2)
ax=plot(v2_con,dt:dt:dt*nt,'color','red','linewidth',3);


set(gca,'Ydir','reverse');
xlabel('v2 [m/s]');
ylabel('t [s]');
%ylim([.8,3])
hold on;
ax2=plot(v2_cur,dt:dt:dt*nt,'color','blue','linewidth',3);


subplot(1,3,3)
ax=plot(v3_con,dt:dt:dt*nt,'color','red','linewidth',3);


set(gca,'Ydir','reverse');
xlabel('v3 [m/s]');
ylabel('t [s]');
%ylim([.8,3])
hold on;
ax2=plot(v3_cur,dt:dt:dt*nt,'color','blue','linewidth',3);


legend([ax,ax2],'Conventional','Curvilinear','Location','southeast');
%%
corrcoef(v1_con,v1_cur)
corrcoef(v2_con,v2_cur)
corrcoef(v3_con,v3_cur)