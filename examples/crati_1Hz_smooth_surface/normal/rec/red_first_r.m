dt=10^-3;
nt=10000;
tt=load('./rec_1.mat');
R1=tt.data(:,1);


tt=load('./rec_2.mat');
R2=tt.data(:,1);


tt=load('./rec_3.mat');
R3=tt.data(:,1);

figure;
subplot(3,1,1)
plot(dt:dt:dt*nt,R1)
subplot(3,1,2)
plot(dt:dt:dt*nt,R2)
subplot(3,1,3)
plot(dt:dt:dt*nt,R3)