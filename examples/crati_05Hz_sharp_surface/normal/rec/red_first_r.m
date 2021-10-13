dt=10^-3;
nt=20000;
tt=load('./rec_1.mat');
R1=tt.data(:,1);


tt=load('./rec_2.mat');
R2=tt.data(:,1);


tt=load('./rec_3.mat');
R3=tt.data(:,1);

figure;
subplot(3,1,1)
plot(dt:dt:dt*nt,R1)
xlabel('t [s]');
ylabel('v1 [m/s]');
subplot(3,1,2)
plot(dt:dt:dt*nt,R2)
xlabel('t [s]');
ylabel('v2 [m/s]');
subplot(3,1,3)
plot(dt:dt:dt*nt,R3)
xlabel('t [s]');
ylabel('v3 [m/s]');