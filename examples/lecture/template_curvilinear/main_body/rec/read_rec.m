clear all;
close all;
dt=10^-3;
nt=1000;
%% load data
tt=load('./rec_1.mat');
V1=tt.data;
tt=load('./rec_2.mat');
V2=tt.data;
tt=load('./rec_3.mat');
V3=tt.data;
nr=size(V32);
%% plot data
figure;
for i=1:nr
    subplot(1,nr,i);
    plot(V1(:,i),dt:dt:dt*nt);
    set(gca,'ydir','reverse');
    xlabel('v1 [m/s]');
    if i==1
        ylabel('t [s]');
    end
end

