tt=load('./rec_1.mat');
R1=tt.data;
figure;
imagesc([1,50],[0,10^-3*1000],R1);
xlabel('receiver');
ylabel('time [s]')