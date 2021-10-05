dt=10^-3;
nt=10000;
tt=load('./rec_2.mat');
R=tt.data;
%%
cqwva(R,dt:dt:dt*nt,1:size(R,2),1,1,3,'black','black','new','max');
title('v2 [m/s]');
xlabel('N receiver');
ylabel('t [s]');
%%