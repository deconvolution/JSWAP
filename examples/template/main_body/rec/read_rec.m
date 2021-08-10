dt=10^-3;
nt=1000;
tt=load('./rec_p.mat');
R=tt.data;
%%
cqwva(R,dt:dt:dt*nt,1:size(R,2),1,1,3,'black','black','new','max');
title('v3 [m/s]');
xlabel('N receiver');
ylabel('t [s]');