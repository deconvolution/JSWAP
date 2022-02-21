clear all;
close all;
dt=.01;
tt=load('./rec/rec_3.mat');
R1=tt.data;
locr=readmatrix('./model/receiver location.csv');
locs=readmatrix('./model/source location.csv');
d=vecnorm(locr-locs,2,2);
[d,I]=sort(d);
index_incre=1;
lvl=1;
clip=10;
line_color=[0,0,0];
face_color='red';
mode='new';
trace_balance='';
cqwva(R1(:,I),dt,d,index_incre,lvl,clip,line_color,face_color,mode,trace_balance)
%%
clear all;
close all;
dt=.01;
pwd
tt=load('./rec/rec_3.mat');
R1=tt.data;
locr=readmatrix('./model/receiver location.csv');
locs=readmatrix('./model/source location.csv');
locr(2,:)=[];
R1(:,2)=[];

d = zeros(52,1);
for i = 1:52
d(i) = norm(locr(i,:)-locs,2);
end
[d,I]=sort(d);
index_incre=1;
lvl=1;
clip=10;
line_color=[0,0,0];
face_color='red';
mode='new';
trace_balance='max';
cqwva(R1(:,I),dt,d,index_incre,lvl,clip,line_color,face_color,mode,trace_balance);