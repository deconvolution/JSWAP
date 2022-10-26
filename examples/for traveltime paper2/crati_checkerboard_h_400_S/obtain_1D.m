clear all;
close all;
tt=load("./inversion_progress/output_vp_1.mat");
vp=tt.data.vp;
vs=vp/1.7321;
dz=200;
nz=size(vp,3);

vp_1D=mean(vp,[1,2]);
vs_1D=mean(vs,[1,2]);
Z_1D=reshape(tt.data.Z(1,1,:),[length(vp_1D),1]);

vp_1D=reshape(vp_1D,[length(vp_1D),1]);
vs_1D=reshape(vs_1D,[length(vs_1D),1]);

tt=load("./m.mat");
vp_old=tt.data.vp;
vs_old=vp_old/1.7321;
dz=200;
nz=size(vp_old,3);

vp_1D_old=mean(vp_old,[1,2]);
vs_1D_old=mean(vs_old,[1,2]);
Z_1D_old=reshape(tt.data.Z(1,1,:),[length(vp_1D_old),1]);

vp_1D_old=reshape(vp_1D_old,[length(vp_1D),1]);
vs_1D_old=reshape(vs_1D_old,[length(vs_1D),1]);
figure;
subplot(1,2,1)
plot(vp_1D,Z_1D);
set(gca,'ydir','reverse');
hold on;
plot(vp_1D_old,Z_1D);

subplot(1,2,2)
plot(vs_1D,Z_1D);
hold on;
plot(vs_1D_old,Z_1D);
set(gca,'ydir','reverse');
%%
data=[Z_1D,vp_1D_old,vs_1D_old];
writematrix(data,"./new_1D_velocity.txt");