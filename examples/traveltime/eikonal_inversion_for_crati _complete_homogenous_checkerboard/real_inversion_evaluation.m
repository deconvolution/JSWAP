clear all;
close all;
tt=load('./output_vp_3.mat');
v=tt.data.vp;
X=tt.data.X;
Y=tt.data.Y;
Z=tt.data.Z;

[nx,ny,nz]=size(v);
n=30;

%% load checkerboard effective
tt=load('./checkerboard_effective.mat');
cr=tt.data;
%%
figure(2)
for i=80
    set(gcf,'position',[80,80,1200,600]);
    subplot(1,2,1)
    tt=cr(:,:,i);
    imAlpha=ones(size(tt));
    imAlpha(isnan(tt))=0;
    imagesc(tt,'AlphaData',imAlpha);
    set(gca,'color',0*[1 1 1]);
    title(Z(1,1,i));
    colorbar;
    subplot(1,2,2)
    tt=v(:,:,i).*cr(:,:,i);
    imAlpha=ones(size(tt));
    imAlpha(isnan(tt))=0;
    imagesc(tt,'AlphaData',imAlpha);
    set(gca,'color',0*[1 1 1]);
    colorbar;
    drawnow;
end
