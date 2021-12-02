close all;
clear all;

tt=load('./inversion_progress/output_vp_50.mat');
v=tt.data.vp;
X=tt.data.X;
Y=tt.data.Y;
Z=tt.data.Z;

[nx,ny,nz]=size(v);
n=30;

mv=ones(nx, ny, nz);

for i=1:floor(nx/n)
    for j=1:floor(ny/n)
        for k=1:floor(nz/n)
            if mod(i + j + k, 2)==0
                tt=1.2;
            else
                tt=.8;
            end
            mv(((1:n)+(i-1)*n), ((1:n)+(j-1)*n), ((1:n)+(k-1)*n))=tt;
        end
    end
end

vt=6800*mv;



figure(1)
for i=5
    set(gcf,'position',[80,80,1200,600]);
    subplot(1,2,1)
    imagesc(vt(:,:,i));
    title(i)
    colorbar;
    subplot(1,2,2)
    imagesc(v(:,:,i));
    colorbar;
    drawnow;
end
%%
cr=ones(nx,ny,nz)*nan;
tt=abs(v-vt)./vt;

for i=1:floor(nx/n)
    for j=1:floor(ny/n)
        for k=1:floor(nz/n)
            
            if min(tt(((1:n)+(i-1)*n), ((1:n)+(j-1)*n), ((1:n)+(k-1)*n)),[],'all')<.1
                cr(((1:n)+(i-1)*n), ((1:n)+(j-1)*n), ((1:n)+(k-1)*n))=1;
            end
        end
    end
end

figure(2)
for i=1:floor(n/2):size(v,3)
    set(gcf,'position',[80,80,1200,600]);
    subplot(1,2,1)
    tt=cr(:,:,i).*v(:,:,i);
    imAlpha=ones(size(tt));
    imAlpha(isnan(tt))=0;
    imagesc(tt,'AlphaData',imAlpha);
    set(gca,'color',0*[1 1 1]);
    title(i);
    colorbar;
    drawnow;
    pause(5)
end
%%
data=cr;
save('./checkerboard_effective.mat','data');
%%
data(isnan(data))=0;
vtkwrite('./checkerboard_effective.vtk','structured_grid',X,Y,Z,'scalars','cr',data);