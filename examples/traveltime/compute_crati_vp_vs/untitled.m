clear all;
tt=load('./output_vp_40.mat');
X=tt.data.X;
Y=tt.data.Y;
Z=tt.data.Z;
vp=tt.data.vp;
%%
tt=load('./output_vs_40.mat');
X=tt.data.X;
Y=tt.data.Y;
Z=tt.data.Z;
vs=tt.data.vp;
vp_vs=vp./vs;
%%
vtkwrite('./vp_vs.vtk','structured_grid',X,Y,Z,'scalars','vp_vs',vp_vs)