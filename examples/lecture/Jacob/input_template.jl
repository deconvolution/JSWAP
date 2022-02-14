## load data
dimensions      = DelimitedFiles.readdlm("./dimension_variables_for_input.csv",  ',', Float64, '\n',skipstart=0,header=false);
grids           = DelimitedFiles.readdlm("./grids_for_input.csv",   ',', Float64, '\n',skipstart=0,header=false);
lamda_unshaped  = DelimitedFiles.readdlm("./lamda_for_input.csv",   ',', Float64, '\n',skipstart=0,header=false);
mu_unshaped     = DelimitedFiles.readdlm("./mu_for_input.csv",   ',', Float64, '\n',skipstart=0,header=false);
rho_unshaped    = DelimitedFiles.readdlm("./rho_for_input.csv",   ',', Float64, '\n',skipstart=0,header=false);
sources_grid    = floor.(Int64,DelimitedFiles.readdlm("./sources_grid.csv", ',', Float64, '\n',skipstart=0,header=false));
sources_mag     = DelimitedFiles.readdlm("./Sources_latlon+magn.csv", ',', Float64, '\n',skipstart=0,header=false);
receiver_grid   = DelimitedFiles.readdlm("./receiver_grid.csv", ',', Float64, '\n',skipstart=0,header=false);
## dimensions
# Time increment
dt=.01;
# dx
dx=dimensions[1,1];
# dy
dy=dimensions[1,2];
# dz
dz=abs(dimensions[1,3]);
# number of time steps
nt=4000;
# nx
nx=100;
# ny
ny=100;
# nz
nz=100;
# 3D true coordinate X, Y and Z
Y,X,Z=JSWAP.meshgrid((1:ny)*dy,(1:nx)*dx,(1:nz)*dz);
## material properties
lambda=zeros(nx,ny,nz);
mu=zeros(nx,ny,nz);
rho=zeros(nx,ny,nz);
lambda   =   zeros(nx, ny, nz);
for ix = 100:100:10000
    j = ix/ny; j = trunc(Int, j);
    lambda[:,:,j] = lamda_unshaped[:,ix-99:ix]
end

mu       =   zeros(nx, ny, nz);
for ix = 100:100:10000
    j = ix/ny; j = trunc(Int, j);
    mu[:,:,j] = mu_unshaped[:,ix-99:ix]
end
rho      =   zeros(nx, ny, nz);
for ix = 100:100:10000
    j = ix/ny; j = trunc(Int, j);
    rho[:,:,j] = rho_unshaped[:,ix-99:ix]
end
inv_Qa=ones(nx,ny,nz)*0.0;
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver grid location x
r1=zeros(Int64,1,53);
r1[:] .=receiver_grid[1,:];

# receiver grid location y
r2=zeros(Int64,1,53);
r2[:] .=receiver_grid[2,:];
# receiver grid location z
r3=zeros(Int64,1,53);
r3[:] .=receiver_grid[3,:];
# source grid location x
s1=zeros(Int64,1,1);
s1[:] .= 20;
# source grid location y
s2=zeros(Int64,1,1);
s2[:] .= 20;
# source grid location z
s3=zeros(Int64,1,1);
s3[:] .= 20;
# moment tensor source
M11=zeros(nt,1);
M22=zeros(nt,1);
M33=zeros(nt,1);
M23=zeros(nt,1);
M13=zeros(nt,1);
M12=zeros(nt,1);
freq=1.5;
M11[:]=-1*rickerWave(freq,dt,nt,2);
M22[:]=-1*rickerWave(freq,dt,nt,2);
M33[:]=0*rickerWave(freq,dt,nt,2);
M23[:]=0*rickerWave(freq,dt,nt,2);
M13[:]=0*rickerWave(freq,dt,nt,2);
M12[:]=0*rickerWave(freq,dt,nt,2);

# receiver true location x
r1t=r1*dx;
# receiver true location y
r2t=r2*dy;
# receiver true location z
r3t=r3*dz;
# source true location x
s1t=s1*dx;
# source true location y
s2t=s2*dy;
# source true location z
s3t=s3*dz;
## PML
# PML layers
lp=10;
# PML power
nPML=2;
# PML theorecital coefficient
Rc=.1;
# set PML active
# xminus,xplus,yminus,yplus,zminus,zplus
PML_active=[1 1 1 1 1 1];
## plot
# path for storage. This must be the master branch of the following pathes.
path=p3;
# plot interval
plot_interval=100;
# wavefield interval
wavefield_interval=nothing;
