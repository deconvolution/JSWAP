# input template MSH Jacob Thesis
using DelimitedFiles, GeophysicalModelGenerator, JSWAP

## load data
dimensions      = DelimitedFiles.readdlm("./dimension_variables_for_input.csv",  ',', Float64, '\n',skipstart=0,header=false);
grids           = DelimitedFiles.readdlm("./grids_for_input.csv",   ',', Float64, '\n',skipstart=0,header=false);
lamda_unshaped  = DelimitedFiles.readdlm("./lamda_for_input.csv",   ',', Float64, '\n',skipstart=0,header=false);
mu_unshaped     = DelimitedFiles.readdlm("./mu_for_input.csv",   ',', Float64, '\n',skipstart=0,header=false);
rho_unshaped    = DelimitedFiles.readdlm("./rho_for_input.csv",   ',', Float64, '\n',skipstart=0,header=false);
sources_grid    = DelimitedFiles.readdlm("./sources_grid.csv", ',', Int, '\n',skipstart=0,header=false);
sources_mag     = DelimitedFiles.readdlm("./Sources_latlon+magn.csv", ',', Float64, '\n',skipstart=0,header=false);
receiver_grid   = DelimitedFiles.readdlm("./receiver_grid.csv", ',', Float64, '\n',skipstart=0,header=false);
## dimensions
# Time increment
input2.dt   = 0.01;
# dx
input2.dx   = dimensions[1,1];
# dy
input2.dy   = dimensions[1,2];
# dz
input2.dz   = dimensions[1,3];
# number of time steps
input2.nt   = 1000;
# nx
input2.nx   = dimensions[2,1]; input2.nx = trunc(Int32, input2.nx);
# ny
input2.ny   = dimensions[2,2]; input2.ny = trunc(Int32, input2.ny);
# nz
input2.nz   = dimensions[2,3]; input2.nz = trunc(Int32, input2.nz);
# 3D true coordinate X, Y and Z
input2.Y,input2.X,input2.Z  =   JSWAP.meshgrid(grids[2,:],grids[1,:],grids[3,:])
## material properties
input2.lambda   =   zeros(input2.nx, input2.ny, input2.nz);
for ix = 100:100:10000
    j = ix/input2.ny; j = trunc(Int, j);
    input2.lambda[:,:,j] = lamda_unshaped[:,ix-99:ix]
end

input2.mu       =   zeros(input2.nx, input2.ny, input2.nz);
for ix = 100:100:10000
    j = ix/input2.ny; j = trunc(Int, j);
    input2.mu[:,:,j] = mu_unshaped[:,ix-99:ix]
end
input2.rho      =   zeros(input2.nx, input2.ny, input2.nz);
for ix = 100:100:10000
    j = ix/input2.ny; j = trunc(Int, j);
    input2.rho[:,:,j] = rho_unshaped[:,ix-99:ix]
end

input2.inv_Qa   =   zeros(size(input2.X));
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver grid location x
input2.r1   =   zeros(Int32,1,53);
input2.r1[:]   =   receiver_grid[1,:];
for i = 1:length(input2.r1)
    if input2.r1[i] > 100
        input2.r1[i] = 100
    end
end
# receiver grid location y
input2.r2   =   zeros(Int32,1,53);
input2.r2[:]   =   receiver_grid[2,:];
# receiver grid location z
input2.r3   =   zeros(Int32,1,53);
input2.r3[:]   =   receiver_grid[3,:];
# source grid location x
input2.s1   =   zeros(Int32,1,1);
input2.s1[:].=   sources_grid[1,29];
# source grid location y
input2.s2   =   zeros(Int32,1,1);
input2.s2[:].=   sources_grid[2,29];
# source grid location z
input2.s3   =   zeros(Int32,1,1);
input2.s3[:].=   sources_grid[3,29];
# moment tensor source
frequency   =   1;
input2.M11  =   zeros(input2.nt,1);
input2.M22  =   zeros(input2.nt,1);
input2.M33  =   zeros(input2.nt,1);
input2.M23  =   zeros(input2.nt,1);
input2.M13  =   zeros(input2.nt,1);
input2.M12  =   zeros(input2.nt,1);
input2.M11[:] = rickerWave(frequency,input2.dt,input2.nt,1)*(-sources_mag[29,6]/sqrt(2));
input2.M22[:] = rickerWave(frequency,input2.dt,input2.nt,1)*(-sources_mag[29,6]/sqrt(2));
input2.M33[:] = rickerWave(frequency,input2.dt,input2.nt,1)*0;
input2.M23[:] = rickerWave(frequency,input2.dt,input2.nt,1)*0;
input2.M13[:] = rickerWave(frequency,input2.dt,input2.nt,1)*0;
input2.M12[:] = rickerWave(frequency,input2.dt,input2.nt,1)*0;
# receiver true location x
input2.r1t= minimum(grids[1]) .+ input2.r1*input2.dx;
# receiver true location y
input2.r2t= minimum(grids[2]) .+ input2.r2*input2.dy;
# receiver true location z
input2.r3t= maximum(grids[3]) .- input2.r3*input2.dz;
# activate receivers
input2.Rm=ones(length(input2.r3),4);
# source true location x
input2.s1t= minimum(grids[1]) .+ input2.s1*input2.dx;
# source true location y
input2.s2t= minimum(grids[2]) .+ input2.s2*input2.dy;
# source true location z
input2.s3t= maximum(grids[3]) .- input2.s3*input2.dz;
## PML
# PML layers
input2.lp=10;
# PML power
input2.nPML=2;
# PML theorecital coefficient
input2.Rc=.01;
# set PML active
# xminus,xplus,yminus,yplus,zminus,zplus
input2.PML_active=[1 1 1 1 1 1];
## plot
# path for storage. This must be the master branch of the following pathes.
input2.path=p3;
# path for wavefield .vtk
input2.path_pic=string(input2.path,"/pic");
# path for model
input2.path_model=string(input2.path,"/model");
# path for wavefield .mat
input2.path_wavefield=string(input2.path,"/wavefield");
# path for recordings
input2.path_rec=string(input2.path,"/rec");
# plot interval
input2.plot_interval=50;
# wavefield interval
input2.wavefield_interval=0;
