In JSWAP, all of the input parameters are assigned to a struct input2. The following input parameters are needed to run the simulation:
# Overview of inputs #
* Dimensions
 * [dx,dy,dz](#dx,dy,dz)
 * [dt](#dt)
 * [nx,ny,nz](#nx,ny,nz)
 * [X,Y,Z](#X,Y,Z)
* Material properties
 * [lambda](#lambda)
 * [mu](#mu)
 * [rho](#rho)
 * [inv_Qa](#inv_Qa)
* Receivers
 * [r1,r2,r3](#r1,r2,r3)
 * [r1t,r2t,r3t](#r1t,r2t,r3t)
 * [Rm](#Rm)
* Sources
 * [s1,s2,s3](#s1,s2,s3)
 * [src1,src2,src3,srcp](#src1,src2,src3,srcp)
 * [s1t,s2t,s3t](#s1t,s2t,s3t)
* PML
 * [lp](#lp)
 * [nPML](#nPML)
 * [Rc](#Rc)
 * [PML_active](#nPML_active)
* Storage
 * [path](#path)
 * [path_pic](#path_pic)
 * [path_model](#path_model)
 * [path_wavefield](#path_wavefield)
 * [path_recordings](#path_rec)
* Interval
 * [plot_interval](#plot_interval)
 * [wavefield_interval](#wavefield_interval)

## dx,dy,dz
Grid spacing in the x, y and z directions.

Example:
```julia
# dx
julia> input2.dx=10.0;
# dy
julia> input2.dy=10.0;
# dz
julia> input2.dz=10.0;
```
## dt
Time interval.

Example:
```julia
# dx
julia> input2.dt=10.0^-3;
```
## nx, ny, nz
Grid points in the x, y and z directions.
```julia
# nx
julia> input2.nx=80;
# ny
julia> input2.ny=80;
# nz
julia> input2.nz=90;
```
## nt
Time interval.

Example:
```julia
# dx
julia> input2.nt=10^3;
```
## X,Y,Z
3D true coordinate, 3D tensor.

Example:
```julia
# 3D true coordinate X, Y and Z
input2.Y,input2.X,input2.Z=JSWAP.meshgrid(1:80,1:80,1:90);
```
## lambda
Lame constant, lambda

Example:
```julia
input2.lambda=ones(80,80,90)*10^9*1.0;
```
## mu
Lame constant, mu

Example:
```julia
input2.mu=ones(80,80,90)*10^9*1.0;
```
## rho
Density.

Example:
```julia
input2.rho=ones(80,80,90)*1000.0;
```
## inv_Qa
Apparent attenuation
Example:
```julia
input2.inv_Qa=ones(80,80,90)*0.0;
```

to be continued

An example of input file is given [here](https://github.com/deconvolution/JSWAP/blob/main/examples/template/input_template.jl).
