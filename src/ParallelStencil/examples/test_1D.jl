## This is to test the first derivative of sin
using MAT,PlotlyJS,LinearAlgebra
const USE_GPU=false;  # Use GPU? If this is set false, then no GPU needs to be available
using ParallelStencil
using ParallelStencil.FiniteDifferences1D
@static if USE_GPU
    @init_parallel_stencil(CUDA,Float64,3);
else
    @init_parallel_stencil(Threads,Float64,3);
end
dx=.1;
x=0:dx:pi;
nx=length(x);
y=cos.(x);
## analytical expression of d^2y/dx^2
dy_dx=-sin.(x);
## dy/dx
@parallel function Dx_1(in,out,dx,xs,xs2)
    @pick(out,xs,xs2)=@dx_1(in)/dx;
    return nothing
end

@parallel function Dx_8(in,out,dx,xs,xs2)
    @pick(out,xs,xs2)=@dx_8(in)/dx;
    return nothing
end

@parallel function Dx_12(in,out,dx,xs,xs2)
    @pick(out,xs,xs2)=@dx_12(in)/dx;
    return nothing
end
##
n=size(x,1);
dy_dx_1=@zeros(n,1);
dy_dx_8=@zeros(n,1);
dy_dx_12=@zeros(n,1);

@parallel Dx_1(y,dy_dx_1,dx,1,1)
@parallel Dx_8(y,dy_dx_8,dx,4,3)
@parallel Dx_12(y,dy_dx_12,dx,6,5)

plot(x[1:end-1],dy_dx_1)
plot(x[1:end-1],dy_dx_8)
plot(x[1:end-1],dy_dx_12)
