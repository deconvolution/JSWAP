## import packages
const USE_GPU=false;  # Use GPU? If this is set false, then no GPU needs to be available
using ParallelStencil
using ParallelStencil.FiniteDifferences3D
@static if USE_GPU
    @init_parallel_stencil(CUDA,Float64,3);
else
    @init_parallel_stencil(Threads,Float64,3);
end
##
@parallel function Dx_12(in,out,xs,xs2,ys,ys2,zs,zs2)
    @pick(out,xs,xs2,ys,ys2,zs,zs2)=@dx_12(in);
    return nothing
end

@parallel function Dy_12(in,out,xs,xs2,ys,ys2,zs,zs2)
    @pick(out,xs,xs2,ys,ys2,zs,zs2)=@dy_12(in);
    return nothing
end

@parallel function Dz_12(in,out,xs,xs2,ys,ys2,zs,zs2)
    @pick(out,xs,xs2,ys,ys2,zs,zs2)=@dz_12(in);
    return nothing
end
##
A=@zeros(15,15,15);
A[8,8,8]=1;
A_x=@zeros(15,15,15);
A_y=@zeros(15,15,15);
A_z=@zeros(15,15,15);
B=copy(A);
@parallel Dx_12(A,A_x,6,5,0,0,0,0);
@parallel Dy_12(A,A_y,0,0,6,5,0,0);
@parallel Dz_12(A,A_z,0,0,0,0,6,5);

@parallel function Dx_1(in,out,xs,xs2,ys,ys2,zs,zs2)
    @pick(out,xs,xs2,ys,ys2,zs,zs2)=@dx_1(in);
    return nothing;
end

@parallel Dx_1(A,B,1,1,0,0,0,0);
