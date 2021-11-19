using ParallelStencil,ParallelStencil.FiniteDifferences2D
const USE_GPU=false;
@static if USE_GPU
    @init_parallel_stencil(CUDA,Float64,3);
else
    @init_parallel_stencil(Threads,Float64,3);
end
##
n=5;
a=@zeros(n,9);
a[1:n,5]=1:n;
b=@zeros(n+4,9);
xs=2;
ys=0;
xs2=2;
ys2=0;
@parallel function S(in,out,xs,ys,xs2,ys2)
    @pick(out,xs,ys,xs2,ys2)=@all(in);
    return nothing
end
@parallel S(a,b,xs,ys,xs2,ys2);
display(b)
