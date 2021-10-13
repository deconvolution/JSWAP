"
3D solver on CPU, including isotropic, anisotropic, forward and adjoint solvers.
"
module CPU_3D
export isotropic_forward_solver,isotropic_adjoint_solver,PML_configuration,ParallelStencil
## Using ParallelStencil
include("../ParallelStencil/ParallelStencil.jl");
## utilities
include("..//utilities/utilities.jl");
## using packages
using Random,MAT,Plots,Dates,TimerOutputs,WriteVTK,ProgressMeter,DataFrames,CSV,
.ParallelStencil,.ParallelStencil.FiniteDifferences3D
## Use CPU for ParallelStencil
const USE_GPU=false
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3);
else
    @init_parallel_stencil(Threads, Float64, 3);
end
## forward isotropic
include("./JSWAP_CPU_3D_forward.jl");
end
