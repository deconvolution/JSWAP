"
3D Eikonal equation solver
"
module eikonal
export acoustic_eikonal_forward
## Using ParallelStencil
include("../ParallelStencil/ParallelStencil.jl");
## utilities
include("..//utilities/utilities.jl");
## using packages
using Random,MAT,Plots,Dates,TimerOutputs,WriteVTK,ProgressMeter,DataFrames,CSV
#.ParallelStencil,.ParallelStencil.FiniteDifferences3D
## Use CPU for ParallelStencil
#=
const USE_GPU=false
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3);
else
    @init_parallel_stencil(Threads, Float64, 3);
end
=#
## acoustic forward
include("./acoustic_eikonal_forward.jl");
end
