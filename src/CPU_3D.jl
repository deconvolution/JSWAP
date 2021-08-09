"
3D solver on CPU, including isotropic, anisotropic, forward and adjoint solvers.
"
module CPU_3D
export isotropic_forward_solver,isotropic_adjoint_solver
## Using ParallelStencil
include("./ParallelStencil/ParallelStencil.jl");
using Random,MAT,Plots,Dates,TimerOutputs,WriteVTK,ProgressMeter,DataFrames,CSV,
.ParallelStencil,.ParallelStencil.FiniteDifferences3D
## Use CPU for ParallelStencil
const USE_GPU=false
@static if USE_GPU
    @init_parallel_stencil(CUDA, Float64, 3);
else
    @init_parallel_stencil(Threads, Float64, 3);
end
## finite-difference method package
include("./finite_difference_method.jl");
## timing
ti=TimerOutput();
## function configure_PML
include("./configure_PML.jl");
## function JSWAP_CPU_3D_isotropic_forward_solver and its dependencies
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_v
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_v.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver
include("./JSWAP_CPU_3D_isotropic_forward_solver.jl");
## function JSWAP_CPU_3D_isotropic_adjoint_solver and its dependencies
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_v
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_v.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver
include("./JSWAP_CPU_3D_isotropic_forward_solver.jl");
end
