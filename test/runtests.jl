using Test,JSWAP,JSWAP.LinearAlgebra

include("./test_JSWAP_CPU_3D_isotropic_forward_solver.jl");

@test pass_test_forward
