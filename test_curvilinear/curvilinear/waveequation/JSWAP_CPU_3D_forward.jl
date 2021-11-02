## timing
ti=TimerOutput();
## finite-difference method package
include("./finite_difference_method.jl");
## tri
## function configure_PML
include("./tri_PML_configuration.jl");
## function JSWAP_CPU_3D_tri_forward_solver and its dependencies
# function JSWAP_CPU_3D_tri_forward_solver_compute_au_for_sigma
include("./JSWAP_CPU_3D_tri_forward_solver_compute_au_for_sigma.jl");
# function JSWAP_CPU_3D_tri_forward_solver_compute_sigma
include("./JSWAP_CPU_3D_tri_forward_solver_compute_sigma.jl");
# function JSWAP_CPU_3D_tri_forward_solver_compute_v
include("./JSWAP_CPU_3D_tri_forward_solver_compute_v.jl");
# function JSWAP_CPU_3D_tri_forward_solver
include("./JSWAP_CPU_3D_tri_forward_solver.jl");
## function JSWAP_CPU_3D_tri_adjoint_solver and its dependencies
# function JSWAP_CPU_3D_tri_forward_solver_compute_au_for_sigma
include("./JSWAP_CPU_3D_tri_forward_solver_compute_au_for_sigma.jl");
# function JSWAP_CPU_3D_tri_forward_solver_compute_sigma
include("./JSWAP_CPU_3D_tri_forward_solver_compute_sigma.jl");
# function JSWAP_CPU_3D_tri_forward_solver_compute_v
include("./JSWAP_CPU_3D_tri_forward_solver_compute_v.jl");
include("./JSWAP_CPU_3D_forward_tri_solver.jl");
## isotropic
## function configure_PML
include("./isotropic_PML_configuration.jl");
## function JSWAP_CPU_3D_isotropic_forward_solver and its dependencies
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_v
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_v.jl");
## function JSWAP_CPU_3D_isotropic_forward_curvilinear_solver and its dependencies
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma_curvilinear.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma_curvilinear.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_v
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_v_curvilinear.jl");
include("./JSWAP_CPU_3D_forward_isotropic_curvilinear_solver.jl");
## function JSWAP_CPU_3D_isotropic_adjoint_solver and its dependencies
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma.jl");
# function JSWAP_CPU_3D_isotropic_forward_solver_compute_v
include("./JSWAP_CPU_3D_isotropic_forward_solver_compute_v.jl");
include("./JSWAP_CPU_3D_forward_isotropic_solver.jl");
## function JSWAP_CPU_3D_forward_solver
include("./JSWAP_CPU_3D_forward_solver.jl");
