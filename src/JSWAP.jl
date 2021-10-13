"
module JSWAP

Julia Sound WAve Propagation
"
module JSWAP
## utilities
include("./utilities//utilities.jl");
## finite-difference solver for isotropic media
include("./waveequation/CPU_3D.jl");
## Eikonal
include("./eikonalequation/eikonal.jl");
end
