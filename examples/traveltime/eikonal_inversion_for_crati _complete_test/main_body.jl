## Initialization of JSWAP
include("./JSWAP_initialization_eikonal.jl");
## create process
path_to_input=["./input_template_start.jl"];
## post processing
include("./post_processing_acoustic_eikonal_inversion.jl");
