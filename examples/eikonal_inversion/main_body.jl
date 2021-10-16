## Initialization of JSWAP
include("./JSWAP_initialization_eikonal.jl");
## Specify where the input file is
path_to_input=["input_template_true.jl"];
Threads.@spawn for I=1:length(path_to_input)
    include(path_to_input[I]);
    T,T_rec=JSWAP.eikonal.acoustic_eikonal_forward_implement(input2);
end
## create process
input2.n_iteration=10;
input2.max_gradient=.5;
input2.fu=50;
path_to_input=["./input_template_start.jl"];
## post processing
include("./post_processing_acoustic_eikonal_inversion.jl");
