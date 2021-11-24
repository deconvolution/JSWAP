## Initialization of JSWAP
include("./JSWAP_initialization_eikonal.jl");
## Specify where the input file is
path_to_input=["input_template.jl"];
for I=1:length(path_to_input)
    include(path_to_input[I]);
    T,T_rec=JSWAP.eikonal.acoustic_eikonal_forward_implement(input2);
end
