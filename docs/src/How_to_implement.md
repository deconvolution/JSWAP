An example of the simulation implementation is given [here](https://github.com/deconvolution/JSWAP/tree/main/examples/template).

The main file is `main_body.jl`. Only one line is needed to change:
```julia
# Specify where the input file is
path_to_input=["./input_template.jl"];
```
This line is used to tell JSWAP where the input file is. Multiple input files can be added. For instance,
```julia
# Specify where the input file is
path_to_input=["./input_template.jl","./input_template2.jl",];
```

The main file will go through the input file and the solver during the running, then outputs will be generated under the path specified in the input file. Only the main file is needed to run.
