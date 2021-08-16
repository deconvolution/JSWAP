An example of the simulation implementation is given [here](https://github.com/deconvolution/JSWAP/tree/main/examples/template).

The main file is `main_body.jl`. Only one line is needed to change:
```
# Specify where the input file is
path_to_input=["./input_template.jl"];
```
This line is used to tell JSWAP where the input file is. Multiple input files can be added. For instance,
```
# Specify where the input file is
path_to_input=["./input_template.jl","./input_template2.jl",];
```

Please navigate to the same path with `main_body.jf` before running it. Otherwise, some folders might not be found.

The main file will go through the input file and the solver during the running, then outputs will be generated under the path specified in the input file. Only the main file is needed to run.
