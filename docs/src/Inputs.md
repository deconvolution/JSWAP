# Inputs
In JSWAP, all of the input parameters are assigned to a struct input2. The following input parameters are needed to run the simulation. Please go through one by one to check the input is correct.

# 1. dt
Explanation: Time interval.

Type: Float64.

Dimension: [].

Example:

```
input2.dt=10.0^-3;
```

# 2. dx dy dz
Explanation: Grid spacing in the x, y and z directions.

Type: Float64.

Dimension: [].

Example:

```
input2.dx=10.0;
input2.dy=10.0;
input2.dz=10.0;
```

# 3. nt
Explanation: Time interval.

Type: Float64.

Dimension: [].

Example:

```
input2.nt=10.0^3;
```

# 4. nx, ny, nz
Explanation: Grid points in the x, y and z directions.

Type: Int32.

Dimension: [].

```
input2.nx=80;
input2.ny=80;
input2.nz=90;
```

# 5. X,Y,Z
Explanation: 3D true coordinate of the model.

Type: Float or Int.

Dimension: [nx, ny, nz]

Example:

```
input2.Y,input2.X,input2.Z=JSWAP.meshgrid(1:80,1:80,1:90);
```

# 6. lambda, mu
Explanation: Lame constants, &lambda; and &mu;.

Type: Float64.

Dimension: [nx, ny, nz]


# 7. rho
Explanation: Density.

Type: Float64.

Dimension: [nx, ny, nz].

Example:

```
input2.rho=ones(80,80,90)*1000.0;
```

# 8. inv_Qa
Explanation: Apparent attenuation.

Type: Float64.

Dimension: [nx, ny, nz].

Example:

```
input2.inv_Qa=ones(80,80,90)*0.0;
```

# 8. r1, r2, r3
Explanation: Grid locations of receivers in the x, y and z directions.

Type: Int32.

Dimension: Each of them is a 2-dimensional matrix, with the first direction to be 1 and the second direction as the number of receiver numbers.

Here is an example of 50 receivers.

```
input2.r1=zeros(Int32,1,50);
input2.r1[:]=30:79;

input2.r2=zeros(Int32,1,50);
input2.r2[:]=30:79;

input2.r3=zeros(Int32,1,50);
input2.r3 .=15;
```

# 9. s1, s2, s3
Explanation: Grid locations of sources in the x, y and z directions.

Type: Int32.

Dimension: Each of them is a 2-dimensional matrix, with the first direction to be 1 and the second direction as the number of receiver numbers.

Here is an example of 2 sources.

```
input2.s1=zeros(Int32,1,2);
input2.s1[:] =[60 60];

input2.s2=zeros(Int32,1,2);
input2.s2[:] =[49 51];

input2.s3=zeros(Int32,1,2);
input2.s3[:] =[49 51];
```

# 10. Source signals
Explanation: Source signal in each component. The source can be either directional source, P-wave source or moment-tensor source. There are 2 options:
* If one chooses directional source and P-wave source, the input is supposed to be input2.src1, input2.src2, input2.src3 and input2.srcp, corresponding to the source signals in the x, y, z directions and the P-wave components.
* If moment tensor source is wanted, M11, M22, M33, M23, M13 and M12 are desired.

Type: Float64.

Dimensions: [nt, number of sources].

# 11. r1t, r2t, r3t
Explanation: Receiver true locations in the x, y and z directions.

Type: Int32 or Float64.

Dimension: [1, number of receivers].

Example of computation true locations based on above information.

```
input2.r1t=input2.r1*input2.dx;
input2.r2t=input2.r2*input2.dy;
input2.r3t=input2.r3*input2.dz;
```

# 12. Rm
Explanation: One can mute some receivers with this option. 0 - mute, 1 - unmute.

Type: Int32 or Float64.

Dimension: [number of receivers, component].

Component: 1 - x component, 2 - y component, 3 - z component, 4 - pressure component.

If one does not want to mute any receivers, then

```
input2.Rm=ones(length(input2.r3),4);
```


# 13. s1t, s2t, s3t
Explanation: Source true locations in the x, y and z directions.

Type: Int32 or Float64.

Dimension: [1, number of sources].

Here is an example of 50 receivers.

```
input2.s1t=input2.s1*input2.dx;
input2.s2t=input2.s2*input2.dy;
input2.s3t=input2.s3*input2.dz;
```

# 14. lp
Explanation: PML layers.

Type: Int32.

Dimension: [].

For a 10 layer PML,

```
input2.lp=10;
```

#  15. nPML
Explanation: power of PML, normally 2.

Type: Int32.

Dimension: [].

Example:

```
input2.nPML=2;
```

# 16. Rc
Explanation: PML theoretical reflection coefficient.

Type: Float64.

Dimension: [].

Some common coefficients:

|PML layers (lp)|Rc     |
|---------------|-------|
|10             |10.0^-1|
|20             |10.0^-2|
|30             |10.0^-3|
|40             |10.0^-4|

Example:

```
input2.Rc=.1;
```

# 17. PML_active
Explanation: One can set if PML at one edge is working. Each element of the 1 by 6 matrix means if the PML is activated on the edge in an order of xminus, xplus, yminus, yplus, zminus, zplus. "0" is deactivated and "1" is activated.

Type: Float64.

Dimension: [1, 6].

Example for only zminus PML is deactivated:

```
input2.PML_active=[1 1 1 1 0 1];
```

# 18. path
Explanation: Path of master branch on where to store the wavefield vtk and wavefield mat file. This must be a parent path of path_pic, path_model and path_wavefield, path_rec.

Type: String.

Dimension: NA.

```
input2.path=p3;
```

# 19. path_pic
Explanation: Path to store wavefield snapshot vtk file. If one assign it "nothing", then no wavefield snapshot will be saved.

Type: String.

Dimension: NA.

```
input2.path_pic=string(input2.path,"/pic");
input2.path_pic=nothing;
```

# 20. path_model
Explanation: Path to store the model, including material parameters and source and receiver locations.

Type: String.

Dimension: NA.

```
input2.path_model=string(input2.path,"/model");
```
# 21. path_wavefield
Explanation: Path to store mat wavefield. No wavefield matfile will be stored if path_wavefield is "nothing".

Type: String.

Dimension: NA.

```
input2.path_wavefield=string(input2.path,"/wavefield");
```

# 22. path_rec
Explanation: Path to store recordings.

Type: String.

Dimension: NA.

```
input2.path_rec=string(input2.path,"/rec");
```

# 23. plot_interval
Explanation: Plot frequency/interval of vtk file. "0" for saving nothing.

Type: Int32.

Dimension: [].

```
input2.plot_interval=100;
```

# 24. wavefield_interval
Explanation: How frequent the wavefield is saved. "0" for saving nothing.

Type: Int32.

Dimension: [].

```
input2.wavefield_interval=0;
```

An example of input file is given [here](https://github.com/deconvolution/JSWAP/blob/main/examples/template/input_template.jl).
