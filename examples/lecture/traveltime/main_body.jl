using JSWAP
## create folder for saving
p2= @__FILE__;
p3=chop(p2,head=0,tail=3);
if isdir(p3)==0
    mkdir(p3);
end
##
include("./input_template.jl");
##
T,R_cal=JSWAP.eikonal.acoustic_eikonal_forward(
nx=nx,
ny=ny,
nz=nz,
h=h,
v=v,
s1=s1,
s2=s2,
s3=s3,
T0=0,
s1t=s1t,
s2t=s2t,
s3t=s3t,
r1=r1,
r2=r2,
r3=r3,
r1t=r1t,
r2t=r2t,
r3t=r3t,
X=X,
Y=Y,
Z=Z,
path=path,
write_t=write_t);
