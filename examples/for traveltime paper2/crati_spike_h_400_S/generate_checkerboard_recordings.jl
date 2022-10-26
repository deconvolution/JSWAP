## read data
using JSWAP, Statistics
tt=JSWAP.readmat(string("./m.mat"), "data");
nx=round(Int64, tt["nx"]);
ny=round(Int64, tt["ny"]);
nz=round(Int64, tt["nz"]);
X=tt["X"];
Y=tt["Y"];
Z=tt["Z"];
h=tt["dx"];
## create checkerboard
w=ones(size(X));
w[62:93,75:118,end-26:end] .=.9;
v=3000*w;

vtkfile=JSWAP.vtk_grid("spike_model",X,Y,Z);
vtkfile["v"]=v;
JSWAP.vtk_save(vtkfile);

zero_Z=findmin(abs.(Z[1, 1, :]));

tt=readdir("./crati_traveltime_input/");

path="./crati_traveltime_checkerboard_input/";
if isdir(path)==0
    mkdir(path)
end

M=[0:30:size(tt, 1); size(tt, 1)];

for m=1:size(M,1)-1
    for I=(M[m]+1):M[m+1]
        tt2=JSWAP.readmat(string("./crati_traveltime_input/",tt[I]),"data");
        if size(tt2["Rs"],1)!=0
            nr=size(tt2["Rs"],1);
            s1=zeros(Int64,1,1);
            s2=zeros(Int64,1,1);
            s3=zeros(Int64,1,1);
            r1=zeros(Int64,1,nr);
            r2=zeros(Int64,1,nr);
            r3=zeros(Int64,1,nr);

            s1[:] .=tt2["S"][1];
            s2[:] .=tt2["S"][2];
            s3[:] .=tt2["S"][3];

            r1[:]=tt2["Rs"][:, 1];
            r2[:]=tt2["Rs"][:, 2];
            r3[:]=tt2["Rs"][:, 3];

            s1t=h*s1;
            s2t=h*s2;
            s3t=h*(zero_Z[2] .-s3);

            r1t=h*r1;
            r2t=h*r2;
            r3t=h*(zero_Z[2] .-r3);
            path="./source_1";
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
            write_t=0);

            data=copy(tt2);
            data["Rs"][:,4]=R_cal;
            file=JSWAP.matopen(string("./crati_traveltime_checkerboard_input/", tt[I]),"w");
            write(file, "data", data);
            close(file);
            global temp;
            temp=data;
        end
    end
end
