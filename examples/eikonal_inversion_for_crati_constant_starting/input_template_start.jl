## inversion paramemters
input2.n_iteration=10^3;
input2.max_gradient=400;
input2.fu=50;

input2.R_true=Vector{Vector{Float64}}();
input2.s1=Vector{Vector{Int32}}();
input2.s2=Vector{Vector{Int32}}();
input2.s3=Vector{Vector{Int32}}();
input2.r1=Vector{Vector{Int32}}();
input2.r2=Vector{Vector{Int32}}();
input2.r3=Vector{Vector{Int32}}();
input2.s1t=Vector{Vector{Float64}}();
input2.s2t=Vector{Vector{Float64}}();
input2.s3t=Vector{Vector{Float64}}();
input2.r1t=Vector{Vector{Float64}}();
input2.r2t=Vector{Vector{Float64}}();
input2.r3t=Vector{Vector{Float64}}();
#input2.path=Vector{String}();
##
tt=readmat("./m.mat","data");
input2.nx=round(Int32,tt["nx"]);
input2.ny=round(Int32,tt["ny"]);
input2.nz=round(Int32,tt["nz"]);
input2.X=tt["X"];
input2.Y=tt["Y"];
input2.Z=tt["Z"];
input2.h=tt["dx"];
input2.v=tt["vp"];
input2.v[:] .=5000;

input2.v[input2.v .==0] .=340;

input2.path="./source_1";
tt=readdir("./data/receiver/crati_traveltime_input/");

for I=1:size(tt,1)
    tt2=JSWAP.readmat(string("./data/receiver/crati_traveltime_input/",tt[I]),"data");
    input2.R_true=push!(input2.R_true,tt2["Rp"][:,4]);
    input2.s1=push!(input2.s1,round.(Int32,tt2["S"][:,1]));
    input2.s2=push!(input2.s2,round.(Int32,tt2["S"][:,2]));
    input2.s3=push!(input2.s3,round.(Int32,tt2["S"][:,3]));
    input2.r1=push!(input2.r1,round.(Int32,tt2["Rp"][:,1]));
    input2.r2=push!(input2.r2,round.(Int32,tt2["Rp"][:,2]));
    input2.r3=push!(input2.r3,round.(Int32,tt2["Rp"][:,3]));
    #=
    input2.s1t=push!(input2.s1t,);
    input2.s2t=push!(input2.s2t,);
    input2.s3t=push!(input2.s3t,);
    input2.r1t=push!(input2.r1t,);
    input2.r2t=push!(input2.r2t,);
    input2.r3t=push!(input2.r3t,);
    =#
end
## receiver and source configuration.
"
The type of r1,r2,r3,s1,s2 and s3 should not be changed.
"
# receiver true location x
input2.r1t=input2.r1*input2.h;
# receiver true location y
input2.r2t=input2.r2*input2.h;
# receiver true location z
input2.r3t=input2.r3*input2.h;
# source true location x
input2.s1t=input2.s1*input2.h;
# source true location y
input2.s2t=input2.s2*input2.h;
# source true location z
input2.s3t=input2.s3*input2.h;
## plot
# path for storage. This must be the master branch of the following pathes.
input2.write_rec=0;
