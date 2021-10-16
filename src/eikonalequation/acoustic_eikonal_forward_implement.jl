function acoustic_eikonal_forward_implement(input2)
    T,T_rec=acoustic_eikonal_forward(input2.nx,input2.ny,input2.nz,
    input2.h,input2.v,input2.s1,input2.s2,input2.s3,input2.s1t,
    input2.s2t,input2.s3t,
    input2.r1,input2.r2,input2.r3,input2.r1t,input2.r2t,input2.r3t,
    input2.X,input2.Y,input2.Z,input2.path,input2.write_rec);
    return T,T_rec
end
