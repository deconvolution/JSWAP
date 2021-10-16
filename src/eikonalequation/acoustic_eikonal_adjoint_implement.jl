function acoustic_eikonal_adjoint_implement(input2)
    lambda=acoustic_eikonal_adjoint(input2.nx,input2.ny,input2.nz,
    input2.h,input2.T,input2.r1,input2.r2,input2.r3,
    input2.s1,input2.s2,input2.s3,input2.R_cal,input2.R_true);
    return lambda
end
