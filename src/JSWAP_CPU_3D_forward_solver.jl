"
Finite-difference solver for isotropic viscoelastic media.
"
@timeit ti "iso_3D" function forward_solver(input2)

if isdefined(input2,:lambda)
    v1,v2,v3,R1,R2,R3,P=JSWAP_CPU_3D_forward_isotropic_solver(input2);
else
    v1,v2,v3,R1,R2,R3,P=JSWAP_CPU_3D_forward_tri_solver(input2);
end


return v1,v2,v3,R1,R2,R3,P
end
