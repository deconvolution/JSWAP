"
d/dx, 12-point, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function Dx_inn(in,out)
    @inn(out)=@d_xi_12(in);
    return nothing
end
"
d/dy, 12-point, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function Dy_inn(in,out)
    @inn(out)=@d_yi_12(in);
    return nothing
end
"
d/dz, 12-point, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function Dz_inn(in,out)
    @inn(out)=@d_zi_12(in);
    return nothing
end
"
shift grid in x direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[6:end-4,iy,iz]=in[:,iy,iz]
"
@parallel_indices (iy,iz) function u_1_plus(in,out)
out[6:end-4,iy,iz]=in[:,iy,iz];
return nothing
end

"
shift grid in x direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[5:end-5,iy,iz]=in[:,iy,iz]
"
@parallel_indices (iy,iz) function u_1_minus(in,out)
out[5:end-5,iy,iz]=in[:,iy,iz];
return nothing
end

"
shift grid in y direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[ix,6:end-4,iz]=in[ix,:,iz]
"
@parallel_indices (ix,iz) function u_2_plus(in,out)
out[ix,6:end-4,iz]=in[ix,:,iz];
return nothing
end

"
shift grid in y direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[ix,5:end-5,iz]=in[ix,:,iz]
"
@parallel_indices (ix,iz) function u_2_minus(in,out)
out[ix,5:end-5,iz]=in[ix,:,iz];
return nothing
end

"
shift grid in z direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[ix,iy,6:end-4]=in[ix,iy,:]
"
@parallel_indices (ix,iy) function u_3_plus(in,out)
out[ix,iy,6:end-4]=in[ix,iy,:];
return nothing
end

"
shift grid in z direction, subfunction of JSWAP_CPU_3D_isotropic_solver.

out[ix,iy,5:end-5]=in[ix,iy,:]
"
@parallel_indices (ix,iy) function u_3_minus(in,out)
out[ix,iy,5:end-5]=in[ix,iy,:];
return nothing
end
##
#=
Forward: use minus
Backward: use plus
=#
