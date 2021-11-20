"
d/dx, 12-point, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function Dx_12(in,out,xs,xs2,ys,ys2,zs,zs2)
    @pick(out,xs,xs2,ys,ys2,zs,zs2)=@dx_12(in);
    return nothing
end

"
d/dy, 12-point, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function Dy_12(in,out,xs,xs2,ys,ys2,zs,zs2)
    @pick(out,xs,xs2,ys,ys2,zs,zs2)=@dy_12(in);
    return nothing
end

"
d/dz, 12-point, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function Dz_12(in,out,xs,xs2,ys,ys2,zs,zs2)
    @pick(out,xs,xs2,ys,ys2,zs,zs2)=@dz_12(in);
    return nothing
end
