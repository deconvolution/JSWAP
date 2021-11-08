"
Computes v, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function JSWAP_CPU_3D_isotropic_forward_solver_compute_v_curvilinear(dt,dx,dy,dz,rho,beta,
    v1,v2,v3,
    sigmas11_1_minus,
    sigmas22_2_minus,
    sigmas33_3_minus,
    sigmas23_2_plus,sigmas23_3_plus,
    sigmas13_1_plus,sigmas13_3_plus,
    sigmas12_1_plus,sigmas12_2_plus,
    p_1_minus,p_2_minus,p_3_minus,
    sigmas11_3_cur,sigmas12_3_cur,sigmas13_3_cur,sigmas22_3_cur,sigmas23_3_cur,p_3_cur,
    Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax)

    @all(v1)=dt./@all(rho) .*((@all(sigmas11_1_minus)-@all(p_1_minus))/dx-
    @all(Z_Kmax) .* @all(Kmax_x) .*(@all(sigmas11_3_cur)-@all(p_3_cur))/dz+
    @all(sigmas12_2_plus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(sigmas12_3_cur)/dz+
    @all(Zmax_Kmax) .*@all(sigmas13_3_plus)/dz)+
    @all(v1)-
    dt*@all(beta) .*@all(v1);

    @all(v2)=dt./@all(rho) .*(@all(sigmas12_1_plus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(sigmas12_3_cur)/dz+
    (@all(sigmas22_2_minus)-@all(p_2_minus))/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*(@all(sigmas22_3_cur)-@all(p_3_cur))/dz+
    @all(Zmax_Kmax) .*@all(sigmas23_3_plus)/dz)+
    @all(v2)-
    dt*@all(beta) .*@all(v2);

    @all(v3)=dt./@all(rho) .*(@all(sigmas13_1_plus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(sigmas13_3_cur)/dz+
    @all(sigmas23_2_plus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(sigmas23_3_cur)/dz+
    @all(Zmax_Kmax) .*(@all(sigmas33_3_minus)-@all(p_3_minus))/dz)+
    @all(v3)-
    dt*@all(beta) .*@all(v3);

    return nothing
end
