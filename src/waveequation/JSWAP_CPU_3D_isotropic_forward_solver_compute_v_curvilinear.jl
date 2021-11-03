"
Computes v, subfunction of JSWAP_CPU_3D_isotropic_solver.
"

@parallel function JSWAP_CPU_3D_isotropic_forward_solver_compute_v_curvilinear(dt,dx,dy,dz,rho,beta,
    v1_bar,v2_bar,v3_bar,
    v1,v2,v3,
    sigmas11_bar_1_plus,sigmas11_bar_3_plus,
    sigmas22_bar_2_plus,sigmas22_bar_3_plus,
    sigmas33_bar_3_plus,
    sigmas23_bar_2_plus,sigmas23_bar_3_plus,
    sigmas13_bar_1_plus,sigmas13_bar_3_plus,
    sigmas12_bar_1_plus,sigmas12_bar_2_plus,sigmas12_bar_3_plus,
    p_bar_1_plus,p_bar_2_plus,p_bar_3_plus,
    Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax)

    @all(v1)=.5*(dt./@all(rho) .*((@all(sigmas11_bar_1_plus)-@all(p_bar_1_plus))/dx-
    @all(Z_Kmax) .* @all(Kmax_x) .*(@all(sigmas11_bar_3_plus)-@all(p_bar_3_plus))/dz+
    @all(sigmas12_bar_2_plus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(sigmas12_bar_3_plus)/dz+
    @all(Zmax_Kmax) .*@all(sigmas13_bar_3_plus)/dz)+
    @all(v1)+@all(v1_bar));

    @all(v2)=.5*(dt./@all(rho) .*(@all(sigmas12_bar_1_plus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(sigmas12_bar_3_plus)/dz+
    (@all(sigmas22_bar_2_plus)-@all(p_bar_2_plus))/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*(@all(sigmas22_bar_3_plus)-@all(p_bar_3_plus))/dz+
    @all(Zmax_Kmax) .*@all(sigmas23_bar_3_plus)/dz)+
    @all(v2)+@all(v2_bar));

    @all(v3)=.5*(dt./@all(rho) .*(@all(sigmas13_bar_1_plus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(sigmas13_bar_3_plus)/dz+
    @all(sigmas23_bar_2_plus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(sigmas23_bar_3_plus)/dz+
    @all(Zmax_Kmax) .*(@all(sigmas33_bar_3_plus)-@all(p_bar_3_plus))/dz)+
    @all(v3)+@all(v3_bar));

    return nothing
end
