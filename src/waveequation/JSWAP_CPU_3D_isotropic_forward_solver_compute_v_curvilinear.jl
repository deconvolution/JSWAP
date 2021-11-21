"
Computes v, subfunction of JSWAP_CPU_3D_isotropic_solver.
"
@parallel function JSWAP_CPU_3D_isotropic_forward_solver_compute_v_curvilinear(dt,dx,dy,dz,rho,beta,
    v1_iph_j_k,v2_i_jph_k,v3_i_j_kph,
    sigmas11_ip1_j_k_1,
    sigmas22_i_jp1_k_2,
    sigmas33_i_j_kp1_3,
    sigmas23_i_jph_kph_2,sigmas23_i_jph_kph_3,
    sigmas13_iph_j_k_1,sigmas13_iph_j_kph_3,
    sigmas12_iph_jph_k_1,sigmas12_iph_jph_k_2,
    p_ip1_j_k_1,p_i_jp1_k_2,p_i_j_kp1_3,
    sigmas11_3_v1,
    p_3_v1,
    sigmas12_3_v1,
    sigmas12_3_v2,
    sigmas22_3_v2,
    p_3_v2,
    sigmas13_3_v3,
    sigmas23_3_v3,
    Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax)

    @all(v1_iph_j_k)=dt./@all(rho) .*((@all(sigmas11_ip1_j_k_1)-@all(p_ip1_j_k_1))/dx-
    @all(Z_Kmax) .* @all(Kmax_x) .*(@all(sigmas11_3_v1)-@all(p_3_v1))/dz+
    @all(sigmas12_iph_jph_k_2)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(sigmas12_3_v1)/dz+
    @all(Zmax_Kmax) .*@all(sigmas13_iph_j_kph_3)/dz)+
    @all(v1_iph_j_k)-
    dt*@all(beta) .*@all(v1_iph_j_k);

    @all(v2_i_jph_k)=dt./@all(rho) .*(@all(sigmas12_iph_jph_k_1)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(sigmas12_3_v2)/dz+
    (@all(sigmas22_i_jp1_k_2)-@all(p_i_jp1_k_2))/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*(@all(sigmas22_3_v2)-@all(p_3_v2))/dz+
    @all(Zmax_Kmax) .*@all(sigmas23_i_jph_kph_3)/dz)+
    @all(v2_i_jph_k)-
    dt*@all(beta) .*@all(v2_i_jph_k);

    @all(v3_i_j_kph)=dt./@all(rho) .*(@all(sigmas13_iph_j_k_1)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(sigmas13_3_v3)/dz+
    @all(sigmas23_i_jph_kph_2)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(sigmas23_3_v3)/dz+
    @all(Zmax_Kmax) .*(@all(sigmas33_i_j_kp1_3)-@all(p_i_j_kp1_3))/dz)+
    @all(v3_i_j_kph)-
    dt*@all(beta) .*@all(v3_i_j_kph);

    return nothing
end
