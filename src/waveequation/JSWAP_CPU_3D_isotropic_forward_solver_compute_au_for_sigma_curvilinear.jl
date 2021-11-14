"
Compuites auxiliary variable used to compute sigma, subfunction of JSWAP_CPU_3D_isotropic_solver.
"

@parallel function JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma_curvilinear(dt,dx,dy,dz,inv_Qa,lambda,mu,
    beta,
    v1_iph_j_k_1,v1_iph_jp1_k_2,v1_ip1_jph_kph_3,
    v2_ip1_jph_k_1,v2_i_jph_k_2,v2_i_jph_kp1_3,
    v3_ip1_jph_kph_3,v3_i_jp1_kph_2,v3_i_j_kph_3,
    v1_3_sigmas11,
    v2_3_sigmas11,
    v3_3_sigmas23,
    v3_3_sigmas13,
    v2_3_sigmas12,
    v1_3_sigmas12,
    sigmas11_i_j_k,
    sigmas22_i_j_k,
    sigmas33_i_j_k,
    sigmas23_i_jph_kph,
    sigmas13_iph_j_kph,
    sigmas12_iph_jph_k,
    ax,ax2,ax3,ax4,ax5,ax6,ax7,
    Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
    ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt,
    Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax)

    @all(Ax)=@all(ax);
    @all(Ax2)=@all(ax2);
    @all(Ax3)=@all(ax3);
    @all(Ax4)=@all(ax4);
    @all(Ax5)=@all(ax5);
    @all(Ax6)=@all(ax6);
    @all(Ax7)=@all(ax7);

    # sigmas11
    @all(ax)=4*@all(mu) .*(@all(v1_iph_j_k_1)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_sigmas11)/dz)+
    (-2*@all(mu)) .*(@all(v2_i_jph_k_2)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_sigmas11)/dz)+
    (-2*@all(mu)) .*@all(Zmax_Kmax) .*@all(v3_i_j_kph_3)/dz;

    # sigmas22

    @all(ax2)=(-2*@all(mu)) .*(@all(v1_iph_j_k_1)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_sigmas11)/dz)+
    (4*@all(mu)) .* (@all(v2_i_jph_k_2)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_sigmas11)/dz)+
    (-2*@all(mu)) .* @all(Zmax_Kmax) .*@all(v3_i_j_kph_3)/dz;

    # sigmas33
    @all(ax3)=(-2*@all(mu)) .*(@all(v1_iph_j_k_1)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_sigmas11)/dz)+
    (-2*@all(mu)) .*(@all(v2_i_jph_k_2)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_sigmas11)/dz)+
    (4*@all(mu)) .*@all(Zmax_Kmax) .*@all(v3_i_j_kph_3)/dz;

    # sigmas12
    @all(ax4)=@all(mu).*(@all(v2_ip1_jph_k_1)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v2_3_sigmas12)/dz)+
    @all(mu).*(@all(v1_iph_jp1_k_2)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v1_3_sigmas12)/dz);

    # sigmas13
    @all(ax5)=@all(mu) .*(@all(v3_ip1_jph_kph_3)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v3_3_sigmas13)/dz)+
    @all(mu) .*@all(Zmax_Kmax).*@all(v1_ip1_jph_kph_3)/dz;

    # sigmas23
    @all(ax6)=@all(mu).*(@all(v3_i_jp1_kph_2)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v3_3_sigmas23)/dz)+
    @all(mu) .*@all(Zmax_Kmax) .*@all(v2_i_jph_kp1_3)/dz;

    # p
    @all(ax7)=(3*@all(lambda)+2*@all(mu)) .*(@all(v1_iph_j_k_1)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_sigmas11)/dz)+
    (3*@all(lambda)+2*@all(mu)) .*(@all(v2_i_jph_k_2)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_sigmas11)/dz)+
    (3*@all(lambda)+2*@all(mu)) .*@all(Zmax_Kmax) .*@all(v3_i_j_kph_3)/dz;

    @all(ax_dt)=(@all(ax)-@all(Ax))/dt;
    @all(ax2_dt)=(@all(ax2)-@all(Ax2))/dt;
    @all(ax3_dt)=(@all(ax3)-@all(Ax3))/dt;
    @all(ax4_dt)=(@all(ax4)-@all(Ax4))/dt;
    @all(ax5_dt)=(@all(ax5)-@all(Ax5))/dt;
    @all(ax6_dt)=(@all(ax6)-@all(Ax6))/dt;
    @all(ax7_dt)=(@all(ax7)-@all(Ax7))/dt;

    return nothing
end
