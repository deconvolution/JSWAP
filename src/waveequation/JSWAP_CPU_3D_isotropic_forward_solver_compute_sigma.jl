"
Comp_i_j_kutes sigma, subfunction of JSWAp_i_j_k_Cp_i_j_kU_3D_isotrop_i_j_kic_solver.
"

@parallel function JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma_curvilinear(dt,dx,dy,dz,inv_Qa,lambda,mu,
    beta,
    sigmas11_i_j_k,
    sigmas22_i_j_k,
    sigmas33_i_j_k,
    sigmas23_i_jph_kph,
    sigmas13_iph_j_kph,
    sigmas12_iph_jph_k,
    p_i_j_k,
    ax,ax2,ax3,ax4,ax5,ax6,ax7,
    Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
    ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt,
    Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax)

    @all(sigmas11_i_j_k)=1/3*dt*(
    @all(ax)+@all(inv_Qa) .*@all(ax_dt))+
    @all(sigmas11_i_j_k)-
    dt*@all(beta).*@all(sigmas11_i_j_k);

    @all(sigmas22_i_j_k)=1/3*dt*(
    @all(ax2)+@all(inv_Qa) .*@all(ax2_dt))+
    @all(sigmas22_i_j_k)-
    dt*@all(beta).*@all(sigmas22_i_j_k);

    @all(sigmas33_i_j_k)=1/3*dt*(
    @all(ax3)+@all(inv_Qa) .*@all(ax3_dt))+
    @all(sigmas33_i_j_k)-
    dt*@all(beta).*@all(sigmas33_i_j_k);

    @all(sigmas12_iph_jph_k)=dt*(
    @all(ax4)+@all(inv_Qa) .*@all(ax4_dt))+
    @all(sigmas12_iph_jph_k)-
    dt*@all(beta).*@all(sigmas12_iph_jph_k);

    @all(sigmas13_iph_j_kph)=dt*(
    @all(ax5)+@all(inv_Qa) .*@all(ax5_dt))+
    @all(sigmas13_iph_j_kph)-
    dt*@all(beta).*@all(sigmas13_iph_j_kph);

    @all(sigmas23_i_jph_kph)=dt*(
    @all(ax6)+@all(inv_Qa) .*@all(ax6_dt))+
    @all(sigmas23_i_jph_kph)-
    dt*@all(beta).*@all(sigmas23_i_jph_kph);

    @all(p_i_j_k)=-1/3*dt*(
    @all(ax7)+@all(inv_Qa) .*@all(ax7_dt))+
    @all(p_i_j_k)-
    dt*@all(beta).*@all(p_i_j_k);

    return nothing
end
