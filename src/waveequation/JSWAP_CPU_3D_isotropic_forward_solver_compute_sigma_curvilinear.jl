"
Computes sigma, subfunction of JSWAP_CPU_3D_isotropic_solver.
"

@parallel function JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma_curvilinear(dt,dx,dy,dz,inv_Qa,lambda,mu,
    beta,
    v1_bar_1_plus,v1_bar_2_plus,v1_bar_3_plus,
    v2_bar_1_plus,v2_bar_2_plus,v2_bar_3_plus,
    v3_bar_1_plus,v3_bar_2_plus,v3_bar_3_plus,
    sigmas11_bar,sigmas22_bar,sigmas33_bar,sigmas13_bar,sigmas23_bar,sigmas12_bar,p_bar,
    sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p,
    ax,ax2,ax3,ax4,ax5,ax6,ax7,
    Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
    ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt,
    Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax)

    @all(sigmas11)=.5*(1/3*dt*(
    @all(ax)+@all(inv_Qa) .*@all(ax_dt))+
    @all(sigmas11)+@all(sigmas11_bar));

    @all(sigmas22)=.5*(1/3*dt*(
    @all(ax2)+@all(inv_Qa) .*@all(ax2_dt))+
    @all(sigmas11)+@all(sigmas22_bar));

    @all(sigmas33)=.5*(1/3*dt*(
    @all(ax3)+@all(inv_Qa) .*@all(ax3_dt))+
    @all(sigmas33)+@all(sigmas33_bar));

    @all(sigmas12)=.5*(dt*(
    @all(ax4)+@all(inv_Qa) .*@all(ax4_dt))+
    @all(sigmas12)+@all(sigmas12_bar));

    @all(sigmas13)=.5*(dt*(
    @all(ax5)+@all(inv_Qa) .*@all(ax5_dt))+
    @all(sigmas13)+@all(sigmas13_bar));

    @all(sigmas23)=.5*(dt*(
    @all(ax6)+@all(inv_Qa) .*@all(ax6_dt))+
    @all(sigmas23)+@all(sigmas23_bar));

    @all(p)=.5*(-1/3*dt*(
    @all(ax7)+@all(inv_Qa) .*@all(ax7_dt))+
    @all(p)+@all(p_bar));

    return nothing
end
