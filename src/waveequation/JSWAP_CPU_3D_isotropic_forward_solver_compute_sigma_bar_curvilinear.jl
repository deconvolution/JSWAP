"
Computes sigma, subfunction of JSWAP_CPU_3D_isotropic_solver.
"

@parallel function JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma_bar_curvilinear(dt,dx,dy,dz,inv_Qa,lambda,mu,
    beta,
    v1_1_plus,v1_2_minus,v1_3_minus,
    v2_1_minus,v2_2_plus,v2_3_minus,
    v3_1_minus,v3_2_minus,v3_3_plus,
    sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p,
    sigmas11_bar,sigmas22_bar,sigmas33_bar,sigmas23_bar,sigmas13_bar,sigmas12_bar,p_bar,
    ax_bar,ax2_bar,ax3_bar,ax4_bar,ax5_bar,ax6_bar,ax7_bar,
    Ax_bar,Ax2_bar,Ax3_bar,Ax4_bar,Ax5_bar,Ax6_bar,Ax7_bar,
    ax_dt_bar,ax2_dt_bar,ax3_dt_bar,ax4_dt_bar,ax5_dt_bar,ax6_dt_bar,ax7_dt_bar,
    Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax)

    @all(sigmas11_bar)=1/3*dt*(
    @all(ax_bar)+@all(inv_Qa) .*@all(ax_dt_bar))+
    @all(sigmas11)-
    dt*@all(beta).*@all(sigmas11);

    @all(sigmas22_bar)=1/3*dt*(
    @all(ax2_bar)+@all(inv_Qa) .*@all(ax2_dt_bar))+
    @all(sigmas22)-
    dt*@all(beta).*@all(sigmas22);

    @all(sigmas33_bar)=1/3*dt*(
    @all(ax3_bar)+@all(inv_Qa) .*@all(ax3_dt_bar))+
    @all(sigmas33)-
    dt*@all(beta).*@all(sigmas33);

    @all(sigmas12_bar)=dt*(
    @all(ax4_bar)+@all(inv_Qa) .*@all(ax4_dt_bar))+
    @all(sigmas12)-
    dt*@all(beta).*@all(sigmas12);

    @all(sigmas13_bar)=dt*(
    @all(ax5_bar)+@all(inv_Qa) .*@all(ax5_dt_bar))+
    @all(sigmas13)-
    dt*@all(beta).*@all(sigmas13);

    @all(sigmas23_bar)=dt*(
    @all(ax6_bar)+@all(inv_Qa) .*@all(ax6_dt_bar))+
    @all(sigmas23)-
    dt*@all(beta).*@all(sigmas23);

    @all(p_bar)=-1/3*dt*(
    @all(ax7_bar)+@all(inv_Qa) .*@all(ax7_dt_bar))+
    @all(p)-
    dt*@all(beta).*@all(p);

    return nothing
end
