"
Computes sigma, subfunction of JSWAP_CPU_3D_isotropic_solver.
"

@parallel function JSWAP_CPU_3D_tri_forward_solver_compute_sigma(dt,dx,dy,dz,inv_Qa,
    beta,
    v1_1_plus,v1_2_minus,v1_3_minus,
    v2_1_minus,v2_2_plus,v2_3_minus,
    v3_1_minus,v3_2_minus,v3_3_plus,
    sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p,
    ax,ax2,ax3,ax4,ax5,ax6,ax7,
    Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
    ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt)

    @all(sigmas11)=1/3*dt*(
    @all(ax)+@all(inv_Qa) .*@all(ax_dt))+
    @all(sigmas11)-
    dt*@all(beta).*@all(sigmas11);

    @all(sigmas22)=1/3*dt*(
    @all(ax2)+@all(inv_Qa) .*@all(ax2_dt))+
    @all(sigmas22)-
    dt*@all(beta).*@all(sigmas22);

    @all(sigmas33)=1/3*dt*(
    @all(ax3)+@all(inv_Qa) .*@all(ax3_dt))+
    @all(sigmas33)-
    dt*@all(beta).*@all(sigmas33);

    @all(sigmas12)=dt*(
    @all(ax4)+@all(inv_Qa) .*@all(ax4_dt))+
    @all(sigmas12)-
    dt*@all(beta).*@all(sigmas12);

    @all(sigmas13)=dt*(
    @all(ax5)+@all(inv_Qa) .*@all(ax5_dt))+
    @all(sigmas13)-
    dt*@all(beta).*@all(sigmas13);

    @all(sigmas23)=dt*(
    @all(ax6)+@all(inv_Qa) .*@all(ax6_dt))+
    @all(sigmas23)-
    dt*@all(beta).*@all(sigmas23);

    @all(p)=-1/3*dt*(
    @all(ax7)+@all(inv_Qa) .*@all(ax7_dt))+
    @all(p)-
    dt*@all(beta).*@all(p);

    return nothing
end
