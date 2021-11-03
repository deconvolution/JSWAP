"
Compuites auxiliary variable used to compute sigma, subfunction of JSWAP_CPU_3D_isotropic_solver.
"

@parallel function JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma_curvilinear(dt,dx,dy,dz,inv_Qa,lambda,mu,
    beta,
    v1_1_plus,v1_2_minus,v1_3_minus,
    v2_1_minus,v2_2_plus,v2_3_minus,
    v3_1_minus,v3_2_minus,v3_3_plus,
    sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p,
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

    @all(ax)=4*@all(mu) .*(@all(v1_1_plus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_minus)/dz)+
    (-2*@all(mu)) .*(@all(v2_2_plus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_minus)/dz)+
    (-2*@all(mu)) .*@all(Zmax_Kmax) .*@all(v3_3_plus)/dz;

    @all(ax2)=(-2*@all(mu)) .*(@all(v1_1_plus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_minus)/dz)+
    (4*@all(mu)) .* (@all(v2_2_plus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_minus)/dz)+
    (-2*@all(mu)) .* @all(Zmax_Kmax) .*@all(v3_3_plus)/dz;

    @all(ax3)=(-2*@all(mu)) .*(@all(v1_1_plus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_minus)/dz)+
    (-2*@all(mu)) .*(@all(v2_2_plus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_minus)/dz)+
    (4*@all(mu)) .*@all(Zmax_Kmax) .*@all(v3_3_plus)/dz;

    @all(ax4)=@all(mu).*(@all(v2_1_minus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v2_3_minus)/dz)+
    @all(mu).*(@all(v1_2_minus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v1_3_minus)/dz);

    @all(ax5)=@all(mu) .*(@all(v3_1_minus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v3_3_plus)/dz)+
    @all(mu) .*@all(Zmax_Kmax).*@all(v1_3_minus)/dz;

    @all(ax6)=@all(mu).*(@all(v3_2_minus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v3_3_plus)/dz)+
    @all(mu) .*@all(Zmax_Kmax) .*@all(v2_3_minus)/dz;

    @all(ax7)=(3*@all(lambda)+2*@all(mu)) .*(@all(v1_1_plus)/dx-
    0*@all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_minus)/dz)+
    (3*@all(lambda)+2*@all(mu)) .*(@all(v2_2_plus)/dy-
    0*@all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_minus)/dz)+
    (3*@all(lambda)+2*@all(mu)) .*@all(Zmax_Kmax) .*@all(v3_3_plus)/dz;

    @all(ax_dt)=(@all(ax)-@all(Ax))/dt;
    @all(ax2_dt)=(@all(ax2)-@all(Ax2))/dt;
    @all(ax3_dt)=(@all(ax3)-@all(Ax3))/dt;
    @all(ax4_dt)=(@all(ax4)-@all(Ax4))/dt;
    @all(ax5_dt)=(@all(ax5)-@all(Ax5))/dt;
    @all(ax6_dt)=(@all(ax6)-@all(Ax6))/dt;
    @all(ax7_dt)=(@all(ax7)-@all(Ax7))/dt;

    return nothing
end