"
Compuites auxiliary variable used to compute sigma, subfunction of JSWAP_CPU_3D_isotropic_solver.
"

@parallel function JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma_bar_curvilinear(dt,dx,dy,dz,inv_Qa,lambda,mu,
    beta,
    v1_1_plus,v1_2_minus,v1_3_minus,
    v2_1_minus,v2_2_plus,v2_3_minus,
    v3_1_minus,v3_2_minus,v3_3_plus,
    sigmas11_bar,sigmas22_bar,sigmas33_bar,sigmas23_bar,sigmas13_bar,sigmas12_bar,p_bar,
    ax_bar,ax2_bar,ax3_bar,ax4_bar,ax5_bar,ax6_bar,ax7_bar,
    Ax_bar,Ax2_bar,Ax3_bar,Ax4_bar,Ax5_bar,Ax6_bar,Ax7_bar,
    ax_dt_bar,ax2_dt_bar,ax3_dt_bar,ax4_dt_bar,ax5_dt_bar,ax6_dt_bar,ax7_dt_bar,
    Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax)

    @all(Ax_bar)=@all(ax_bar);
    @all(Ax2_bar)=@all(ax2_bar);
    @all(Ax3_bar)=@all(ax3_bar);
    @all(Ax4_bar)=@all(ax4_bar);
    @all(Ax5_bar)=@all(ax5_bar);
    @all(Ax6_bar)=@all(ax6_bar);
    @all(Ax7_bar)=@all(ax7_bar);

    @all(ax_bar)=4*@all(mu) .*(@all(v1_1_plus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_minus)/dz)+
    (-2*@all(mu)) .*(@all(v2_2_plus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_minus)/dz)+
    (-2*@all(mu)) .*@all(Zmax_Kmax) .*@all(v3_3_plus)/dz;

    @all(ax2_bar)=(-2*@all(mu)) .*(@all(v1_1_plus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_minus)/dz)+
    (4*@all(mu)) .* (@all(v2_2_plus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_minus)/dz)+
    (-2*@all(mu)) .* @all(Zmax_Kmax) .*@all(v3_3_plus)/dz;

    @all(ax3_bar)=(-2*@all(mu)) .*(@all(v1_1_plus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_minus)/dz)+
    (-2*@all(mu)) .*(@all(v2_2_plus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_minus)/dz)+
    (4*@all(mu)) .*@all(Zmax_Kmax) .*@all(v3_3_plus)/dz;

    @all(ax4_bar)=@all(mu).*(@all(v2_1_minus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v2_3_minus)/dz)+
    @all(mu).*(@all(v1_2_minus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v1_3_minus)/dz);

    @all(ax5_bar)=@all(mu) .*(@all(v3_1_minus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v3_3_plus)/dz)+
    @all(mu) .*@all(Zmax_Kmax).*@all(v1_3_minus)/dz;

    @all(ax6_bar)=@all(mu).*(@all(v3_2_minus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v3_3_plus)/dz)+
    @all(mu) .*@all(Zmax_Kmax) .*@all(v2_3_minus)/dz;

    @all(ax7_bar)=(3*@all(lambda)+2*@all(mu)) .*(@all(v1_1_plus)/dx-
    @all(Z_Kmax) .*@all(Kmax_x) .*@all(v1_3_minus)/dz)+
    (3*@all(lambda)+2*@all(mu)) .*(@all(v2_2_plus)/dy-
    @all(Z_Kmax) .*@all(Kmax_y) .*@all(v2_3_minus)/dz)+
    (3*@all(lambda)+2*@all(mu)) .*@all(Zmax_Kmax) .*@all(v3_3_plus)/dz;

    @all(ax_dt_bar)=(@all(ax_bar)-@all(Ax_bar))/dt;
    @all(ax2_dt_bar)=(@all(ax2_bar)-@all(Ax2_bar))/dt;
    @all(ax3_dt_bar)=(@all(ax3_bar)-@all(Ax3_bar))/dt;
    @all(ax4_dt_bar)=(@all(ax4_bar)-@all(Ax4_bar))/dt;
    @all(ax5_dt_bar)=(@all(ax5_bar)-@all(Ax5_bar))/dt;
    @all(ax6_dt_bar)=(@all(ax6_bar)-@all(Ax6_bar))/dt;
    @all(ax7_dt_bar)=(@all(ax7_bar)-@all(Ax7_bar))/dt;

    return nothing
end
