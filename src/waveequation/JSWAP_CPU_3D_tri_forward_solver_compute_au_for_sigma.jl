"
Compuites auxiliary variable used to compute sigma, subfunction of JSWAP_CPU_3D_isotropic_solver.
"

@parallel function JSWAP_CPU_3D_tri_forward_solver_compute_au_for_sigma(dt,dx,dy,dz,inv_Qa,
    C11,C12,C13,C14,C15,C16,
    C22,C23,C24,C25,C26,
    C33,C34,C35,C36,
    C44,C45,C46,
    C55,C56,
    C66,
    beta,
    v1_1_plus,v1_2_minus,v1_3_minus,
    v2_1_minus,v2_2_plus,v2_3_minus,
    v3_1_minus,v3_2_minus,v3_3_plus,
    sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p,
    ax,ax2,ax3,ax4,ax5,ax6,ax7,
    Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
    ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt)

    @all(Ax)=@all(ax);
    @all(Ax2)=@all(ax2);
    @all(Ax3)=@all(ax3);
    @all(Ax4)=@all(ax4);
    @all(Ax5)=@all(ax5);
    @all(Ax6)=@all(ax6);
    @all(Ax7)=@all(ax7);

    @all(ax)=(2*@all(C11)-@all(C12)-@all(C13)) .*@all(v1_1_plus)/dx+
    (2*@all(C16)-@all(C26)-@all(C36)) .*@all(v1_2_minus)/dy+
    (2*@all(C15)-@all(C25)-@all(C35)) .*@all(v1_3_minus)/dz+
    (2*@all(C16)-@all(C26)-@all(C36)) .*@all(v2_1_minus)/dx+
    (2*@all(C12)-@all(C22)-@all(C23)) .*@all(v2_2_plus)/dy+
    (2*@all(C14)-@all(C24)-@all(C34)) .*@all(v2_3_minus)/dz+
    (2*@all(C15)-@all(C25)-@all(C35)) .*@all(v3_1_minus)/dx+
    (2*@all(C14)-@all(C24)-@all(C34)) .*@all(v3_2_minus)/dy+
    (2*@all(C13)-@all(C23)-@all(C33)) .*@all(v3_3_plus)/dz;

    @all(ax2)=(-@all(C11)+2*@all(C12)-@all(C13)) .*@all(v1_1_plus)/dx+
    (-@all(C16)+2*@all(C26)-@all(C36)) .*@all(v1_2_minus)/dy+
    (-@all(C15)+2*@all(C25)-@all(C35)) .*@all(v1_3_minus)/dz+
    (-@all(C16)+2*@all(C26)-@all(C36)) .*@all(v2_1_minus)/dx+
    (-@all(C12)+2*@all(C22)-@all(C23)) .*@all(v2_2_plus)/dy+
    (-@all(C14)+2*@all(C24)-@all(C34)) .*@all(v2_3_minus)/dz+
    (-@all(C15)+2*@all(C25)-@all(C35)) .*@all(v3_1_minus)/dx+
    (-@all(C14)+2*@all(C24)-@all(C34)) .*@all(v3_2_minus)/dy+
    (-@all(C13)+2*@all(C23)-@all(C33)) .*@all(v3_3_plus)/dz;

    @all(ax3)=(-1*@all(C11)-@all(C12)+2*@all(C13)) .*@all(v1_1_plus)/dx+
    (-1*@all(C16)-@all(C26)+2*@all(C36)) .*@all(v1_2_minus)/dy+
    (-1*@all(C15)-@all(C25)+2*@all(C35)) .*@all(v1_3_minus)/dz+
    (-1*@all(C16)-@all(C26)+2*@all(C36)) .*@all(v2_1_minus)/dx+
    (-1*@all(C12)-@all(C22)+2*@all(C23)) .*@all(v2_2_plus)/dy+
    (-1*@all(C14)-@all(C24)+2*@all(C34)) .*@all(v2_3_minus)/dz+
    (-1*@all(C15)-@all(C25)+2*@all(C35)) .*@all(v3_1_minus)/dx+
    (-1*@all(C14)-@all(C24)+2*@all(C34)) .*@all(v3_2_minus)/dy+
    (-1*@all(C13)-@all(C23)+2*@all(C33)) .*@all(v3_3_plus)/dz;

    @all(ax4)=@all(C16).*(@all(v1_1_plus)/dx)+
    @all(C56).*(@all(v3_1_minus)/dx)+
    @all(C66).*(@all(v2_1_minus)/dx)+
    @all(C26).*(@all(v2_2_plus)/dy)+
    @all(C46).*(@all(v3_2_minus)/dy)+
    @all(C66).*(@all(v1_2_minus)/dy)+
    @all(C36).*(@all(v3_3_plus)/dz)+
    @all(C46).*(@all(v2_3_minus)/dz)+
    @all(C56).*(@all(v1_3_minus)/dz);

    @all(ax5)=@all(C15) .*(@all(v1_1_plus)/dx)+
    @all(C56) .*(@all(v2_1_minus)/dx)+
    @all(C55) .*(@all(v3_1_minus)/dx)+
    @all(C25) .*(@all(v2_2_plus)/dy)+
    @all(C45) .*(@all(v3_2_minus)/dy)+
    @all(C56) .*(@all(v1_2_minus)/dy)+
    @all(C35) .*(@all(v3_3_plus)/dz)+
    @all(C45) .*(@all(v2_3_minus)/dz)+
    @all(C55) .*(@all(v1_3_minus)/dz);

    @all(ax6)=@all(C14).*(@all(v1_1_plus)/dx)+
    @all(C46).*(@all(v2_1_minus)/dx)+
    @all(C45).*(@all(v3_1_minus)/dx)+
    @all(C24).*(@all(v2_2_plus)/dy)+
    @all(C46).*(@all(v1_2_minus)/dy)+
    @all(C44).*(@all(v3_2_minus)/dy)+
    @all(C34).*(@all(v3_3_plus)/dz)+
    @all(C45).*(@all(v1_3_minus)/dz)+
    @all(C44).*(@all(v2_3_minus)/dz);

    @all(ax7)=(@all(C11)+@all(C12)+@all(C13)) .*@all(v1_1_plus)/dx+
    (@all(C16)+@all(C26)+@all(C36)) .*@all(v1_2_minus)/dy+
    (@all(C15)+@all(C25)+@all(C35)) .*@all(v1_3_minus)/dz+
    (@all(C16)+@all(C26)+@all(C36)) .*@all(v2_1_minus)/dx+
    (@all(C12)+@all(C22)+@all(C23)) .*@all(v2_2_plus)/dy+
    (@all(C14)+@all(C24)+@all(C34)) .*@all(v2_3_minus)/dz+
    (@all(C15)+@all(C25)+@all(C35)) .*@all(v3_1_minus)/dx+
    (@all(C14)+@all(C24)+@all(C34)) .*@all(v3_2_minus)/dy+
    (@all(C13)+@all(C23)+@all(C33)) .*@all(v3_3_plus)/dz;

    @all(ax_dt)=(@all(ax)-@all(Ax))/dt;
    @all(ax2_dt)=(@all(ax2)-@all(Ax2))/dt;
    @all(ax3_dt)=(@all(ax3)-@all(Ax3))/dt;
    @all(ax4_dt)=(@all(ax4)-@all(Ax4))/dt;
    @all(ax5_dt)=(@all(ax5)-@all(Ax5))/dt;
    @all(ax6_dt)=(@all(ax6)-@all(Ax6))/dt;
    @all(ax7_dt)=(@all(ax7)-@all(Ax7))/dt;
    return nothing
end
