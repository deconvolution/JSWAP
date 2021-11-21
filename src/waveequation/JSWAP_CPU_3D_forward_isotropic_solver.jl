"
CPU 3D isotropic forward solver
"
function JSWAP_CPU_3D_forward_isotropic_curvilinear_solver(input2)
    #global data
    global data
    # zero stress condition at the boundaries
    input2.lambda[1:5,:,:] .=0;
    input2.lambda[end-4:end,:,:] .=0;
    input2.lambda[:,1:5,:] .=0;
    input2.lambda[:,end-4:end,:] .=0;
    input2.lambda[:,:,1:5] .=0;
    input2.lambda[:,:,end-4:end] .=0;

    input2.mu[1:5,:,:] .=0;
    input2.mu[end-4:end,:,:] .=0;
    input2.mu[:,1:5,:] .=0;
    input2.mu[:,end-4:end,:] .=0;
    input2.mu[:,:,1:5] .=0;
    input2.mu[:,:,end-4:end] .=0;



    IND_accuracy=9;
    IND_accuracy_c=3;
    d0=Dates.now();

    # source number
    ns=length(input2.s3);

    # create main folder
    if input2.path!=nothing
        if isdir(input2.path)==0
            mkdir(input2.path);
        end
    end

    # create folder for picture
    n_picture=1;
    n_wavefield=1;
    if input2.path_pic!=nothing
        if isdir(input2.path_pic)==0
            mkdir(input2.path_pic);
        end
        # initialize pvd
        pvd=paraview_collection(string(input2.path,"/time_info"));
    end

    # create folder for model
    if input2.path_model!=nothing && input2.path!=nothing
        if isdir(input2.path_model)==0
            mkdir(input2.path_model);
        end
        vtkfile = vtk_grid(string(input2.path_model,"/material_properties"),input2.X,input2.Y,input2.K);
        vtkfile["lambda"]=input2.lambda;
        vtkfile["mu"]=input2.mu;
        vtkfile["rho"]=input2.rho;
        vtk_save(vtkfile);
        CSV.write(string(input2.path_model,"/receiver location.csv"),DataFrame([input2.r1t' input2.r2t' input2.r3t'],:auto));
        CSV.write(string(input2.path_model,"/source location.csv"),DataFrame([input2.s1t' input2.s2t' input2.s3t'],:auto));
    end

    # create folder for wavefield
    if input2.path_wavefield!=nothing && input2.path!=nothing
        if isdir(input2.path_wavefield)==0
            mkdir(input2.path_wavefield)
        end
    end

    # create folder for rec
    if input2.path_rec!=nothing && input2.path!=nothing
        if isdir(input2.path_rec)==0
            mkdir(input2.path_rec)
        end
    end

    beta=isotropic_PML_configuration(input2.nx,input2.ny,input2.nz,input2.dx,input2.dy,input2.dz,
    input2.lambda,input2.mu,input2.rho,input2.nPML,input2.Rc,input2.lp,input2.PML_active);

    # receiver configuration
    R1=@zeros(input2.nt,length(input2.r3));
    R2=copy(R1);
    R3=copy(R1);
    P=@zeros(input2.nt,length(input2.r3));

    # decompose moment tensor
    if isdefined(input2,:M33)
        Mp=-1/3*(input2.M11+input2.M22+input2.M33);
        Ms11=input2.M11+Mp;
        Ms22=input2.M22+Mp;
        Ms33=input2.M33+Mp;
        Ms23=input2.M23;
        Ms13=input2.M13;
        Ms12=input2.M12;

        Mp_t=@zeros(size(input2.M33));
        Ms11_t=copy(Mp_t);
        Ms22_t=copy(Mp_t);
        Ms33_t=copy(Mp_t);
        Ms23_t=copy(Mp_t);
        Ms13_t=copy(Mp_t);
        Ms12_t=copy(Mp_t);

        Mp_t[1:end-1,:]=diff(Mp,dims=1)/input2.dt;
        Ms11_t[1:end-1,:]=diff(Ms11,dims=1)/input2.dt;
        Ms22_t[1:end-1,:]=diff(Ms22,dims=1)/input2.dt;
        Ms33_t[1:end-1,:]=diff(Ms33,dims=1)/input2.dt;
        Ms23_t[1:end-1,:]=diff(Ms23,dims=1)/input2.dt;
        Ms13_t[1:end-1,:]=diff(Ms13,dims=1)/input2.dt;
        Ms12_t[1:end-1,:]=diff(Ms12,dims=1)/input2.dt;
    end
    dx=input2.dx;
    dy=input2.dy;
    dz=input2.dz;
    # wave vector
    v1_iph_j_k=@zeros(input2.nx,input2.ny,input2.nz);
    v2_i_jph_k=copy(v1_iph_j_k);
    v3_i_j_kph=copy(v1_iph_j_k);

    sigmas11_i_j_k=copy(v1_iph_j_k);
    sigmas22_i_j_k=copy(v1_iph_j_k);
    sigmas33_i_j_k=copy(v1_iph_j_k);
    sigmas23_i_jph_kph=copy(v1_iph_j_k);
    sigmas13_iph_j_kph=copy(v1_iph_j_k);
    sigmas12_iph_jph_k=copy(v1_iph_j_k);
    p_i_j_k=copy(v1_iph_j_k);

    # derivatives
    v1_iph_j_k_1=copy(v1_iph_j_k);
    v1_iph_jp1_k_2=copy(v1_iph_j_k);
    v1_ip1_jph_kph_3=copy(v1_iph_j_k);
    v2_ip1_jph_k_1=copy(v1_iph_j_k);
    v2_i_jph_k_2=copy(v1_iph_j_k);
    v2_i_jph_kp1_3=copy(v1_iph_j_k);
    v3_ip1_jph_kph_3=copy(v1_iph_j_k);
    v3_i_jp1_kph_2=copy(v1_iph_j_k);
    v3_i_j_kph_3=copy(v1_iph_j_k);

    sigmas11_ip1_j_k_1=copy(v1_iph_j_k);
    sigmas22_i_jp1_k_2=copy(v1_iph_j_k);
    sigmas33_i_j_kp1_3=copy(v1_iph_j_k);
    sigmas23_i_jph_kph_2=copy(v1_iph_j_k);
    sigmas23_i_jph_kph_3=copy(v1_iph_j_k);
    sigmas13_iph_j_k_1=copy(v1_iph_j_k);
    sigmas13_iph_j_kph_3=copy(v1_iph_j_k);
    sigmas12_iph_jph_k_1=copy(v1_iph_j_k);
    sigmas12_iph_jph_k_2=copy(v1_iph_j_k);

    sigmas11_3_minus=copy(v1_iph_j_k);
    sigmas22_3_minus=copy(v1_iph_j_k);
    sigmas12_3_plus=copy(v1_iph_j_k);
    sigmas12_3_minus=copy(v1_iph_j_k);
    sigmas23_3_minus=copy(v1_iph_j_k);
    sigmas13_3_minus=copy(v1_iph_j_k);

    p_ip1_j_k_1=copy(v1_iph_j_k);
    p_i_jp1_k_2=copy(v1_iph_j_k);
    p_i_j_kp1_3=copy(v1_iph_j_k);
    # for curvilinear
    #=
    v1_3_sigmas11=copy(v1_iph_j_k);
    v2_3_sigmas11=copy(v1_iph_j_k);
    v3_3_sigmas23=copy(v1_iph_j_k);
    v3_3_sigmas13=copy(v1_iph_j_k);
    v2_3_sigmas12=copy(v1_iph_j_k);
    v1_3_sigmas12=copy(v1_iph_j_k);

    sigmas11_1_v1=copy(v1_iph_j_k);
    p_1_v1=copy(v1_iph_j_k);
    sigmas12_3_v1=copy(v1_iph_j_k);
    sigmas12_3_v2=copy(v1_iph_j_k);
    sigmas22_3_v2=copy(v1_iph_j_k);
    p_3_v2=copy(v1_iph_j_k);
    sigmas13_3_v3=copy(v1_iph_j_k);
    sigmas23_3_v3=copy(v1_iph_j_k);
    =#
    ax=copy(v1_iph_j_k);
    ax2=copy(v1_iph_j_k);
    ax3=copy(v1_iph_j_k);
    ax4=copy(v1_iph_j_k);
    ax5=copy(v1_iph_j_k);
    ax6=copy(v1_iph_j_k);
    ax7=copy(v1_iph_j_k);
    Ax=copy(v1_iph_j_k);
    Ax2=copy(v1_iph_j_k);
    Ax3=copy(v1_iph_j_k);
    Ax4=copy(v1_iph_j_k);
    Ax5=copy(v1_iph_j_k);
    Ax6=copy(v1_iph_j_k);
    Ax7=copy(v1_iph_j_k);
    ax_dt=copy(v1_iph_j_k);
    ax2_dt=copy(v1_iph_j_k);
    ax3_dt=copy(v1_iph_j_k);
    ax4_dt=copy(v1_iph_j_k);
    ax5_dt=copy(v1_iph_j_k);
    ax6_dt=copy(v1_iph_j_k);
    ax7_dt=copy(v1_iph_j_k);

    l=1;
    # save wavefield
    if input2.path_wavefield!=nothing && input2.wavefield_interval!=0 && input2.path!=nothing
        if mod(l,input2.wavefield_interval)==0
            data=zeros(input2.nx,input2.ny,input2.nz);
            write2mat(string(input2.path_wavefield,"/v1_",n_wavefield,".mat"),data);
            data=zeros(input2.nx,input2.ny,input2.nz);
            write2mat(string(input2.path_wavefield,"/v3_",n_wavefield,".mat"),data);
            data=zeros(input2.nx,input2.ny,input2.nz);
            write2mat(string(input2.path_wavefield,"/sigmas11_",n_wavefield,".mat"),data);
            data=zeros(input2.nx,input2.ny,input2.nz);
            write2mat(string(input2.path_wavefield,"/sigmas33_",n_wavefield,".mat"),data);
            data=zeros(input2.nx,input2.ny,input2.nz);
            write2mat(string(input2.path_wavefield,"/sigmas13_",n_wavefield,".mat"),data);
            data=zeros(input2.nx,input2.ny,input2.nz);
            write2mat(string(input2.path_wavefield,"/p_",n_wavefield,".mat"),data);
            n_wavefield=n_wavefield+1;
        end
    end

    pro_bar=Progress(input2.nt,1,"forward_simulation...",50);
    #
    Kmax=(input2.Kmax);
    Kmax_x0=@zeros(input2.nx,input2.ny);
    Kmax_y0=copy(Kmax_x0);

    Kmax_x0[1:end-1,:]=(Kmax[2:end,:]-Kmax[1:end-1,:])/input2.dx;
    Kmax_y0[:,1:end-1]=(Kmax[:,2:end]-Kmax[:,1:end-1])/input2.dy;
    Kmax_x=repeat(Kmax_x0,1,1,input2.nz);
    Kmax_y=repeat(Kmax_y0,1,1,input2.nz);
    # simulation coordinate
    Y,X,Z=meshgrid((1:input2.ny)*input2.dy,(1:input2.nx)*input2.dx,
    (1:input2.nz)*input2.dz);
    Z_Kmax=Z./repeat(Kmax,1,1,input2.nz);
    Zmax_Kmax=repeat(maximum(Z)./Kmax,1,1,input2.nz);
    K=input2.K;

    for l=1:input2.nt-1
        # plus: 6:5
        # minus: 5:6
        #@parallel Dx_inn(v1_iph_j_k,dtt1);
        #@parallel (2:input2.ny-1,2:input2.nz-1) u_1_plus(dtt1,v1_iph_j_k_1);
        @parallel Dx_12(v1_iph_j_k,v1_iph_j_k_1,6,5,0,0,0,0);

        #@parallel Dy_inn(v1_iph_j_k,dtt2);
        #@parallel (2:input2.nx-1,2:input2.nz-1) u_2_minus(dtt2,v1_iph_jp1_k_2);
        @parallel Dy_12(v1_iph_j_k,v1_iph_jp1_k_2,0,0,5,6,0,0);

        #@parallel Dz_inn(v1_iph_j_k,dtt3);
        #@parallel (2:input2.nx-1,2:input2.ny-1) u_3_minus(dtt3,v1_ip1_jph_kph_3);
        @parallel Dz_12(v1_iph_j_k,v1_ip1_jph_kph_3,0,0,0,0,5,6);

        #@parallel Dx_inn(v2_i_jph_k,dtt1);
        #@parallel (2:input2.ny-1,2:input2.nz-1) u_1_minus(dtt1,v2_ip1_jph_k_1);
        @parallel Dx_12(v2_i_jph_k,v2_ip1_jph_k_1,5,6,0,0,0,0);

        # @parallel Dy_inn(v2_i_jph_k,dtt2);
        # @parallel (2:input2.nx-1,2:input2.nz-1) u_2_plus(dtt2,v2_i_jph_k_2);
        @parallel Dy_12(v2_i_jph_k,v2_i_jph_k_2,0,0,6,5,0,0);

        #@parallel Dz_inn(v2_i_jph_k,dtt3);
        #@parallel (2:input2.nx-1,2:input2.ny-1) u_3_minus(dtt3,v2_i_jph_kp1_3);
        @parallel Dz_12(v2_i_jph_k,v2_i_jph_kp1_3,0,0,0,0,5,6);

        #@parallel Dx_inn(v3_i_j_kph,dtt1);
        #@parallel (2:input2.ny-1,2:input2.nz-1) u_1_minus(dtt1,v3_ip1_jph_kph_3);
        @parallel Dx_12(v3_i_j_kph,v3_ip1_jph_kph_3,5,6,0,0,0,0);

        #@parallel Dy_inn(v3_i_j_kph,dtt2);
        #@parallel (2:input2.nx-1,2:input2.nz-1) u_2_minus(dtt2,v3_i_jp1_kph_2);
        @parallel Dy_12(v3_i_j_kph,v3_i_jp1_kph_2,0,0,5,6,0,0);

        #@parallel Dz_inn(v3_i_j_kph,dtt3);
        #@parallel (2:input2.nx-1,2:input2.ny-1) u_3_plus(dtt3,v3_i_j_kph_3);
        @parallel Dz_12(v3_i_j_kph,v3_i_j_kph_3,0,0,0,0,6,5);

        @timeit ti "compute_sigma" @parallel JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma_curvilinear(input2.dt,input2.dx,input2.dy,input2.dz,input2.inv_Qa,input2.lambda,input2.mu,
        beta,
        v1_iph_j_k_1,v1_iph_jp1_k_2,v1_ip1_jph_kph_3,
        v2_ip1_jph_k_1,v2_i_jph_k_2,v2_i_jph_kp1_3,
        v3_ip1_jph_kph_3,v3_i_jp1_kph_2,v3_i_j_kph_3,
        sigmas11_i_j_k,
        sigmas22_i_j_k,
        sigmas33_i_j_k,
        sigmas23_i_jph_kph,
        sigmas13_iph_j_kph,
        sigmas12_iph_jph_k,
        ax,ax2,ax3,ax4,ax5,ax6,ax7,
        Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
        ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt,
        Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax);

        @timeit ti "compute_sigma" @parallel JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma_curvilinear(input2.dt,input2.dx,input2.dy,input2.dz,input2.inv_Qa,input2.lambda,input2.mu,
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
        Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax);

        # moment tensor source
        if isdefined(input2,:M33)
            if ns==1 && l<=size(Ms33_t,1)
                sigmas11_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas11_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms11_t[l];
                sigmas22_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas22_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms22_t[l];
                sigmas33_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas33_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms33_t[l];
                sigmas23_i_jph_kph[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas23_i_jph_kph[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms23_t[l];
                sigmas13_iph_j_kph[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas13_iph_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms13_t[l];
                sigmas12_iph_jph_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas12_iph_jph_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms12_t[l];
                p_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=p_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Mp_t[l];
            end

            if ns>=2 && l<=size(Ms33_t,1)
                sigmas11_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas11_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms11_t[l,:]';
                sigmas22_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas22_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms22_t[l,:]';
                sigmas33_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas33_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms33_t[l,:]';
                sigmas23_i_jph_kph[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas23_i_jph_kph[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms23_t[l,:]';
                sigmas13_iph_j_kph[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas13_iph_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms13_t[l,:]';
                sigmas12_iph_jph_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas12_iph_jph_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms12_t[l,:]';
                p_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=p_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Mp_t[l,:]';
            end
        end

        #@parallel Dx_inn(sigmas11_i_j_k,dtt1);
        #@parallel (2:input2.ny-1,2:input2.nz-1) u_1_minus(dtt1,sigmas11_ip1_j_k_1);
        @parallel Dx_12(sigmas11_i_j_k,sigmas11_ip1_j_k_1,5,6,0,0,0,0);

        #@parallel Dy_inn(sigmas22_i_j_k,dtt2);
        #@parallel (2:input2.nx-1,2:input2.nz-1) u_2_minus(dtt2,sigmas22_i_jp1_k_2);
        @parallel Dy_12(sigmas22_i_j_k,sigmas22_i_jp1_k_2,0,0,5,6,0,0);

        #@parallel Dz_inn(sigmas33_i_j_k,dtt3);
        #@parallel (2:input2.nx-1,2:input2.ny-1) u_3_minus(dtt3,sigmas33_i_j_kp1_3);
        @parallel Dz_12(sigmas33_i_j_k,sigmas33_i_j_kp1_3,0,0,0,0,5,6);

        #@parallel Dy_inn(sigmas23_i_jph_kph,dtt2);
        #@parallel (2:input2.nx-1,2:input2.nz-1) u_2_plus(dtt2,sigmas23_i_jph_kph_2);
        @parallel Dy_12(sigmas23_i_jph_kph,sigmas23_i_jph_kph_2,0,0,6,5,0,0);
        #@parallel Dz_inn(sigmas23_i_jph_kph,dtt3);
        #@parallel (2:input2.nx-1,2:input2.ny-1) u_3_plus(dtt3,sigmas23_i_jph_kph_3);
        @parallel Dz_12(sigmas23_i_jph_kph,sigmas23_i_jph_kph_3,0,0,0,0,6,5);

        #@parallel Dx_inn(sigmas13_iph_j_kph,dtt1);
        #@parallel (2:input2.ny-1,2:input2.nz-1) u_1_plus(dtt1,sigmas13_iph_j_k_1);
        @parallel Dx_12(sigmas13_iph_j_kph,sigmas13_iph_j_k_1,6,5,0,0,0,0);
        #@parallel Dz_inn(sigmas13_iph_j_kph,dtt3);
        #@parallel (2:input2.nx-1,2:input2.ny-1) u_3_plus(dtt3,sigmas13_iph_j_kph_3);
        @parallel Dz_12(sigmas13_iph_j_kph,sigmas13_iph_j_kph_3,0,0,0,0,6,5);

        #@parallel Dx_inn(sigmas12_iph_jph_k,dtt1);
        #@parallel (2:input2.ny-1,2:input2.nz-1) u_1_plus(dtt1,sigmas12_iph_jph_k_1);
        @parallel Dx_12(sigmas12_iph_jph_k,sigmas12_iph_jph_k_1,6,5,0,0,0,0);
        #@parallel Dy_inn(sigmas12_iph_jph_k,dtt2);
        #@parallel (2:input2.nx-1,2:input2.nz-1) u_2_plus(dtt2,sigmas12_iph_jph_k_2);
        @parallel Dy_12(sigmas12_iph_jph_k,sigmas12_iph_jph_k_2,0,0,6,5,0,0);

        #@parallel Dx_inn(p_i_j_k,dtt1);
        #@parallel (2:input2.ny-1,2:input2.nz-1) u_1_minus(dtt1,p_ip1_j_k_1);
        @parallel Dx_12(p_i_j_k,p_ip1_j_k_1,5,6,0,0,0,0);
        #@parallel Dy_inn(p_i_j_k,dtt2);
        #@parallel (2:input2.nx-1,2:input2.nz-1) u_2_minus(dtt2,p_i_jp1_k_2);
        @parallel Dy_12(p_i_j_k,p_i_jp1_k_2,0,0,5,6,0,0);
        #@parallel Dz_inn(p_i_j_k,dtt3);
        #@parallel (2:input2.nx-1,2:input2.ny-1) u_3_minus(dtt3,p_i_j_kp1_3);
        @parallel Dz_12(p_i_j_k,p_i_j_kp1_3,0,0,0,0,5,6);


        @timeit ti "compute_v" @parallel JSWAP_CPU_3D_isotropic_forward_solver_compute_v_curvilinear(input2.dt,input2.dx,input2.dy,input2.dz,input2.rho,beta,
        v1_iph_j_k,v2_i_jph_k,v3_i_j_kph,
        sigmas11_ip1_j_k_1,
        sigmas22_i_jp1_k_2,
        sigmas33_i_j_kp1_3,
        sigmas23_i_jph_kph_2,sigmas23_i_jph_kph_3,
        sigmas13_iph_j_k_1,sigmas13_iph_j_kph_3,
        sigmas12_iph_jph_k_1,sigmas12_iph_jph_k_2,
        p_ip1_j_k_1,p_i_jp1_k_2,p_i_j_kp1_3,
        Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax);

        if isdefined(input2,:src3)
            if ns==1 && l<=size(input2.src3,1)
                v1_iph_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v1_iph_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src1[l];
                v2_i_jph_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v2_i_jph_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src2[l];
                v3_i_j_kph[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v3_i_j_kph[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src3[l];
                p_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=p_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+@ones(1,1) .*input2.srcp[l];
            end
            if ns>=2 && l<=size(input2.src3,1)
                v1_iph_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v1_iph_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src1[l,:]';
                v2_i_jph_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v2_i_jph_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src2[l,:]';
                v3_i_j_kph[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v3_i_j_kph[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src3[l,:]';
                p_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=p_i_j_k[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+@ones(1,ns) .*input2.srcp[l,:]';
            end
        end

        # assign recordings
        @timeit ti "receiver" R1[l+1,:]=reshape(v1_iph_j_k[CartesianIndex.(input2.r1,input2.r2,input2.r3)],length(input2.r3),);
        @timeit ti "receiver" R2[l+1,:]=reshape(v2_i_jph_k[CartesianIndex.(input2.r1,input2.r2,input2.r3)],length(input2.r3),);
        @timeit ti "receiver" R3[l+1,:]=reshape(v3_i_j_kph[CartesianIndex.(input2.r1,input2.r2,input2.r3)],length(input2.r3),);
        @timeit ti "receiver" P[l+1,:]=reshape(p_i_j_k[CartesianIndex.(input2.r1,input2.r2,input2.r3)],length(input2.r3),);
        # save wavefield
        if input2.path_wavefield!=nothing && input2.wavefield_interval!=0 && input2.path!=nothing
            if mod(l,input2.wavefield_interval)==0
                data=v1_iph_j_k;
                write2mat(string(input2.path_wavefield,"/v1_",n_wavefield,".mat"),data);
                data=v2_i_jph_k;
                write2mat(string(input2.path_wavefield,"/v2_",n_wavefield,".mat"),data);
                data=v3_i_j_kph;
                write2mat(string(input2.path_wavefield,"/v3_",n_wavefield,".mat"),data);
                data=sigmas11_i_j_k;
                write2mat(string(input2.path_wavefield,"/sigmas11_",n_wavefield,".mat"),data);
                data=sigmas22_i_j_k;
                write2mat(string(input2.path_wavefield,"/sigmas22_",n_wavefield,".mat"),data);
                data=sigmas33_i_j_k;
                write2mat(string(input2.path_wavefield,"/sigmas33_",n_wavefield,".mat"),data);
                data=sigmas23_i_jph_kph;
                write2mat(string(input2.path_wavefield,"/sigmas23_",n_wavefield,".mat"),data);
                data=sigmas13_iph_j_kph;
                write2mat(string(input2.path_wavefield,"/sigmas13_",n_wavefield,".mat"),data);
                data=sigmas12_iph_jph_k;
                write2mat(string(input2.path_wavefield,"/sigmas12_",n_wavefield,".mat"),data);
                data=p_i_j_k;
                write2mat(string(input2.path_wavefield,"/p_",n_wavefield,".mat"),data);
                n_wavefield=n_wavefield+1;
            end
        end

        # plot
        if input2.path_pic!=nothing && input2.plot_interval!=0 && input2.path!=nothing
            if mod(l,input2.plot_interval)==0 || l==input2.nt-1
                vtkfile = vtk_grid(string(input2.path_pic,"/wavefield_pic_",n_picture),input2.X,input2.Y,K);
                vtkfile["v1"]=v1_iph_j_k;
                vtkfile["v2"]=v2_i_jph_k;
                vtkfile["v3"]=v3_i_j_kph;
                vtkfile["p"]=p_i_j_k;
                vtkfile["sigmas33"]=sigmas33_i_j_k;
                vtkfile["lambda"]=input2.lambda;
                vtkfile["mu"]=input2.mu;
                vtkfile["rho"]=input2.rho;
                pvd[input2.dt*(l+1)]=vtkfile;
                n_picture=n_picture+1;
            end
        end

        next!(pro_bar);
    end


    R1=R1 .*repeat(input2.Rm[:,1]',input2.nt,1);
    R2=R2 .*repeat(input2.Rm[:,2]',input2.nt,1);
    R3=R3 .*repeat(input2.Rm[:,3]',input2.nt,1);
    P=P .*repeat(input2.Rm[:,4]',input2.nt,1);

    if input2.path_rec!=nothing && input2.path!=nothing
        data=R1;
        write2mat(string(input2.path_rec,"/rec_1.mat"),data);
        data=R2;
        write2mat(string(input2.path_rec,"/rec_2.mat"),data);
        data=R3;
        write2mat(string(input2.path_rec,"/rec_3.mat"),data);
        data=P;
        write2mat(string(input2.path_rec,"/rec_p.mat"),data);
    end

    if input2.path_pic!=nothing && input2.plot_interval!=0 && input2.path!=nothing
        vtk_save(pvd);
    end
    return v1_iph_j_k,v2_i_jph_k,v3_i_j_kph,R1,R2,R3,P
end
