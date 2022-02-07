"
CPU 3D isotropic curvilinear forward solver.

Input:
  nt is the number of time steps.
  nx is the number of grids in the x-direciton.
  ny is the number of grids in the y-direction.
  nz is the number of grids in the z-direction.
  dt is the time interval for the simulation.
  dx is the grid spacing in the x-direction.
  dy is the grid spacing in the y-direction.
  dz is the grid spacing in the z-direction.
  X is the x true coordinate.
  Y is the y true coordinate.
  K is the z true coordinate.
  Kmax is the topography.
  lambda and mu are the Lame constants.
  rho is the density.
  inv_Qa is the apparent attenuation.
  s1 is the grid location of source in the x-direction.
  s2 is the grid location of source in the y-direction.
  s3 is the grid location of source in the z-direction.
  s1t is the true location of source in the x-direction.
  s2t is the true location of source in the y-direction.
  s3t is the true location of source in the z-direction.
  r1 is the grid location of receiver in the x-direction.
  r2 is the grid location of receiver in the y-direction.
  r3 is the grid location of receiver in the z-direction.
  r1t is the true location of receiver in the x-direction.
  r2t 2 is the true location of receiver in the y-direction.
  r3t is the true location of receiver in the z-direction.
  path is the path where one wants to save the output.
  wavefield_interval is the time interval for saving wavefield, nothing for not saving.
  plot_interval is the time interval for plot the wavefield, nothing for not plotting.
  M11, M22, M33, M23, M13 and M12 are the moment tensors.
  src1, src2, src3 and srcp are the source components for the x-, y-, z- and pressure components.
  lp is the number of PML layers.
  Rc is the theoretical reflectivity for PML.
  nPML is the power of the PML problem, normally 2.
  PML_active is whether or not to activate PML in each direciton.

Output:
  <path>/model/material_properties.vts is the pareview file to visualize material properties.
  <path>/model/receiver location.csv is the receiver location in a .csv file.
  <path>/model/source location.csv is the source location in a .csv file.
  <path>/pic/ contains .vts files for wavefield visualization in paraview.
  <path>/rec/ contains .mat files for seismograms.
  <path>/wavefield contains .mat files for wavefield.
  <path>/time_info.pvd is a file that can be loaded to paraview for wavefield movie.

Return:
  v1_iph_j_k is the v1 component at the last time step.
  v2_i_jph_k is the v2 component at the last time step.
  v3_i_j_kph is the v3 component at the last time step.
  R1 is the recording of v1.
  R2 is the recording of v2.
  R3 is the recording of v3.
  P is the recording of pressure.
"
function JSWAP_CPU_3D_forward_isotropic_curvilinear_solver(;nt,
    nx,
    ny,
    nz,
    dt,
    dx,
    dy,
    dz,
    X,
    Y,
    K,
    Kmax,
    lambda,
    mu,
    rho,
    inv_Qa,
    s1,
    s2,
    s3,
    s1t,
    s2t,
    s3t,
    r1,
    r2,
    r3,
    r1t,
    r2t,
    r3t,
    path,
    wavefield_interval=nothing,
    plot_interval=nothing,
    M11=nothing,
    M22=nothing,
    M33=nothing,
    M23=nothing,
    M13=nothing,
    M12=nothing,
    src1=nothing,
    src2=nothing,
    src3=nothing,
    srcp=nothing,
    lp,
    Rc,
    nPML,
    PML_active)

    #global data
    global data
    # zero stress condition at the boundaries
    lambda[1:5,:,:] .=0;
    lambda[end-4:end,:,:] .=0;
    lambda[:,1:5,:] .=0;
    lambda[:,end-4:end,:] .=0;
    lambda[:,:,1:5] .=0;
    lambda[:,:,end-4:end] .=0;

    mu[1:5,:,:] .=0;
    mu[end-4:end,:,:] .=0;
    mu[:,1:5,:] .=0;
    mu[:,end-4:end,:] .=0;
    mu[:,:,1:5] .=0;
    mu[:,:,end-4:end] .=0;

    nx=floor(Int64,nx);
    ny=floor(Int64,ny);
    nz=floor(Int64,nz);

    d0=Dates.now();

    # source number
    ns=length(s3);

    # create main folder
    if path!=nothing
        if isdir(path)==0
            mkdir(path);
        end
    end

    # create folder for picture
    n_picture=1;
    n_wavefield=1;
    path_pic=string(path,"/pic");
    if path!=nothing
        if isdir(path_pic)==0
            mkdir(path_pic);
        end
        # initialize pvd
        pvd=paraview_collection(string(path,"/time_info"));
    end

    # create folder for model
    path_model=string(path,"/model");
    if path!=nothing
        if isdir(path_model)==0
            mkdir(path_model);
        end
        vtkfile = vtk_grid(string(path_model,"/material_properties"),X,Y,K);
        vtkfile["lambda"]=lambda;
        vtkfile["mu"]=mu;
        vtkfile["rho"]=rho;
        vtk_save(vtkfile);
        CSV.write(string(path_model,"/receiver location.csv"),DataFrame([r1t' r2t' r3t'],:auto));
        CSV.write(string(path_model,"/source location.csv"),DataFrame([s1t' s2t' s3t'],:auto));
    end

    # create folder for wavefield
    path_wavefield=string(path,"/wavefield");
    if path!=nothing
        if isdir(path_wavefield)==0
            mkdir(path_wavefield)
        end
    end

    # create folder for rec
    path_rec=string(path,"/rec");
    if path!=nothing
        if isdir(path_rec)==0
            mkdir(path_rec)
        end
    end

    beta=isotropic_PML_configuration(nx,ny,nz,dx,dy,dz,
    lambda,mu,rho,nPML,Rc,lp,PML_active);

    # receiver configuration
    R1=@zeros(nt,length(r3));
    R2=copy(R1);
    R3=copy(R1);
    P=@zeros(nt,length(r3));

    # decompose moment tensor
    if M11!=nothing
        Mp=-1/3*(M11+M22+M33);
        Ms11=M11+Mp;
        Ms22=M22+Mp;
        Ms33=M33+Mp;
        Ms23=M23;
        Ms13=M13;
        Ms12=M12;

        Mp_t=@zeros(size(M33));
        Ms11_t=copy(Mp_t);
        Ms22_t=copy(Mp_t);
        Ms33_t=copy(Mp_t);
        Ms23_t=copy(Mp_t);
        Ms13_t=copy(Mp_t);
        Ms12_t=copy(Mp_t);

        Mp_t[1:end-1,:]=diff(Mp,dims=1)/dt;
        Ms11_t[1:end-1,:]=diff(Ms11,dims=1)/dt;
        Ms22_t[1:end-1,:]=diff(Ms22,dims=1)/dt;
        Ms33_t[1:end-1,:]=diff(Ms33,dims=1)/dt;
        Ms23_t[1:end-1,:]=diff(Ms23,dims=1)/dt;
        Ms13_t[1:end-1,:]=diff(Ms13,dims=1)/dt;
        Ms12_t[1:end-1,:]=diff(Ms12,dims=1)/dt;
    end
    dx=dx;
    dy=dy;
    dz=dz;
    # wave vector
    v1_iph_j_k=@zeros(nx,ny,nz);
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
    v1_3_sigmas11=copy(v1_iph_j_k);
    v2_3_sigmas11=copy(v1_iph_j_k);
    v3_3_sigmas23=copy(v1_iph_j_k);
    v3_3_sigmas13=copy(v1_iph_j_k);
    v2_3_sigmas12=copy(v1_iph_j_k);
    v1_3_sigmas12=copy(v1_iph_j_k);

    sigmas11_3_v1=copy(v1_iph_j_k);
    p_3_v1=copy(v1_iph_j_k);
    sigmas12_3_v1=copy(v1_iph_j_k);
    sigmas12_3_v2=copy(v1_iph_j_k);
    sigmas22_3_v2=copy(v1_iph_j_k);
    p_3_v2=copy(v1_iph_j_k);
    sigmas13_3_v3=copy(v1_iph_j_k);
    sigmas23_3_v3=copy(v1_iph_j_k);

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
    if wavefield_interval!=nothing && path!=nothing
        if mod(l,wavefield_interval)==0
            data=zeros(nx,ny,nz);
            write2mat(string(path_wavefield,"/v1_",n_wavefield,".mat"),data);
            data=zeros(nx,ny,nz);
            write2mat(string(path_wavefield,"/v3_",n_wavefield,".mat"),data);
            data=zeros(nx,ny,nz);
            write2mat(string(path_wavefield,"/sigmas11_",n_wavefield,".mat"),data);
            data=zeros(nx,ny,nz);
            write2mat(string(path_wavefield,"/sigmas33_",n_wavefield,".mat"),data);
            data=zeros(nx,ny,nz);
            write2mat(string(path_wavefield,"/sigmas13_",n_wavefield,".mat"),data);
            data=zeros(nx,ny,nz);
            write2mat(string(path_wavefield,"/p_",n_wavefield,".mat"),data);
            n_wavefield=n_wavefield+1;
        end
    end

    pro_bar=Progress(nt,1,"forward_simulation...",50);
    #
    Kmax_x0=@zeros(nx,ny);
    Kmax_y0=copy(Kmax_x0);

    Kmax_x0[1:end-1,:]=(Kmax[2:end,:]-Kmax[1:end-1,:])/dx;
    Kmax_y0[:,1:end-1]=(Kmax[:,2:end]-Kmax[:,1:end-1])/dy;
    Kmax_x=repeat(Kmax_x0,1,1,nz);
    Kmax_y=repeat(Kmax_y0,1,1,nz);
    # simulation coordinate
    ~,~,Z=meshgrid((1:ny)*dy,(1:nx)*dx,
    (1:nz)*dz);
    Z_Kmax=Z./repeat(Kmax,1,1,nz);
    Zmax_Kmax=repeat(maximum(Z)./Kmax,1,1,nz);

    for l=1:nt-1
        # plus: 6:5
        # minus: 5:6
        #@parallel Dx_inn(v1_iph_j_k,dtt1);
        #@parallel (2:ny-1,2:nz-1) u_1_plus(dtt1,v1_iph_j_k_1);
        @parallel Dx_12(v1_iph_j_k,v1_iph_j_k_1,6,5,0,0,0,0);

        #@parallel Dy_inn(v1_iph_j_k,dtt2);
        #@parallel (2:nx-1,2:nz-1) u_2_minus(dtt2,v1_iph_jp1_k_2);
        @parallel Dy_12(v1_iph_j_k,v1_iph_jp1_k_2,0,0,5,6,0,0);

        #@parallel Dz_inn(v1_iph_j_k,dtt3);
        #@parallel (2:nx-1,2:ny-1) u_3_minus(dtt3,v1_ip1_jph_kph_3);
        @parallel Dz_12(v1_iph_j_k,v1_ip1_jph_kph_3,0,0,0,0,5,6);

        #@parallel Dx_inn(v2_i_jph_k,dtt1);
        #@parallel (2:ny-1,2:nz-1) u_1_minus(dtt1,v2_ip1_jph_k_1);
        @parallel Dx_12(v2_i_jph_k,v2_ip1_jph_k_1,5,6,0,0,0,0);

        # @parallel Dy_inn(v2_i_jph_k,dtt2);
        # @parallel (2:nx-1,2:nz-1) u_2_plus(dtt2,v2_i_jph_k_2);
        @parallel Dy_12(v2_i_jph_k,v2_i_jph_k_2,0,0,6,5,0,0);

        #@parallel Dz_inn(v2_i_jph_k,dtt3);
        #@parallel (2:nx-1,2:ny-1) u_3_minus(dtt3,v2_i_jph_kp1_3);
        @parallel Dz_12(v2_i_jph_k,v2_i_jph_kp1_3,0,0,0,0,5,6);

        #@parallel Dx_inn(v3_i_j_kph,dtt1);
        #@parallel (2:ny-1,2:nz-1) u_1_minus(dtt1,v3_ip1_jph_kph_3);
        @parallel Dx_12(v3_i_j_kph,v3_ip1_jph_kph_3,5,6,0,0,0,0);

        #@parallel Dy_inn(v3_i_j_kph,dtt2);
        #@parallel (2:nx-1,2:nz-1) u_2_minus(dtt2,v3_i_jp1_kph_2);
        @parallel Dy_12(v3_i_j_kph,v3_i_jp1_kph_2,0,0,5,6,0,0);

        #@parallel Dz_inn(v3_i_j_kph,dtt3);
        #@parallel (2:nx-1,2:ny-1) u_3_plus(dtt3,v3_i_j_kph_3);
        @parallel Dz_12(v3_i_j_kph,v3_i_j_kph_3,0,0,0,0,6,5);

        # for curvilinear
        #v1_3_sigmas11[k1,k2,k3]=(v1_iph_j_k[k1,k2,(k3 .+1)]+v1_iph_j_k[(k1 .-1),k2,(k3 .+1)]-
        #v1_iph_j_k[k1,k2,(k3 .-1)]-v1_iph_j_k[(k1 .-1),k2,(k3 .-1)])/4;
        @parallel CUR1(v1_iph_j_k,v1_3_sigmas11,1,1,0,0,2,1);

        #v2_3_sigmas11[k1,k2,k3]=(v2_i_jph_k[k1,k2,(k3 .+1)]+v2_i_jph_k[k1,(k2 .-1),(k3 .+1)]-
        #v2_i_jph_k[k1,k2,(k3 .-1)]-v2_i_jph_k[k1,(k2 .-1),(k3 .-1)])/4;
        @parallel CUR3(v2_i_jph_k,v2_3_sigmas11,0,0,1,1,2,1);

        #v3_3_sigmas23[k1,k2,k3]=(v3_i_j_kph[k1,(k2 .+1),(k3 .+1)]+v3_i_j_kph[k1,k2,(k3 .+1)]-
        #v3_i_j_kph[k1,(k2 .+1),(k3 .-1)]-v3_i_j_kph[k1,k2,(k3 .-1)])/4;
        @parallel CUR3(v3_i_j_kph,v3_3_sigmas23,0,0,2,0,2,1);

        #v3_3_sigmas13[k1,k2,k3]=(v3_i_j_kph[(k1 .+1),k2,(k3 .+1)]+v3_i_j_kph[k1,k2,(k3 .+1)]-
        #v3_i_j_kph[(k1 .+1),k2,(k3 .-1)]-v3_i_j_kph[k1,k2,(k3 .-1)])/4;
        @parallel CUR1(v3_i_j_kph,v3_3_sigmas13,2,0,0,0,2,1);

        #v2_3_sigmas12[k1,k2,k3]=(v2_i_jph_k[(k1 .+1),k2,(k3 .+1)]+v2_i_jph_k[k1,k2,(k3 .+1)]-
        #v2_i_jph_k[(k1 .+1),k2,(k3 .-1)]-v2_i_jph_k[k1,k2,(k3 .-1)])/4;
        @parallel CUR1(v2_i_jph_k,v2_3_sigmas12,2,0,0,0,2,1);

        #v1_3_sigmas12[k1,k2,k3]=(v1_iph_j_k[k1,(k2 .+1),k3]+v1_iph_j_k[k1,k2,(k3 .+1)]-
        #v1_iph_j_k[k1,k2,(k3 .-1)]-v1_iph_j_k[k1,(k2 .+1),(k3 .-1)])/4;
        @parallel CUR3(v1_iph_j_k,v1_3_sigmas12,0,0,2,0,2,1);

        @timeit ti "compute_sigma" @parallel JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma_curvilinear(dt,dx,dy,dz,inv_Qa,lambda,mu,
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
        Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax);

        @timeit ti "compute_sigma" @parallel JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma_curvilinear(dt,dx,dy,dz,inv_Qa,lambda,mu,
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
        if M11!=nothing
            if ns==1 && l<=size(Ms33_t,1)
                sigmas11_i_j_k[CartesianIndex.(s1,s2,s3)]=sigmas11_i_j_k[CartesianIndex.(s1,s2,s3)]-@ones(1,1)*dt/dx/dy/dz*Ms11_t[l];
                sigmas22_i_j_k[CartesianIndex.(s1,s2,s3)]=sigmas22_i_j_k[CartesianIndex.(s1,s2,s3)]-@ones(1,1)*dt/dx/dy/dz*Ms22_t[l];
                sigmas33_i_j_k[CartesianIndex.(s1,s2,s3)]=sigmas33_i_j_k[CartesianIndex.(s1,s2,s3)]-@ones(1,1)*dt/dx/dy/dz*Ms33_t[l];
                sigmas23_i_jph_kph[CartesianIndex.(s1,s2,s3)]=sigmas23_i_jph_kph[CartesianIndex.(s1,s2,s3)]-@ones(1,1)*dt/dx/dy/dz*Ms23_t[l];
                sigmas13_iph_j_kph[CartesianIndex.(s1,s2,s3)]=sigmas13_iph_j_k[CartesianIndex.(s1,s2,s3)]-@ones(1,1)*dt/dx/dy/dz*Ms13_t[l];
                sigmas12_iph_jph_k[CartesianIndex.(s1,s2,s3)]=sigmas12_iph_jph_k[CartesianIndex.(s1,s2,s3)]-@ones(1,1)*dt/dx/dy/dz*Ms12_t[l];
                p_i_j_k[CartesianIndex.(s1,s2,s3)]=p_i_j_k[CartesianIndex.(s1,s2,s3)]-@ones(1,1)*dt/dx/dy/dz*Mp_t[l];
            end

            if ns>=2 && l<=size(Ms33_t,1)
                sigmas11_i_j_k[CartesianIndex.(s1,s2,s3)]=sigmas11_i_j_k[CartesianIndex.(s1,s2,s3)]-dt/dx/dy/dz*Ms11_t[l,:]';
                sigmas22_i_j_k[CartesianIndex.(s1,s2,s3)]=sigmas22_i_j_k[CartesianIndex.(s1,s2,s3)]-dt/dx/dy/dz*Ms22_t[l,:]';
                sigmas33_i_j_k[CartesianIndex.(s1,s2,s3)]=sigmas33_i_j_k[CartesianIndex.(s1,s2,s3)]-dt/dx/dy/dz*Ms33_t[l,:]';
                sigmas23_i_jph_kph[CartesianIndex.(s1,s2,s3)]=sigmas23_i_jph_kph[CartesianIndex.(s1,s2,s3)]-dt/dx/dy/dz*Ms23_t[l,:]';
                sigmas13_iph_j_kph[CartesianIndex.(s1,s2,s3)]=sigmas13_iph_j_k[CartesianIndex.(s1,s2,s3)]-dt/dx/dy/dz*Ms13_t[l,:]';
                sigmas12_iph_jph_k[CartesianIndex.(s1,s2,s3)]=sigmas12_iph_jph_k[CartesianIndex.(s1,s2,s3)]-dt/dx/dy/dz*Ms12_t[l,:]';
                p_i_j_k[CartesianIndex.(s1,s2,s3)]=p_i_j_k[CartesianIndex.(s1,s2,s3)]-dt/dx/dy/dz*Mp_t[l,:]';
            end
        end

        #@parallel Dx_inn(sigmas11_i_j_k,dtt1);
        #@parallel (2:ny-1,2:nz-1) u_1_minus(dtt1,sigmas11_ip1_j_k_1);
        @parallel Dx_12(sigmas11_i_j_k,sigmas11_ip1_j_k_1,5,6,0,0,0,0);

        #@parallel Dy_inn(sigmas22_i_j_k,dtt2);
        #@parallel (2:nx-1,2:nz-1) u_2_minus(dtt2,sigmas22_i_jp1_k_2);
        @parallel Dy_12(sigmas22_i_j_k,sigmas22_i_jp1_k_2,0,0,5,6,0,0);

        #@parallel Dz_inn(sigmas33_i_j_k,dtt3);
        #@parallel (2:nx-1,2:ny-1) u_3_minus(dtt3,sigmas33_i_j_kp1_3);
        @parallel Dz_12(sigmas33_i_j_k,sigmas33_i_j_kp1_3,0,0,0,0,5,6);

        #@parallel Dy_inn(sigmas23_i_jph_kph,dtt2);
        #@parallel (2:nx-1,2:nz-1) u_2_plus(dtt2,sigmas23_i_jph_kph_2);
        @parallel Dy_12(sigmas23_i_jph_kph,sigmas23_i_jph_kph_2,0,0,6,5,0,0);
        #@parallel Dz_inn(sigmas23_i_jph_kph,dtt3);
        #@parallel (2:nx-1,2:ny-1) u_3_plus(dtt3,sigmas23_i_jph_kph_3);
        @parallel Dz_12(sigmas23_i_jph_kph,sigmas23_i_jph_kph_3,0,0,0,0,6,5);

        #@parallel Dx_inn(sigmas13_iph_j_kph,dtt1);
        #@parallel (2:ny-1,2:nz-1) u_1_plus(dtt1,sigmas13_iph_j_k_1);
        @parallel Dx_12(sigmas13_iph_j_kph,sigmas13_iph_j_k_1,6,5,0,0,0,0);
        #@parallel Dz_inn(sigmas13_iph_j_kph,dtt3);
        #@parallel (2:nx-1,2:ny-1) u_3_plus(dtt3,sigmas13_iph_j_kph_3);
        @parallel Dz_12(sigmas13_iph_j_kph,sigmas13_iph_j_kph_3,0,0,0,0,6,5);

        #@parallel Dx_inn(sigmas12_iph_jph_k,dtt1);
        #@parallel (2:ny-1,2:nz-1) u_1_plus(dtt1,sigmas12_iph_jph_k_1);
        @parallel Dx_12(sigmas12_iph_jph_k,sigmas12_iph_jph_k_1,6,5,0,0,0,0);
        #@parallel Dy_inn(sigmas12_iph_jph_k,dtt2);
        #@parallel (2:nx-1,2:nz-1) u_2_plus(dtt2,sigmas12_iph_jph_k_2);
        @parallel Dy_12(sigmas12_iph_jph_k,sigmas12_iph_jph_k_2,0,0,6,5,0,0);

        #@parallel Dx_inn(p_i_j_k,dtt1);
        #@parallel (2:ny-1,2:nz-1) u_1_minus(dtt1,p_ip1_j_k_1);
        @parallel Dx_12(p_i_j_k,p_ip1_j_k_1,5,6,0,0,0,0);
        #@parallel Dy_inn(p_i_j_k,dtt2);
        #@parallel (2:nx-1,2:nz-1) u_2_minus(dtt2,p_i_jp1_k_2);
        @parallel Dy_12(p_i_j_k,p_i_jp1_k_2,0,0,5,6,0,0);
        #@parallel Dz_inn(p_i_j_k,dtt3);
        #@parallel (2:nx-1,2:ny-1) u_3_minus(dtt3,p_i_j_kp1_3);
        @parallel Dz_12(p_i_j_k,p_i_j_kp1_3,0,0,0,0,5,6);

        # compute curvilinear derivatives
        #sigmas11_3_v1[k1,k2,k3]=(sigmas11_i_j_k[(k1 .+1),k2,(k3 .+1)]+sigmas11_i_j_k[k1,k2,(k3 .+1)]-
        #sigmas11_i_j_k[(k1 .+1),k2,(k3 .-1)]-sigmas11_i_j_k[k1,k2,(k3 .-1)])/4;
        @parallel CUR1(sigmas11_i_j_k,sigmas11_3_v1,2,0,0,0,2,1);

        #p_3_v1[k1,k2,k3]=(p_i_j_k[(k1 .+1),k2,(k3 .+1)]+p_i_j_k[k1,k2,(k3 .+1)]-
        #p_i_j_k[(k1 .+1),k2,(k3 .-1)]-p_i_j_k[k1,k2,(k3 .-1)])/4;
        @parallel CUR1(p_i_j_k,p_3_v1,2,0,0,0,2,1);

        #sigmas12_3_v1[k1,k2,k3]=(sigmas12_iph_jph_k[k1,k2,(k3 .+1)]+sigmas12_iph_jph_k[k1,(k2 .-1),(k3 .+1)]-
        #sigmas12_iph_jph_k[k1,k2,(k3 .-1)]-sigmas12_iph_jph_k[k1,(k2 .-1),(k3 .-1)])/4;
        @parallel CUR3(sigmas12_iph_jph_k,sigmas12_3_v1,0,0,1,1,2,1);

        #sigmas12_3_v2[k1,k2,k3]=(sigmas12_iph_jph_k[k1,k2,(k3 .+1)]+sigmas12_iph_jph_k[(k1 .-1),k2,(k3 .+1)]-
        #sigmas12_iph_jph_k[k1,k2,(k3 .-1)]-sigmas12_iph_jph_k[(k1 .-1),k2,(k3 .-1)])/4;
        @parallel CUR1(sigmas12_iph_jph_k,sigmas12_3_v2,1,1,0,0,2,1);

        #sigmas22_3_v2[k1,k2,k3]=(sigmas22_i_j_k[k1,(k2 .+1),(k3 .+1)]+sigmas22_i_j_k[k1,k2,(k3 .+1)]-
        #sigmas22_i_j_k[k1,(k2 .+1),(k3 .-1)]-sigmas22_i_j_k[k1,k2,(k3 .-1)])/4;
        @parallel CUR3(sigmas22_i_j_k,sigmas22_3_v2,0,0,2,0,2,1);

        #p_3_v2[k1,k2,k3]=(p_i_j_k[k1,(k2 .+1),(k3 .+1)]+p_i_j_k[k1,k2,(k3 .+1)]-
        #p_i_j_k[k1,(k2 .+1),(k3 .-1)]-p_i_j_k[k1,k2,(k3 .-1)])/4;
        @parallel CUR3(p_i_j_k,p_3_v2,0,0,2,0,2,1);

        #sigmas13_3_v3[k1,k2,k3]=(sigmas13_iph_j_kph[k1,k2,(k3 .+1)]+sigmas13_iph_j_kph[(k1 .-1),k2,(k3 .+1)]-
        #sigmas13_iph_j_kph[k1,k2,(k3 .-1)]-sigmas13_iph_j_kph[(k1 .-1),k2,(k3 .-1)])/4;
        @parallel CUR1(sigmas13_iph_j_kph,sigmas13_3_v3,1,1,0,0,2,1);

        #sigmas23_3_v3[k1,k2,k3]=(sigmas23_i_jph_kph[k1,k2,(k3 .+1)]+sigmas23_i_jph_kph[k1,(k2 .-1),(k3 .+1)]-
        #sigmas23_i_jph_kph[k1,k2,(k3 .-1)]-sigmas23_i_jph_kph[k1,(k2 .-1),(k3 .-1)])/4;
        @parallel CUR3(sigmas23_i_jph_kph,sigmas23_3_v3,0,0,1,1,2,1);

        @timeit ti "compute_v" @parallel JSWAP_CPU_3D_isotropic_forward_solver_compute_v_curvilinear(dt,dx,dy,dz,rho,beta,
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
        Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax);

        if src1!=nothing
            if ns==1 && l<=size(src3,1)
                v1_iph_j_k[CartesianIndex.(s1,s2,s3)]=v1_iph_j_k[CartesianIndex.(s1,s2,s3)]+1 ./rho[CartesianIndex.(s1,s2,s3)] .*src1[l];
                v2_i_jph_k[CartesianIndex.(s1,s2,s3)]=v2_i_jph_k[CartesianIndex.(s1,s2,s3)]+1 ./rho[CartesianIndex.(s1,s2,s3)] .*src2[l];
                v3_i_j_kph[CartesianIndex.(s1,s2,s3)]=v3_i_j_kph[CartesianIndex.(s1,s2,s3)]+1 ./rho[CartesianIndex.(s1,s2,s3)] .*src3[l];
                p_i_j_k[CartesianIndex.(s1,s2,s3)]=p_i_j_k[CartesianIndex.(s1,s2,s3)]+@ones(1,1) .*srcp[l];
            end
            if ns>=2 && l<=size(src3,1)
                v1_iph_j_k[CartesianIndex.(s1,s2,s3)]=v1_iph_j_k[CartesianIndex.(s1,s2,s3)]+1 ./rho[CartesianIndex.(s1,s2,s3)] .*src1[l,:]';
                v2_i_jph_k[CartesianIndex.(s1,s2,s3)]=v2_i_jph_k[CartesianIndex.(s1,s2,s3)]+1 ./rho[CartesianIndex.(s1,s2,s3)] .*src2[l,:]';
                v3_i_j_kph[CartesianIndex.(s1,s2,s3)]=v3_i_j_kph[CartesianIndex.(s1,s2,s3)]+1 ./rho[CartesianIndex.(s1,s2,s3)] .*src3[l,:]';
                p_i_j_k[CartesianIndex.(s1,s2,s3)]=p_i_j_k[CartesianIndex.(s1,s2,s3)]+@ones(1,ns) .*srcp[l,:]';
            end
        end

        # assign recordings
        @timeit ti "receiver" R1[l+1,:]=reshape(v1_iph_j_k[CartesianIndex.(r1,r2,r3)],length(r3),);
        @timeit ti "receiver" R2[l+1,:]=reshape(v2_i_jph_k[CartesianIndex.(r1,r2,r3)],length(r3),);
        @timeit ti "receiver" R3[l+1,:]=reshape(v3_i_j_kph[CartesianIndex.(r1,r2,r3)],length(r3),);
        @timeit ti "receiver" P[l+1,:]=reshape(p_i_j_k[CartesianIndex.(r1,r2,r3)],length(r3),);
        # save wavefield
        if wavefield_interval!=nothing && path!=nothing
            if mod(l,wavefield_interval)==0
                data=v1_iph_j_k;
                write2mat(string(path_wavefield,"/v1_",n_wavefield,".mat"),data);
                data=v2_i_jph_k;
                write2mat(string(path_wavefield,"/v2_",n_wavefield,".mat"),data);
                data=v3_i_j_kph;
                write2mat(string(path_wavefield,"/v3_",n_wavefield,".mat"),data);
                data=sigmas11_i_j_k;
                write2mat(string(path_wavefield,"/sigmas11_",n_wavefield,".mat"),data);
                data=sigmas22_i_j_k;
                write2mat(string(path_wavefield,"/sigmas22_",n_wavefield,".mat"),data);
                data=sigmas33_i_j_k;
                write2mat(string(path_wavefield,"/sigmas33_",n_wavefield,".mat"),data);
                data=sigmas23_i_jph_kph;
                write2mat(string(path_wavefield,"/sigmas23_",n_wavefield,".mat"),data);
                data=sigmas13_iph_j_kph;
                write2mat(string(path_wavefield,"/sigmas13_",n_wavefield,".mat"),data);
                data=sigmas12_iph_jph_k;
                write2mat(string(path_wavefield,"/sigmas12_",n_wavefield,".mat"),data);
                data=p_i_j_k;
                write2mat(string(path_wavefield,"/p_",n_wavefield,".mat"),data);
                n_wavefield=n_wavefield+1;
            end
        end

        # plot
        if plot_interval!=nothing && path!=nothing
            if mod(l,plot_interval)==0 || l==nt-1
                vtkfile = vtk_grid(string(path_pic,"/wavefield_pic_",n_picture),X,Y,K);
                vtkfile["v1"]=v1_iph_j_k;
                vtkfile["v2"]=v2_i_jph_k;
                vtkfile["v3"]=v3_i_j_kph;
                vtkfile["p"]=p_i_j_k;
                vtkfile["sigmas33"]=sigmas33_i_j_k;
                vtkfile["lambda"]=lambda;
                vtkfile["mu"]=mu;
                vtkfile["rho"]=rho;
                pvd[dt*(l+1)]=vtkfile;
                n_picture=n_picture+1;
            end
        end

        next!(pro_bar);
    end

    if path!=nothing
        data=R1;
        write2mat(string(path_rec,"/rec_1.mat"),data);
        data=R2;
        write2mat(string(path_rec,"/rec_2.mat"),data);
        data=R3;
        write2mat(string(path_rec,"/rec_3.mat"),data);
        data=P;
        write2mat(string(path_rec,"/rec_p.mat"),data);
    end

    if plot_interval!=nothing && path!=nothing
        vtk_save(pvd);
    end
    return v1_iph_j_k,v2_i_jph_k,v3_i_j_kph,R1,R2,R3,P
end
