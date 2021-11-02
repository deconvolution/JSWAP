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
        vtkfile = vtk_grid(string(input2.path_model,"/material_properties"),input2.X,input2.Y,input2.Z);
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

    # wave vector
    v1=@zeros(input2.nx,input2.ny,input2.nz);
    v2=copy(v1);
    v3=copy(v1);

    sigmas11=copy(v1);
    sigmas22=copy(v1);
    sigmas33=copy(v1);
    sigmas23=copy(v1);
    sigmas13=copy(v1);
    sigmas12=copy(v1);
    p=copy(v1);

    v1_1_plus=copy(v1);
    v1_2_minus=copy(v1);
    v1_3_minus=copy(v1);
    v2_1_minus=copy(v1);
    v2_2_plus=copy(v1);
    v2_3_minus=copy(v1);
    v3_1_minus=copy(v1);
    v3_2_minus=copy(v1);
    v3_3_plus=copy(v1);

    sigmas11_1_minus=copy(v1);
    sigmas22_2_minus=copy(v1);
    sigmas33_3_minus=copy(v1);
    sigmas23_2_plus=copy(v1);
    sigmas23_3_plus=copy(v1);
    sigmas13_1_plus=copy(v1);
    sigmas13_3_plus=copy(v1);
    sigmas12_1_plus=copy(v1);
    sigmas12_2_plus=copy(v1);

    sigmas11_3_plus=copy(v1);
    sigmas22_3_plus=copy(v1);
    sigmas12_3_minus=copy(v1);

    p_1_minus=copy(v1);
    p_2_minus=copy(v1);
    p_3_minus=copy(v1);

    dtt1=@zeros(input2.nx-IND_accuracy,input2.ny,input2.nz);
    dtt2=@zeros(input2.nx,input2.ny-IND_accuracy,input2.nz);
    dtt3=@zeros(input2.nx,input2.ny,input2.nz-IND_accuracy);

    ax=copy(v1);
    ax2=copy(v1);
    ax3=copy(v1);
    ax4=copy(v1);
    ax5=copy(v1);
    ax6=copy(v1);
    ax7=copy(v1);
    Ax=copy(v1);
    Ax2=copy(v1);
    Ax3=copy(v1);
    Ax4=copy(v1);
    Ax5=copy(v1);
    Ax6=copy(v1);
    Ax7=copy(v1);
    ax_dt=copy(v1);
    ax2_dt=copy(v1);
    ax3_dt=copy(v1);
    ax4_dt=copy(v1);
    ax5_dt=copy(v1);
    ax6_dt=copy(v1);
    ax7_dt=copy(v1);

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

    Z_Kmax=input2.Z ./repeat(Kmax,1,1,input2.nz);
    Zmax_Kmax=repeat(maximum(input2.Z)./Kmax,1,1,input2.nz);

    for l=1:input2.nt-1

        @parallel Dx_inn(v1,dtt1);
        @parallel (2:input2.ny-1,2:input2.nz-1) u_1_plus(dtt1,v1_1_plus);

        @parallel Dy_inn(v1,dtt2);
        @parallel (2:input2.nx-1,2:input2.nz-1) u_2_minus(dtt2,v1_2_minus);

        @parallel Dz_inn(v1,dtt3);
        @parallel (2:input2.nx-1,2:input2.ny-1) u_3_minus(dtt3,v1_3_minus);

        @parallel Dx_inn(v2,dtt1);
        @parallel (2:input2.ny-1,2:input2.nz-1) u_1_minus(dtt1,v2_1_minus);

        @parallel Dy_inn(v2,dtt2);
        @parallel (2:input2.nx-1,2:input2.nz-1) u_2_plus(dtt2,v2_2_plus);

        @parallel Dz_inn(v2,dtt3);
        @parallel (2:input2.nx-1,2:input2.ny-1) u_3_minus(dtt3,v2_3_minus);

        @parallel Dx_inn(v3,dtt1);
        @parallel (2:input2.ny-1,2:input2.nz-1) u_1_minus(dtt1,v3_1_minus);

        @parallel Dy_inn(v3,dtt2);
        @parallel (2:input2.nx-1,2:input2.nz-1) u_2_minus(dtt2,v3_2_minus);

        @parallel Dz_inn(v3,dtt3);
        @parallel (2:input2.nx-1,2:input2.ny-1) u_3_plus(dtt3,v3_3_plus);

        @timeit ti "compute_sigma" @parallel JSWAP_CPU_3D_isotropic_forward_solver_compute_au_for_sigma_curvilinear(input2.dt,input2.dx,input2.dy,input2.dz,input2.inv_Qa,input2.lambda,input2.mu,
            beta,
            v1_1_plus,v1_2_minus,v1_3_minus,
            v2_1_minus,v2_2_plus,v2_3_minus,
            v3_1_minus,v3_2_minus,v3_3_plus,
            sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p,
            ax,ax2,ax3,ax4,ax5,ax6,ax7,
            Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
            ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt,
            Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax);

        @timeit ti "compute_sigma" @parallel JSWAP_CPU_3D_isotropic_forward_solver_compute_sigma_curvilinear(input2.dt,input2.dx,input2.dy,input2.dz,input2.inv_Qa,input2.lambda,input2.mu,
            beta,
            v1_1_plus,v1_2_minus,v1_3_minus,
            v2_1_minus,v2_2_plus,v2_3_minus,
            v3_1_minus,v3_2_minus,v3_3_plus,
            sigmas11,sigmas22,sigmas33,sigmas23,sigmas13,sigmas12,p,
            ax,ax2,ax3,ax4,ax5,ax6,ax7,
            Ax,Ax2,Ax3,Ax4,Ax5,Ax6,Ax7,
            ax_dt,ax2_dt,ax3_dt,ax4_dt,ax5_dt,ax6_dt,ax7_dt,
            Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax);

        # moment tensor source
        if isdefined(input2,:M33)
            if ns==1 && l<=size(Ms33_t,1)
                sigmas11[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas11[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms11_t[l];
                sigmas22[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas22[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms22_t[l];
                sigmas33[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas33[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms33_t[l];
                sigmas23[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas23[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms23_t[l];
                sigmas13[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas13[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms13_t[l];
                sigmas12[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas12[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Ms12_t[l];
                p[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=p[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-@ones(1,1)*input2.dt/input2.dx/input2.dy/input2.dz*Mp_t[l];
            end

            if ns>=2 && l<=size(Ms33_t,1)
                sigmas11[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas11[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms11_t[l,:]';
                sigmas22[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas22[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms22_t[l,:]';
                sigmas33[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas33[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms33_t[l,:]';
                sigmas23[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas23[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms23_t[l,:]';
                sigmas13[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas13[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms13_t[l,:]';
                sigmas12[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=sigmas12[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Ms12_t[l,:]';
                p[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=p[CartesianIndex.(input2.s1,input2.s2,input2.s3)]-input2.dt/input2.dx/input2.dy/input2.dz*Mp_t[l,:]';
            end
        end

        @parallel Dx_inn(sigmas11,dtt1);
        @parallel (2:input2.ny-1,2:input2.nz-1) u_1_minus(dtt1,sigmas11_1_minus);

        @parallel Dy_inn(sigmas22,dtt2);
        @parallel (2:input2.nx-1,2:input2.nz-1) u_2_minus(dtt2,sigmas22_2_minus);

        @parallel Dz_inn(sigmas33,dtt3);
        @parallel (2:input2.nx-1,2:input2.ny-1) u_3_minus(dtt3,sigmas33_3_minus);

        @parallel Dy_inn(sigmas23,dtt2);
        @parallel (2:input2.nx-1,2:input2.nz-1) u_2_plus(dtt2,sigmas23_2_plus);
        @parallel Dz_inn(sigmas23,dtt3);
        @parallel (2:input2.nx-1,2:input2.ny-1) u_3_plus(dtt3,sigmas23_3_plus);

        @parallel Dx_inn(sigmas13,dtt1);
        @parallel (2:input2.ny-1,2:input2.nz-1) u_1_plus(dtt1,sigmas13_1_plus);
        @parallel Dz_inn(sigmas13,dtt3);
        @parallel (2:input2.nx-1,2:input2.ny-1) u_3_plus(dtt3,sigmas13_3_plus);

        @parallel Dx_inn(sigmas12,dtt1);
        @parallel (2:input2.ny-1,2:input2.nz-1) u_1_plus(dtt1,sigmas12_1_plus);
        @parallel Dy_inn(sigmas12,dtt2);
        @parallel (2:input2.nx-1,2:input2.nz-1) u_2_plus(dtt2,sigmas12_2_plus);

        @parallel Dx_inn(p,dtt1);
        @parallel (2:input2.ny-1,2:input2.nz-1) u_1_minus(dtt1,p_1_minus);
        @parallel Dy_inn(p,dtt2);
        @parallel (2:input2.nx-1,2:input2.nz-1) u_2_minus(dtt2,p_2_minus);
        @parallel Dz_inn(p,dtt3);
        @parallel (2:input2.nx-1,2:input2.ny-1) u_3_minus(dtt3,p_3_minus);

        @timeit ti "compute_v" @parallel JSWAP_CPU_3D_isotropic_forward_solver_compute_v_curvilinear(input2.dt,input2.dx,input2.dy,input2.dz,input2.rho,beta,
            v1,v2,v3,
            sigmas11_1_minus,sigmas11_3_plus,
            sigmas22_2_minus,sigmas22_3_plus,
            sigmas33_3_minus,
            sigmas23_2_plus,sigmas23_3_plus,
            sigmas13_1_plus,sigmas13_3_plus,
            sigmas12_1_plus,sigmas12_2_plus,sigmas12_3_minus,
            p_1_minus,p_2_minus,p_3_minus,
            Kmax_x,Kmax_y,Z_Kmax,Zmax_Kmax);

        if isdefined(input2,:src3)
            if ns==1 && l<=size(input2.src3,1)
                v1[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v1[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src1[l];
                v2[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v2[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src2[l];
                v3[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v3[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src3[l];
                p[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=p[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+@ones(1,1) .*input2.srcp[l];
            end
            if ns>=2 && l<=size(input2.src3,1)
                v1[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v1[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src1[l,:]';
                v2[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v2[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src2[l,:]';
                v3[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=v3[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+1 ./input2.rho[CartesianIndex.(input2.s1,input2.s2,input2.s3)] .*input2.src3[l,:]';
                p[CartesianIndex.(input2.s1,input2.s2,input2.s3)]=p[CartesianIndex.(input2.s1,input2.s2,input2.s3)]+@ones(1,ns) .*input2.srcp[l,:]';
            end
        end

        # assign recordings
        @timeit ti "receiver" R1[l+1,:]=reshape(v1[CartesianIndex.(input2.r1,input2.r2,input2.r3)],length(input2.r3),);
        @timeit ti "receiver" R2[l+1,:]=reshape(v2[CartesianIndex.(input2.r1,input2.r2,input2.r3)],length(input2.r3),);
        @timeit ti "receiver" R3[l+1,:]=reshape(v3[CartesianIndex.(input2.r1,input2.r2,input2.r3)],length(input2.r3),);
        @timeit ti "receiver" P[l+1,:]=reshape(p[CartesianIndex.(input2.r1,input2.r2,input2.r3)],length(input2.r3),);
        # save wavefield
        if input2.path_wavefield!=nothing && input2.wavefield_interval!=0 && input2.path!=nothing
            if mod(l,input2.wavefield_interval)==0
                data=v1;
                write2mat(string(input2.path_wavefield,"/v1_",n_wavefield,".mat"),data);
                data=v2;
                write2mat(string(input2.path_wavefield,"/v2_",n_wavefield,".mat"),data);
                data=v3;
                write2mat(string(input2.path_wavefield,"/v3_",n_wavefield,".mat"),data);
                data=sigmas11;
                write2mat(string(input2.path_wavefield,"/sigmas11_",n_wavefield,".mat"),data);
                data=sigmas22;
                write2mat(string(input2.path_wavefield,"/sigmas22_",n_wavefield,".mat"),data);
                data=sigmas33;
                write2mat(string(input2.path_wavefield,"/sigmas33_",n_wavefield,".mat"),data);
                data=sigmas23;
                write2mat(string(input2.path_wavefield,"/sigmas23_",n_wavefield,".mat"),data);
                data=sigmas13;
                write2mat(string(input2.path_wavefield,"/sigmas13_",n_wavefield,".mat"),data);
                data=sigmas12;
                write2mat(string(input2.path_wavefield,"/sigmas12_",n_wavefield,".mat"),data);
                data=p;
                write2mat(string(input2.path_wavefield,"/p_",n_wavefield,".mat"),data);
                n_wavefield=n_wavefield+1;
            end
        end

        # plot
        if input2.path_pic!=nothing && input2.plot_interval!=0 && input2.path!=nothing
            if mod(l,input2.plot_interval)==0 || l==input2.nt-1
                vtkfile = vtk_grid(string(input2.path_pic,"/wavefield_pic_",n_picture),input2.X,input2.Y,input2.Z);
                vtkfile["v1"]=v1;
                vtkfile["v2"]=v2;
                vtkfile["v3"]=v3;
                vtkfile["p"]=p;
                vtkfile["sigmas33"]=sigmas33;
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
    return v1,v2,v3,R1,R2,R3,P
end
