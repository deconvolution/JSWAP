function acoustic_eikonal_forward(nx,ny,nz,h,v,s1,s2,s3,s1t,s2t,s3t,r1,r2,r3,r1t,r2t,r3t,path)
    ## create path
    if path!=nothing
        if isdir(path)==0
            mkdir(path);
        end
    end
    ## create CSV for source and receivers
    CSV.write(string(path,"/receiver location.csv"),DataFrame([r1t' r2t' r3t'],:auto));
    CSV.write(string(path,"/source location.csv"),DataFrame([s1t' s2t' s3t'],:auto));

    T=ones(nx,ny,nz)*3.1415926*10^12;
    T[s1,s2,s3]=0;

    ## compute distance to the source
    Y3D,X3D,Z3D=JSWAP.meshgrid(1:ny,1:nx,1:nz);
    Y,X,Z=JSWAP.meshgrid(2:ny-1,2:nx-1,2:nz-1);
    Z=reshape(Z,(nx-2)*(ny-2)*(nz-2),1);
    Y=reshape(Y,(nx-2)*(ny-2)*(nz-2),1);
    X=reshape(X,(nx-2)*(ny-2)*(nz-2),1);
    dis_s=(X .-s1).^2+(Y .-s2).^2+(Z .-s3).^2;

    tt=mapslices(sortperm,dis_s,dims=(1));

    X=X[tt];
    Y=Y[tt];
    Z=Z[tt];
    ## fast sweeping
    u=zeros(1,1);
    a=zeros(3,1);
    for l=1
        for i=1:size(X,1)
            a[1]=min(T[X[i]-1,Y[i],Z[i]],T[X[i]+1,Y[i],Z[i]]);
            a[2]=min(T[X[i],Y[i]-1,Z[i]],T[X[i],Y[i]+1,Z[i]]);
            a[3]=min(T[X[i],Y[i],Z[i]-1],T[X[i],Y[i],Z[i]+1]);
            a[:]=sort(a,dims=1);
            u[1]=a[1]+h/v[X[i],Y[i],Z[i]];
            if u[1]>a[2]
                u[1]=(a[1]+a[2]+sqrt(2*h^2/v[X[i],Y[i],Z[i]]^2-a[1]^2+2*a[1]*a[2]-a[2]^2))/2;
                if u[1]>a[3]
                    u[1]=(a[1]+a[2]+a[3]+sqrt(-2*a[1]^2+2*a[1]*a[2]+2*a[1]*a[3]-2*a[2]^2+
                    2*a[2]*a[3]-2*a[3]^2+3*h^2/v[X[i],Y[i],Z[i]]^2))/3;
                end
            end
            T[X[i],Y[i],Z[i]]=min(T[X[i],Y[i],Z[i]],u[1]);
        end
    end
    ## approximate boundaries
    T[1,:,:]=T[2,:,:];
    T[end,:,:]=T[end-1,:,:];
    T[:,1,:]=T[:,2,:];
    T[:,end,:]=T[:,end-1,:];
    T[:,:,1]=T[:,:,2];
    T[:,:,end]=T[:,:,end-1];
    ## save to vtk
    vtkfile=vtk_grid(string(path,"/traveltime_visualization"),X3D,Y3D,Z3D);
    vtkfile["v"]=v;
    vtkfile["T"]=T;
    vtk_save(vtkfile);
    ## assign receivers
    T_obs=T[CartesianIndex.(r1,r2,r3)];
    ## save receivers
    data=T_obs;
    write2mat(string(path,"recording.mat"),data);
    return T,T_obs
end
