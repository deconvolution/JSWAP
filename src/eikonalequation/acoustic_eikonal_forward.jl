function acoustic_eikonal_forward(;nx,ny,nz,h,v,s1,s2,s3,T0,s1t,s2t,s3t,r1,r2,r3,
    r1t,r2t,r3t,X,Y,Z,path,write_t=0)
    global data;
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
    Y1col,X1col,Z1col=meshgrid(2:ny-1,2:nx-1,2:nz-1);
    Z1col=reshape(Z1col,(nx-2)*(ny-2)*(nz-2),1);
    Y1col=reshape(Y1col,(nx-2)*(ny-2)*(nz-2),1);
    X1col=reshape(X1col,(nx-2)*(ny-2)*(nz-2),1);
    dis_s=(X1col .-s1).^2+(Y1col .-s2).^2+(Z1col .-s3).^2;

    tt=mapslices(sortperm,dis_s,dims=(1));

    X1col=X1col[tt];
    Y1col=Y1col[tt];
    Z1col=Z1col[tt];
    ## fast sweeping
    u=zeros(1,1);
    a=zeros(3,1);
    for l=1
        for i=1:size(X1col,1)
            a[1]=min(T[X1col[i]-1,Y1col[i],Z1col[i]],T[X1col[i]+1,Y1col[i],Z1col[i]]);
            a[2]=min(T[X1col[i],Y1col[i]-1,Z1col[i]],T[X1col[i],Y1col[i]+1,Z1col[i]]);
            a[3]=min(T[X1col[i],Y1col[i],Z1col[i]-1],T[X1col[i],Y1col[i],Z1col[i]+1]);
            a[:]=sort(a,dims=1);
            u[1]=a[1]+h/v[X1col[i],Y1col[i],Z1col[i]];
            if u[1]>a[2]
                u[1]=(a[1]+a[2]+sqrt(2*h^2/v[X1col[i],Y1col[i],Z1col[i]]^2-a[1]^2+2*a[1]*a[2]-a[2]^2))/2;
                if u[1]>a[3]
                    u[1]=(a[1]+a[2]+a[3]+sqrt(-2*a[1]^2+2*a[1]*a[2]+2*a[1]*a[3]-2*a[2]^2+
                    2*a[2]*a[3]-2*a[3]^2+3*h^2/v[X1col[i],Y1col[i],Z1col[i]]^2))/3;
                end
            end
            T[X1col[i],Y1col[i],Z1col[i]]=min(T[X1col[i],Y1col[i],Z1col[i]],u[1]);
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
    if write_t==1
        vtkfile=vtk_grid(string(path,"/traveltime_visualization"),X,Y,Z);
        vtkfile["v"]=v;
        vtkfile["T"]=T;
        vtk_save(vtkfile);
    end
    ## assign receivers
    T_obs=T[CartesianIndex.(r1,r2,r3)];
    ## save receivers
    data=T_obs;
    if write_t==1
        write2mat(string(path,"/recording.mat"),data);
    end
    return T,T_obs
end
