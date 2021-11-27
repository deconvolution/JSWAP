function acoustic_eikonal_adjoint(;nx,ny,nz,h,T,r1,r2,r3,s1,s2,s3,R_cal,R_true)
    ## calculate misfit
    lambda=zeros(nx,ny,nz);
    R_diff=R_cal-R_true;

    ## calculate dirivative to T
    Tx=@zeros(nx,ny,nz);
    Ty=copy(Tx);
    Tz=copy(Tx);

    @parallel Dx_1(T,Tx,1,1,0,0,0,0);
    @parallel Dy_1(T,Ty,0,0,1,1,0,0);
    @parallel Dz_1(T,Tz,0,0,0,0,1,1);

    #Tx[1:end-1,:,:]=T[2:end,:,:]-T[1:end-1,:,:];
    #Ty[:,1:end-1,:]=T[:,2:end,:]-T[:,1:end-1,:];
    #Tz[:,:,1:end-1]=T[:,:,2:end]-T[:,:,1:end-1];

    a1_plus=(Tx+abs.(Tx))/h/2;
    a1_minus=(Tx-abs.(Tx))/h/2;
    a2_plus=(Ty+abs.(Ty))/h/2;
    a2_minus=(Ty-abs.(Ty))/h/2;
    a3_plus=(Tz+abs.(Tz))/h/2;
    a3_minus=(Tz-abs.(Tz))/h/2;

    ##
    A2=zeros(nx,ny,nz);
    A2[2:end,2:end,2:end]=1/h*(a1_minus[2:end,2:end,2:end]-a1_plus[1:end-1,2:end,2:end]+
    a2_minus[2:end,2:end,2:end]-a2_plus[2:end,1:end-1,2:end]+
    a3_minus[2:end,2:end,2:end]-a3_plus[2:end,2:end,1:end-1]);

    A2[1,:,:]=A2[2,:,:];
    A2[end,:,:]=A2[end-1,:,:];
    A2[:,1,:]=A2[:,2,:];
    A2[:,end,:]=A2[:,end-1,:];
    A2[:,:,1]=A2[:,:,2];
    A2[:,:,end]=A2[:,:,end-1];

    A2[s1,s2,s3]=1/6*(A2[s1-1,s2,s3]+A2[s1+1,s2,s3]+
    A2[s1,s2-1,s3]+A2[s1,s2+1,s3]+
    A2[s1,s2,s3-1]+A2[s1,s2,s3+1]);

    ## boundary condition of lambda
    lambda[CartesianIndex.(r1,r2,r3)]=R_diff./Tz[CartesianIndex.(r1,r2,r3)];
    ##
    lambda_old=copy(lambda);
    for l=1
        lambda_old[:]=lambda[:];
        for i=2:nx-1
            for j=2:ny-1
                for k=2:nz-1
                    A=A2[i,j,k];
                    B=1/h*(a1_minus[i-1,j,k]*lambda[i-1,j,k]-a1_plus[i,j,k]*lambda[i+1,j,k]+
                    a2_minus[i,j-1,k]*lambda[i,j-1,k]-a2_plus[i,j,k]*lambda[i,j+1,k]+
                    a3_minus[i,j,k-1]*lambda[i,j,k-1]-a3_plus[i,j,k]*lambda[i,j,k+1]);
                    lambda[i,j,k]=B/A;
                    if abs(lambda[i,j,k])<abs(lambda_old[i,j,k])
                        lambda[i,j,k]=lambda_old[i,j,k];
                    end
                end
            end
        end

        for i=nx-1:-1:2
            for j=2:ny-1
                for k=2:nz-1
                    A=A2[i,j,k];
                    B=1/h*(a1_minus[i-1,j,k]*lambda[i-1,j,k]-a1_plus[i,j,k]*lambda[i+1,j,k]+
                    a2_minus[i,j-1,k]*lambda[i,j-1,k]-a2_plus[i,j,k]*lambda[i,j+1,k]+
                    a3_minus[i,j,k-1]*lambda[i,j,k-1]-a3_plus[i,j,k]*lambda[i,j,k+1]);
                    lambda[i,j,k]=B/A;
                    if abs(lambda[i,j,k])<abs(lambda_old[i,j,k])
                        lambda[i,j,k]=lambda_old[i,j,k];
                    end
                end
            end
        end

        for i=2:nx-1
            for j=ny-1:-1:2
                for k=2:nz-1
                    A=A2[i,j,k];
                    B=1/h*(a1_minus[i-1,j,k]*lambda[i-1,j,k]-a1_plus[i,j,k]*lambda[i+1,j,k]+
                    a2_minus[i,j-1,k]*lambda[i,j-1,k]-a2_plus[i,j,k]*lambda[i,j+1,k]+
                    a3_minus[i,j,k-1]*lambda[i,j,k-1]-a3_plus[i,j,k]*lambda[i,j,k+1]);
                    lambda[i,j,k]=B/A;
                    if abs(lambda[i,j,k])<abs(lambda_old[i,j,k])
                        lambda[i,j,k]=lambda_old[i,j,k];
                    end
                end
            end
        end

        for i=2:nx-1
            for j=2:ny-1
                for k=nz-1:-1:2
                    A=A2[i,j,k];
                    B=1/h*(a1_minus[i-1,j,k]*lambda[i-1,j,k]-a1_plus[i,j,k]*lambda[i+1,j,k]+
                    a2_minus[i,j-1,k]*lambda[i,j-1,k]-a2_plus[i,j,k]*lambda[i,j+1,k]+
                    a3_minus[i,j,k-1]*lambda[i,j,k-1]-a3_plus[i,j,k]*lambda[i,j,k+1]);
                    lambda[i,j,k]=B/A;
                    if abs(lambda[i,j,k])<abs(lambda_old[i,j,k])
                        lambda[i,j,k]=lambda_old[i,j,k];
                    end
                end
            end
        end

        for i=2:nx-1
            for j=ny-1:-1:2
                for k=nz-1:-1:2
                    A=A2[i,j,k];
                    B=1/h*(a1_minus[i-1,j,k]*lambda[i-1,j,k]-a1_plus[i,j,k]*lambda[i+1,j,k]+
                    a2_minus[i,j-1,k]*lambda[i,j-1,k]-a2_plus[i,j,k]*lambda[i,j+1,k]+
                    a3_minus[i,j,k-1]*lambda[i,j,k-1]-a3_plus[i,j,k]*lambda[i,j,k+1]);
                    lambda[i,j,k]=B/A;
                    if abs(lambda[i,j,k])<abs(lambda_old[i,j,k])
                        lambda[i,j,k]=lambda_old[i,j,k];
                    end
                end
            end
        end

        for i=nx-1:-1:2
            for j=2:ny-1
                for k=nz-1:-1:2
                    A=A2[i,j,k];
                    B=1/h*(a1_minus[i-1,j,k]*lambda[i-1,j,k]-a1_plus[i,j,k]*lambda[i+1,j,k]+
                    a2_minus[i,j-1,k]*lambda[i,j-1,k]-a2_plus[i,j,k]*lambda[i,j+1,k]+
                    a3_minus[i,j,k-1]*lambda[i,j,k-1]-a3_plus[i,j,k]*lambda[i,j,k+1]);
                    lambda[i,j,k]=B/A;
                    if abs(lambda[i,j,k])<abs(lambda_old[i,j,k])
                        lambda[i,j,k]=lambda_old[i,j,k];
                    end
                end
            end
        end

        for i=nx-1:-1:2
            for j=ny-1:-1:2
                for k=2:nz-1
                    A=A2[i,j,k];
                    B=1/h*(a1_minus[i-1,j,k]*lambda[i-1,j,k]-a1_plus[i,j,k]*lambda[i+1,j,k]+
                    a2_minus[i,j-1,k]*lambda[i,j-1,k]-a2_plus[i,j,k]*lambda[i,j+1,k]+
                    a3_minus[i,j,k-1]*lambda[i,j,k-1]-a3_plus[i,j,k]*lambda[i,j,k+1]);
                    lambda[i,j,k]=B/A;
                    if abs(lambda[i,j,k])<abs(lambda_old[i,j,k])
                        lambda[i,j,k]=lambda_old[i,j,k];
                    end
                end
            end
        end

        for i=nx-1:-1:2
            for j=ny-1:-1:2
                for k=nz-1:-1:2
                    A=A2[i,j,k];
                    B=1/h*(a1_minus[i-1,j,k]*lambda[i-1,j,k]-a1_plus[i,j,k]*lambda[i+1,j,k]+
                    a2_minus[i,j-1,k]*lambda[i,j-1,k]-a2_plus[i,j,k]*lambda[i,j+1,k]+
                    a3_minus[i,j,k-1]*lambda[i,j,k-1]-a3_plus[i,j,k]*lambda[i,j,k+1]);
                    lambda[i,j,k]=B/A;
                    if abs(lambda[i,j,k])<abs(lambda_old[i,j,k])
                        lambda[i,j,k]=lambda_old[i,j,k];
                    end
                end
            end
        end
    end
    return lambda
end
