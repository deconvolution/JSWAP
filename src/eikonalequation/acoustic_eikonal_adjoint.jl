function acoustic_eikonal_adjoint(nx,ny,nz,h,T,r1,r2,r3,s1,s2,s3,R_cal,R_true)
    ## calculate misfit
    lambda=zeros(nx,ny,nz);
    R_diff=R_cal-R_true;

    ## calculate dirivative to T
    Tx=zeros(nx,ny,nz);
    Ty=copy(Tx);
    Tz=copy(Tx);

    Tx[1:end-1,:,:]=(T[2:end,:,:]-T[1:end-1,:,:])/h;
    Ty[:,1:end-1,:]=(T[:,2:end,:]-T[:,1:end-1,:])/h;
    Tz[:,:,1:end-1]=(T[:,:,2:end]-T[:,:,1:end-1])/h;

    a1_plus=(Tx+abs.(Tx))/2;
    a1_minus=(Tx-abs.(Tx))/2;
    a2_plus=(Ty+abs.(Ty))/2;
    a2_minus=(Ty-abs.(Ty))/2;
    a3_plus=(Tz+abs.(Tz))/2;
    a3_minus=(Tz-abs.(Tz))/2;
    ##
    A2=zeros(nx,ny,nz);
    for i=2:nx-1
        for j=2:ny-1
            for k=2:nz-1
                A2[i,j,k]=1/h*(a1_minus[i,j,k]-a1_plus[i-1,j,k]+
                a2_minus[i,j,k]-a2_plus[i,j-1,k]+
                a3_minus[i,j,k]-a3_plus[i,j,k-1]);
            end
        end
    end
    A2[1,:,:]=A2[2,:,:];
    A2[end,:,:]=A2[end-1,:,:];
    A2[:,1,:]=A2[:,2,:];
    A2[:,end,:]=A2[:,end-1,:];
    A2[:,:,1]=A2[:,:,2];
    A2[:,:,end]=A2[:,:,end-1];

    A2[s1,s2,s3]=1/6*(A2[s1[1]-1,s2,s3]+A2[s1[1]+1,s2,s3]+
    A2[s1,s2[1]-1,s3]+A2[s1,s2[1]+1,s3]+
    A2[s1,s2,s3[1]-1]+A2[s1,s2,s3[1]+1]);

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
