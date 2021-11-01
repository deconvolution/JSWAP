function ricker=rickerWave(freq,dt,ns,M)
%% calculate scale
E=10.^(5.24+1.44.*M);
s=sqrt(E.*freq/.299);
%%
t=0:dt:dt*ns;
t0=1./freq;
t=t-t0;
ricker=s*(1-2*pi^2*freq.^2*t.^2).*exp(-pi^2*freq^2*t.^2);
ricker=ricker(2:end)';
end