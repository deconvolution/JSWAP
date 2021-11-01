dt=10^-3;
nt=1000;
st=[600,600,600];
rt=[900,750,825];

tt=load('./rec/rec_1.mat');
v1=tt.data;
tt=load('./rec/rec_2.mat');
v2=tt.data;
tt=load('./rec/rec_3.mat');
v3=tt.data;
r=norm(rt-st);

v=zeros(nt,3);
v(:,1)=v1;
v(:,2)=v2;
v(:,3)=v3;
%%
freq=10;
ns=nt;
M=2;
s=rickerWave(freq,dt,ns,M);
m=[1,0,0;
    0,-1,0;
    0,0,0];
%% propagation direction
d=rt-st;
n=[1,0,0;
    0,1,0;
    0,0,1];
gamma=zeros(3,1);
for i=1:3
    alp=d*n(i,:)'/norm(d);
    gamma(i)=(alp);
end
%%
Vp=1732;
Vs=1000;
rho=1000;
s_t=zeros(size(s));
s_t(1:end-1)=diff(s)/dt;
s_t_shift_p=zeros(size(s));
s_t_shift_s=s_t_shift_p;
s_t_shift_p(fix(r/Vp/dt):end)=s_t(1:end-fix(r/Vp/dt)+1);
s_t_shift_s(fix(r/Vs/dt):end)=s_t(1:end-fix(r/Vs/dt)+1);
up0=1/4/pi/r/rho/Vp^3;
us0=1/4/pi/r/rho/Vs^3;
%%
shift=8;
v_cal=zeros(nt,3);
figure
for i=1:3
    up=zeros(nt,3);
    us=zeros(nt,3);
    tt=0;
    tt2=0;
    up(:,i)=up0*s_t_shift_p;
    us(:,i)=us0*s_t_shift_s;
    for j=1:3
        for k=1:3
            tt=tt+gamma(i)*gamma(j)*gamma(k)*m(j,k);
            tt2=tt2+(delta(i,j)-gamma(i)*gamma(j))*gamma(k)*m(j,k);
        end
    end
    up(:,i)=up(:,i)*tt;
    us(:,i)=us(:,i)*tt2;
    
    vp=zeros(size(up));
    vs=zeros(size(us));
    vp(1:end-1,:)=diff(up,1,1)/dt;
    vs(1:end-1,:)=diff(us,1,1)/dt;
    vp(fix(r/Vp/dt)-1,:)=0;
    vs(fix(r/Vs/dt)-1,:)=0;
    
    
    subplot(1,3,i)
    v_cal(:,i)=(vp(:,i)+vs(:,i));
    v_cal(shift+1:end,i)=v_cal(1:end-shift,i,:);
    
    ax=plot(v_cal(:,i),dt:dt:dt*nt,'color','blue','linewidth',3);
    set(gca,'ydir','reverse');
    hold on;
    
    ax2=plot(v(:,i),dt:dt:dt*nt,'--','color','red','linewidth',3);
    set(gca,'ydir','reverse');
    xlabel({['v',num2str(i),'[m/s]']});
    if i==1
        ylabel('t [s]');
    end
    hold off;
    if i==3
        legend([ax,ax2],'Reference','FD','Location','southeast');
    end
end
%% calculate
for i=1:3
    corrcoef(v_cal(:,i),v(:,i))
end