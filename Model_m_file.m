clc
clear all

Times=54000;
dt=0.005;
t=0:dt:Times;
N=length(t);

[rho_s,rho_g,H_Hmax,P_eqa,P_eqd,T,P,f_r,f_in,Q,T_out]=deal(zeros(N,1));%同时初始化多个变量
rho_s0=8400;%density of empty metal hydride, kg/m3
P_0=101325;%Pa
a=13.7;
b=3704;
beta=0.11414;
m_MH=47;%kg
Phi=0.33;
Phi_0=0.008584;
R_u=8.314472;%J/molK
MW_H2=2.01598e-3;%kg/mol
C_a=59.187;
C_d=9.57;
E_a=21179.6;
E_d=16420;
C_pg=14890;
epsilon=0.5;
C_ps=419;
Delta_Ha=30478;
Delta_Hd=30800;
U=25;%W/(m2k)
f_w=0.1667;%kg/s
C_pw=1860;
T_win=293.15;%K
A_s=1.347;
rho_w=1000;
V_tank=0.0165;
V_H2=7.849;%m3
rho_H2=0.0897;%kg/m3
P_ref=1500000;

rho_s(1)=8412.4;
rho_g(1)=0.0831;
T(1)=293.15;

 for i=1:N
    V_MH = m_MH/(rho_s0*epsilon);
    V1=V_tank/V_MH-1+epsilon;
    V2=1-epsilon;
    rho_ss = rho_s0*(1+V_H2/m_MH*rho_H2);
    
    H_Hmax(i)=(rho_s(i)-rho_s0)/(rho_ss-rho_s0);
    P_eqa(i) = P_0*exp(a-b/T(i)+(Phi+Phi_0)*tan(pi*(H_Hmax(i)-0.5))+beta/2);
    P_eqd(i) = P_0*exp(a-b/T(i)+(Phi-Phi_0)*tan(pi*(H_Hmax(i)-0.5))-beta/2);
    
    P(i)=rho_g(i)*R_u/MW_H2*T(i);
    alpha=U*A_s/(f_w*C_pw);
    Q(i) = f_w*C_pw*(T_win-T(i))*(1-exp(-alpha));
    T_out(i)=T(i)+(T_win-T(i))*exp(-alpha);
    
    if P(i)<P_ref
        f_in(i)=0.0012;
    else
        f_in(i)=0;
    end

    if P(i)>P_eqa(i)
        f_r(i) = C_a*exp(-E_a/(R_u*T(i)))*log(P(i)/P_eqa(i))*(rho_ss-rho_s(i));
        T(i+1) = T(i)+ (-C_pg*(f_in(i)-f_r(i))*T(i)-C_ps*f_r(i)*T(i)+f_r(i)*Delta_Ha/MW_H2+ Q(i)/V_MH)/(V1*C_pg*rho_g(i)+V2*C_ps*rho_s(i))*dt;
    elseif P(i)<P_eqd(i)
        f_r(i) = C_d*exp(-E_d/(R_u*T(i)))*((P(i)-P_eqd(i))/P_eqd(i))*(rho_s(i)-rho_s0);
        T(i+1) = T(i)+ (-C_pg*(f_in(i)-f_r(i))*T(i)-C_ps*f_r(i)*T(i)+f_r(i)*Delta_Hd/MW_H2+ Q(i)/V_MH)/(V1*C_pg*rho_g(i)+V2*C_ps*rho_s(i))*dt;
    else
        f_r(i)=0;
        T(i+1) = T(i)+ (-C_pg*f_in(i)*T(i)+ Q(i)/V_MH)/(V1*C_pg*rho_g(i)+V2*C_ps*rho_s(i))*dt;
    end
    
    rho_g(i+1) = rho_g(i)+(f_in(i)-f_r(i))/V1*dt;
    rho_s(i+1) = rho_s(i)+f_r(i)/V2*dt;
end