function T=geotherm_iso(m, p, qin, qout, t_final)
labda_s=2.6; % W/(m.K)
labda_w=0.6; % W/(m.K)
rho_w=1000; % kg/m^3
rho_s=2650; % kg/m^3
cp_w=4184; % J/(kg.K)
cp_s=900; % J/(kg.K)
poros=0.2;
perm=1e-12; % m^2
mu_w=1e-3; % Pa.s
a=poros*rho_w*cp_w+(1-poros)*rho_s*cp_s;
labda=createCellVariable(m,labda_w^poros*labda_s^(1-poros));
T_init=80+273.15; % K
T_inj=40+273.15;  % K
T_0=25+273.15;
BC=createBC(m);
BC.left.a(:)=0.0;
BC.left.b(:)=1.0;
BC.left.c(:)=T_init;
BC.right.a(:)=0.0;
BC.right.b(:)=1.0;
BC.right.c(:)=T_init;
[Mbc,RHSbc]=boundaryCondition(BC);
T0=createCellVariable(m, T_init, BC)
u=-perm/mu_w*gradientTerm(p); % m/s
Mconv=convectionUpwindTerm(rho_w*cp_w*u);
RHSconv=divergenceTerm(rho_w*cp_w*T_0*u);
Mcond=diffusionTerm(harmonicMean(labda));
RHSs1=constantSourceTerm(rho_w*cp_w*(T_inj-T_0)*qin);
Ms=linearSourceTerm(-rho_w*cp_w*qout);
RHSs2=constantSourceTerm(-rho_w*cp_w*T_0*qout);
% dt=100*24*3600; % s
% t_final=200*dt;
dt = t_final/100
for t=dt:dt:t_final
    [Mt, RHSt]=transientTerm(T0, dt, a);
    M=Mt+Mconv-Mcond+Mbc+Ms;
    RHS=RHSbc+RHSt+RHSs1+RHSs2+RHSconv;
    T=solvePDE(m, M, RHS);
    T0=T;
    visualizeCells(T); drawnow;
end
end
