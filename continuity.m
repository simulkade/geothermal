function p_val=continuity(m, q)
mu_val=1e-3; % (Pa.s)
poros=0.2;
perm_val=1.0e-12; % (m^2)
% physical system
p_out=250e5; % (Pa)
% assign values to the domain
k=createCellVariable(m, perm_val);
labda_face=harmonicMean(k/mu_val);
phi=createCellVariable(m,poros);
% Define the boundaries
BCp = createBC(m); % Neumann BC for pressure
BCc = createBC(m); % Neumann BC for concentration
% change the right boandary to constant pressure (Dirichlet)
BCp.right.a(:)=0.0;
BCp.right.b(:)=1.0;
BCp.right.c(:)=p_out;
% left boundary
BCp.left.a(:)=0.0;
BCp.left.b(:)=1.0;
BCp.left.c(:)=p_out;
Mdiffp=diffusionTerm(-labda_face);
[Mbcp, RHSbcp] = boundaryCondition(BCp);
RHSsp=constantSourceTerm(q);
Mp= Mdiffp+Mbcp;
RHSp=RHSsp+RHSbcp;
p_val=solvePDE(m, Mp, RHSp);
end
