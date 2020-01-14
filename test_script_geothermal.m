Nx=100;
Ny=50;
Nz=1;
Lx=2000;
Ly=500;
Lz=100;
v_cell=(Lx*Ly*Lz)/(Nx*Ny*Nz)
q_in=150/3600; %[m^3/s]
m=createMesh2D(Nx,Ny, 2000, 400);
q=createCellVariable(m, 0.0);
qin=q;
qout=q;
q.value(25,25)=q_in/v_cell;
qin.value(25,25)=q_in/v_cell;
q.value(75,25)=-q_in/v_cell;
qout.value(75,25)=-q_in/v_cell;
p=continuity(m, q);
visualizeCells(p)
T=geotherm_iso(m, p, qin, qout)
shading interp