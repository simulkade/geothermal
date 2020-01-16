Nx=100;
Ny=50;
Nz=1;
Lx=300;
Ly=200;
Lz=50;
v_cell=(Lx*Ly*Lz)/(Nx*Ny*Nz)
q_in=200/3600; %[m^3/s]
m=createMesh2D(Nx,Ny, Lx, Ly);
q=createCellVariable(m, 0.0);
qin=q;
qout=q;
q.value(25,25)=q_in/v_cell;
qin.value(25,25)=q_in/v_cell;
q.value(75,25)=-q_in/v_cell;
qout.value(75,25)=-q_in/v_cell;
p=continuity(m, q);
visualizeCells(p)
t_final = 5*30*24*3600; % [s]
T=geotherm_iso(m, p, qin, qout, t_final)
shading interp