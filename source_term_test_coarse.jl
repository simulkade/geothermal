# geothermal energy simulation; see http://fvt.simulkade.com/posts/2017-04-07-heat-transport-in-porous-media.html
# written by Ali A. Eftekhari, while on paternity leave!
# 14 June 2017; Lyngby public library
# output:
# pressure gradient (pressure difference between the injection and the production well)
# outlet temperature: temperature at the middle of the right boundary
# vs simulation time
#

# You can install CoolProp by:
# Pkg.clone("https://github.com/vimalaad/CoolProp.jl")
# Pkg.build("CoolProp")
using CoolProp, DataFrames, PyPlot, Polynomials, MUMPS, JFVM
import GR
# a table of density values
depth_reservoir = 2500         # [m]
p_res = depth_reservoir/10*1e5 # [Pa]
T_res = 80 + 273.15            # [K]
T_inj = 40 + 273.15            # [K]
T     = collect(linspace(T_inj, T_res, 10))
rho_water = zeros(length(T))
mu_water  = zeros(length(T))
for i in eachindex(T)
    rho_water[i] = PropsSI("D", "T", T[i], "P", p_res, "water")           # [kg/m^3]
    mu_water[i]  = PropsSI("viscosity", "T", T[i], "P", p_res, "water")   # [Pa.s]
end

mu_fit  = polyfit(T, mu_water, 2)
rho_fit = polyfit(T, rho_water, 2)
# figure()
# subplot(2,1,1)
# scatter(T, mu_water, label = "μ data")
# plot(T, polyval(mu_fit, T))
# xlabel("T [K]")
# ylabel("Viscosity [Pa.s]")
# legend()
# axis([minimum(T), maximum(T), minimum(mu_water), maximum(mu_water)])
#
# subplot(2,1,2)
# scatter(T, rho_water, label = "ρ data")
# plot(T, polyval(rho_fit, T))
# xlabel("T [K]")
# ylabel("Density [kg/m^3]")
# legend()

flow_m3_h = 100        # [m3/h]
flow_inj  = flow_m3_h/3600 # m3/s
perm_D  = 0.02        # Darcy
perm    = perm_D*1e-12 # [m^2]
L_res    = 500            # [m] will be changed later in a loop
L_left = 500           # [m] distance between the left boundary and injection well
L_right = 500          # [m] distance between the right boundary and injection well
L = L_res + L_left + L_right # [m] domain length
W     = 500            # [m]
thick_res = 100        # [m]
poros = 0.20           # [-]
T0    = 25.0 + 273.15  # [K] reference temperature
V_dp  = 0.2            # Dykstra-Parsons coef
clx   = 0.1            # Correlation length in x direction
cly   = 0.1            # Correlation length in y direction
T_inj = 40 + 273.15    # [K]
theta_inj = T_inj - T0 # [K]
T_res = 80 + 273.15    # [K]
theta_res = T_res - T0 # [K]
λ_water  = PropsSI("conductivity", "T", T_res, "P", p_res, "water") # W/m/K
λ_rock   = 0.16        # [W/m/K] rock conductivity (for sandstone from my thesis)
λ_eff    = λ_water^poros*λ_rock^(1-poros) # effective conductivity
cp_water = PropsSI("CPMASS", "T", T_res, "P", p_res, "water") # J/kg/K
cp_rock  = 837.0       # [J/kg/K] rock heat capacity (for sandstone)
rho_rock = 2650.0      # [kg/m^3] rock density (for sandstone)


well_diameter = 7*0.0254 # [m]; converted from a well diameter of 7 inch
# nx_left = 20
# nx_right = 20
# nx_well           = 10
# ny_well           = 10
# nx                = 40
# ny                = 15
# x_range =
#   collect([linspace(0, L_left-(nx_well+1)*well_diameter, nx_left);
#            linspace(L_left-nx_well*well_diameter, L_left+nx_well*well_diameter, 2*nx_well);
#            linspace(L_left+(nx_well+1)*well_diameter, L_left+L_res-(nx_well+1)*well_diameter, nx);
#            linspace(L_left+L_res-nx_well*well_diameter, L_left+L_res+nx_well*well_diameter, 2*nx_well);
#            linspace(L_left+L_res+(nx_well+1)*well_diameter, L, nx_right)])
# y_range = collect([linspace(0, 0.5*W-(ny+1)*well_diameter, ny);
#                    linspace(0.5*W-ny*well_diameter, 0.5*W+ny*well_diameter, 2*ny_well);
#                    linspace(0.5*W+(ny+1)*well_diameter, W, ny)])
# ind_inj  = [nx_left+nx_well+1, ny+ny_well]
# ind_prod = [nx_left+nx_well+nx_well+nx+nx_well, ny+ny_well]
mesh_refine = 8
Nx = mesh_refine*150
Ny = mesh_refine*50
m = createMesh2D(Nx, Ny, L, W)
ind_inj = mesh_refine*[50, 25]
ind_prod = mesh_refine*[100,25]
# m   = createMesh2D(x_range, y_range)  # 2D Cartesian grid
v_cell = thick_res*cellVolume(m)
r_eq = 0.2*m.cellsize.x[ind_inj[1]+1] # Peaceman well model
# Nx  = length(x_range)-1
# Ny  = length(y_range)-1

perm_val   = perm # permfieldlogrnde(Nx, Ny, perm, V_dp, clx, cly)
perm_field = createCellVariable(m, perm_val)

# left_range  = ny+ny_well:ny+ny_well
# right_range = ny+ny_well:ny+ny_well
# u_inj = flow_inj/(m.cellsize.y[ny+ny_well+1]*thick_res)
BCp = createBC(m)                 # pressure boundary
BCp.left.a[:] = 0.0
BCp.left.b[:] = 1.0
BCp.left.c[:] = p_res
BCp.right.a[:] = 0.0
BCp.right.b[:] = 1.0
BCp.right.c[:] = p_res

# BCp.left.a[left_range]   =
#     perm_field.value[2, left_range]/polyval(mu_fit, T_inj)
# BCp.left.b[left_range]   = 0.0
# BCp.left.c[left_range]   = -u_inj
# BCp.right.a[right_range] = 0.0
# BCp.right.b[right_range] = 1.0
# BCp.right.c[right_range] = p_res

BCt = createBC(m)                 # temperature boundary
# Danckwertz with assumptions (density of water is almost constant)
rho_water_inj  = polyval(rho_fit, T_inj)
# BCt.left.a[left_range]  = λ_eff
# BCt.left.b[left_range]  = -rho_water_inj*cp_water*u_inj
# BCt.left.c[left_range]  = -rho_water_inj*cp_water*u_inj*theta_inj

# figure(figsize=(8,2))
# visualizeCells(perm_field)
# title("Permeability [m^2]", fontsize = 10)
# colorbar()

theta_init = createCellVariable(m, theta_res, BCt)
theta_val  = createCellVariable(m, theta_res, BCt)
p_init     = createCellVariable(m, p_res, BCp)
p_val      = createCellVariable(m, p_res, BCp)
rho_init   = createCellVariable(m,
    polyval(rho_fit, theta_init.value+T0))
rho_val    = copyCell(rho_init)
mu_init    = createCellVariable(m,
    polyval(mu_fit, theta_init.value+T0))
mu_val     = copyCell(mu_init)

dt_init    = (L*W*thick_res*poros)/flow_inj/50 # (L/Nx)/(u_inj/poros)   # [s] time step
final_time = 100*dt_init              # [s]
t_step     = collect([linspace(0, 5*dt_init, 30);
                      linspace(5.5*dt_init, final_time, 70)])
n_steps = length(t_step) # number of time steps
# variables to be stored:
dp_res = zeros(n_steps)
T_out  = zeros(n_steps)
p_inj  = zeros(n_steps)
p_out  = zeros(n_steps)
T_out[1] = T_res

# discretization
M_BCp, RHS_BCp = boundaryConditionTerm(BCp)
M_BCt, RHS_BCt = boundaryConditionTerm(BCt)
M_conductivity = diffusionTerm(-harmonicMean(
        createCellVariable(m, λ_eff)))
q_inj = createCellVariable(m, 0)
q_inj.value[ind_inj+1 ...] = flow_inj/v_cell.value[ind_inj+1 ...]
RHSinjt = constantSourceTerm(rho_water_inj*cp_water*(T_inj-T0)*q_inj)

q_prod = createCellVariable(m, 0)
q_prod.value[ind_prod+1 ...] = -flow_inj/v_cell.value[ind_prod+1 ...]

q_cont = createCellVariable(m, 0)

  # solve pressure equation
  water_mobil = harmonicMean(perm_field./mu_val)
  rho_face = arithmeticMean(rho_val)
  M_diff_p    = diffusionTerm(-rho_face.*water_mobil)
  q_cont.value[ind_inj+1 ...] = rho_val.value[ind_inj+1 ...]*flow_inj/v_cell.value[ind_inj+1 ...]
  q_cont.value[ind_prod+1 ...] = -rho_val.value[ind_prod+1 ...]*flow_inj/v_cell.value[ind_prod+1 ...]
  RHScont = constantSourceTerm(q_cont)
  p_val       = solveLinearPDE(m, M_diff_p + M_BCp,
      RHS_BCp + RHScont)

  # velocity vector
  u = -water_mobil.*gradientTerm(p_val)

  WI_inj = 2*π*perm_val*thick_res/(mu_val.value[ind_inj+1 ...]*log(r_eq/(well_diameter/2)))
  WI_prod = 2*π*perm_val*thick_res/(mu_val.value[ind_prod+1 ...]*log(r_eq/(well_diameter/2)))

figure()
visualizeCells(p_val)
figure()
plot(m.cellcenters.x[ind_inj[1]:ind_prod[1]], p_val.value[ind_inj[1]+1:ind_prod[1]+1, ind_inj[2]+1], ".")
p_inj  = p_val.value[ind_inj+1 ...]+flow_inj/WI_inj
p_out  = p_val.value[ind_prod+1 ...]-flow_inj/WI_prod
dp_res = p_val.value[ind_inj+1 ...]-p_val.value[ind_prod+1 ...]
dp_peaceman = p_inj - p_out
println("p_inj= $p_inj, p_out= $p_out, dp = $dp_res, dp_peaceman= $dp_peaceman")
