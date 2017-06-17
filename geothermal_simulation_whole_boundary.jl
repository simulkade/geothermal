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
using CoolProp, DataFrames, PyPlot, Polynomials, JFVM
import GR
# a table of density values
depth_reservoir = 1000         # [m]
p_res = depth_reservoir/10*1e5 # [Pa]
T_res = 80 + 273.15            # [K]
T_inj = 25 + 273.15            # [K]
T     = collect(linspace(T_inj, T_res, 10))
rho_water = zeros(length(T))
mu_water  = zeros(length(T))
for i in eachindex(T)
    rho_water[i] = PropsSI("D", "T", T[i], "P", p_res, "water")           # [kg/m^3]
    mu_water[i]  = PropsSI("viscosity", "T", T[i], "P", p_res, "water")   # [Pa.s]
end

mu_fit  = polyfit(T, mu_water, 2)
rho_fit = polyfit(T, rho_water, 2)
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

T0    = 25.0 + 273.15  # [K] reference temperature
L     = 500            # [m] will be changed later in a loop
W     = 500            # [m]
thick_res = 100        # [m]
poros = 0.35           # [-]
perm_D  = 0.02        # Darcy
perm    = perm_D*1e-12 # [m^2]
V_dp  = 0.2            # Dykstra-Parsons coef
clx   = 0.1            # Correlation length in x direction
cly   = 0.1            # Correlation length in y direction
T_inj = 40 + 273.15    # [K]
theta_inj = T_inj - T0 # [K]
T_res = 80 + 273.15    # [K]
theta_res = T_res - T0 # [K]
flow_m3_h = 50        # [m3/h]
flow_inj  = flow_m3_h/3600 # m3/s
λ_water  = PropsSI("conductivity", "T", T_res, "P", p_res, "water") # W/m/K
λ_rock   = 0.16        # [W/m/K] rock conductivity (for sandstone from my thesis)
λ_eff    = λ_water^poros*λ_rock^(1-poros) # effective conductivity
cp_water = PropsSI("CPMASS", "T", T_res, "P", p_res, "water") # J/kg/K
cp_rock  = 837.0       # [J/kg/K] rock heat capacity (for sandstone)
rho_rock = 2650.0      # [kg/m^3] rock density (for sandstone)


well_diameter = 0.1
nx_well           = 10
ny_well           = 10
nx                = 100
ny                = 15
x_range = collect([linspace(0, nx_well*well_diameter, nx_well);
                   linspace((nx_well+1)*well_diameter, L-(nx_well+1)*well_diameter, nx);
                   linspace(L-nx_well*well_diameter, L, nx_well)])
y_range = collect([linspace(0, 0.5*W-(ny+1)*well_diameter, ny);
                   linspace(0.5*W-ny*well_diameter, 0.5*W+ny*well_diameter, 2*ny_well);
                   linspace(0.5*W+(ny+1)*well_diameter, W, ny)])
m   = createMesh2D(x_range, y_range)  # 2D Cartesian grid
Nx  = length(x_range)-1
Ny  = length(y_range)-1

perm_val   = perm # permfieldlogrnde(Nx, Ny, perm, V_dp, clx, cly)
perm_field = createCellVariable(m, perm_val)

left_range  = 1:Ny # ny+ny_well:ny+ny_well
right_range = 1:Ny #ny+ny_well:ny+ny_well
u_inj = flow_inj/(W*thick_res) # (m.cellsize.y[ny+ny_well+1]*thick_res)
BCp = createBC(m)                 # pressure boundary
BCp.left.a[left_range]   =
    perm_field.value[2, left_range]/polyval(mu_fit, T_inj)
BCp.left.b[left_range]   = 0.0
BCp.left.c[left_range]   = -u_inj
BCp.right.a[right_range] = 0.0
BCp.right.b[right_range] = 1.0
BCp.right.c[right_range] = p_res

BCt = createBC(m)                 # temperature boundary
# Danckwertz with assumptions (density of water is almost constant)
rho_water_inj  = polyval(rho_fit, T_inj)
BCt.left.a[left_range]  = λ_eff
BCt.left.b[left_range]  = -rho_water_inj*cp_water*u_inj
BCt.left.c[left_range]  = -rho_water_inj*cp_water*u_inj*theta_inj

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
final_time = 150*dt_init              # [s]
t_step     = collect([linspace(0, 5*dt_init, 30);
                      linspace(5.5*dt_init, final_time, 100)])
n_steps = length(t_step) # number of time steps
# variables to be stored:
dp_res = zeros(n_steps)
T_out  = zeros(n_steps)
T_out[1] = T_res

# discretization
M_BCp, RHS_BCp = boundaryConditionTerm(BCp)
M_BCt, RHS_BCt = boundaryConditionTerm(BCt)
M_conductivity = diffusionTerm(-harmonicMean(
        createCellVariable(m, λ_eff)))

for t_ind in 2:length(t_step)
  dt = t_step[t_ind]-t_step[t_ind-1]
  t  = t_step[t_ind]
  println(t_ind)
  for i in 1:3 # internal loop
    # solve pressure equation
    RHS_ddt_p   = constantSourceTerm(poros/dt*(rho_val - rho_init))
    water_mobil = harmonicMean(perm_field./mu_val)
    rho_face = arithmeticMean(rho_val)
    M_diff_p    = diffusionTerm(-rho_face.*water_mobil)
    p_val       = solveLinearPDE(m, M_diff_p + M_BCp,
        RHS_BCp - RHS_ddt_p)

    # velocity vector
    u = -water_mobil.*gradientTerm(p_val)

    # solve heat equation
    α                  = poros*rho_val*cp_water+(1-poros)*rho_rock*cp_rock
    RHS_ddt_t          = constantSourceTerm(cp_water*poros*theta_val.*(rho_val - rho_init)/dt)
    M_trans, RHS_trans = transientTerm(theta_init, dt, α)
    M_conv             = convectionUpwindTerm(cp_water*rho_face.*u)
    theta_val          = solveLinearPDE(m,
        M_BCt + M_conv + M_trans + M_conductivity,
        RHS_BCt + RHS_trans - RHS_ddt_t)

    # update density and viscosity values
    rho_val.value[:] = polyval(rho_fit, theta_val.value + T0)
    mu_val.value[:]  = polyval(mu_fit, theta_val.value + T0)
  end # end of inner loop
  # GR.imshow(internalCells(theta_val)+T0)
  T_out[t_ind]  = 0.5*(theta_val.value[end, 1+ny+ny_well]+theta_val.value[end-1, 1+ny+ny_well])+T0
  dp_res[t_ind] = 0.5*(p_val.value[1, 1+ny+ny_well]+p_val.value[2, 1+ny+ny_well])-p_res
                #  0.5*(p_val.value[end, 1+ny+ny_well]+p_val.value[end-1, 1+ny+ny_well])
  rho_init   = copyCell(rho_val)
  p_init     = copyCell(p_val) # not necessary
  theta_init = copyCell(theta_val)
end
dp_res[1] = dp_res[2] # pressure grad at time zero; assume that process has already started
df = DataFrame(T_s = t_step, dp_Pa = dp_res, T_K = T_out)
# DataFrames.writetable("Q-$flow_m3_h-L-$flow_m3_h-k$perm_D.csv", df)
figure(figsize=(8,2))
visualizeCells(theta_val+T0)
title("temperature profile")
colorbar()
figure()
plot(t_step/(365*24*3600), T_out)
xlabel("time [year]")
ylabel("Water temperature [K]")
figure()
plot(t_step/(365*24*3600), dp_res/1e5)
xlabel("time [year]")
ylabel("Pressure difference [bar]")
