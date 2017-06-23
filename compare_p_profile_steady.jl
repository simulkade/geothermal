# compare pressure profiles
using DataFrames
df1 = readtable("p_profile_local_refine")
df2 = readtable("p_profile_coarse_peaceman")
figure()
plot(df1[:x_m], df1[:p_pa], "--o", label="locally refined grid")
plot(df2[:x_m], df2[:p_pa], label="coarse, Peaceman well model")
xlabel("x [m]")
ylabel("p [Pa]")
legend()
