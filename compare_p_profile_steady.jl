# compare pressure profiles
using DataFrames
df1 = readtable("p_profile_local_refine")
df2 = readtable("p_profile_coarse_peaceman")
figure(figsize=(5,2.5))
plot(df1[:x_m], df1[:p_pa]./1e5, "--o", label="locally refined grid")
plot(df2[:x_m], df2[:p_pa]./1e5, label="coarse, Peaceman well model")
xlabel("x [m]")
ylabel("p [bar]")
legend()
savefig("dp_t.png", bbox_inches="tight")

