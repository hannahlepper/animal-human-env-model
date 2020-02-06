using JLD2

@load "D:/results_workspace.jld2" dat
@load "D:/results_workspace_P.jld2" P

dat_E_1, dat_E_2, dat_E_3, dat_E_4, dat_E_5, dat_E_6, dat_E_7,
       dat_B_1, dat_B_2, dat_B_3, dat_B_4, dat_B_5, dat_B_6, dat_B_7,
       dat_Bd_1, dat_Bd_2, dat_Bd_3, dat_Bd_4, dat_Bd_5, dat_Bd_6, dat_Bd_7,
       dat_A_1, dat_A_2, dat_A_3, dat_A_4, dat_A_5, dat_A_6, dat_A_7,
       dat_H_1, dat_H_2, dat_H_3, dat_H_4, dat_H_5, dat_H_6, dat_H_7 = dat
p_E, p_E_2, p_E_3, p_E_4, p_E_5, p_E_6, p_E_7,
     p_B, p_B_2, p_B_3, p_B_4, p_B_5, p_B_6, p_B_7,
     p_Bd, p_Bd_2, p_Bd_3, p_Bd_4, p_Bd_5, p_Bd_6, p_Bd_7,
     p_A, p_A_2, p_A_3, p_A_4, p_A_5, p_A_6, p_A_7,
     p_H, p_H_2, p_H_3, p_H_4, p_H_5, p_H_6, p_H_7 = P

#GET PRESENCE/ABSENCE OF TARGETS REACHED

#Get measures of impact
#1. for fixed LA 0.1 -> 0.0
function get_impact(dat_pre_int, dat_post_int)
    if !(dat_pre_int==0.) #can't have a 0 in denominator
        1 - dat_post_int/dat_pre_int #will get %decrease and %increase here
    else
        1 - dat_post_int/(dat_pre_int + 0.0000001) #will get %decrease and %increase here
    end
end

function map_impacts(dat_pre_col, dat_post_col)
    map(x -> get_impact(dat_pre_col[x], dat_post_col[x]), 1:length(dat_pre_col))
end

#Bins for parameters
lower_bin = [0.:0.05:1;]
bin_N = size(lower_bin)[1]

#1.fixed LA intervention, varying bEH
impact_1 = [map_impacts(dat[i,1][:,2], dat[i,2][:,2]) for i in 1:5]

#2. varyng LA intervention, varying bEH
impact_2 = [map_impacts(dat[i,3][:,2], dat[i,2][:,2]) for i in 1:5]

#3. varying bEH intervention, fixed LA
impact_3 = [map_impacts(dat[i,1][:,2], dat[i,5][:,2]) for i in 1:5]

#4. varying LA intervention, fixed bEH
impact_4 = [map_impacts(dat[i,4][:,2], dat[i,6][:,2]) for i in 1:5]

#5. varying LA intervention, varying bEH, low bHA
impact_5 = [map_impacts(dat[i,9][:,2], dat[1, 8][:,2]) for i in 1:5]

#Did simulations reach target RH of 0.65 - 0.75?
n_target = hcat([Int.(0.65 .< dat[x,1][:,1] .< 0.75) for x in 1:5])

#Did simulations have a low impact of less than 2%?
n_low_impact = hcat([Int.(impact_1[x] .< 0.02) for x in 1:5])


#Using RCall to get % reaching target
using RCall
R"""
source("RGetPercAndPlot.R")
"""

@rput lower_bin
@rput P
@rput n_target n_low_impact
@rput impact_1 impact_2 impact_3 impact_4 impact_5


R"str(P)"
R"str(impact_1)"

@rget get_mean_var
get_mean_var(lower_bin, P[1][:,[11,15]], impact_1[1])

bEHmuE_df_E_impact = get_mean_var(lower_bin, p_E[, c(11, 15)], impact_E)
impact_1
# Conclusion 1: realistic RHs are attainable for environmental transmission scenarios
R"""
bEH_df = lapply(P[[1]]get_perc_target()
bEHmuE_df_E = get_perc_target(lower_bin, P[[1]][[]][, c(11, 15)], n_target_E)
bEHmuH_df_E = get_perc_target(lower_bin, p_E[, c(11, 13)], n_target_E)
bEHLH_df_E = get_perc_target(lower_bin, p_E[, c(11, 1)], n_target_E)

bEHmuE_df_B = get_perc_target(lower_bin, p_B[, c(11, 15)], n_target_B)
bEHmuH_df_B = get_perc_target(lower_bin, p_B[, c(11, 13)], n_target_B)
bEHLH_df_B = get_perc_target(lower_bin, p_B[, c(11, 1)], n_target_B)

bEHmuE_df_Bd = get_perc_target(lower_bin, p_Bd[, c(11, 15)], n_target_Bd)
bEHmuH_df_Bd = get_perc_target(lower_bin, p_Bd[, c(11, 13)], n_target_Bd)
bEHLH_df_Bd = get_perc_target(lower_bin, p_Bd[, c(11, 1)], n_target_Bd)

bEHmuE_df_A = get_perc_target(lower_bin, p_A[, c(11, 15)], n_target_A)
bEHmuH_df_A = get_perc_target(lower_bin, p_A[, c(11, 13)], n_target_A)
bEHLH_df_A = get_perc_target(lower_bin, p_A[, c(11, 1)], n_target_A)

bEHmuE_df_H = get_perc_target(lower_bin, p_H[, c(11, 15)], n_target_H)
bEHmuH_df_H = get_perc_target(lower_bin, p_H[, c(11, 13)], n_target_H)
bEHLH_df_H = get_perc_target(lower_bin, p_H[, c(11, 1)], n_target_H)
"""

#plotting!
R"""
p1 = plot_heatmap(bEHmuE_df_E, c('bEH', 'muE'), '% Target Achieved, ts = E', '', limits = c(0, 0.75))
p2 = plot_heatmap(bEHmuH_df_E, c('bEH', 'muH'), '% Target Achieved, ts = E', '', limits = c(0, 0.75))
p3 = plot_heatmap(bEHLH_df_E, c('bEH', 'LH'), '% Target Achieved, ts = E', '', limits = c(0, 0.75))

p4 = plot_heatmap(bEHmuE_df_B, c('bEH', 'muE'), '% Target Achieved, ts = B', '', limits = c(0, 0.75))
p5 = plot_heatmap(bEHmuH_df_B, c('bEH', 'muH'), '% Target Achieved, ts = B', '', limits = c(0, 0.75))
p6 = plot_heatmap(bEHLH_df_B, c('bEH', 'LH'), '% Target Achieved, ts = B', '', limits = c(0, 0.75))

p4bd = plot_heatmap(bEHmuE_df_Bd, c('bEH', 'muE'), '% Target Achieved, ts = Bd', '', limits = c(0, 0.75))
p5bd = plot_heatmap(bEHmuH_df_Bd, c('bEH', 'muH'), '% Target Achieved, ts = Bd', '', limits = c(0, 0.75))
p6bd = plot_heatmap(bEHLH_df_Bd, c('bEH', 'LH'), '% Target Achieved, ts = Bd', '', limits = c(0, 0.75))

p7 = plot_heatmap(bEHmuE_df_H, c('bEH', 'muE'), '% Target Achieved, ts = H', '', limits = c(0, 0.75))
p8 = plot_heatmap(bEHmuH_df_H, c('bEH', 'muH'), '% Target Achieved, ts = H', '', limits = c(0, 0.75))
p9 = plot_heatmap(bEHLH_df_H, c('bEH', 'LH'), '% Target Achieved, ts = H', '', limits = c(0, 0.75))

p10 = plot_heatmap(bEHmuE_df_A, c('bEH', 'muE'), '% Target Achieved, ts = A', '', limits = c(0, 0.75))
p11 = plot_heatmap(bEHmuH_df_A, c('bEH', 'muH'), '% Target Achieved, ts = A', '', limits = c(0, 0.75))
p12 = plot_heatmap(bEHLH_df_A, c('bEH', 'LH'), '% Target Achieved, ts = A', '', limits = c(0, 0.75))
"""

R"""
#svg('M:/Project folders/Model env compartment/Plots/ptaplot.svg', height=20, width = 20)
ptaplot <- grid.arrange(p1,p2,p3,p4,p5,p6,p4bd,p5bd,p6bd,p7,p8,p9,p10,p11,p12, nrow = 5, ncol = 3)
#dev.off()
"""

#Conclusions 2: impact of reducing LA is low for parameter combinations of interest
R"""
bEHmuE_df_E_low_impact = get_perc_target(lower_bin, p_E[, c(11, 15)], n_lowimpact_E)
bEHmuH_df_E_low_impact = get_perc_target(lower_bin, p_E[, c(11, 13)], n_lowimpact_E)
bEHLH_df_E_low_impact = get_perc_target(lower_bin, p_E[, c(11, 1)], n_lowimpact_E)

bEHmuE_df_B_low_impact = get_perc_target(lower_bin, p_B[, c(11, 15)], n_lowimpact_B)
bEHmuH_df_B_low_impact = get_perc_target(lower_bin, p_B[, c(11, 13)], n_lowimpact_B)
bEHLH_df_B_low_impact = get_perc_target(lower_bin, p_B[, c(11, 1)], n_lowimpact_B)

bEHmuE_df_Bd_low_impact = get_perc_target(lower_bin, p_Bd[, c(11, 15)], n_lowimpact_Bd)
bEHmuH_df_Bd_low_impact = get_perc_target(lower_bin, p_Bd[, c(11, 13)], n_lowimpact_Bd)
bEHLH_df_Bd_low_impact = get_perc_target(lower_bin, p_Bd[, c(11, 1)], n_lowimpact_Bd)

bEHmuE_df_A_low_impact = get_perc_target(lower_bin, p_A[, c(11, 15)], n_lowimpact_A)
bEHmuH_df_A_low_impact = get_perc_target(lower_bin, p_A[, c(11, 13)], n_lowimpact_A)
bEHLH_df_A_low_impact = get_perc_target(lower_bin, p_A[, c(11, 1)], n_lowimpact_A)

bEHmuE_df_H_low_impact = get_perc_target(lower_bin, p_H[, c(11, 15)], n_lowimpact_H)
bEHmuH_df_H_low_impact = get_perc_target(lower_bin, p_H[, c(11, 13)], n_lowimpact_H)
bEHLH_df_H_low_impact = get_perc_target(lower_bin, p_H[, c(11, 1)], n_lowimpact_H)
"""

R"""
p13 = plot_heatmap(bEHmuE_df_E_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = E', '', limits = c(0,1))
p14 = plot_heatmap(bEHmuH_df_E_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = E', '', limits = c(0,1))
p15 = plot_heatmap(bEHLH_df_E_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = E', '', limits = c(0,1))

p16 = plot_heatmap(bEHmuE_df_B_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = B', '', limits = c(0,1))
p17 = plot_heatmap(bEHmuH_df_B_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = B', '', limits = c(0,1))
p18 = plot_heatmap(bEHLH_df_B_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = B', '', limits = c(0,1))

p16bd = plot_heatmap(bEHmuE_df_Bd_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = Bd', '', limits = c(0,1))
p17bd = plot_heatmap(bEHmuH_df_Bd_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = Bd', '', limits = c(0,1))
p18bd = plot_heatmap(bEHLH_df_Bd_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = Bd', '', limits = c(0,1))

p19 = plot_heatmap(bEHmuE_df_H_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = H', '', limits = c(0,1))
p20 = plot_heatmap(bEHmuH_df_H_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = H', '', limits = c(0,1))
p21 = plot_heatmap(bEHLH_df_H_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = H', '', limits = c(0,1))

p22 = plot_heatmap(bEHmuE_df_A_low_impact, c('bEH', 'muE'), '% Low (<2) impact, ts = A', '', limits = c(0,1))
p23 = plot_heatmap(bEHmuH_df_A_low_impact, c('bEH', 'muH'), '% Low (<2) impact, ts = A', '', limits = c(0,1))
p24 = plot_heatmap(bEHLH_df_A_low_impact, c('bEH', 'LH'), '% Low (<2) impact, ts = A', '', limits = c(0,1))
"""

R"""
#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/liplot.svg', height=20, width = 20)
liplot = grid.arrange(p13,p14,p15,p16,p17,p18,p16bd,p17bd,p18bd,p19,p20,p21,p22,p23,p24,nrow = 5, ncol = 3)
#dev.off()
"""


#Conclusion 3: increasing βEH in many cases reduces the impacts of increased LA
R"""
bEHmuE_df_E_impact = get_mean_var(lower_bin, p_E[, c(11, 15)], impact_E)
bEHmuH_df_E_impact = get_mean_var(lower_bin, p_E[, c(11, 13)], impact_E)
bEHLH_df_E_impact = get_mean_var(lower_bin, p_E[, c(11, 1)], impact_E)

bEHmuE_df_B_impact = get_mean_var(lower_bin, p_B[, c(11, 15)], impact_B)
bEHmuH_df_B_impact = get_mean_var(lower_bin, p_B[, c(11, 13)], impact_B)
bEHLH_df_B_impact = get_mean_var(lower_bin, p_B[, c(11, 1)], impact_B)

bEHmuE_df_Bd_impact = get_mean_var(lower_bin, p_Bd[, c(11, 15)], impact_Bd)
bEHmuH_df_Bd_impact = get_mean_var(lower_bin, p_Bd[, c(11, 13)], impact_Bd)
bEHLH_df_Bd_impact = get_mean_var(lower_bin, p_Bd[, c(11, 1)], impact_Bd)

bEHmuE_df_A_impact = get_mean_var(lower_bin, p_A[, c(11, 15)], impact_A)
bEHmuH_df_A_impact = get_mean_var(lower_bin, p_A[, c(11, 13)], impact_A)
bEHLH_df_A_impact = get_mean_var(lower_bin, p_A[, c(11, 1)], impact_A)

bEHmuE_df_H_impact = get_mean_var(lower_bin, p_H[, c(11, 15)], impact_H)
bEHmuH_df_H_impact = get_mean_var(lower_bin, p_H[, c(11, 13)], impact_H)
bEHLH_df_H_impact = get_mean_var(lower_bin, p_H[, c(11, 1)], impact_H)
"""

R"""
p25_1 = plot_heatmap(bEHmuE_df_E_impact, c('bEH', 'muE'), 'Mean impact, ts = E', '', limits = c(0, 0.4))
p25 = plot_heatmap(bEHmuH_df_E_impact, c('bEH', 'muH'), 'Mean impact, ts = E', '', limits = c(0, 0.4))
p26 = plot_heatmap(bEHLH_df_E_impact, c('bEH', 'LH'), 'Mean impact, ts = E', '', limits = c(0, 0.4))

p27 = plot_heatmap(bEHmuE_df_B_impact, c('bEH', 'muE'), 'Mean impact, ts = B', '', limits = c(0, 0.4))
p28 = plot_heatmap(bEHmuH_df_B_impact, c('bEH', 'muH'), 'Mean impact, ts = B', '', limits = c(0, 0.4))
p29 = plot_heatmap(bEHLH_df_B_impact, c('bEH', 'LH'), 'Mean impact, ts = B', '', limits = c(0, 0.4))

p27bd = plot_heatmap(bEHmuE_df_Bd_impact, c('bEH', 'muE'), 'Mean impact, ts = Bd', '', limits = c(0, 0.4))
p28bd = plot_heatmap(bEHmuH_df_Bd_impact, c('bEH', 'muH'), 'Mean impact, ts = Bd', '', limits = c(0, 0.4))
p29bd = plot_heatmap(bEHLH_df_Bd_impact, c('bEH', 'LH'), 'Mean impact, ts = Bd', '', limits = c(0, 0.4))

p30 = plot_heatmap(bEHmuE_df_H_impact, c('bEH', 'muE'), 'Mean impact, ts = H', '', limits = c(0, 0.4))
p31 = plot_heatmap(bEHmuH_df_H_impact, c('bEH', 'muH'), 'Mean impact, ts = H', '', limits = c(0, 0.4))
p32 = plot_heatmap(bEHLH_df_H_impact, c('bEH', 'LH'), 'Mean impact, ts = H', '', limits = c(0, 0.4))

p33 = plot_heatmap(bEHmuE_df_A_impact, c('bEH', 'muE'), 'Mean impact, ts = A', '', limits = c(0, 0.4))
p34 = plot_heatmap(bEHmuH_df_A_impact, c('bEH', 'muH'), 'Mean impact, ts = A', '', limits = c(0, 0.4))
p35 = plot_heatmap(bEHLH_df_A_impact, c('bEH', 'LH'), 'Mean impact, ts = A', '', limits = c(0, 0.4))
"""

R"""
#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/miplot.svg', height=20, width = 20)
miplot = grid.arrange(p25_1,p25,p26,p27,p28,p29,p27bd,p28bd,p29bd,p30,p31,p32,p33,p34,p35,nrow=5,ncol=3)
#dev.off()
"""

R"""
p36 = plot_heatmap_var(bEHmuH_df_E_impact, c('bEH', 'muH'), 'Impact variance, ts = E', '', limits = c(0, 0.1))
p37 = plot_heatmap_var(bEHmuH_df_E_impact, c('bEH', 'muH'), 'Impact variance, ts = E', '', limits = c(0, 0.1))
p38 = plot_heatmap_var(bEHLH_df_E_impact, c('bEH', 'LH'), 'Impact variance, ts = E', '', limits = c(0, 0.1))

p39 = plot_heatmap_var(bEHmuE_df_B_impact, c('bEH', 'muE'), 'Impact variance, ts = B', '', limits = c(0, 0.1))
p40 = plot_heatmap_var(bEHmuH_df_B_impact, c('bEH', 'muH'), 'Impact variance, ts = B', '', limits = c(0, 0.1))
p41 = plot_heatmap_var(bEHLH_df_B_impact, c('bEH', 'LH'), 'Impact variance, ts = B', '', limits = c(0, 0.1))

p39bd = plot_heatmap_var(bEHmuE_df_Bd_impact, c('bEH', 'muE'), 'Impact variance, ts = Bd', '', limits = c(0, 0.1))
p40bd = plot_heatmap_var(bEHmuH_df_Bd_impact, c('bEH', 'muH'), 'Impact variance, ts = Bd', '', limits = c(0, 0.1))
p41bd = plot_heatmap_var(bEHLH_df_Bd_impact, c('bEH', 'LH'), 'Impact variance, ts = Bd', '', limits = c(0, 0.1))

p42 = plot_heatmap_var(bEHmuE_df_H_impact, c('bEH', 'muE'), 'Impact variance, ts = H', '', limits = c(0, 0.1))
p43 = plot_heatmap_var(bEHmuH_df_H_impact, c('bEH', 'muH'), 'Impact variance, ts = H', '', limits = c(0, 0.1))
p44 = plot_heatmap_var(bEHLH_df_H_impact, c('bEH', 'LH'), 'Impact variance, ts = H', '', limits = c(0, 0.1))

p45 = plot_heatmap_var(bEHmuE_df_A_impact, c('bEH', 'muE'), 'Impact variance, ts = A', '', limits = c(0, 0.1))
p46 = plot_heatmap_var(bEHmuH_df_A_impact, c('bEH', 'muH'), 'Impact variance, ts = A', '', limits = c(0, 0.1))
p47 = plot_heatmap_var(bEHLH_df_A_impact, c('bEH', 'LH'), 'Impact variance, ts = A', '', limits = c(0, 0.1))

#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/viplot.svg', height=12, width = 20)
#viplot = grid.arrange(p36,p37,p38,p39,p40,p41,p39bd,p40bd,p41bd,p42,p43,p44,p45,p46,p47,nrow=5,ncol=3)
#dev.off()
"""

R"""
impact_bEHLA_E = get_mean_var(lower_bin, p_E_3[, c(11, 2)], impact_E_2)
impact_bEHLA_H = get_mean_var(lower_bin, p_H_3[, c(11, 2)], impact_H_2)
impact_bEHLA_A = get_mean_var(lower_bin, p_A_3[, c(11, 2)], impact_A_2)
impact_bEHLA_B = get_mean_var(lower_bin, p_B_3[, c(11, 2)], impact_B_2)
impact_bEHLA_Bd = get_mean_var(lower_bin, p_Bd_3[, c(11, 2)], impact_Bd_2)

p48 = plot_heatmap(impact_bEHLA_E, c('bEH', 'LA'), 'Mean impact, ts = E', '', limits = c(0.,0.126))#max(impact_bEHLA_E$mean)))
p49 = plot_heatmap(impact_bEHLA_B, c('bEH', 'LA'), 'Mean impact, ts = B', '', limits = c(0.,0.126))#max(impact_bEHLA_B$mean)))
p49bd = plot_heatmap(impact_bEHLA_Bd, c('bEH', 'LA'), 'Mean impact, ts = Bd', '', limits = c(0.,0.126))#max(impact_bEHLA_Bd$mean)))
p50 = plot_heatmap(impact_bEHLA_H, c('bEH', 'LA'), 'Mean impact, ts = H', '', limits = c(0.,0.126))#max(impact_bEHLA_H$mean)))
p51 = plot_heatmap(impact_bEHLA_A, c('bEH', 'LA'), 'Mean impact, ts = A', '', limits = c(0.,0.126))#max(impact_bEHLA_A$mean)))

#p52 = plot_heatmap_var(impact_bEHLA_E, c('bEH', 'LA'), 'Impact variance, ts = E', '', limits = c(0.,0.126))#max(impact_bEHLA_E$var)))
#p53 = plot_heatmap_var(impact_bEHLA_B, c('bEH', 'LA'), 'Impact variance, ts = B', '', limits = c(0.,0.126))#max(impact_bEHLA_B$var)))
#p53bd = plot_heatmap_var(impact_bEHLA_Bd, c('bEH', 'LA'), 'Impact variance, ts = B'd, '', limits = c(0.,0.126))#max(impact_bEHLA_Bd$var)))
#p54 = plot_heatmap_var(impact_bEHLA_H, c('bEH', 'LA'), 'Impact variance, ts = H', '', limits = c(0.,0.126))#max(impact_bEHLA_H$var)))
#p55 = plot_heatmap_var(impact_bEHLA_A, c('bEH', 'LA'), 'Impact variance, ts = A', '', limits = c(0.,0.126))#max(impact_bEHLA_A$var)))
"""

R"""
#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/miLAplot.svg', height=20, width = 5)
miLAplot = grid.arrange(p48, p49, p49bd, p50, p51, nrow=5)
#dev.off()

#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/viLAplot.svg', height=12, width = 7)
#viLAplot = grid.arrange(p52, p53, p54, p55, nrow=4)
#dev.off()
"""
R"""dev.off()"""

@rput p_B_9 impact_B_8
R"""
impact_lowbHA = get_mean_var(lower_bin, p_B_9[, c(11, 2)], impact_B_8)
p_lowbHA = plot_heatmap(impact_lowbHA, c('bEH', 'LA'), 'mean impact low bHA',limits = c(0., 126))
"""

#Plot just the Human dominated and environment dominated scenarios

R"""
df_E <- data.frame(param_val = c(p_E[1:10000,11],  p_E_4[1:10000,2]),
                   param = rep(c("beH", "LA"), each = 10000),
                   impact = c(impact_bEH_E[1:10000], impact_E_3[1:10000]))
df_A <- data.frame(param_val = c(p_A[1:10000,11],  p_A_4[1:10000,2]),
                  param = rep(c("beH", "LA"), each = 10000),
                  impact = c(impact_bEH_A[1:10000], impact_A_3[1:10000]))
df_B <- data.frame(param_val = c(p_B[1:10000,11],  p_B_4[1:10000,2]),
                 param = rep(c("beH", "LA"), each = 10000),
                  impact = c(impact_bEH_B[1:10000], impact_B_3[1:10000]))
df_Bd <- data.frame(param_val = c(p_Bd[1:10000,11],  p_Bd_4[1:10000,2]),
              param = rep(c("beH", "LA"), each = 10000),
              impact = c(impact_bEH_Bd[1:10000], impact_Bd_3[1:10000]))
df_H <- data.frame(param_val = c(p_H[1:10000,11],  p_H_4[1:10000,2]),
                   param = rep(c("beH", "LA"), each = 10000),
                   impact = c(impact_bEH_H[1:10000], impact_H_3[1:10000]))

df <- rbind.data.frame(df_E, df_A, df_B, df_Bd, df_H) %>%
    cbind.data.frame(., TS = rep(c("E", "A", "B", "Bd", "H"), each = 20000))
df <- mutate(df, pv_cut = cut(df$param_val, c(0, 0.095, 0.105, 0.495, 0.505, 1.0)))
df_low <- df[df$pv_cut == "(0.095,0.105]",]
df_high <- df[df$pv_cut == "(0.495,0.505]",]
print(head(df_short))

loess_plot <- function(df, ts) {
    ggplot(df, aes(param_val, impact, col = param)) +
                    geom_point(shape = 1) +
                    labs(x = "Parameter value",
                         y = "Impact",
                         col = "Parameter targeted \nby intervention",
                        title = paste("Transmission scenario: ", ts, sep = "")) +
                    geom_smooth(method = "loess", colour = "black") +
                    facet_grid(.~param) +
                    theme_bw()
}

p49 <- loess_plot(df_E, "Environment dominated")
p50 <- loess_plot(df_H, "Human dominated")
p51 <- loess_plot(df_A, "Animal dominated")
p52 <- loess_plot(df_B, "Baseline")
p52bd <- loess_plot(df_Bd, "Balanced")

#svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/bEHLAimpactcomparison.svg', height=12, width = 7)
#bEHLAimpactcomparison = grid.arrange(p49, p50, p51, p52, p52bd, nrow = 5)
#dev.off()

ggplot(df_low, aes(TS, impact, col = param)) + geom_boxplot()
ggplot(df_high, aes(TS, impact, col = param)) + geom_boxplot()
"""

#What is the range of impacts?
#Get parameter sets with +- 10% of bEH value of that transmission scenario
#bEH_E = findall((0.95 * 0.14) .< p_E[:,11] .< (1.05 * 0.14))[1:1:200]
#dat_E_bEH_target = dat_E[bEH_E,2]

#bEH_B = findall((0.95 * 0.01) .< p_B[:,11] .< (1.05 * 0.01))[1:1:200]
#dat_B_bEH_target = dat_B[bEH_B,2]

#bEH_A = findall((0.95 * 0.001) .< p_A[:,11] .< (1.05 * 0.001))[1:1:200]
#dat_A_bEH_target = dat_A[bEH_A,2]

#bEH_H = findall((0.95 * 0.001) .< p_H[:,11] .< (1.05 * 0.001))[1:1:200]
#dat_H_bEH_target = dat_H[bEH_H,2]

#Get together RH measures for these indices
@rput dat_E_7 dat_B_7 dat_Bd_7 dat_H_7 dat_A_7

R"""
RH <- data.frame(ts = c(rep('E', 10000),
                        rep('B', 10000),
                        rep('Bd', 10000),
                        rep('A', 10000),
                        rep('H', 10000)),
                    bEH = rep(c(0.14, 0.01, 0.07432092, 0.0001, 0.0001), each = 10000),
                 RH = c(dat_E_7[sample(1:190000,10000, replace = FALSE),2],
                        dat_B_7[sample(1:190000,10000, replace = FALSE),2],
                        dat_Bd_7[sample(1:190000,10000, replace = FALSE),2],
                        dat_A_7[sample(1:190000,10000, replace = FALSE),2],
                        dat_H_7[sample(1:190000,10000, replace = FALSE),2]))

svg('M:/Project\ folders\\/Model\ env\ compartment\\/Plots\\/TransmissionScenarioRHVar.svg', height=7, width = 7)
ggplot(RH, aes(y=ts, x=RH, fill = factor(bEH))) +
    geom_density_ridges(scale=1) +
    #geom_violin(outlier.shape = NA, kernel = "rectangular") +
    #geom_boxplot(width = 0.1, outlier.shape = NA) +
    geom_vline(xintercept = 0.71) +
    theme(panel.background = element_rect(fill = "white", colour = NA),
          panel.border = element_rect(fill = NA, colour= "black"),
          panel.grid = element_line(colour = NA),
          panel.grid.minor = element_line(colour = NA),
          strip.background = element_rect(fill = NA, colour = NA),
          axis.text = element_text(size = 15))
dev.off()
"""

R"dev.off()"

#compare impacts of 0.1 -> 0.0 across TS
