set terminal  pngcairo size 1400, 1000

#Variable Initialization
#Sets analysis data file to use, sets title, sets output

#PATH = "C:\\Users\\georgeau\\Box Sync\\August\\ZuckermanLab\\Proof-maker\\RUNS\\DEBUG_NSTEPS1e4_ALPHA1_SEED456789"
PATH = "C:\\Users\\georgeau\\Box Sync\\August\\ZuckermanLab\\Proof-maker\\RUNS\\proof_maker_alpha0_s987654\\analysis_n5000"


DATAFILE1 = PATH."\\analysis-vary_dmu_N-dmu_init__-6__to__dmu_fin__0"
#DATAFILE2 = PATH."\\evolver_rates2.dat"
#DATAFILE3 = PATH."\\evolver_rates3.dat"
#DATAFILE4 = PATH."\\evolver_rates_s9999999.dat"

NEGATIVE_N = 2
NEGATIVE_S = 0
NEGATIVE_W = -6

set arrow from NEGATIVE_N, graph 0 to NEGATIVE_N, graph 1 heads lw 2 lc rgb "black"
set arrow from NEGATIVE_S, graph 0 to NEGATIVE_S, graph 1 heads lw 2 lc rgb "orange"
set arrow from NEGATIVE_W, graph 0 to NEGATIVE_W, graph 1 heads lw 2 lc rgb "red"

#TITLE = "Debugging proof-maker: MC Energy vs MC n \n{/*0.85 MC Energy = -sflow*|sflow/wflow|^a^l^p^h^a} \n{/*0.85 seed=456789, alpha=1, n=1e4, nbeta=3e2}"
TITLE = "Proof-Maker Model: Substrate, Toxin, and Sodium Flows \n{/*0.85 MC Energy = -sflow*|sflow/wflow|^a^l^p^h^a} \n{/*0.85 alpha=0, seed=987654, n steps=1e4}"
OUTPUT = PATH."\\analysis_graph_a0_s987654.png"

set output OUTPUT
set title TITLE
set xlabel "dMu Sodium (kT)"
set ylabel "Flow Out->In (s^-^1)"
#set xrange [-1.0:1.0]
set yrange [-2:2]
set zeroaxis
set grid

#plot DATAFILE1 u ($1):($2) w lines t 'MC Model Default' lt -1 lw 2 lc rgb "black"

plot DATAFILE1 u ($1):($3/$2) w lines t 'Substrate:Sodium' lt -1 lw 2 lc rgb "purple",\
DATAFILE1 u ($1):($4/$2) w lines t 'Toxin:Sodium' lt -1 lw 2 lc rgb "green",\
DATAFILE1 u ($1):($3/$4) w lines t 'Substrate:Toxin' lt -1 lw 2 lc rgb "blue",\
1/0 t 'Negative Sodium Flow' lw 2 lc rgb "black", 1/0 t 'Negative Substrate Flow' lw 2 lc rgb "orange",\
1/0 t 'Negative Toxin Flow' lw 2 lc rgb "red"

##DATAFILE4 u ($1):($2) w lines t 'S=9999999' lt -1 lw 2 lc rgb "red"
