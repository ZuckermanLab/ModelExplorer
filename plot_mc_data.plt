
set terminal  pngcairo size 1400, 1000

#Variable Initialization
#Sets analysis data file to use, sets title, sets output

#PATH = "C:\\Users\\georgeau\\Box Sync\\August\\ZuckermanLab\\Proof-maker\\RUNS\\DEBUG_NSTEPS1e4_ALPHA1_SEED456789"
PATH = "C:\\Users\\georgeau\\Box Sync\\August\\ZuckermanLab\\Proof-maker\\RUNS\\proof_maker_alpha0_5_s654321_1e5"

DATAFILE1 = PATH."\\evolver_rates.dat"
#DATAFILE2 = PATH."\\evolver_rates2.dat"
#DATAFILE3 = PATH."\\evolver_rates3.dat"
#DATAFILE4 = PATH."\\evolver_rates_s9999999.dat"

#TITLE = "Proof-maker: MC Energy vs MC n \n{/*0.85 MC Energy = -sflow*|sflow/wflow|^a^l^p^h^a} \n{/*0.85 seed=456789, alpha=1, n=1e4}"
TITLE = "Proof-Maker: MC Energy Landscape \n{/*0.85 MC Energy = -sflow*|sflow/wflow|^a^l^p^h^a} \n{/*0.85 alpha=0.5, seed=654321, n steps=1e5}"
OUTPUT = PATH."\\mc_graph_a0_5.png"

set output OUTPUT
set title TITLE
set xlabel "MC n"
set ylabel "MC Energy"
#set xrange [-1.0:1.0]
#set yrange [-1e-04:1e-04]
set rmargin 5
set zeroaxis
set grid

plot DATAFILE1 u ($1):($2) w lines t 'MC Run 1' lt -1 lw 2 lc rgb "magenta"

#plot DATAFILE1 u ($1):($2) w lines t 'MC Run1' lt -1 lw 2 lc rgb "purple",\
#DATAFILE2 u ($1):($2) w lines t 'MC Run2' lt -1 lw 2 lc rgb "green",\
#DATAFILE3 u ($1):($2) w lines t 'MC Run3' lt -1 lw 2 lc rgb "blue"#,\
##DATAFILE4 u ($1):($2) w lines t 'S=9999999' lt -1 lw 2 lc rgb "red"
