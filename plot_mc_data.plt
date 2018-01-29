
set terminal  pngcairo size 1400, 1000


PATH = "C:\\Users\\georgeau\\Desktop\\Zuckerman_Lab\\Proof_maker_runs\\RUNS\\proof_maker_alpha0_s901234_n1e5_dmuW-2_na_first_off"

DATAFILE1 = PATH."\\evolver_rates.dat"
#DATAFILE2 = PATH."\\evolver_rates2.dat"
#DATAFILE3 = PATH."\\evolver_rates3.dat"
#DATAFILE4 = PATH."\\evolver_rates_s9999999.dat"

TITLE = "Proof-Maker: MC Energy Landscape \n{/*0.85 MC Energy = -sflow*|sflow/wflow|^a^l^p^h^a} \n{/*0.85 alpha = 0, seed = 901234, n steps = 1e5, dmu_w = -2, Na first constraint off}"
OUTPUT = PATH."\\mc_graph_a0_dmuw-2_nafirstoff.png"

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
