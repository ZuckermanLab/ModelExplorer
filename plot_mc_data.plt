
set terminal  pngcairo size 1400, 1000


PATH = "C:\\Users\\georgeau\\Desktop\\runs\\testing_dgsw2_1e5_a2_s456789"

DATAFILE1 = PATH."\\evolver_rates.dat"

TITLE = "Model-Space Explorer: MC Energy Trajectory \n{/*0.85 MC Energy = -sflow*|sflow/wflow|^a^l^p^h^a} \n{/*0.85 alpha = 2, seed = 456789, n steps = 1e5, dmu w = 2, Na first constraint off, dg sw = 2}"
OUTPUT = PATH."\\mc_graph.png"

set output OUTPUT
set title TITLE
set xlabel "MC n"
set ylabel "MC Energy"
#set xrange [-1.0:1.0]
#set yrange [-1e-04:1e-04]
set rmargin 5
set zeroaxis
set grid

plot DATAFILE1 u ($1):($2) w lines t 'MC Run 1' lt -1 lw 2 lc rgb "purple"
