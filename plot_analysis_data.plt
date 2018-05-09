set terminal pngcairo enhanced size 1600, 1200 color dashed

PATH = "C:\\Users\\georgeau\\Desktop\\debug\\testing_dgsw2\\n9000"

DATAFILE1 = PATH."\\analysis-vary_dmu_N-dmu_init__-10__to__dmu_fin__0"
#DATAFILE2 = PATH."\\evolver_rates2.dat"
#DATAFILE3 = PATH."\\evolver_rates3.dat"
#DATAFILE4 = PATH."\\evolver_rates_s9999999.dat"

NEGATIVE_N = -0.6
NEGATIVE_S = -3.6
NEGATIVE_W = 1

#for relative flows
#set arrow from NEGATIVE_N, graph 0 to NEGATIVE_N, graph 1 heads lt 2 lw 2 lc rgb "pink"
#set arrow from NEGATIVE_S, graph 0 to NEGATIVE_S, graph 1 heads lw 2 lc rgb "purple"
#set arrow from NEGATIVE_W, graph 0 to NEGATIVE_W, graph 1 heads lw 2 lc rgb "orange"

TITLE = "{/:Bold Transporter: Substrate, Toxin, and Sodium Flows} \n {/*0.85 dmu-W = 2, dg-SW = 2, Na-binds-first = FALSE, Vary dMu-N} \n {/*0.85 alpha = 0, seed = 456789, n-steps = 1e4, n = 9000}"
OUTPUT = PATH."\\analysis_graph.png"


set output OUTPUT
set title TITLE
set title font "Arial,26"
#set title format "{/:Bold}"

set border 31 lw 4
set lmargin screen 0.1
set bmargin screen 0.075

set xlabel "{/:Bold d{/Symbol m} (Chemical Potential) Na^+ [kT] }"
set xlabel font "Arial,22"



set ylabel "{/:Bold Flow (Into Cell) [Ms^-^1] }" offset -2
set ylabel font "Arial,22"

set key font "Arial,20"
set key right top
set ytics format "{/Arial:Bold {/=18 %0.1t{/Symbol \264}10^{%T}}}"
set xtics format "{/Arial:Bold {/=18 %h}}"

#set yrange [-2:2]
set zeroaxis
set grid

## for relative flows
#plot DATAFILE1 u ($1):($3/$2) w lines t "{/:Bold Substrate:Sodium }" lt -1 lw 4 lc rgb "green",\
#DATAFILE1 u ($1):($4/$2) w lines t "{/:Bold Toxin:Sodium }" lt -1 lw 4 lc rgb "red",\
#DATAFILE1 u ($1):($4/$3) w lines t "{/:Bold Toxin:Substrate }" dashtype 2 lw 4 lc rgb "blue"#,\
##1/0 t 'Negative Sodium Flow' lw 2 lc rgb "pink", 1/0 t 'Negative Substrate Flow' lw 2 lc rgb "purple",\
##1/0 t 'Negative Toxin Flow (not in graph)' lw 2 lc rgb "orange"

## for absolute flows
plot DATAFILE1 u ($1):($2) w lines t "{/:Bold Sodium }" dashtype 2 lw 4 lc rgb "#40000000",\
DATAFILE1 u ($1):($3) w lines t '{/:Bold Substrate }' dashtype 1 lw 4 lc rgb "#4000FF00",\
DATAFILE1 u ($1):($4) w lines t '{/:Bold Toxin}' dashtype 4 lw 4 lc rgb "#40FF0000"
