set terminal pngcairo enhanced size 1600, 1800 color dashed

PATH = ".\\runs\\na_first\\na_first_a1_dgsw1_s123456\\99999"
DATAFILE1 = PATH."\\analysis-vary_dmu_N-dmu_init__-10__to__dmu_fin__0"
TITLE = "{/:Bold Transporter: Substrate, Toxin, and Sodium Flows} \n {/*0.85 dmu-W = 2, dg-SW = 1, Na-binds-first = TRUE, Vary dMu-N} \n {/*0.85 alpha = 1.0 seed = 123456, n-steps = 1e5, n = 99999}"
OUTPUT = PATH."\\analysis_graph.png"

set output OUTPUT
set tmargin 5
set lmargin 18
set multiplot layout 2, 1 title TITLE font "Arial:Bold,26"
set ytics format "{/Arial {/=18 %0.1t{/Symbol \264}10^{%T}}}"
set xtics format "{/Arial {/=18 %h}}"

## for absolute flows
set title "Absolute Flows"
set title font "Arial:Bold,22"
set xlabel "d{/Symbol m} (Chemical Potential) Na^+ [kT]"
set xlabel font "Arial,22"
set ylabel "Flow (Into Cell) [Ms^-^1] " offset -2
set ylabel font "Arial,22"
set key font "Arial,20"
set key right top
set zeroaxis
set grid
#set yrange [-2.5e-5:2.5e-5]
set autoscale y

plot DATAFILE1 u ($1):($2) w lines t "Sodium" dashtype 2 lw 4 lc rgb "#40000000",\
DATAFILE1 u ($1):($3) w lines t 'Substrate' dashtype 1 lw 4 lc rgb "#4000FF00",\
DATAFILE1 u ($1):($4) w lines t 'Toxin' dashtype 4 lw 4 lc rgb "#40FF0000"


## for relative flows
set title "Relative Flows"
set title font "Arial:Bold,22"
set xlabel "d{/Symbol m} (Chemical Potential) Na^+ [kT]"
set xlabel font "Arial,22"
set ylabel "Relative Flow (Into Cell) [Ms^-^1]" offset -2
set ylabel font "Arial,22"
set key font "Arial,20"
set key right bot
set zeroaxis
set grid
set yrange [0:25]
#set autoscale y

plot DATAFILE1 u ($1):(abs($3/$2)) w lines t "|S:N|" dashtype 2 lw 4 lc rgb "#40000000",\
DATAFILE1 u ($1):(abs($3/$4)) w lines t "|S/W|" dashtype 1 lw 4 lc rgb "#4000FF00"
unset multiplot
