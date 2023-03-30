# Example command line: 
# $ gnuplot -e "nTargets=11" -e "logXY=1" -e "noKey=1" ~/programs/gnuplotScripts/bestFitTargetData.gnu
# 
# Example plot: p for [i=2:12] 'targetData.dat' u 1:i w p lc i, for [i=2:12] 'bestFitData.dat' u 1:i w l lc i
# 
set terminal qt size 600,400 font 'Arial, 14'
set key font 'Arial, 10'
if (exists("noKey")) { unset key; }
if (exists("logXY")) { set log xy; set xrange [0.5:]; }
if (exists("logX")) { set log x; set xrange [0.5:]; }
if (exists("logY")) { set log y; }
if (exists("yLabel")) { set ylabel yLabel; } else { set xlabel 'x' }
if (exists("xLabel")) { set xlabel xLabel; } else { set ylabel 'f(x)' }
if (exists("xMin")) { set xrange [xMin:]; }
if (exists("xMax")) { set xrange [:xMax]; }
if (exists("yMin")) { set yrange [yMin:]; }
if (exists("yMax")) { set yrange [:yMax]; }
p for [i=2:nTargets+1] 'targetData.dat' u 1:i w p lc i, for [i=2:nTargets+1] 'bestFitData.dat' u 1:i w l lc i
pause -1
