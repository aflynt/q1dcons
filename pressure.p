#Gnuplot script file for plotting data in file "p.dat"
# this file is called pressure.p
set autoscale
unset log
unset label
set xtic auto
set ytic auto
set title "Pressure Distribution"
set xlabel "x/L"
set ylabel "p/p0"
#set key at 0.01,100
#set label "Yield Point" at 0.003,260
#set arrow from 0.0028,250 to 0.003,260
set xr[0.0:3.0]
set yr[0.7:1.0]
plot "p0500.dat" using 1:2 title ' 500 dt' with linespoints, \
     "p1000.dat" using 1:2 title '1000 dt' with linespoints, \
     "p5000.dat" using 1:2 title '5000 dt' with linespoints
