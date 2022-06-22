# set title "step through volume"
# rlow=-2
# rhigh=2
# set xrange [rlow:rhigh]; set yrange [rlow:rhigh]; set zrange [rlow:rhigh]
# set border 0
# set xyplane at 0
# set xzeroaxis lt -1; set yzeroaxis lt -1; set zzeroaxis lt -1;
# set view 68, 28, 1.6, 0.9
# set sample 50,50; set isosample 50,50
# set cbrange [1:3]

# set palette cubehelix negative
# set log cb; unset cbtics;
# set cblabel "point density"

# sp "3d_test.txt" u 1:2:3:4 w points pt 7 ps 0.5 lc palette

# set pm3d
# set view map
# set contour
# set cntrparam levels 15
# unset key
# set cbrange[0:1.6]
plot "parallel_2d_error_40.txt" w l

# Two dimension plot (pm3d)
# set view map
# set cntrparam levels 100
# unset key
# set samples 20, 20
# set isosamples 21, 21
# set zrange [ 0 : 1.6 ] noreverse writeback
# set cbrange [0:1.4]
# set title "The contour plot"
# set palette rgb 7,5,15
# splot "cos_test.txt" w pm3d

# Two dimension plot (pm3d)
# set pm3d at s
# set cntrparam levels 15
# set cntrparam levels incremental -1, 0.05, 1
# set contour base
# set view 0,0
# unset key
# unset surface
# set samples 20, 20
# set isosamples 21, 21
# set zrange [ 0 : 100]
# set cbrange [0:60]
# set title "The contour plot"
# set palette rgb 7,5,15
# set arrow from 0,0 to 0.1,0.1
# splot "cos_test.txt" w l, "origin.txt" u 1:2:(0) w points pt 5