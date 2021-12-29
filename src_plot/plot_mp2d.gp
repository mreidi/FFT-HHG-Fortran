reset

input_file="itemp"
output_file="otemp"
label_x="lbxtemp"
label_y="lbytemp"

#________________________________________________________________________________
mpl_top    = 0.3 #inch  outer top margin, title goes here
mpl_bot    = 0.5 #inch  outer bottom margin, x label goes here
mpl_left   = 1.5 #inch  outer left margin, y label goes here
mpl_right  = 1.5 #inch  outer right margin, y2 label goes here
mpl_height = 3.0 #inch  height of individual plots
mpl_width  = 3.0 #inch  width of individual plots
mpl_dx     = 0.2 #inch  inter-plot horizontal spacing
mpl_dy     = 0.2 #inch  inter-plot vertical spacing
mpl_ny     = 1   #number of rows
mpl_nx     = 1   #number of columns

# calculate full dimensions
xsize = mpl_left+mpl_right+(mpl_width*mpl_nx)+(mpl_nx-1)*mpl_dx
ysize = mpl_top+mpl_bot+(mpl_ny*mpl_height)+(mpl_ny-1)*mpl_dy

# placement functions
#   rows are numbered from bottom to top
bot(n) = (mpl_bot+(n-1)*mpl_height+(n-1)*mpl_dy)/ysize
top(n)  = 1-((mpl_top+(mpl_ny-n)*(mpl_height+mpl_dy))/ysize)
#   columns are numbered from left to right
left(n) = (mpl_left+(n-1)*mpl_width+(n-1)*mpl_dx)/xsize
right(n)  = 1-((mpl_right+(mpl_nx-n)*(mpl_width+mpl_dx))/xsize)
#__________________________________________________________________________________

#set terminal postscript eps enhanced color size xsize,ysize "Helvetica" 28
#set terminal postscript eps enhanced color size 12in,6in "Helvetica" 28
#set terminal epslatex color
set terminal pdf enhanced color size 12in,6in font "Helvetica{,28}"
#_________________________________________________________________________________

#set colorsequence default
#set colorsequence podo
#set colorsequence classic

#MATLAB

# set style line 2  lc rgb '#000090' lt 1 lw 1.5 # 
# set style line 3  lc rgb '#000fff' lt 1 lw 1.5 #  
# set style line 4  lc rgb '#0090ff' lt 1 lw 1.5 #  
# set style line 5  lc rgb '#0fffee' lt 1 lw 1.5 #  
# set style line 6  lc rgb '#90ff70' lt 1 lw 1.5 #  
# set style line 7  lc rgb '#ffee00' lt 1 lw 1.5 #  
# set style line 8  lc rgb '#ff7000' lt 1 lw 1.5 #   
# set style line 9  lc rgb '#7f0000' lt 1 lw 1.5 #  

# lc : line color
# dt : dash type
# lw : line width
# ls : line style

set style line 2   lc rgb '#000000' dt 1 lw 4.0 #  black
set style line 3   lc rgb '#8A2BE2' dt 1 lw 4.0 #  blueviolet
set style line 4   lc rgb '#DC143C' dt 1 lw 4.0 #  crimson
set style line 5   lc rgb '#006400' dt 1 lw 4.0 #  darkgreen
set style line 6   lc rgb '#1E90FF' dt 1 lw 4.0 #  dodgerblue
set style line 7   lc rgb '#FF7F50' dt 1 lw 4.0 #  coral
set style line 8   lc rgb '#00FFFF' dt 1 lw 4.0 #  cyan
set style line 9   lc rgb '#B8860B' dt 1 lw 4.0 #  darkgoldenrod
set style line 10  lc rgb '#FF1493' dt 1 lw 4.0 #  deeppink
set style line 11  lc rgb '#2F4F4F' dt 1 lw 4.0 #  darkslategray
set style line 12  lc rgb '#B22222' dt 1 lw 4.0 #  firebrick
set style line 13  lc rgb '#FFD700' dt 1 lw 4.0 #  gold
set style line 14  lc rgb '#FF4500' dt 1 lw 4.0 #  orangered
set style line 15  lc rgb '#4169E1' dt 1 lw 4.0 #  royalblue
set style line 16  lc rgb '#D2B48C' dt 1 lw 4.0 #  tan
set style line 17  lc rgb '#9ACD32' dt 1 lw 4.0 #  yellowgreen
set style line 18  lc rgb '#008080' dt 1 lw 4.0 #  teal
set style line 19  lc rgb '#4B0082' dt 1 lw 4.0 #  indigo
set style line 20  lc rgb '#FF8C00' dt 1 lw 4.0 #  darkorange
set style line 21  lc rgb '#F4A460' dt 1 lw 4.0 #  sandybrown

#_____________________________________________________________________________

set output output_file


# stats input_file using 1 noout
# xmin=STATS_min
# xmax=STATS_max
# set xrange[xmin:xmax]

set xrange[*:*]

# stats input_file using nc noout
# ymin=STATS_min
# ymax=STATS_max
# set yrange[ymin:ymax]

set yrange[*:*]


# to get number of columns from the second row
stats input_file every ::2 noout # no output
nc = STATS_columns

# nc=system("head -1 input_file | wc -w") : cannot get input_file string from gnuplot

set xlabel label_x
set ylabel label_y
#set key off
set key right center
set key outside
set key autotitle columnheader

# plot for [i=2:nc] input_file u 1:i smooth csplines ls i
plot input_file u i:nc with line 

#___________________________________________________________________
