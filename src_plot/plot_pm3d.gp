reset

input_file="itemp"
output_file="otemp"
nc=column_number
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
#set terminal postscript eps enhanced color size 8in,6in "Helvetica" 28
#set terminal epslatex color
set terminal pdf enhanced color size 10in,6in font "Helvetica{,28}"
#__________________________________________________________________________________

set output output_file


stats input_file using 1 noout
xmin=STATS_min
xmax=STATS_max
set xrange[xmin:xmax]

stats input_file using 2 noout
ymin=STATS_min
ymax=STATS_max
set yrange[ymin:ymax]

# stats input_file using nc noout
# cbmin=STATS_min
# cbmax=STATS_max
# set cbrange[cbmin:cbmax]

set cbrange[*:*]    # cbrange[a:b] for a fixed range of colorbar
    
set title "title_temp"
set xlabel label_x
set ylabel label_y
set key off

#set palette defined 
#set palette rgbformulae 7,5,15  # "pm3d default"
# set palette rgbformulae 33,13,10  # "rainbow"
#set palette rgbformulae 3,11,6  # "Green_Red-Violet"
#set palette rgbformulae 23,28,3  # "Ocean"
#set palette rgbformulae 30,31,32  # "Color printable on gray"
set palette rgbformulae 21,22,23  # "Hot"

#set pm3d map
set view map
splot input_file u 1:2:nc with points palette pt 9

#splot input_file1 every ::1::10001 u 1:2:3 with points palette
#________________________________________________________________________________
