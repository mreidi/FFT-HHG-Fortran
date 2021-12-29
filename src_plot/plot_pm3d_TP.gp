reset

input_file="itemp"
output_file="otemp"
label_x="lbxtemp"
label_y="lbytemp"

#__________________________________________________________________________________
set terminal pdf enhanced color size 10in,6in font "Helvetica{,28}"
#__________________________________________________________________________________
#set palette rgbformulae 7,5,15  # "pm3d default"
# set palette rgbformulae 33,13,10  # "rainbow"
#set palette rgbformulae 3,11,6  # "Green_Red-Violet"
#set palette rgbformulae 23,28,3  # "Ocean"
#set palette rgbformulae 30,31,32  # "Color printable on gray"
# set palette rgbformulae 21,22,23  # "Hot"
set palette defined 
#__________________________________________________________________________________


set output output_file

nc =system("head -1 itemp | wc -w")
nr=system("wc -l itemp")

# stats input_file every ::2::nr using 1 noout
# xmin=STATS_min
# xmax=STATS_max

set xrange[xmin:xmax]
set yrange[ymin:ymax]

# for an automatic range of colorbar cbrange[*:*] 
set cbrange[cbmin:cbmax] 

set title "title_temp"

set xlabel label_x
set ylabel label_y
set key autotitle columnheader

set multiplot
set pm3d map
# plot input_file nonuniform matrix every 1 using 2:1:(log10($3))   w image t ""
plot input_file nonuniform matrix every 1 using 2:1:3  w image t ""
unset multiplot

# The data organization for non-uniform matrix input is
# 
# <N+1>  <x0>   <x1>   <x2>  ...  <xN>
#  <y0> <z0,0> <z0,1> <z0,2> ... <z0,N>
#  <y1> <z1,0> <z1,1> <z1,2> ... <z1,N>
#   :      :      :      :   ...    :
# which is then converted into triplets:
# 
# <x0> <y0> <z0,0>
# <x0> <y1> <z0,1>
# <x0> <y2> <z0,2>
#  :    :     :
# <x0> <yN> <z0,N>
# <x1> <y0> <z1,0>
# <x1> <y1> <z1,1>
#  :    :     :

# --------------------------------------------------------------------
# every {<column_incr>}
#         {:{<row_incr>}
#             {:{<start_column>}
#                 {:{<start_row>}
#                     {:{<end_column>}
#                         {:<end_row>}}}}}
# Examples:
# every :::N::N     # plot all values in row N of matrix
# every ::3::7:     # plot columns 3 to 7 of all rows
# every ::3:2:7:4   # submatrix of columns 3-7 rows 2-4
#________________________________________________________________________________
