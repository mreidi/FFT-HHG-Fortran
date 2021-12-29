reset

input_file="itemp"
output_file="otemp"
label_x="lbxtemp"
label_y="lbytemp"

nf=nu_frames   # number of frames

#__________________________________________________________________________________
set terminal gif animate delay 1.0 
#__________________________________________________________________________________
#set palette defined 
#set palette rgbformulae 7,5,15  # "pm3d default"
# set palette rgbformulae 33,13,10  # "rainbow"
#set palette rgbformulae 3,11,6  # "Green_Red-Violet"
#set palette rgbformulae 23,28,3  # "Ocean"
#set palette rgbformulae 30,31,32  # "Color printable on gray"
set palette rgbformulae 21,22,23  # "Hot"
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

set cbrange[0:1]    # cbrange[*:*] for an automatic range of colorbar

set title "title_temp"
set xlabel label_x
set ylabel label_y
set key autotitle columnheader

stats input_file every ::2 noout # no output
nc = STATS_columns

nr=system("wc -l itemp")  # (STATS_records/nc)nc	

set view map
do for [i=3:nc:nf] {
    
#     stats input_file using i noout
#     cbmin=STATS_min
#     cbmax=STATS_max
#     set cbrange[cbmin:cbmax]
    
    set multiplot
    splot input_file every ::2::nr u 1:2:i with points palette pt 9
    unset multiplot
}

#________________________________________________________________________________
