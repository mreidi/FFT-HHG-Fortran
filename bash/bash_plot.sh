#!/bin/sh


    inp_path="inp/"
    out_path="out/"
    src_plot_path="src_plot/"
    plot_path="plot/"

    input="1_case.inp"
        
    i_read=0 
    while read line
    do 
        ((i_read=i_read+1))
        
        if [ $i_read -eq 4 ]
        then 
            HHG_case=$(echo "$line" | cut -c1-1)
        elif [ $i_read -eq 5 ]
        then 
            amp_spec_case=$(echo "$line" | cut -c1-1)
        elif [ $i_read -eq 6 ]
        then 
            phase_spec_case=$(echo "$line" | cut -c1-1)
        elif [ $i_read -eq 7 ]
        then 
            att_pulse_case=$(echo "$line" | cut -c1-1)
        elif [ $i_read -eq 8 ]
        then 
            Time_profile_case=$(echo "$line" | cut -c1-1)
        
        elif [ $i_read -eq 9 ]
        then 
            harm_scan_case=$(echo "$line" | cut -c1-1)
        fi 
    done < "$inp_path""$input"
    
    input="2_parameter.inp"
    i_read=0 
    while read line
    do 
        ((i_read=i_read+1))
        
        if [ $i_read -eq 3 ]
        then 
            func_col_ini=$(echo "$line" | cut -c1-3)
        elif [ $i_read -eq 4 ]
        then 
            func_col_fin=$(echo "$line" | cut -c1-3)
        fi 
    done < "$inp_path""$input"
    
    
#------------------------------- Batch Plot --------------------------------

    
    if [ $func_col_fin -gt $func_col_ini ]
    then
        j=0
        for i in `ls "$out_path"Batch*.gnumeric | sort -V` 
        do  
            ((j=j+1))
            
            plot_file="$src_plot_path"plot_mp2d.gp
            pl_temp="$plot_path"plot_mp2d_$j.gp
            
            cp -f $plot_file $pl_temp
            
            if [ $j -eq 1 ]
            then 
                lbx='Harmonic order'  
                lby='Log_{10}(A({/Symbol w}))'  
            elif [ $j -eq 2 ]
            then 
                lbx='Harmonic order'  
                lby='{/Symbol f}({/Symbol w}) (Degree)'
            elif [ $j -eq 3 ]
            then 
                lbx='Harmonic order'  
                lby='Log_{10}(S({/Symbol w}))'
            elif [ $j -eq 4 ]
            then 
                lbx='Time (Laser cycle)'   
                lby='SAP'
            fi
                        
            sed -i "s=otemp=$i=g"  $pl_temp
            sed -i "s="$out_path"="$plot_path"=g"  $pl_temp
            sed -i "s=.gnumeric=.pdf=g" $pl_temp
            sed -i "s=itemp=$i=g"  $pl_temp
            
            sed -i "s=lbxtemp=$lbx=g"  $pl_temp
            sed -i "s=lbytemp=$lby=g"  $pl_temp
            
            gnuplot $pl_temp
            
            rm -rf $pl_temp
        done
    fi

#-------------------------- HHG: Power spectrum S(w) -----------------------------

    x_min="*"
    x_max="*"
    y_min="*"
    y_max="*"

    
    if [ $HHG_case -eq 1 ]
    then
        
        ename="Power_Spectrum"
            
        lbx='Harmonic order'  
        lby='Log_{10}(S({/Symbol w}))'
        
        t_col=1
        HHG_col=2
        
        j=0
    
        for i in "$out_path"Out_Spectra_*.gnumeric
        do  
            ((j=j+1))
                        
            plot_file="$src_plot_path"plot_2d.gp
            pl_temp="$plot_path"plot_2d_$j.gp
            
            cp -f $plot_file $pl_temp
            
            echo 

            sed -i "s=otemp=$i=g"  $pl_temp
            sed -i "s="$out_path"="$plot_path"=g"  $pl_temp
            sed -i "s=Out_Spectra=$ename=g"  $pl_temp
            sed -i "s=.gnumeric=.pdf=g" $pl_temp
            sed -i "s=itemp=$i=g"  $pl_temp
            sed -i "s=xmin=$x_min=g"  $pl_temp
            sed -i "s=xmax=$x_max=g"  $pl_temp
            sed -i "s=ymin=$y_min=g"  $pl_temp
            sed -i "s=ymax=$y_max=g"  $pl_temp
            sed -i "s=t_column=$t_col=g"  $pl_temp
            sed -i "s=f_column=$HHG_col=g"  $pl_temp
            
            sed -i "s=lbxtemp=$lbx=g"  $pl_temp
            sed -i "s=lbytemp=$lby=g"  $pl_temp
            
            gnuplot $pl_temp
            
#             rm -rf $pl_temp
        done
    fi
#-------------------------- HHG: Amplitude spectrum A(w) -----------------------------

    x_min="*"
    x_max="*"
    y_min="*"
    y_max="*"
    
    
    if [ $amp_spec_case -eq 1 ]
    then
        
        ename="Amplitude_Spectrum"
        
        lbx='Harmonic order'  
        lby='Log_{10}(A({/Symbol w}))'
        
        t_col=1
        HHG_col=3
        
        j=0
   
        for i in "$out_path"Out_Spectra_*.gnumeric
        do  
            ((j=j+1))
                        
            plot_file="$src_plot_path"plot_2d.gp
            pl_temp="$plot_path"plot_2d_$j.gp
            
            cp -f $plot_file $pl_temp

            sed -i "s=otemp=$i=g"  $pl_temp
            sed -i "s="$out_path"="$plot_path"=g"  $pl_temp
            sed -i "s=Out_Spectra=$ename=g"  $pl_temp
            sed -i "s=.gnumeric=.pdf=g" $pl_temp
            sed -i "s=itemp=$i=g"  $pl_temp
            sed -i "s=xmin=$x_min=g"  $pl_temp
            sed -i "s=xmax=$x_max=g"  $pl_temp
            sed -i "s=ymin=$y_min=g"  $pl_temp
            sed -i "s=ymax=$y_max=g"  $pl_temp
            sed -i "s=t_column=$t_col=g"  $pl_temp
            sed -i "s=f_column=$HHG_col=g"  $pl_temp
            
            sed -i "s=lbxtemp=$lbx=g"  $pl_temp
            sed -i "s=lbytemp=$lby=g"  $pl_temp
            
            gnuplot $pl_temp
            
            rm -rf $pl_temp
        done
    fi

#-------------------------- HHG: Phase spectrum Phi(w) -----------------------------

    x_min="*"
    x_max="*"
    y_min="*"
    y_max="*"

    if [ $phase_spec_case -eq 1 ]
    then
        
            
        lbx='Harmonic order'  
        lby='{/Symbol f}({/Symbol w}) (Degree)'
        
        t_col=1
        HHG_col=4
        
        j=0
        
        ename="Phase_Spectrum"
    
        for i in "$out_path"Out_Spectra_*.gnumeric
        do  
            ((j=j+1))
                        
            plot_file="$src_plot_path"plot_2d.gp
            pl_temp="$plot_path"plot_2d_$j.gp
            
            cp -f $plot_file $pl_temp

            sed -i "s=otemp=$i=g"  $pl_temp
            sed -i "s="$out_path"="$plot_path"=g"  $pl_temp
            sed -i "s=Out_Spectra=$ename=g"  $pl_temp
            sed -i "s=.gnumeric=.pdf=g" $pl_temp
            sed -i "s=itemp=$i=g"  $pl_temp
            sed -i "s=xmin=$x_min=g"  $pl_temp
            sed -i "s=xmax=$x_max=g"  $pl_temp
            sed -i "s=ymin=$y_min=g"  $pl_temp
            sed -i "s=ymax=$y_max=g"  $pl_temp
            sed -i "s=t_column=$t_col=g"  $pl_temp
            sed -i "s=f_column=$HHG_col=g"  $pl_temp
            
            sed -i "s=lbxtemp=$lbx=g"  $pl_temp
            sed -i "s=lbytemp=$lby=g"  $pl_temp
            
            gnuplot $pl_temp
            
            rm -rf $pl_temp
        done
    fi

#------------------------------- Harm_Scan --------------------------------
    
    if [ $harm_scan_case -gt 0 ]
    then 
        
        lbytemp='Log_{10}(S{/Symbol w}))'
        
        j=1
        for i in `ls "$out_path"Batch_Scan_*.gnumeric | sort -V` 
        do  
            
            plot_file="$src_plot_path"plot_mp2d.gp
            pl_temp="$plot_path"plot_mp2d_$j.gp
            
            cp -f $plot_file $pl_temp
            
            if [ $harm_scan_case -eq 1 ]
            then 
                lbxtemp='Gap (eV)'
            elif [ $harm_scan_case -eq 2 ]
            then 
                lbxtemp='ecc'
            elif [ $harm_scan_case -eq 3 ]
            then 
                if [ $j -eq 1 ]
                then
                    lbxtemp='E_0 (V/Å)'
                elif [ $j -eq 2 ]
                then
                    lbxtemp='Intensity (W/{cm}^2)'
                fi
            elif [ $harm_scan_case -eq 4 ]
            then 
                lbxtemp='CEP (Degree)'
            elif [ $harm_scan_case -eq 5 ]
            then 
                lbxtemp='T1'
            elif [ $harm_scan_case -eq 6 ]
            then 
                lbxtemp='T2'
            elif [ $harm_scan_case -eq 7 ]
            then 
                lbxtemp='Pu-Pr-Delay'
            fi

            sed -i "s=otemp=$i=g"  $pl_temp
            sed -i "s="$out_path"="$plot_path"=g"  $pl_temp
            sed -i "s=.gnumeric=.pdf=g" $pl_temp
            sed -i "s=itemp=$i=g"  $pl_temp
            
            sed -i "s=lbxtemp=$lbxtemp=g"  $pl_temp
            sed -i "s=lbytemp=$lbytemp=g"  $pl_temp
            
            gnuplot $pl_temp
            
            rm -rf $pl_temp
            
            ((j=j+1))
        done
    fi
    
#-------------------------- Atto -----------------------------

    x_min="*"
    x_max="*"
    y_min="*"
    y_max="*"

    
    if [ $att_pulse_case -eq 1 ]
    then
        
        t_col=3
        f_col=4
        
        lbx='Time (Laser cycle)'   
        lby='SAP'
        
        j=0
        
        for i in "$out_path"Out_atto_*.gnumeric
        do  
            ((j=j+1))
            
            plot_file="$src_plot_path"plot_2d.gp
            pl_temp="$plot_path"plot_2d_$j.gp
            
            cp -f $plot_file $pl_temp
      
            sed -i "s=otemp=$i=g"  $pl_temp
            sed -i "s="$out_path"="$plot_path"=g"  $pl_temp
            sed -i "s=.gnumeric=.pdf=g" $pl_temp
            sed -i "s=itemp=$i=g"  $pl_temp
            sed -i "s=xmin=$x_min=g"  $pl_temp
            sed -i "s=xmax=$x_max=g"  $pl_temp
            sed -i "s=ymin=$y_min=g"  $pl_temp
            sed -i "s=ymax=$y_max=g"  $pl_temp
            sed -i "s=t_column=$t_col=g"  $pl_temp
            sed -i "s=f_column=$f_col=g"  $pl_temp
            
            sed -i "s=lbxtemp=$lbx=g"  $pl_temp
            sed -i "s=lbytemp=$lby=g"  $pl_temp
            
            gnuplot $pl_temp
            
            rm -rf $pl_temp
        done
    fi    

#------------------------------ Time frequency profile -------------------------------
    
    x_min="*"
    x_max="*"
    y_min="*"
    y_max="*"
    cb_min="*"
    cb_max="*"

    
    if [ $Time_profile_case -gt 0 ]
    then 
        
        lbx='Time (Laser cycle)'  
        lby='Harmonic order'    
        
        
        if [ $Time_profile_case -eq 1 ]
        then 
            tit='Log_{10}(S({/Symbol w}))'
        elif [ $Time_profile_case -eq 2 ]
        then 
            tit='Log_{10}(A({/Symbol w}))'
        elif [ $Time_profile_case -eq 3 ]
        then 
            tit='{/Symbol f}({/Symbol w}) (Degree)'
        fi
        
        j=0
        
        for i in "$out_path"Out_TP_Gabor_*.gnumeric
        do  
            ((j=j+1))
            
            plot_file="$src_plot_path"plot_pm3d_TP.gp
            pl_temp="$plot_path"plot_pm3d_TP_$j.gp
            
            cp -f $plot_file $pl_temp

            sed -i "s=otemp=$i=g"  $pl_temp
            sed -i "s="$out_path"="$plot_path"=g"  $pl_temp
            sed -i "s=.gnumeric=.pdf=g" $pl_temp
            sed -i "s=itemp=$i=g"  $pl_temp
            sed -i "s=xmin=$x_min=g"  $pl_temp
            sed -i "s=xmax=$x_max=g"  $pl_temp
            sed -i "s=ymin=$y_min=g"  $pl_temp
            sed -i "s=ymax=$y_max=g"  $pl_temp
            sed -i "s=cbmin=$cb_min=g"  $pl_temp
            sed -i "s=cbmax=$cb_max=g"  $pl_temp
            
            sed -i "s=lbxtemp=$lbx=g"  $pl_temp
            sed -i "s=lbytemp=$lby=g"  $pl_temp
            sed -i "s=title_temp=$tit=g"  $pl_temp
            
            gnuplot $pl_temp
            
            rm -rf $pl_temp
        done
    fi
