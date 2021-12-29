#!/bin/sh

#----------------------------------------- Batch ------------------------------------------

    inp_path="inp/"
    out_path="out/"
    
    input="1_case.inp"
    i_read=0 
    j_batch=0
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
    
    if [ $func_col_fin -gt $func_col_ini ]
    then
    
        if [ $HHG_case -eq 1 ]
        then 
        
            fp_name=Out_Spectra_
            rep_name=Power_Spectrum
            x_scan=1
            y_scan=2     
            cname='log10(S(w))'
            
            j=0
                    
            for file in `ls "$out_path""$fp_name"* | sort -V`
            do 
                ((j=j+1))
                
                cp  $file "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=$cname=$file=g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=Out_==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=.gnumeric==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=#==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=_=-=g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=out/==g" "$out_path"N_"$rep_name"_$j.gnumeric
                
                if [ $j -eq 1 ]
                then
                    awk '{print $'$x_scan',$'$y_scan'}' "$out_path"N_"$rep_name"_$j.gnumeric > "$out_path"C_"$rep_name"_$j.gnumeric
                else 
                    awk '{print $'$y_scan'}' "$out_path"N_"$rep_name"_$j.gnumeric > "$out_path"C_"$rep_name"_$j.gnumeric
                fi 

            done 

            paste `ls "$out_path"C_"$rep_name"* | sort -V`| column -s $'\t' -t > "$out_path"Batch_"$rep_name".gnumeric
            
        fi
        
        if [ $amp_spec_case -eq 1 ]
        then 
            
            fp_name=Out_Spectra_
            rep_name=Amplitude_Spectrum
            x_scan=1
            y_scan=3   
            cname='log10(A(w))'
            
            j=0
                    
            for file in `ls "$out_path""$fp_name"* | sort -V`
            do 
                ((j=j+1))
                
                cp  $file "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=$cname=$file=g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=Out_==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=.gnumeric==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=#==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=_=-=g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=out/==g" "$out_path"N_"$rep_name"_$j.gnumeric
                
                if [ $j -eq 1 ]
                then
                    awk '{print $'$x_scan',$'$y_scan'}' "$out_path"N_"$rep_name"_$j.gnumeric > "$out_path"C_"$rep_name"_$j.gnumeric
                else 
                    awk '{print $'$y_scan'}' "$out_path"N_"$rep_name"_$j.gnumeric > "$out_path"C_"$rep_name"_$j.gnumeric
                fi 

            done 

            paste `ls "$out_path"C_"$rep_name"* | sort -V`| column -s $'\t' -t > "$out_path"Batch_"$rep_name".gnumeric
            
        fi
        
        if [ $phase_spec_case -eq 1 ]
        then 
            fp_name=Out_Spectra_
            rep_name=Phase_Spectrum
            x_scan=1
            y_scan=4   
            cname='phi(w)'
            
            j=0
                    
            for file in `ls "$out_path""$fp_name"* | sort -V`
            do 
                ((j=j+1))
                
                cp  $file "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=$cname=$file=g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=Out_==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=.gnumeric==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=#==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=_=-=g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=out/==g" "$out_path"N_"$rep_name"_$j.gnumeric
                
                if [ $j -eq 1 ]
                then
                    awk '{print $'$x_scan',$'$y_scan'}' "$out_path"N_"$rep_name"_$j.gnumeric > "$out_path"C_"$rep_name"_$j.gnumeric
                else 
                    awk '{print $'$y_scan'}' "$out_path"N_"$rep_name"_$j.gnumeric > "$out_path"C_"$rep_name"_$j.gnumeric
                fi 

            done 

            paste `ls "$out_path"C_"$rep_name"* | sort -V`| column -s $'\t' -t > "$out_path"Batch_"$rep_name".gnumeric
            
        fi
        
        header=$(head -n 1 "$inp_path"input)
        
        for i in "$out_path"Batch*Spectrum*
        do 
            sed -i.bak "1 s/^.*$/$header/" $i
        done
        
        if [ $att_pulse_case -eq 1 ]
        then 
            fp_name=Out_atto_
            rep_name=SAP
            x_scan=3
            y_scan=4  
            cname='SAP'
            
            j=0
                    
            for file in `ls "$out_path""$fp_name"* | sort -V`
            do 
                ((j=j+1))
                
                cp  $file "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=$cname=$file=g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=Out_==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=.gnumeric==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=#==g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=_=-=g" "$out_path"N_"$rep_name"_$j.gnumeric
                sed -i "s=out/==g" "$out_path"N_"$rep_name"_$j.gnumeric
                
                if [ $j -eq 1 ]
                then
                    awk '{print $'$x_scan',$'$y_scan'}' "$out_path"N_"$rep_name"_$j.gnumeric > "$out_path"C_"$rep_name"_$j.gnumeric
                else 
                    awk '{print $'$y_scan'}' "$out_path"N_"$rep_name"_$j.gnumeric > "$out_path"C_"$rep_name"_$j.gnumeric
                fi 

            done 

            paste `ls "$out_path"C_"$rep_name"* | sort -V`| column -s $'\t' -t > "$out_path"Batch_"$rep_name".gnumeric
            
        fi

        rm -rf "$out_path"N_*
        rm -rf "$out_path"C_*
        
    fi
        
#----------------------------------------- Batch: Harm_Scan ------------------------------------------

    input="1_case.inp"
    i_read=0 
    while read line
    do 
        ((i_read=i_read+1))
        
        if [ $i_read -eq 9 ]
        then 
            harm_scan_case=$(echo "$line" | cut -c1-1)
        fi 
    done < "$inp_path""$input"
    
    if [ $harm_scan_case -gt 0 ]
    then
    
        fp_name=Out_Harmonic_
        
        c_scan=3    
        if [ $harm_scan_case -eq 1 ]
        then 
            rep_name='Delta_gap'
        elif [ $harm_scan_case -eq 2 ]
        then 
            rep_name='ecc'
        elif [ $harm_scan_case -eq 3 ]
        then 
            rep_name='E0'
            rep_name_inten='Intensity'
            c_scan=4
        elif [ $harm_scan_case -eq 4 ]
        then 
            rep_name='CEP'
        elif [ $harm_scan_case -eq 5 ]
        then 
            rep_name='T1'
        elif [ $harm_scan_case -eq 6 ]
        then 
            rep_name='T2'
        elif [ $harm_scan_case -eq 7 ]
        then 
            rep_name='Pu-Pr-Delay'
        fi
        
        cname='Log10(S(w))'

        j=0
                    
        for file in `ls "$out_path""$fp_name"* | sort -V`
        do 
            ((j=j+1))
            
            cp  $file "$out_path"N_"$rep_name"_$j.gnumeric
            sed -i "s=$cname=$file=g" "$out_path"N_"$rep_name"_$j.gnumeric
            sed -i "s=Out_==g" "$out_path"N_"$rep_name"_$j.gnumeric
            sed -i "s=.gnumeric==g" "$out_path"N_"$rep_name"_$j.gnumeric
            sed -i "s=#==g" "$out_path"N_"$rep_name"_$j.gnumeric
            sed -i "s=_=-=g" "$out_path"N_"$rep_name"_$j.gnumeric
            sed -i "s=out/==g" "$out_path"N_"$rep_name"_$j.gnumeric
            
            if [ $j -eq 1 ]
            then
                awk '{print $2,$'$c_scan'}' "$out_path"N_"$rep_name"_$j.gnumeric > "$out_path"C_"$rep_name"_$j.gnumeric
            else 
                awk '{print $'$c_scan'}' "$out_path"N_"$rep_name"_$j.gnumeric > "$out_path"C_"$rep_name"_$j.gnumeric
            fi 

        done 

        paste `ls "$out_path"C_"$rep_name"* | sort -V`| column -s $'\t' -t > "$out_path"Batch_Scan_"$rep_name".gnumeric
        
        
        if [ $harm_scan_case -eq 3 ]
        then 
        
            j=0
        
            for file in `ls "$out_path""$fp_name"* | sort -V`
            do
                ((j=j+1))
                
                cp  $file "$out_path"N_"$rep_name_inten"_$j.gnumeric
                sed -i "s=$cname=$file=g" "$out_path"N_"$rep_name_inten"_$j.gnumeric
                sed -i "s=Out_==g" "$out_path"N_"$rep_name_inten"_$j.gnumeric
                sed -i "s=.gnumeric==g" "$out_path"N_"$rep_name_inten"_$j.gnumeric
                sed -i "s=#==g" "$out_path"N_"$rep_name_inten"_$j.gnumeric
                sed -i "s=_=-=g" "$out_path"N_"$rep_name_inten"_$j.gnumeric
                sed -i "s=out/==g" "$out_path"N_"$rep_name_inten"_$j.gnumeric
                
                if [ $j -eq 1 ]
                then
                    awk '{print $3,$'$c_scan'}' "$out_path"N_"$rep_name_inten"_$j.gnumeric > "$out_path"C_"$rep_name_inten"_$j.gnumeric
                else 
                    awk '{print $'$c_scan'}' "$out_path"N_"$rep_name_inten"_$j.gnumeric > "$out_path"C_"$rep_name_inten"_$j.gnumeric
                fi 
                
            done 
            
            paste `ls "$out_path"C_"$rep_name_inten"* | sort -V`| column -s $'\t' -t > "$out_path"Batch_Scan_"$rep_name_inten".gnumeric
        
        fi
        
        rm -rf "$out_path"N_*
        rm -rf "$out_path"C_*
    fi
        
 
