!-------------------------------------------------
!    12 Aug 2021
!    Fast Fourier Transform - HHG - Time profile
!    Author : Mohammadreza Eidi
!-------------------------------------------------
 
    program FFT
    
    use MKL_DFTI
    use, intrinsic :: iso_c_binding 
    implicit none
    include 'fftw3.f03'

    integer, parameter :: dp = selected_real_kind(15,307)
    real(dp), parameter :: pi=(4.0_dp) * datan (1.0_dp)
    
    integer :: n_column,FFT_win_case,freq_resol_case,HHG_case,Time_profile_case,att_pulse_case,FFT_pack_case,n_FFT,n_Harm,i_Harm,iTime_profile,nTime_profile,label_inp_case,time_case,harm_scan_case,amp_spec_case,phase_spec_case,j_harm,norm_case,wrap_phase_case,zero_pad_case
    integer :: i_column,i_r,it,n_time,n_t,io,i,time_col_index,func_col_ini,func_col_fin,func_col_index,i_SAP,n_SAP_ini,n_SAP_fin,i_harm_scan,i_scan_param,harm_scan_i,harm_scan_f,harm_scan_s,n_zero_pad,n_time_z
    real(dp) :: dt,t_tot,t_tot_z,dNu,dw,dTime_profile,att_harm_i,att_harm_f,Gabor_time_win,landa,t_1c,Nu_0,w_0,Harm_max,t_Time_profile,a0,a1,a2,a3,a4,p_sin,alpha_blackman,sigma_gaus,harm_scan_tol,scan_param_i,d_scan_param,scan_param,scan_param_intensity,Harm_sum,abs_FFT_func_sum,phase_shift,delta_phase
    real(dp), dimension(:) , allocatable :: t,SAP,Harm,FFT_win,func
    real(dp), dimension(:,:) , allocatable :: inp
    complex(dp), dimension(:) , allocatable :: FFT_func, phase_spec
    complex(dp) :: SAP_Harm
    
    ! if (FFT_pack_case==1) then
    integer*8 :: plan
    ! if (FFT_pack_case==2)
    type(dfti_descriptor), pointer :: HHG_hand, TP_hand
    integer :: status
    
    character (len=4) :: inp_path,out_path
    character (len=20) :: func_col_index_char,i_harm_scan_char
    character (len=20) , dimension (:) , allocatable :: label_c
    
    inp_path='inp/'
	out_path='out/'
	
    open (file=inp_path//'1_case.inp' ,unit=1)	
    read(1,*) label_inp_case
    read(1,*) FFT_pack_case
    read(1,*) FFT_win_case
    read(1,*) HHG_case
    read(1,*) amp_spec_case
    read(1,*) phase_spec_case
    read(1,*) att_pulse_case
    read(1,*) Time_profile_case
    read(1,*) harm_scan_case
    read(1,*) time_case
    read(1,*) norm_case
    read(1,*) wrap_phase_case
    read(1,*) zero_pad_case

    open (file=inp_path//'2_parameter.inp' ,unit=2)	
    
    read(2,*) n_column
    
    allocate(label_c(n_column))
    
    read(2,*) time_col_index
    read(2,*) func_col_ini
    read(2,*) func_col_fin
    
    read(2,*) landa
    
    if ((time_case==1).or.(time_case==3)) then
        t_1c=landa*0.1378999
    else if (time_case==2) then
        t_1c=landa*0.1378999*2.4188843D-2
    end if
    
    Nu_0=1/t_1c
    w_0=2.0d0*pi*Nu_0
    
    read(2,*) Harm_max
    
    read(2,*) dTime_profile
    
    read(2,*) Gabor_time_win
    Gabor_time_win=1.0/(Gabor_time_win*w_0)
    
    read(2,*) att_harm_i
    read(2,*) att_harm_f
    read(2,*) p_sin
    read(2,*) sigma_gaus
    read(2,*) n_zero_pad
    
    open (file=inp_path//'3_harm_scan.inp' ,unit=3)
    read(3,*) harm_scan_i
    read(3,*) harm_scan_f
    read(3,*) harm_scan_s
    read(3,*) harm_scan_tol
    read(3,*) scan_param_i
    read(3,*) d_scan_param
    
    open (file=inp_path//'input' ,unit=4)
    
    n_time = 0     ! number of Rows 
    do
        read(4,*,iostat=io)
        if (io/=0) then 
            exit
        end if
        n_time=n_time+1
        
    end do

    rewind(4)
    
    if (label_inp_case==1) then
        n_time=n_time-1
        read (4,*) label_c(:)
    end if
    
    n_zero_pad=n_zero_pad*n_time
    n_time_z=n_zero_pad+n_time
    write(*,*) 'number of rows = ', n_time_z
        
    n_FFT=n_time_z

    allocate (inp(n_time_z,n_column))
    
    do i_r=1,n_time
        read (4,*) inp(i_r,:)
    end do
    
    do i_r=1+n_time,n_time_z
        inp(i_r,:) = 0.0d0
    end do
    
    allocate (t(n_time_z))
    t(:)=inp(:,time_col_index)
    if (time_case==3) then
        t(:)=t(:)*41.341374575751   ! fs > a.u.
    end if
    dt=t(2)-t(1)
    t_tot=n_time*dt
    t_tot_z=n_time_z*dt
    dNu=1/(dt*n_time_z)
    dw=dNu*2.0d0*pi
    
    n_Harm=nint(Harm_max*Nu_0/dNu)+1
    allocate (Harm(n_Harm))
    do i_Harm=1,n_Harm
        Harm(i_Harm)=(i_Harm-1)*dNu/Nu_0
    end do
    
    if (FFT_win_case>0) then
        allocate (FFT_win(n_time))
        a0=0.0d0
        a1=0.0d0
        a2=0.0d0
        a3=0.0d0
        a4=0.0d0
    end if
    
    select case (FFT_win_case)
        case (1)
            write(*,*)  'Applying Hann window'
            a0=0.5d0
            a1=0.5d0
        case (2)
            write(*,*)  'Applying Hamming window'
            a0=25.0d0/46.0d0
            a1=1-a0
        case (3)
            write(*,*)  'Applying blackman window'
            alpha_blackman=0.16d0
            a0=(1.0d0-alpha_blackman)/2.0d0
            a1=0.5d0
            a2=alpha_blackman/2.0d0
        case (4)
            write(*,*)  'Applying nuttall window'
            a0=0.355768d0
            a1=0.487396d0
            a2=0.144232d0
            a3=0.012604d0
        case (5)
            write(*,*)  'Applying blackman_nuttall window'
            a0=0.3635819d0
            a1=0.4891775d0
            a2=0.1365995d0
            a3=0.0106411d0
        case (6)
            write(*,*)  'Applying blackman_harris window'
            a0=0.35875d0
            a1=0.48829d0
            a2=0.14128d0
            a3=0.01168d0
        case (7)
            write(*,*)  'Applying Flat top window'
            a0=0.21557895d0
            a1=0.41663158d0
            a2=0.277263158d0
            a3=0.083578947d0
            a4=0.006947368d0
        case (8)
            write(*,*)  'Applying Triangular (Bartlett or Fejér) window'   
            FFT_win(:)=1.0d0-dabs(t(:)-(t_tot/2.0d0)/(t_tot/2.0d0))
        case (9)
            write(*,*)  'Applying Welch window'   
            FFT_win(:)=1.0d0-(t(:)-(t_tot/2.0d0)/(t_tot/2.0d0))**2.0d0
        case (10)
            write(*,*)  'Applying Modified Bartlett-Hanning window'   
            FFT_win(:)=0.62d0-0.48d0*dabs((t(:)/t_tot)-0.5d0)+0.38d0*dcos(2.0d0*pi*((t(:)/t_tot)-0.5d0))
        case (11)
            write(*,*)  'Applying Sine window'   
            FFT_win(:)=dsin(pi*t(:)/t_tot)
        case (12)
            write(*,*)  'Applying power of Sine window'   
            FFT_win(:)=dsin(pi*t(:)/t_tot)**p_sin
        case (13)
            write(*,*)  'Applying Bohman window'  
            FFT_win(:)=0.0d0
            do it=1,n_time
                if ((t(it)>=-1.0d0).and.(t(it)<=1.0d0)) then
                    FFT_win(it)=(1.0d0-dabs(t(it)))*dcos(pi*dabs(t(it)))+(1.0d0/pi)*dsin(pi*dabs(t(it)))
                end if
            end do
        case (14)
            write(*,*)  'Applying Gaussian window'   
            FFT_win(:)=dexp(-0.5d0*(t(:)-(t_tot/2.0d0)/(sigma_gaus*t_tot/2.0d0))**2.0d0)
    end select
    
    if ((FFT_win_case>0).and.(FFT_win_case<8)) then
        FFT_win(1:n_time)=a0-a1*dcos(2.0d0*pi*t(1:n_time)/t_tot)+a2*dcos(4.0d0*pi*t(1:n_time)/t_tot)-a3*dcos(6.0d0*pi*t(1:n_time)/t_tot)+a4*dcos(2.0d0*pi*t(1:n_time)/t_tot)
    end if
    
    if (Time_profile_case>0) then
        dTime_profile=dTime_profile*dt
        nTime_profile=nint(t_tot/dTime_profile)+1
    end if
    
    open (file=out_path//'param.out' ,unit=10)	
    
    write (10,*) 'n_time= ', n_time
    write (10,*) 'n_time_z= ', n_time_z
    write (10,*) 'n_FFT= ', n_FFT
    write (10,*) 'T_1c= ', t_1c
    write (10,*) 'T_tot= ', t_tot
    write (10,*) 'T_tot_z= ', t_tot_z
    write (10,*) 'Nu_0= ', Nu_0
    write (10,*) 'w_0= ', w_0
    write (10,*) 'dt= ', dt
    write (10,*) 'dNu= ', dNu
    write (10,*) 'dw= ', dw
    write (10,*) 'dTime_profile= ', dTime_profile
    write (10,*) 'Gabor_time_win= ', Gabor_time_win
    write (10,*) 'Harm_max= ', Harm_max
    write (10,*) 'n_Harm= ', n_harm
    write (10,*) 'nTime_profile= ', nTime_profile
    
    allocate (func(n_time_z),FFT_func(n_FFT), phase_spec(n_Harm))
    func=0.0d0
    if (att_pulse_case==1) then 
        allocate(SAP(n_time_z))
    end if 
    
    if (HHG_case==1) then 
    
        if (harm_scan_case>0) then
            
            do i_harm_scan=harm_scan_i,harm_scan_f,harm_scan_s
                            
                write(i_harm_scan_char,11) i_harm_scan
                open (file=out_path//'Out_Harmonic_'//trim(adjustl(i_harm_scan_char))//'.gnumeric' ,unit=10+i_harm_scan+3*func_col_fin)
                
                write (10+i_harm_scan+3*func_col_fin,10,advance="no") '#Harmonic_order'
                
                select case (harm_scan_case)
                    case(1) 
                        write (10+i_harm_scan+3*func_col_fin,10,advance="no") '#Delta_gap'
                        write (10+i_harm_scan+3*func_col_fin,10) '#Log10(S(w))'
                    case(2) 
                        write (10+i_harm_scan+3*func_col_fin,10,advance="no") '#ecc'
                        write (10+i_harm_scan+3*func_col_fin,10) '#Log10(S(w))'
                    case(3) 
                        write (10+i_harm_scan+3*func_col_fin,10,advance="no") '#E0'
                        write (10+i_harm_scan+3*func_col_fin,10,advance="no") '#Intensity'
                        write (10+i_harm_scan+3*func_col_fin,10) '#Log10(S(w))'
                        
                        write (10+i_harm_scan+3*func_col_fin,4) real(i_harm_scan),0.0d0,0.0d0,0d0
                    case(4) 
                        write (10+i_harm_scan+3*func_col_fin,10,advance="no") '#CEP'
                        write (10+i_harm_scan+3*func_col_fin,10) '#Log10(S(w))'
                    case(5) 
                        write (10+i_harm_scan+3*func_col_fin,10,advance="no") '#T1'
                        write (10+i_harm_scan+3*func_col_fin,10) '#Log10(S(w))'
                    case(6) 
                        write (10+i_harm_scan+3*func_col_fin,10,advance="no") '#T2'
                        write (10+i_harm_scan+3*func_col_fin,10) '#Log10(S(w))'
                    case(7) 
                        write (10+i_harm_scan+3*func_col_fin,10,advance="no") '#Pu-Pr-Delay'
                        write (10+i_harm_scan+3*func_col_fin,10) '#Log10(S(w))'
                end select
                 
            end do
            
            i_scan_param=1
            
        end if
        
    end if
    
    do func_col_index=func_col_ini,func_col_fin
        
        if ((HHG_case==1).and.(harm_scan_case>0)) then 
            scan_param=scan_param_i+(i_scan_param-1)*d_scan_param
            
            if (harm_scan_case==3) then
                scan_param_intensity=((10.0d0**10.0d0)*scan_param)**2.0d0/377.0d0
            end if
            
        end if
    
        write(func_col_index_char,11) func_col_index
        
        if ((HHG_case==1).or.(att_pulse_case==1)) then
        
            write(*,*) 'computing HHG of column= ', func_col_index
            
            if (FFT_win_case>0) then
                func(1:n_time)=inp(1:n_time,func_col_index)*FFT_win(1:n_time)
            else 
                func(:)=inp(:,func_col_index)
            end if
            
            if (FFT_pack_case==1) then
                call dfftw_plan_dft_r2c_1d(plan,n_time_z,func,FFT_func,FFTW_ESTIMATE)
                call dfftw_execute_dft_r2c(plan,func,FFT_func)
                call dfftw_destroy_plan(plan)
            else if (FFT_pack_case==2) then   
                status = DftiCreateDescriptor(HHG_hand,DFTI_DOUBLE,DFTI_REAL,1,n_time_z)
                Status = DftiSetValue(HHG_hand,DFTI_PLACEMENT,DFTI_NOT_INPLACE)
                status = DftiCommitDescriptor(HHG_hand)
                status = DftiComputeForward(HHG_hand,func,FFT_func)
                status = DftiFreeDescriptor(HHG_hand)
            end if
            
            do i_Harm=1,n_Harm
                if (norm_case==1) then
                    FFT_func(i_Harm)=FFT_func(i_Harm)/sqrt(float(n_time_z))
                else if (norm_case==2) then
                    FFT_func(i_Harm)=FFT_func(i_Harm)/n_time_z
                end if    
            end do
            
            if (HHG_case==1) then 
            
                open (file=out_path//'Out_Spectra_'//trim(adjustl(func_col_index_char))//'.gnumeric' ,unit=10+func_col_index)	
                
                write (10+func_col_index,10,advance="no") '#Harmonic_order'
                write (10+func_col_index,10,advance="no") '#log10(S(w))'
                write (10+func_col_index,10,advance="no") '#log10(A(w))'
                write (10+func_col_index,10) '#phi(w)'

                phase_spec(:)=datand(dimag(FFT_func(:))/dreal(FFT_func(:)))

                if (wrap_phase_case==1) then
                    phase_shift=0.0d0
                    do i_Harm=1,n_Harm-1
                        delta_phase=phase_spec(i_Harm+1)-(phase_spec(i_Harm)-phase_shift)
                        if (dabs(delta_phase)>180.0d0) then
                            phase_shift=phase_shift-sign(360.0d0,delta_phase)
                        end if
                        phase_spec(i_Harm)=phase_spec(i_Harm)+phase_shift
                    end do
                    phase_spec(n_Harm)=phase_spec(n_Harm)+phase_shift
                end if

                do i_Harm=1,n_Harm
                    write (10+func_col_index,1,advance="no") Harm(i_Harm)
                    write (10+func_col_index,1,advance="no") log10(cdabs(FFT_func(i_Harm))**2.0d0)
                    write (10+func_col_index,1,advance="no") log10(sqrt(cdabs(FFT_func(i_Harm))**2.0d0))
                    write (10+func_col_index,1)  phase_spec(i_Harm)
                end do
                
                if (harm_scan_case>0) then
         
                    do i_harm_scan=harm_scan_i,harm_scan_f,harm_scan_s
                        
                        j_harm=0
                        abs_FFT_func_sum=0.0_dp
                        Harm_sum=0.0_dp

                        do i_Harm=1,n_Harm
                            if ((Harm(i_Harm)>(real(i_harm_scan)-harm_scan_tol)).and.(Harm(i_Harm)<(real(i_harm_scan)+harm_scan_tol))) then
                                j_harm=j_harm+1
                                Harm_sum=Harm_sum+Harm(i_Harm)
                                abs_FFT_func_sum=abs_FFT_func_sum+cdabs(FFT_func(i_Harm))
                            else if (Harm(i_Harm)>=(real(i_harm_scan)+harm_scan_tol)) then
                                exit
                            end if
                        end do
                        Harm_sum=Harm_sum/j_harm
                        abs_FFT_func_sum=abs_FFT_func_sum/j_harm
                        if (harm_scan_case==3) then
                            write (10+i_harm_scan+3*func_col_fin,4) Harm_sum, scan_param , scan_param_intensity , log10(abs_FFT_func_sum**2.0d0)
                        else 
                            write (10+i_harm_scan+3*func_col_fin,3) Harm_sum, scan_param , log10(abs_FFT_func_sum**2.0d0) 
                        end if
  
                    end do
                    
                end if 
   
            end if
            
            if (att_pulse_case==1) then 
            
                write(*,*) 'Computing isolated attosecond pulse filtering between Harmonic number= ', att_harm_i,' and Harmonic number= ', att_harm_f
                
                open (file=out_path//'Out_atto_'//trim(adjustl(func_col_index_char))//'.gnumeric' ,unit=10+func_col_index++func_col_fin)	
                
                n_SAP_ini=nint(att_harm_i*Nu_0/dNu)+1
                n_SAP_fin=nint(att_harm_f*Nu_0/dNu)+1

                do it=1,n_time
                    SAP_Harm=(0.0d0,0.0d0)
                    do i_SAP=n_SAP_ini,n_SAP_fin
                        SAP_Harm=SAP_Harm+dw*FFT_func(i_SAP)*cdexp((0.0d0,1.0d0)*i_SAP*dw*t(it))
                    end do
                    SAP(it)=cdabs(SAP_Harm)**2.0d0
                end do
                
                write (10+func_col_index+func_col_fin,10,advance="no") 'T(a.u.)'
                write (10+func_col_index+func_col_fin,10,advance="no") 'T(fs)'
                write (10+func_col_index+func_col_fin,10,advance="no") 'T(cycle)'
                write (10+func_col_index+func_col_fin,10) 'SAP'
                
                do it=1,n_time
                    write (10+func_col_index+func_col_fin,3,advance="no") t(it),t(it)*2.4188843D-2,t(it)/t_1c
                    write (10+func_col_index+func_col_fin,1) SAP(it)
                end do
            end if
        end if
        
        if (Time_profile_case>0) then
        
            write(*,*) 'Using Gabor transform computing time frequnecy profile of column= ', func_col_index
            
            open (file=out_path//'Out_TP_Gabor_'//trim(adjustl(func_col_index_char))//'.gnumeric' ,unit=10+func_col_index+2*func_col_fin)	
            write (10+func_col_index+2*func_col_fin,10,advance="no") 'T(cycle)/w(w0)'
            
            do i_Harm=1,n_Harm
                if (i_Harm<n_harm) then
                    write (10+func_col_index+2*func_col_fin,1,advance="no") Harm(i_Harm)
                else 
                    write (10+func_col_index+2*func_col_fin,1) Harm(i_Harm)
                end if
            end do

            do iTime_profile=1,nTime_profile
            
                t_Time_profile=(iTime_profile-1)*dTime_profile

                if (FFT_win_case>0) then
                    func(1:n_time)=inp(1:n_time,func_col_index)*FFT_win(1:n_time)*dexp(-((t_Time_profile-t(1:n_time))**2.0d0)/(2.0d0*Gabor_time_win**2.0d0))
                else 
                    ! Gabor transform from: DOI: 10.1002/wcms.78 Eq. 14
                    func(:)=inp(:,func_col_index)*dexp(-((t_Time_profile-t(:))**2.0d0)/(2.0d0*Gabor_time_win**2.0d0))
                end if
          
                if (FFT_pack_case==1) then
                    call dfftw_plan_dft_r2c_1d(plan,n_time_z,func,FFT_func,FFTW_ESTIMATE)
                    call dfftw_execute_dft_r2c(plan,func,FFT_func)
                    call dfftw_destroy_plan(plan)
                else if (FFT_pack_case==2) then   
                    status = DftiCreateDescriptor(HHG_hand,DFTI_DOUBLE,DFTI_REAL,1,n_time_z)
                    Status = DftiSetValue(HHG_hand,DFTI_PLACEMENT,DFTI_NOT_INPLACE)
                    status = DftiCommitDescriptor(HHG_hand)
                    status = DftiComputeForward(HHG_hand,func,FFT_func)
                    status = DftiFreeDescriptor(HHG_hand)
                end if
                
                write (10+func_col_index+2*func_col_fin,1,advance="no") t_Time_profile/t_1c

                do i_Harm=1,n_Harm
                    if (norm_case==1) then
                        FFT_func(i_Harm)=FFT_func(i_Harm)/sqrt(float(n_time_z))
                    else if (norm_case==2) then
                        FFT_func(i_Harm)=FFT_func(i_Harm)/n_time_z
                    end if 
                end do

                phase_spec(:)=datand(dimag(FFT_func(:))/dreal(FFT_func(:)))

                if (wrap_phase_case==1) then
                    phase_shift=0.0d0
                    do i_Harm=1,n_Harm-1
                        delta_phase=phase_spec(i_Harm+1)-(phase_spec(i_Harm)-phase_shift)
                        if (dabs(delta_phase)>180.0d0) then
                            phase_shift=phase_shift-sign(360.0d0,delta_phase)
                        end if
                        phase_spec(i_Harm)=phase_spec(i_Harm)+phase_shift
                    end do
                    phase_spec(n_Harm)=phase_spec(n_Harm)+phase_shift
                end if

                do i_Harm=1,n_Harm

                    if (Time_profile_case==1) then
                        if (i_Harm<n_harm) then
                            write (10+func_col_index+2*func_col_fin,1,advance="no") log10(cdabs(FFT_func(i_Harm))**2.0d0)
                        else 
                            write (10+func_col_index+2*func_col_fin,1) log10(cdabs(FFT_func(i_Harm))**2.0d0)
                        end if
                    else if (Time_profile_case==2) then
                        if (i_Harm<n_harm) then
                            write (10+func_col_index+2*func_col_fin,1,advance="no") log10(sqrt(cdabs(FFT_func(i_Harm))**2.0d0))
                        else 
                            write (10+func_col_index+2*func_col_fin,1) log10(sqrt(cdabs(FFT_func(i_Harm))**2.0d0))
                        end if
                    else if (Time_profile_case==3) then
                        if (i_Harm<n_harm) then
                            write (10+func_col_index+2*func_col_fin,1,advance="no") phase_spec(i_Harm)
                        else 
                            write (10+func_col_index+2*func_col_fin,1) phase_spec(i_Harm)
                        end if
                    end if
                end do
            end do
        end if
        
        if ((HHG_case==1).and.(harm_scan_case>0)) then 
            i_scan_param=i_scan_param+1
        end if
        
    end do
    
10   format(1(a,9x)) 	 
11   format(1(i10,1x))
1    format(1(ES14.7,1x))
2    format(2(ES14.7,1x))
3    format(3(ES14.7,1x))
4    format(4(ES14.7,1x))
    
    
    end program FFT
