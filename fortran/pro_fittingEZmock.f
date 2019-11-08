!Fitting procedure for EZmock
      program fitting ezmock
      implicit none 
!width of a step
      integer,parameter :: max_main_step = 4
      integer,parameter :: num_step1 = 5
      integer,parameter :: num_step2 = 5
      integer,parameter :: num_step3 = 5
      integer,parameter :: num_step4 = 3

      real :: delta_modify_pdf = 0.5
      real :: delta_density_cut = 0.1
      real :: delta_density_sat = 0.1
      real :: delta_grow2z0 = 0.1
      real :: delta_modify_pk = 10.
      real :: delta_antidamping = 0.01
      real :: delta_zdist_rate = 0.5
      real :: delta_zdist_fog = 10

      character*200 :: log_file
      character*200 :: next_par_file!='./next.par'
      character*200 :: pk_mono_file,pk_zdist_mono_file,
     &   pk_zdist_quad_file
      character*200 :: cf_mono_file,cf_zdist_mono_file,
     &   cf_zdist_quad_file
      character*200 :: bik_file !!----------------------- The format should be {theta, b(k), ...} !!!
      character*200 :: sim_pk_mono_file,sim_pk_zdist_mono_file,
     &  sim_pk_zdist_quad_file
      character*200 :: sim_cf_mono_file,sim_cf_zdist_mono_file,
     &  sim_cf_zdist_quad_file
      character*200 :: sim_bik_file
      character*200 :: sim_file, mock_file, sim_file2, mock_file2

      real :: chi_old,chi_new,chi_new1,chi_new2,
     &  kk,pkk,k_min,k_max,k_min2,k_max2,p_old, p_new, p_new2
      integer :: main_step, sub_step !the key flag to determine what to do
      integer :: errorflag, count_sim, count_mock
      real    :: total_pk_sim, total_pk_mock, delta
      integer :: param_changed = 0, i
      logical :: compute_bik
      logical :: compute_simPDF
      character*100 :: format_string
!-- all the input for ezmock. Will read from input file
      character*200     :: datafile_xyz
      character*200     :: datafile_prefix
      real    :: boxsize  =2500
      integer :: grid_num =480!! 512 using 4 cores per node, 1024 use 16 core (1 node)
      real :: grow2z0 = 0.5 !0.744657**2 !0.744657^2 for Planck Multidark at 0.57
      real :: scatter = 10
      real :: scatter2 = -0.
      real :: modify_pk = 170  !fiducial 350
      real :: modify_pdf = 0.04
      real :: antidamping = 0.1 ! not used if >=1
      real :: density_cut = 2.5
      real :: density_sat = 4.6 !saturation or upper bound
      real :: zdist_rate = 78 !something similar to linear growth factor 
      real :: zdist_fog = 100
      character*200  :: whitenoise_file
      logical        :: use_whitenoise_file = .true. 
      character*200     :: pkfile
      character*200     :: pknwfile
      character*200     :: datafile_path
      logical :: compute_CF
      logical :: compute_CF_zdist
      integer       :: skiplines = 6
      character*200 :: twod_corr_suffix = '.bin5.corr'
      real          :: max_r = 250
      integer       :: bin_size = 5
      real          :: redshift = 0.5618
      real          :: om = 0.307115 !planck

      real :: old_modify_pdf, 
     & old_density_cut, old_density_sat, old_grow2z0, old_modify_pk,
     & old_antidamping, old_zdist_rate, old_zdist_fog
      real :: new_modify_pdf, 
     & new_density_cut, new_density_sat, new_grow2z0, new_modify_pk,
     & new_antidamping, new_zdist_rate, new_zdist_fog

!read input file
       namelist /EZmock_v0_input/
     &  datafile_xyz, datafile_prefix, boxsize, grid_num,
     & grow2z0, scatter, scatter2, modify_pk, modify_pdf,
     & antidamping, density_cut, density_sat, zdist_rate,
     & zdist_fog, whitenoise_file, use_whitenoise_file,
     & pkfile, pknwfile, datafile_path,
     & compute_CF, compute_CF_zdist,
     & skiplines, twod_corr_suffix, max_r, bin_size,
     & redshift, om

!input                                                                                                                 
      read(*,EZmock_v0_input)

      write(*,*) "print out input parameters from input file:",
     &  datafile_xyz, datafile_prefix, boxsize, grid_num,
     & grow2z0, scatter, scatter2, modify_pk, modify_pdf,
     & antidamping, density_cut, density_sat, zdist_rate,
     & zdist_fog, whitenoise_file, use_whitenoise_file,
     & pkfile, pknwfile, datafile_path,
     & compute_CF, compute_CF_zdist,
     & skiplines, twod_corr_suffix, max_r, bin_size,
     & redshift, om

       new_modify_pdf =  modify_pdf
       new_density_cut = density_cut 
       new_density_sat = density_sat
       new_grow2z0 = grow2z0
       new_modify_pk = modify_pk
       new_antidamping = antidamping 
       new_zdist_rate = zdist_rate 
       new_zdist_fog = zdist_fog

!initialize
      errorflag = 0
      next_par_file=Trim(datafile_path)//Trim(datafile_prefix)//
     &  '.next.par'      
      log_file = Trim(datafile_path)//Trim(datafile_prefix)//'.log'      
      pk_mono_file = Trim(datafile_path)//
     & Trim(datafile_prefix)//'.dat.CICassign.NGPbin.pk.mono'
      pk_zdist_mono_file = Trim(datafile_path)// 
     &  Trim(datafile_prefix)//'.dat.CICassign.NGPbin.pk.zdist.mono'
      pk_zdist_quad_file = Trim(datafile_path)//
     &  Trim(datafile_prefix)//'.dat.CICassign.NGPbin.pk.zdist.quad'
      cf_mono_file = Trim(datafile_prefix)//'.dat.bin5.corr.mono'
      cf_zdist_mono_file =Trim(datafile_path)//
     &  Trim(datafile_prefix)//'.dat.bin5.corr.zdist.mono'
      cf_zdist_quad_file =Trim(datafile_path)//
     &  Trim(datafile_prefix)//'.dat.bin5.corr.zdist.quad'
      bik_file = Trim(datafile_path)//
     &  Trim(datafile_prefix)//'.bik'

      sim_pk_mono_file =Trim(datafile_path)//
     & 'Sim_'//Trim(datafile_prefix)//
     &  '.CICassign.NGPbin.pk'
      sim_pk_zdist_mono_file =Trim(datafile_path)// 'Sim_'//
     &  Trim(datafile_prefix)//'.CICassign.NGPbin.pk.zdist.mono'
      sim_pk_zdist_quad_file =Trim(datafile_path)// 'Sim_'//
     &  Trim(datafile_prefix)//'.CICassign.NGPbin.pk.zdist.quad'
      sim_cf_mono_file =Trim(datafile_path)//
     &  Trim(datafile_prefix)//'.bin5.corr.mono'
      sim_cf_zdist_mono_file =Trim(datafile_path)//'Sim_'//
     &  Trim(datafile_prefix)//'.bin5.corr.zdist.mono'
      sim_cf_zdist_quad_file =Trim(datafile_path)//'Sim_'//
     &  Trim(datafile_prefix)//'.bin5.corr.zdist.quad'
      sim_bik_file = Trim(datafile_path)//
     &  'Sim_'//Trim(datafile_prefix)//'.bik'

      format_string = '(2I4,8E12.5,I2,2E12.5)'
      main_step=0
      sub_step=0
      open(10, file=Trim(log_file))
 100  read(10, format_string, end=200) main_step, sub_step,
     & old_modify_pdf, 
     & old_density_cut, old_density_sat, old_grow2z0, old_modify_pk,
     & old_antidamping, old_zdist_rate, old_zdist_fog,
     & param_changed, chi_old, chi_new

      goto 100
 200  close(10)

      write(*,*) "print out the params read from log file:"
      write(*,*)  main_step, sub_step, old_modify_pdf, 
     & old_density_cut, old_density_sat, old_grow2z0, old_modify_pk,
     & old_antidamping, old_zdist_rate, old_zdist_fog,
     & param_changed, chi_old, chi_new

! while sub_step = 0, main_step add 1.
      if(sub_step .eq. 0) then
           main_step = main_step + 1
      end if

      if(main_step .eq. 1) then !find pdf parameter
        sim_file = Trim(sim_pk_mono_file)
        mock_file = Trim(pk_mono_file)
        k_min = 0.8
        k_max = 0.9
        sim_file2 = Trim(sim_file)
        mock_file2 = Trim(mock_file)
        k_min2 = k_min
        k_max2 = k_max
        param_changed = 1

      else if(main_step .eq. 2) then !find density_sat with small box
        sim_file = Trim(sim_pk_mono_file)
        mock_file = Trim(pk_mono_file)
        k_min = 0.05
        k_max = 0.1
        sim_file2 = Trim(sim_file)
        mock_file2 = Trim(mock_file)
        k_min2 = k_min
        k_max2 = k_max
        param_changed = 2

      else if(main_step .eq. 3) then !find pdf parameter
        sim_file = Trim(sim_pk_mono_file)
        mock_file = Trim(pk_mono_file)
        k_min = 0.4
        k_max = 0.5
        sim_file2 = Trim(sim_file)
        mock_file2 = Trim(mock_file)
        k_min2 = 0.05
        k_max2 = 0.1
        param_changed = 1

      else if(main_step .eq. 4) then !find density_sat with big box
        sim_file = Trim(sim_pk_mono_file)
        mock_file = Trim(pk_mono_file)
        k_min = 0.05
        k_max = 0.1
        sim_file2 = Trim(sim_file)
        mock_file2 = Trim(mock_file)
        k_min2 = k_min
        k_max2 = k_max
        param_changed = 2

      else if(main_step .eq. max_main_step+1) then !same as previous step
        sim_file = Trim(sim_pk_mono_file)
        mock_file = Trim(pk_mono_file)
        k_min = 0.05
        k_max = 0.1
        sim_file2 = Trim(sim_file)
        mock_file2 = Trim(mock_file)
        k_min2 = k_min
        k_max2 = k_max
        param_changed = 2

      end if

! compute chi
      write(*,*) "start to compute chi"
      count_sim=0
      total_pk_sim = 0
      open(20,file=Trim(sim_file), status='old')
      read(20,*)
 300  read(20,*,end=400) kk, pkk
      count_sim = count_sim + 1
      if(kk.ge.k_min .and. kk.le.k_max) then
        pkk = pkk*kk
        total_pk_sim = total_pk_sim + pkk
      end if
      goto 300
 400  close(20)

      count_mock=0
      total_pk_mock = 0
      open(20,file=Trim(mock_file), status='old')
      read(20,*)
 500  read(20,*,end=600) kk, pkk
      count_mock = count_mock + 1
      if(kk.ge.k_min .and. kk.le.k_max) then
        pkk = pkk*kk
        total_pk_mock = total_pk_mock + pkk
      end if
      goto 500
 600  close(20)
      if(count_sim .ne. count_mock) then
         write(*,*) "bins of sim and mock are not matching, check!!"
         errorflag = 1
      end if

      chi_new1 = total_pk_mock - total_pk_sim

      count_sim=0
      total_pk_sim = 0
      open(20,file=Trim(sim_file2), status='old')
      read(20,*)
 302  read(20,*,end=402) kk, pkk

      count_sim = count_sim + 1
      if(kk.ge.k_min2 .and. kk.le.k_max2) then
        pkk = pkk*kk
        total_pk_sim = total_pk_sim + pkk
      end if
      goto 302
 402  close(20)


      count_mock=0
      total_pk_mock = 0
      open(20,file=Trim(mock_file2), status='old')
      read(20,*)
 502  read(20,*,end=602) kk, pkk
      count_mock = count_mock + 1
      if(kk.ge.k_min2 .and. kk.le.k_max2) then
        pkk = pkk*kk
        total_pk_mock = total_pk_mock + pkk
      end if
      goto 502
 602  close(20)
      if(count_sim .ne. count_mock) then
         write(*,*) "bins of sim and mock are not matching, check!!"
         errorflag = 1
      end if

      chi_new2 = total_pk_mock - total_pk_sim

! determine next move

! main step 1: adjust modify_pdf
      if(main_step.eq.1) then

        if(mod(sub_step,2) .eq. 0) then
          if(chi_new1 .ge. 0) then
            new_modify_pdf = modify_pdf + delta_modify_pdf
          else
            new_modify_pdf = modify_pdf - delta_modify_pdf
          end if
        else
          p_old = old_modify_pdf
          p_new = modify_pdf
          p_new2 = (p_old*chi_new1-p_new*chi_old)/(chi_new1-chi_old) 
          delta = delta_modify_pdf
          if(p_new2-p_new .gt. delta*3) p_new2 = 
     &                                    p_new + delta*3
          if(p_old-p_new2 .gt. delta*3) p_new2 = 
     &                                    p_old - delta*3
          new_modify_pdf = p_new2
        end if

        if(sub_step .eq. num_step1*2-1) then
           sub_step = 0
        else
           sub_step = sub_step + 1
        end if
        boxsize = 1250
        grid_num= 480
        compute_bik = .false.
        compute_simPDF = .false.

 !------------------------main step 2: adjust density_sat in small box
      else if(main_step .eq. 2) then

        if(mod(sub_step,2) .eq. 0) then
          if(chi_new1 .ge. 0) then
            new_density_sat = density_sat - delta_density_sat
          else
            new_density_sat = density_sat + delta_density_sat
          end if
        else
          p_old = old_density_sat
          p_new = density_sat
          p_new2 = (p_old*chi_new1-p_new*chi_old)/(chi_new1-chi_old) 
          delta = delta_density_sat
          if(p_new2-p_new .gt. delta*5) p_new2 = 
     &                                    p_new + delta*5
          if(p_old-p_new2 .gt. delta*5) p_new2 = 
     &                                    p_old - delta*5
          new_density_sat = p_new2
        end if

        if( sub_step .eq. num_step2*2 - 1) then
          sub_step = 0
        else
          sub_step = sub_step + 1
        end if
        compute_bik = .false.
        boxsize = 1250
        grid_num= 480
        compute_simPDF = .false.

!---- main step 3: adjust modify_pdf
      else if(main_step.eq.3) then 

        if(mod(sub_step,2) .eq. 0) then
          if(chi_new1 .ge. 0) then
            new_modify_pdf = modify_pdf + delta_modify_pdf
          else
            new_modify_pdf = modify_pdf - delta_modify_pdf
          end if
        else if(mod(sub_step,2) .eq. 1) then
          p_old = old_modify_pdf
          p_new = modify_pdf
          p_new2 = (p_old*chi_new1-p_new*chi_old)/(chi_new1-chi_old) 
          delta = delta_modify_pdf
          if(p_new2-p_new .gt. delta*3) p_new2 = 
     &                                    p_new + delta*3
          if(p_old-p_new2 .gt. delta*3) p_new2 = 
     &                                    p_old - delta*3
          new_modify_pdf = p_new2

        else if(mod(sub_step,2) .eq. 2) then
          if(chi_new2 .ge. 0) then
            new_modify_pdf = modify_pdf + delta_modify_pdf
          else
            new_modify_pdf = modify_pdf - delta_modify_pdf
          end if
        else
          p_old = old_modify_pdf
          p_new = modify_pdf
          p_new2 = (p_old*chi_new2-p_new*chi_new)/(chi_new2-chi_new) 
          delta = delta_modify_pdf
          if(p_new2-p_new .gt. delta*3) p_new2 = 
     &                                    p_new + delta*3
          if(p_old-p_new2 .gt. delta*3) p_new2 = 
     &                                    p_old - delta*3
          new_modify_pdf = p_new2
        end if

        if(sub_step .eq. num_step3*4-1) then
           sub_step = 0
           compute_simPDF = .true.
           boxsize = 2500
           grid_num= 960
        else
           compute_simPDF = .false.
           boxsize = 1250
           grid_num= 480
           sub_step = sub_step + 1
        end if
         compute_bik = .false.

! main step 4: big box for density_sat
      else if(main_step .eq. 4) then

        if(mod(sub_step,2) .eq. 0) then
          if(chi_new1 .ge. 0) then
            new_density_sat = density_sat - delta_density_sat
          else
            new_density_sat = density_sat + delta_density_sat
          end if
        else
          p_old = old_density_sat
          p_new = density_sat
          p_new2 = (p_old*chi_new1-p_new*chi_old)/(chi_new1-chi_old) 
          delta = delta_density_sat
          if(p_new2-p_new .gt. delta*3) p_new2 = 
     &                                    p_new + delta*3
          if(p_old-p_new2 .gt. delta*3) p_new2 = 
     &                                    p_old - delta*3
          new_density_sat = p_new2
        end if

        if( sub_step .eq. num_step4*2 - 1) then
          sub_step = 0
        else
          sub_step = sub_step + 1
        end if
        boxsize = 2500
        grid_num= 960
        compute_bik = .false.
        compute_simPDF = .false.

      end if

      if(main_step .eq. max_main_step+1 .and. sub_step .eq. 0) then
          errorflag = 1
          write(*,*) "errorflag = 1!!!!!"
          main_step = 0
          sub_step  = 0
      end if

! write log file
      
      open(10, file=Trim(log_file), access='append')
        write(10,format_string)  main_step, sub_step,
     & modify_pdf, 
     & density_cut, density_sat,  grow2z0, modify_pk,
     & antidamping, zdist_rate, zdist_fog,
     & param_changed, chi_new1, chi_new2
      close(10)

      open(30,file=Trim(next_par_file))

      write(30,*) errorflag !0 for ok, 1 for stop or problems
      write(30,*) boxsize
      write(30,*) grid_num
      write(30,*) new_grow2z0
      write(30,*) new_modify_pk
      write(30,*) new_modify_pdf
      write(30,*) new_antidamping
      write(30,*) new_density_cut
      write(30,*) new_density_sat
      write(30,*) new_zdist_rate 
      write(30,*) new_zdist_fog

      if(compute_CF .eq. .true.) then
        write(30,*) ".true."
      else
        write(30,*) ".false."
      end if
      if(compute_bik .eq. .true.) then
        write(30,*) ".true."
      else
        write(30,*) ".false."
      end if
      if(compute_simPDF .eq. .true.) then
        write(30,*) ".true."
      else
        write(30,*) ".false."
      end if

      close(30)      
      
!--------------------------------

      end program      
