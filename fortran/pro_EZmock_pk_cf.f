! 20180413 add scatter2 to boost low density for modifying bik
! use both particles and CIC to assign halos (when running out of particles for given grid point) -- 20171201
!     EZmock2017 0107 use CIC displacement to assign velocity
!     note: cleaner memory usage, optimize CIC, assign galaxies to particles (20161228)
!
      
      program EZmock
! make clean
! make
      implicit none
!output file prefix
      character*500     :: datafile_prefix
      real    :: boxsize 
      integer :: grid_num
      real :: grow2z0
      real :: scatter
      real :: scatter2
      real :: modify_pk
      real :: modify_pdf
      real :: antidamping ! not used if >=1
      real :: density_cut
      real :: density_sat !saturation or upper bound
      real :: zdist_rate !something similar to linear growth factor 
      real :: zdist_fog !for FOG
      character*500  :: whitenoise_file
      logical        :: use_whitenoise_file
      character*500     :: pkfile
      character*500     :: pknwfile
      character*500     :: datafile_path
      character*500 :: output_file
      integer       ::  data_num, i_dilute
      real :: expect_sum_pdf, expect_A_pdf
      integer :: iseed
      
!!
      character*500 :: output_suffix = ".dat"
      CHARACTER*5  :: rank_string='00000'
      integer      :: i, num_loop,j,k,l, ip,grid_pdf
      real :: temp,ran3_omp

      integer,parameter :: max_data=1000000000 !for keeping particles to assign galaxies
      real :: ran3
      integer :: omp_get_max_threads
      real,parameter    :: data_density = 0.1 !for output raw particle data !per grid cell
!-- 2PCF
      real :: dilute_factor
      logical :: compute_CF
      logical :: compute_CF_zdist
      integer       :: skiplines
      character*500 :: twod_corr_suffix = '.bin5.corr'
      real          :: max_r ! = 250
      integer       :: bin_size ! = 5
      real          :: redshift ! = 0.5618
      real          :: om ! = 0.307115 !planck
      real :: inverse_aH
!!
      character*500 :: twod_nofr_file
      integer :: num_r
      integer :: num_mu = 100
      integer       :: count
      real          :: r1,r2,r3,r4
      real,allocatable :: x(:), y(:), w(:), wt(:), vz(:)
      integer       ::   rand_num, ii, jj,kk,rr,mm
      real,allocatable :: DDnofr(:,:), DRnofr(:,:), RRnofr(:,:)
      real          :: totalwt_data, totalwt_rand, corr
      real          :: x_min,x_max,y_min,y_max,z_min,z_max
      real          :: ddbin, d2rbin,rrbin, rrbin_de, z, hubble
      integer :: time, time1, time2, thread_id=1, iseedxx(512)

!read input file
       namelist /EZmock_v0_input/
     & datafile_prefix,
     &  boxsize, grid_num,
     & grow2z0, scatter, scatter2, modify_pk, modify_pdf,
     & antidamping, density_cut, density_sat, zdist_rate,
     & zdist_fog, whitenoise_file, use_whitenoise_file,
     & pkfile, pknwfile, datafile_path, output_suffix,
     & compute_CF, compute_CF_zdist,
     & skiplines, twod_corr_suffix, max_r, bin_size,
     & redshift, om,  expect_sum_pdf, expect_A_pdf, iseed,dilute_factor

!input                                                                                                                 
      read(*,EZmock_v0_input)

      write(*,*) "pr_EZmock_pk_cf read input:",
     & datafile_prefix,
     & boxsize, grid_num,
     & grow2z0, scatter, scatter2, modify_pk, modify_pdf,
     & antidamping, density_cut, density_sat, zdist_rate,
     & zdist_fog, whitenoise_file, use_whitenoise_file,
     & pkfile, pknwfile, datafile_path, output_suffix,
     & compute_CF, compute_CF_zdist,
     & skiplines, twod_corr_suffix, max_r, bin_size,
     &     redshift, om,  expect_sum_pdf, expect_A_pdf, iseed,
     & dilute_factor 

      time1=time()

      grid_pdf=grid_num

      z = redshift
      Hubble = 100.*(sqrt((1.+z)*(1.+z)*(1.+z)*om+1.-om))! h units
      inverse_aH=(1.+z)/Hubble
      write(*,*) "1/aH (h-unit) =", inverse_aH

      output_file = Trim(datafile_path)//Trim(datafile_prefix)
     &   //Trim(output_suffix)

        call zeldovich_ic(
     &   skiplines,pkfile,pknwfile,
     &  grow2z0, boxsize,grid_num,grid_pdf,
     &  data_density, iseed,data_num, output_file,max_data,scatter,
     &  scatter2,modify_pk,modify_pdf,antidamping,density_cut,
     &  density_sat, zdist_rate, zdist_fog, inverse_aH,
     &  whitenoise_file, use_whitenoise_file,
     &     expect_sum_pdf,expect_A_pdf)
        
        write(*,*)
     &      "finish generating catalog in cubic box. Size:",boxsize,
     &  "output file:", Trim(output_file)

!---  insert 2PCF code here
! initialize
      num_r=int(max_r)

      if(compute_CF .eq. .true.) then
        write(*,*) "start to compute 2PCF in real space"

! read data file
      open(unit=30,file=Trim(output_file),status='old')
      do i=1,skiplines
        read(30,*)
      end do
      data_num = 0
 300  read(30,*,end=301) r1,r2,r3!,r4
       data_num = data_num+1
      goto 300 
 301  close(30)

      allocate(x(data_num))
      allocate(y(data_num))
      allocate(w(data_num))
      allocate(wt(data_num))
      allocate(vz(data_num))
      allocate(DDnofr(num_r,num_mu))
      allocate(RRnofr(num_r,num_mu))

      z = redshift
      open(unit=30,file=Trim(output_file),status='old')
      do i=1,skiplines
        read(30,*)
      end do

      do i=1,data_num
          read(30,*) x(i),y(i),w(i)
          wt(i)=1
      end do
      close(30)
      write(*,*) "finish reading data to arrays"

!dilute sample to make 2PCF faster
      i_dilute = 1
      write(*,*) "dilute sample with factor",dilute_factor
      do i=1, data_num
         if(ran3_omp(iseedxx(thread_id),thread_id).lt.dilute_factor)then
         x(i_dilute) = x(i)
         y(i_dilute) = y(i)
         w(i_dilute) = w(i)
         i_dilute = i_dilute + 1
         end if 
      end do
      i_dilute = i_dilute-1

!      call fast_code_CF_box_linear(x,y,w,wt,data_num,
      call fast_code_CF_box_linear(x,y,w,wt,i_dilute,
     &  boxsize,num_r,num_mu,ddnofr,rrnofr)

  
!output:  

!monopole
      open(unit =11, file=Trim(output_file)//
     &         Trim(twod_corr_suffix)//".mono")
      do rr = 1, num_r,bin_size
          ddbin = 0
          d2rbin = 0
          rrbin = 0
          do i = rr, rr+bin_size-1
            do mm = 1, num_mu
              ddbin  = ddbin + ddnofr(i,mm)
              rrbin  = rrbin + rrnofr(i,mm)
            end do
          end do
          ddbin=ddbin/i_dilute/i_dilute
          write(11,*) 0.5*bin_size+rr-1, ddbin/rrbin-1
      end do   
      close(11)

      end if

      if(compute_CF_zdist .eq. .true.) then
         write(*,*) "start to compute 2PCF in redshift space"
! redshift space
      open(unit=30,file=Trim(output_file),status='old')
      do i=1,skiplines
        read(30,*)
      end do

      do i=1,data_num
          read(30,*) x(i),y(i),w(i),r1,r2,vz(i)
          wt(i)=1
          w(i)=w(i)+vz(i)*inverse_aH

          if(w(i).gt.boxsize) w(i) = w(i) - boxsize
          if(w(i).lt.0.) w(i) = w(i) + boxsize;
      end do

      close(30)
      write(*,*) "finish reading data to arrays"

!     dilute sample to make 2PCF faster
      write(*,*) "dilute sample with factor", dilute_factor
      i_dilute = 1
      do i=1, data_num
       if(ran3_omp(iseedxx(thread_id),thread_id) .lt. dilute_factor)then
         x(i_dilute) = x(i)
         y(i_dilute) = y(i)
         w(i_dilute) = w(i)
         i_dilute = i_dilute + 1
         end if 
      end do
      i_dilute = i_dilute-1

      call fast_code_CF_box_linear(x,y,w,wt,i_dilute,
     &  boxsize,num_r,num_mu,ddnofr,rrnofr)

  
!output:  

!monopole
      open(unit =11, file=Trim(output_file)//
     &         Trim(twod_corr_suffix)//".zdist.mono")
      do rr = 1, num_r,bin_size
          ddbin = 0
          d2rbin = 0
          rrbin = 0
          do i = rr, rr+bin_size-1
            do mm = 1, num_mu
              ddbin  = ddbin + ddnofr(i,mm)
              rrbin  = rrbin + rrnofr(i,mm)
            end do
          end do
          ddbin=ddbin/i_dilute/i_dilute
          write(11,*) 0.5*bin_size+rr-1, ddbin/rrbin-1
      end do   
      close(11)

!quadrupole
      open(unit =11, file=Trim(output_file)//
     &         Trim(twod_corr_suffix)//".zdist.quad")
      do rr = 1, num_r,bin_size
        ddbin = 0
        d2rbin = 0
        rrbin = 0
        rrbin_de = 0
        do i = rr, rr+bin_size-1
          do mm = 1, num_mu
        ddbin = ddbin + ddnofr(i,mm)*(3.*(1./num_mu*(-0.5+mm))**2-1)/2*5
        rrbin = rrbin + rrnofr(i,mm)*(3.*(1./num_mu*(-0.5+mm))**2-1)/2*5
        rrbin_de = rrbin_de + rrnofr(i,mm)
          end do
        end do
        ddbin=ddbin/i_dilute/i_dilute
        write(11,*) 0.5*bin_size+rr-1, (ddbin-rrbin)/rrbin_de
      end do   
      close(11)

      deallocate(x)
      deallocate(y)
      deallocate(w)
      deallocate(wt)
      deallocate(ddnofr)
      deallocate(rrnofr)
!-----
      time2 =  time()
      write(*,*) "time in sec:", time2-time1
      end if

      end program

!-----------------------------------------------------------------

      subroutine zeldovich_ic(  !datafile_xyz,
     &  skip_line,pkfile,
     &     pknwfile,grow2z0,boxsize,grid_num,grid_pdf,density,iseed,
     &     data_num, output_file,max_data,scatter,scatter2,modify_pk,
     &     modify_pdf, antidamping, density_cut,density_sat,zdist_rate,
     &     zdist_fog, inverse_aH,whitenoise_file, use_whitenoise_file,
     &     expect_sum_pdf, expect_A_pdf)
      implicit none
!      include 'fftw_f77.i'
!      FFTW3
      include 'fftw3.f'
!      character*500,intent(in) :: datafile_xyz
      integer,intent(in)          :: skip_line
      character*500,intent(in) :: pkfile
      character*500,intent(in) :: pknwfile
      real,intent(in)          :: grow2z0
      real,intent(in)          :: boxsize
      integer,intent(in)       :: grid_num
      real,intent(in)          :: density
      integer,intent(in)       :: iseed
      integer,intent(out)      :: data_num
      character*500,intent(in) :: output_file
      integer,intent(in)       :: grid_pdf
      integer,intent(in)       :: max_data
      real,intent(in)          :: scatter
      real,intent(in)          :: scatter2
      real,intent(in)          :: modify_pk
      real,intent(in)          :: modify_pdf
      real,intent(in)          :: antidamping
      real,intent(in)          :: density_cut
      real,intent(in)          :: density_sat
      real,intent(in)          :: zdist_rate
      real,intent(in)          :: zdist_fog
      real,intent(in)          :: inverse_aH
      character*500,intent(in) :: whitenoise_file
      logical,intent(in)       :: use_whitenoise_file
!     add pdf manually
      real,intent(in)          :: expect_sum_pdf
      real,intent(in)          :: expect_A_pdf

!      logical,parameter :: zdist_use_modified_pk = .true.

!      logical,parameter :: zeldovich = .true.
!compute pk
      real    :: k_size = 0.005
      logical :: CIC_assignment = .true.
      character*500 :: pkfile_path =
     & "../pk0421/"
      character*500 :: pkfile_prefix =
     & "dk0.005"
!
      logical :: CIC_binning = .false. !I have shown that it's useless
      integer,parameter :: r_num_max = 10000, k_num_max = 10000
      real,parameter    :: kmin = 0.0001, kmax = 100
      real,parameter    :: pi = 3.1415927
!      integer,parameter :: skip_line=2
      character*500     :: sim_pdf_file
      character*500     :: mock_pkfile
      integer :: n_cut
      real    :: nfactor = 1, sigma_G2
      real    :: r_array(r_num_max),corr_array(r_num_max)
      real    :: corr_arrayo(r_num_max)
      real    :: corr_array2(2000)
      real    :: k_array(k_num_max), pk_array(k_num_max), y2(k_num_max)
      real    :: pk_arrayo(k_num_max), y2o(k_num_max) 
      real    :: pk_nw_array(k_num_max)
      integer :: i, r_num, k_num, j
      real    :: volume
      integer :: nx, ny, nz, rx, ry, rz, nx_conj, ny_conj, nz_conj
      real    :: kx, ky, kz, k_mag, k_mag2, pk_out, temp, ak, bk,pk_outo
      real    :: kx1, ky1, kz1, re_in, im_in
      integer :: np, data_index(100)
      integer :: nn(3), ndim=3, isign=-1
      real    :: gridsize, ran3_omp, xx, yy, zz, ran3, temp3,temp4
      integer :: sirko
      integer :: totalcount = 0
      real    :: rescaleG, rescaleLN, temp2, sigma2_densityG, tempo
      integer ::  readtime, xi, yi, zi,xii,yii,zii !, iseed2=-1234
      integer :: subx, suby, subz
c... for fftw
      integer*8  :: plan,planx,plany,planz,planm,palno,planw
      integer*8  :: plan2,planx2,plany2,planz2,planw2
!big arrays
      complex,allocatable, dimension(:,:,:) ::
     & vxarr_in, vyarr_in, vzarr_in
      real, allocatable, dimension(:,:,:) ::
     &     temp_array, temp_array2, temp_array3!,temp_array4
      integer, allocatable, dimension(:,:,:) :: temp_int_array
      integer,parameter :: max_num = 200
      integer :: pdf(0:max_num)
      real :: x, y, z, cell_size, count_real, vx, vy, vz
      integer :: k, ix, iy, iz, count, count2
!
      integer :: sim_pdf(0:max_num) !keep pdf from data (simulaiton) -- the one used for mapping
      integer :: sim_pdf2(0:max_num) !keep pdf from data (simulaiton) -- the one computed for cic
      real    :: pdf_ratio(max_num) !probability to keep each halo in a cell with given number of halos
      integer,parameter :: max_vx_pdf=2000
      real    :: sim_vx_pdf(-max_vx_pdf:max_vx_pdf)
      real    :: za_vx_pdf(-max_vx_pdf:max_vx_pdf)
      integer :: za_sum, za_temp, za_pdf_index, sim_pdf_max
      integer :: temp_int
      real    :: pdf2(0:max_num) 
      integer :: pdf_convert(0:max_num)
      real :: pdf_convert_ratio(0:max_num)
      real :: np_ratio, exp
      real,allocatable :: dataxyz(:,:)
      real,allocatable :: dataxyz_omp(:,:,:) !(6,data,thread_num)
      real :: dataxyz_temp(6)
      logical :: dataxyz_keep(max_data)
      integer :: time, time1, time2

      integer,parameter :: k_index_max = 10000
      real              :: pk(k_index_max),pk_quad(k_index_max)
      integer ::  num_data, mock_data_num, i1,i2,i3,np1,np2,np3
      integer :: mock_data_num2
      real    :: mu,kny, pk_temp, dr, norm1
      integer :: k_index,nozero_k_index_max, box_size
      real    :: k_num_sep(k_index_max), quad_modes(k_index_max)
      integer :: num_to_make
      
!openmp
      integer :: fftw_init_threads, omp_get_max_threads, iret
      integer :: omp_get_thread_num, thread_id, iseedxx(512) !iseed for threads
      integer :: rxx, ryy, rzz, rxi,ryi,rzi, zzi, zzii
      integer :: x_div_omp, y_div_omp, z_div_omp, total_num_thread
      real :: density_total, density_cut_total
!adjust pdf
      real :: A_pdf, B_pdf, sum_pdf,sum_pdf_A


!intitalize for doing cic with openmp (total_num_threads = x_div_omp*y_div_omp*z_div_omp)
      total_num_thread = omp_get_max_threads()
      x_div_omp = 2
      y_div_omp = 2

!      allocate(vxarr_in(grid_num,grid_num,grid_num))

      if(mod(total_num_thread, x_div_omp*y_div_omp) .ne. 0) then
         write(*,*) "total threads is not integer times 2*2"
         stop  
      end if

      z_div_omp = total_num_thread/x_div_omp/y_div_omp
!      if(mod(grid_num,x_div_omp).ne.0 .or.
!     & mod(grid_num,y_div_omp).ne.0 .or.
!     & mod(grid_num,z_div_omp).ne.0) then
      if(mod(grid_num,x_div_omp*2).ne.0 .or.
     & mod(grid_num,y_div_omp*2).ne.0 .or.
     & mod(grid_num,z_div_omp*2).ne.0) then
          write(*,*) 
     & "grid size is not proper divided by x_div_omp(or y or z) *2"
       stop
      end if

      write(*,*)
      time1 = time() 
      write(*,*) "start count time in sec:", 0 

      allocate(vxarr_in(grid_num,grid_num,grid_num)) !vxarr_in
!initial iseed for  openmp
C$OMP PARALLEL PRIVATE(thread_id)
      thread_id =  omp_get_thread_num()+1
      iseedxx(thread_id) =  iseed*1000 + thread_id !number of threads must smaller than 1000
C$OMP END PARALLEL

!! initialize
!C$OMP PARALLEL DO PRIVATE(rx,ry,rz)
!      do rx=1,grid_num
!        do ry=1,grid_num
!          do rz=1,grid_num
!            arr_in(rx,ry,rz)=cmplx(0,0)
!           end do
!        end do
!      end do

      volume = 1 
      data_num = 0
      gridsize = 1.*boxsize/grid_num
      n_cut = grid_num/2 !no cut

      cell_size = boxsize/grid_pdf
      do i=0, max_num
        pdf(i) = 0
        pdf2(i) = 0
        sim_pdf(i) =0
        sim_pdf2(i) =0
        pdf_convert(i)=0
        pdf_convert_ratio(i)=0
      end do

      sim_pdf_file = Trim(output_file)//'.pdf.cic.sim'

      write(*,*) "finish initialization"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

      sum_pdf = expect_sum_pdf
      A_pdf = expect_A_pdf
      
      sum_pdf_A = 0
      do i=1, 20
        sum_pdf_A = sum_pdf_A + A_pdf**i*i
      end do
      B_pdf = sum_pdf/sum_pdf_A
      write(*,*) "new A, B of PDF are", A_pdf, B_pdf
      write(*,*) "new sum of iA^i, A, B", sum_pdf_A, A_pdf, B_pdf

      sum_pdf=0
      temp_int = 2 !initial with any number > 1
      i=1
!      do while (temp_int .ge. 1.)
!        temp_int = int(B_pdf*A_pdf**i)
!        write(*,*) i, temp_int
!        sim_pdf2(i) = temp_int
!        sum_pdf = sum_pdf + temp_int*i
!        i = i+1
!      end do
      do i=1, 20
        temp_int = int(B_pdf*A_pdf**i)
        write(*,*) i, temp_int
        sim_pdf2(i) = temp_int
        sum_pdf = sum_pdf + temp_int*i
      end do

      write(*,*) "sum of actual AB", sum_pdf
!   rewrite pdf of sim
      write(*,*) "new pdf:"
      open(10, file=Trim(sim_pdf_file)//'.new')
      do i = 0, max_num
         write(10,*) i, sim_pdf2(i)
         if(sim_pdf2(i) .gt. 0) write(*,*) i, sim_pdf2(i)
      end do
      close(10)

!read input pk
      open(unit=10, file=Trim(pkfile), status='old')
      do i = 1, k_num_max
         read(10,*,end=300) k_array(i), pk_arrayo(i)
         pk_arrayo(i) = pk_arrayo(i)*grow2z0
      end do
 300  close(10)
      k_num = i-1

      open(unit=10, file=Trim(pknwfile), status='old')
      do i = 1, k_num
         read(10,*) k_array(i), pk_nw_array(i)
         pk_nw_array(i) = pk_nw_array(i)*grow2z0
      end do
      close(10)

      write(*,*) "read pk func with ", k_num, "points"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

      write(*,*) "modify input pk here"
      do i=1,k_num

        if(antidamping .lt. 1) then !turn off BAO modification if antidamping > 1
         if(k_array(i).lt. 0.5)
     &    pk_arrayo(i) = (pk_arrayo(i)-pk_nw_array(i))*
     &          exp(k_array(i)**2/antidamping) + pk_nw_array(i)
        end if
        pk_arrayo(i) = pk_arrayo(i)*
     &    (1.
     &   + (modify_pk*k_array(i))**(1.)
     &     )
      end do

! density on grid array
! to initial the pk function
      call spline(k_array,pk_arrayo,k_num,1.0e30,1.0e30,y2o)

!generate white noise on vxarr_in or read the white noise file
      if(use_whitenoise_file .eq. .false.) then

        thread_id = 1
        do rx = 1, grid_num
          do ry = 1, grid_num
            do rz = 1, grid_num,2  !might use every 2 since gaussion generator provide 2
              call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)
              vxarr_in(rx,ry,rz) = cmplx(ak,0)
              vxarr_in(rx,ry,rz+1) = cmplx(bk,0)
            end do
          end do
        end do
        write(*,*) "finish generating white noise"
        time2 =  time()
        write(*,*) "time in sec:", time2-time1

      else
        write(*,*) "reading white noise file"
        open(10,file=Trim(whitenoise_file),form='unformatted')
        read (10) np1,np2,np3, temp_int
        print*, 'np1,np2,np3,iseed'
        print*, np1,np2,np3,temp_int
!        if (np1 .ne. grid_num) then
!          write(*,*)"the grid size does not match white noise file:",np1 
!          stop
!        end if
        allocate(temp_array(grid_num,grid_num,grid_num)) !vxarr_in, temp_array
        do i3=1,grid_num
         read(10) ((temp_array(i1,i2,i3),i1=1,grid_num),i2=1,grid_num)
        enddo
        close(10)


!read ascii file
!        open(10,file=Trim(whitenoise_file))
!        do i3=1,grid_num
!          do i2=1,grid_num
!         read(10,*) (temp_array(i1,i2,i3),i1=1,grid_num)
!          end do
!        enddo
        close(10)
        write(*,*) "finish reading white noise file"
        time2 =  time()
        write(*,*) "time in sec:", time2-time1

!move white noise to vxarr_in
C$omp parallel do private(rx,ry,rz)
        do rx = 1, grid_num
          do ry = 1, grid_num
            do rz = 1, grid_num
              vxarr_in(rx,ry,rz) = cmplx(temp_array(rx,ry,rz),0)
            end do
          end do
        end do
        deallocate(temp_array) !left: vxarr_in(real part saves white noise)
      end if
     
      
!!for openmp
      call sfftw_init_threads(iret)
      if(iret .eq. 0) then
          write(*,*) "something wrong with using fftw and openmp"
          stop
      end if
      write(*,*) "number of threads=",omp_get_max_threads()

      call sfftw_plan_with_nthreads(omp_get_max_threads())

      call sfftw_plan_dft_3d(plan, grid_num, grid_num, grid_num,
     &  vxarr_in, vxarr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call sfftw_execute_dft(plan,vxarr_in,vxarr_in)
      call sfftw_destroy_plan(plan)


      allocate(vyarr_in(grid_num,grid_num,grid_num))
      allocate(vzarr_in(grid_num,grid_num,grid_num)) !left: vxarr_in(white,0), vyarr_in, vzarr_in
     
      temp = sqrt(1.0*grid_num**3/2) !to normalize white noise with given grid
!      vxarr_in(1,1,1) = cmplx(Real(arr_in(1,1,1))/temp,0)
C$omp parallel private(nx,ny,nz,kx,ky,kz,k_mag,k_mag2,pk_outo,tempo,
C$omp+ak,bk,nx_conj,ny_conj,nz_conj,thread_id)
!!reduction(+:sigma_G2)
      thread_id = omp_get_thread_num()+1
!      sigma_G2 = 0
C$omp do
      do nx = 0, grid_num-1
        do ny = 0, grid_num-1
          do nz = 0, grid_num-1
            if((nx .eq. 0) .and. (ny .eq. 0) .and. (nz .eq. 0)) then
              vxarr_in(1,1,1) = cmplx(0,0)
              vyarr_in(1,1,1) = cmplx(0,0)
              vzarr_in(1,1,1) = cmplx(0,0)
            else
              kx = 2.0*pi/boxsize*nx
              ky = 2.0*pi/boxsize*ny
              kz = 2.0*pi/boxsize*nz
              if (nx .gt. grid_num/2) kx = 2.0*pi/boxsize*(nx-grid_num)
              if (ny .gt. grid_num/2) ky = 2.0*pi/boxsize*(ny-grid_num)
              if (nz .gt. grid_num/2) kz = 2.0*pi/boxsize*(nz-grid_num)
              k_mag2 = kx**2 + ky**2 + kz**2
              k_mag = sqrt(k_mag2)
              call splint(k_array,pk_arrayo,y2o,k_num,k_mag,pk_outo)
              tempo = sqrt(pk_outo/boxsize**3/2)
!              sigma_G2 = sigma_G2 + pk_outo*2 !*2 for nz < 0

              ak = Real(vxarr_in(nx+1,ny+1,nz+1))/temp
              bk = AIMAG(vxarr_in(nx+1,ny+1,nz+1))/temp
              vxarr_in(nx+1,ny+1,nz+1)=
     &              cmplx(tempo*bk*kx/k_mag2,-tempo*ak*kx/k_mag2)
              vyarr_in(nx+1,ny+1,nz+1)=
     &              cmplx(tempo*bk*ky/k_mag2,-tempo*ak*ky/k_mag2)
              vzarr_in(nx+1,ny+1,nz+1)=
     &              cmplx(tempo*bk*kz/k_mag2,-tempo*ak*kz/k_mag2)

            end if
          end do !kz
        end do !ky
      end do !kx
C$omp end do
C$omp end parallel

!      sigma_G2 = sigma_G2/boxsize**3
      write(*,*) "finish assigning pk array"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

      write(*,*) "start FFT"
!!for openmp
!      call sfftw_init_threads(iret)
!      if(iret .eq. 0) then
!          write(*,*) "something wrong with using fftw and openmp"
!          stop
!      end if
!      write(*,*) "number of threads=",omp_get_max_threads()
!      call sfftw_plan_with_nthreads(omp_get_max_threads())

!      FFTW3
      call sfftw_plan_dft_3d(planx, grid_num, grid_num, grid_num,
     &  vxarr_in, vxarr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call sfftw_execute_dft(planx,vxarr_in,vxarr_in)
      call sfftw_destroy_plan(planx)

      time2 =  time()
      write(*,*) "1st FFT, time in sec:", time2-time1

!      FFTW3
      call sfftw_plan_dft_3d(plany, grid_num, grid_num, grid_num,
     &  vyarr_in, vyarr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call sfftw_execute_dft(plany,vyarr_in,vyarr_in)
      call sfftw_destroy_plan(plany)

      time2 =  time()
      write(*,*) "2nd FFT, time in sec:", time2-time1

!      FFTW3
      call sfftw_plan_dft_3d(planz, grid_num, grid_num, grid_num,
     &  vzarr_in, vzarr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call sfftw_execute_dft(planz,vzarr_in,vzarr_in)
      call sfftw_destroy_plan(planz)

      time2 =  time()
      write(*,*) "3rd FFT, time in sec:", time2-time1

!--- displacement field using CIC ----------------------- raw zeldovich particles
! using image part of array for CIC counting
!      write(*,*) "count initial CIC for za particles"
      open(17,file=Trim(output_file)//'.raw')

C$omp parallel do private(rx,ry,rz)
      do rx = 1,grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
              vxarr_in(rx,ry,rz) = cmplx(real(vxarr_in(rx,ry,rz)),0)
            end do
         end do
      end do
C$omp parallel do private(rx,ry,rz)
      do rx = 1,grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
              vyarr_in(rx,ry,rz) = cmplx(real(vyarr_in(rx,ry,rz)),0)
            end do
         end do
      end do
C$omp parallel do private(rx,ry,rz)
      do rx = 1,grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
              vzarr_in(rx,ry,rz) = cmplx(real(vzarr_in(rx,ry,rz)),0)
            end do
         end do
      end do

       write(*,*) "start to construct ZA density array"
        time2 =  time()
        write(*,*) "time in sec:", time2-time1

        allocate(temp_array2(grid_num,grid_num,grid_num))
        allocate(temp_array3(grid_num,grid_num,grid_num))

C$omp parallel do private(rx,ry,rz)
      do rx = 1,grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
              temp_array2(rx,ry,rz) = 0
            end do
         end do
      end do
C$omp parallel do private(rx,ry,rz)
      do rx = 1,grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
              temp_array3(rx,ry,rz) = 0
            end do
         end do
      end do

!     
!assign ZA particles to grid using NGP & CIC
!! no openmp here because of writing ZA particles in the loop
!!C$omp parallel do private(rx,ry,rz,xx,yy,zz,ix,iy,iz,
!!C$omp+ x,y,z,xi,yi,zi,xii,yii,zii,i,rxx,ryy,rzz,rxi,ryi,rzi)
      do i=0, total_num_thread-1
         rxx = i/(y_div_omp*z_div_omp)
         ryy = mod(i,y_div_omp*z_div_omp)/z_div_omp
         rzz = mod(i,z_div_omp)

         do  subx= 0, 1
          do suby = 0, 1
           do subz = 0, 1
         
        do rxi = 1, grid_num/x_div_omp/2
          do ryi = 1, grid_num/y_div_omp/2
            do rzi = 1, grid_num/z_div_omp/2
                
              rx = rxx*grid_num/x_div_omp+rxi+subx*grid_num/x_div_omp/2
              ry = ryy*grid_num/y_div_omp+ryi+suby*grid_num/y_div_omp/2
              rz = rzz*grid_num/z_div_omp+rzi+subz*grid_num/z_div_omp/2

               xx = real(vxarr_in(rx,ry,rz)) + gridsize*rx
               if (xx .ge. 0) then
                   xx = xx - int(xx/boxsize)*boxsize
               else
                   xx = xx - (int(xx/boxsize)-1)*boxsize
               end if
               yy = real(vyarr_in(rx,ry,rz)) + gridsize*ry
               if (yy .ge. 0) then
                   yy = yy - int(yy/boxsize)*boxsize
               else
                   yy = yy - (int(yy/boxsize)-1)*boxsize
               end if
               zz = real(vzarr_in(rx,ry,rz)) + gridsize*rz
               if (zz .ge. 0) then
                   zz = zz - int(zz/boxsize)*boxsize
               else
                   zz = zz - (int(zz/boxsize)-1)*boxsize
               end if

!test generate raw za particles
               if(ran3(iseed) .lt. density) then
                  write(17,*) xx,yy,zz
               end if
!

! for NGP               
               ix = int(xx/cell_size+0.5)+1
               if (ix .gt. grid_pdf) ix=1
               iy = int(yy/cell_size+0.5)+1
               if (iy .gt. grid_pdf) iy=1
               iz = int(zz/cell_size+0.5)+1
               if (iz .gt. grid_pdf) iz=1
! for CIC
               xi = int(xx/cell_size)+1
               if (xi .gt. grid_pdf) xi=grid_pdf
               yi = int(yy/cell_size)+1
               if (yi .gt. grid_pdf) yi=grid_pdf
               zi = int(zz/cell_size)+1
               if (zi .gt. grid_pdf) zi=grid_pdf
               xii = xi + 1
               if(xii.eq.grid_pdf+1) xii=1
               yii = yi + 1
               if(yii.eq.grid_pdf+1) yii=1
               zii = zi + 1
               if(zii.eq.grid_pdf+1) zii=1
               x = xx
               y = yy
               z = zz
               
! avoid boundary collision
       if( real(vxarr_in(rx,ry,rz)).gt.(boxsize/x_div_omp) .or.
     &  real(vyarr_in(rx,ry,rz)).gt.(boxsize/y_div_omp) .or.
     &   real(vzarr_in(rx,ry,rz)).gt.(boxsize/z_div_omp) ) then
  
          !$omp critical
!!NGP          
!               vyarr_in(ix,iy,iz) = cmplx(Real(vyarr_in(ix,iy,iz)),
!     &         AIMAG(vyarr_in(ix,iy,iz)) + 1)
!CIC
               vxarr_in(xi,yi,zi) = cmplx(Real(vxarr_in(xi,yi,zi)),
     &    AIMAG(vxarr_in(xi,yi,zi))+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi)))
               vxarr_in(xii,yi,zi) = cmplx(Real(vxarr_in(xii,yi,zi)),
     &   AIMAG(vxarr_in(xii,yi,zi))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi)))
               vxarr_in(xi,yii,zi) = cmplx(Real(vxarr_in(xi,yii,zi)),
     &    AIMAG(vxarr_in(xi,yii,zi))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi)))
               vxarr_in(xi,yi,zii) = cmplx(Real(vxarr_in(xi,yi,zii)),
     &    AIMAG(vxarr_in(xi,yi,zii))+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1)))
             vxarr_in(xii,yii,zi) = cmplx(Real(vxarr_in(xii,yii,zi)),
     &    AIMAG(vxarr_in(xii,yii,zi))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi)))
             vxarr_in(xi,yii,zii) = cmplx(Real(vxarr_in(xi,yii,zii)),
     &    AIMAG(vxarr_in(xi,yii,zii))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1)))
             vxarr_in(xii,yi,zii) = cmplx(Real(vxarr_in(xii,yi,zii)),
     &    AIMAG(vxarr_in(xii,yi,zii))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1)))
             vxarr_in(xii,yii,zii) =cmplx(Real(vxarr_in(xii,yii,zii)),
     &    AIMAG(vxarr_in(xii,yii,zii))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &            (z-cell_size*(zi-1)))

!CIC of z-axis velocity field
               vyarr_in(xi,yi,zi) = cmplx(Real(vyarr_in(xi,yi,zi)),
     &    AIMAG(vyarr_in(xi,yi,zi))+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vzarr_in(rx,ry,rz)))
               vyarr_in(xii,yi,zi) = cmplx(Real(vyarr_in(xii,yi,zi)),
     &   AIMAG(vyarr_in(xii,yi,zi))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vzarr_in(rx,ry,rz)))
               vyarr_in(xi,yii,zi) = cmplx(Real(vyarr_in(xi,yii,zi)),
     &    AIMAG(vyarr_in(xi,yii,zi))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vzarr_in(rx,ry,rz)))
               vyarr_in(xi,yi,zii) = cmplx(Real(vyarr_in(xi,yi,zii)),
     &    AIMAG(vyarr_in(xi,yi,zii))+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vzarr_in(rx,ry,rz)))
             vyarr_in(xii,yii,zi) = cmplx(Real(vyarr_in(xii,yii,zi)),
     &    AIMAG(vyarr_in(xii,yii,zi))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vzarr_in(rx,ry,rz)))
             vyarr_in(xi,yii,zii) = cmplx(Real(vyarr_in(xi,yii,zii)),
     &    AIMAG(vyarr_in(xi,yii,zii))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vzarr_in(rx,ry,rz)))
             vyarr_in(xii,yi,zii) = cmplx(Real(vyarr_in(xii,yi,zii)),
     &    AIMAG(vyarr_in(xii,yi,zii))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vzarr_in(rx,ry,rz)))
             vyarr_in(xii,yii,zii) =cmplx(Real(vyarr_in(xii,yii,zii)),
     &    AIMAG(vyarr_in(xii,yii,zii))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vzarr_in(rx,ry,rz)))

!CIC of x-axis velocity field
               temp_array2(xi,yi,zi) = 
     &    temp_array2(xi,yi,zi)+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xii,yi,zi) =
     &   temp_array2(xii,yi,zi)+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xi,yii,zi) =
     &    temp_array2(xi,yii,zi)+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xi,yi,zii) =
     &    temp_array2(xi,yi,zii)+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xii,yii,zi) =
     &    temp_array2(xii,yii,zi)+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xi,yii,zii) =
     &    temp_array2(xi,yii,zii)+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xii,yi,zii) =
     &    temp_array2(xii,yi,zii)+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xii,yii,zii) =
     &    temp_array2(xii,yii,zii)+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vxarr_in(rx,ry,rz))

!CIC of y-axis velocity field
               temp_array3(xi,yi,zi) = 
     &    temp_array3(xi,yi,zi)+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xii,yi,zi) =
     &   temp_array3(xii,yi,zi)+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xi,yii,zi) =
     &    temp_array3(xi,yii,zi)+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xi,yi,zii) =
     &    temp_array3(xi,yi,zii)+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xii,yii,zi) =
     &    temp_array3(xii,yii,zi)+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xi,yii,zii) =
     &    temp_array3(xi,yii,zii)+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xii,yi,zii) =
     &    temp_array3(xii,yi,zii)+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xii,yii,zii) =
     &    temp_array3(xii,yii,zii)+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vyarr_in(rx,ry,rz))
             
               
          !$omp end critical                        
          else
!!NGP             
!               vyarr_in(ix,iy,iz) = cmplx(Real(vyarr_in(ix,iy,iz)),
!     &    AIMAG(vyarr_in(ix,iy,iz)) + 1)
!CIC
               vxarr_in(xi,yi,zi) = cmplx(Real(vxarr_in(xi,yi,zi)),
     &    AIMAG(vxarr_in(xi,yi,zi))+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi)))
               vxarr_in(xii,yi,zi) = cmplx(Real(vxarr_in(xii,yi,zi)),
     &   AIMAG(vxarr_in(xii,yi,zi))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi)))
               vxarr_in(xi,yii,zi) = cmplx(Real(vxarr_in(xi,yii,zi)),
     &    AIMAG(vxarr_in(xi,yii,zi))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi)))
               vxarr_in(xi,yi,zii) = cmplx(Real(vxarr_in(xi,yi,zii)),
     &    AIMAG(vxarr_in(xi,yi,zii))+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1)))
             vxarr_in(xii,yii,zi) = cmplx(Real(vxarr_in(xii,yii,zi)),
     &    AIMAG(vxarr_in(xii,yii,zi))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi)))
             vxarr_in(xi,yii,zii) = cmplx(Real(vxarr_in(xi,yii,zii)),
     &    AIMAG(vxarr_in(xi,yii,zii))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1)))
             vxarr_in(xii,yi,zii) = cmplx(Real(vxarr_in(xii,yi,zii)),
     &    AIMAG(vxarr_in(xii,yi,zii))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1)))
             vxarr_in(xii,yii,zii) =cmplx(Real(vxarr_in(xii,yii,zii)),
     &    AIMAG(vxarr_in(xii,yii,zii))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1)))


!CIC of z-axis velocity field
               vyarr_in(xi,yi,zi) = cmplx(Real(vyarr_in(xi,yi,zi)),
     &    AIMAG(vyarr_in(xi,yi,zi))+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vzarr_in(rx,ry,rz)))
               vyarr_in(xii,yi,zi) = cmplx(Real(vyarr_in(xii,yi,zi)),
     &   AIMAG(vyarr_in(xii,yi,zi))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vzarr_in(rx,ry,rz)))
               vyarr_in(xi,yii,zi) = cmplx(Real(vyarr_in(xi,yii,zi)),
     &    AIMAG(vyarr_in(xi,yii,zi))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vzarr_in(rx,ry,rz)))
               vyarr_in(xi,yi,zii) = cmplx(Real(vyarr_in(xi,yi,zii)),
     &    AIMAG(vyarr_in(xi,yi,zii))+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vzarr_in(rx,ry,rz)))
             vyarr_in(xii,yii,zi) = cmplx(Real(vyarr_in(xii,yii,zi)),
     &    AIMAG(vyarr_in(xii,yii,zi))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vzarr_in(rx,ry,rz)))
             vyarr_in(xi,yii,zii) = cmplx(Real(vyarr_in(xi,yii,zii)),
     &    AIMAG(vyarr_in(xi,yii,zii))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vzarr_in(rx,ry,rz)))
             vyarr_in(xii,yi,zii) = cmplx(Real(vyarr_in(xii,yi,zii)),
     &    AIMAG(vyarr_in(xii,yi,zii))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vzarr_in(rx,ry,rz)))
             vyarr_in(xii,yii,zii) =cmplx(Real(vyarr_in(xii,yii,zii)),
     &    AIMAG(vyarr_in(xii,yii,zii))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vzarr_in(rx,ry,rz)))

!CIC of x-axis velocity field
               temp_array2(xi,yi,zi) = 
     &    temp_array2(xi,yi,zi)+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xii,yi,zi) =
     &   temp_array2(xii,yi,zi)+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xi,yii,zi) =
     &    temp_array2(xi,yii,zi)+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xi,yi,zii) =
     &    temp_array2(xi,yi,zii)+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xii,yii,zi) =
     &    temp_array2(xii,yii,zi)+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xi,yii,zii) =
     &    temp_array2(xi,yii,zii)+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xii,yi,zii) =
     &    temp_array2(xii,yi,zii)+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vxarr_in(rx,ry,rz))
               temp_array2(xii,yii,zii) =
     &    temp_array2(xii,yii,zii)+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vxarr_in(rx,ry,rz))

!CIC of y-axis velocity field
               temp_array3(xi,yi,zi) = 
     &    temp_array3(xi,yi,zi)+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xii,yi,zi) =
     &   temp_array3(xii,yi,zi)+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xi,yii,zi) =
     &    temp_array3(xi,yii,zi)+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xi,yi,zii) =
     &    temp_array3(xi,yi,zii)+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xii,yii,zi) =
     &    temp_array3(xii,yii,zi)+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xi,yii,zii) =
     &    temp_array3(xi,yii,zii)+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xii,yi,zii) =
     &    temp_array3(xii,yi,zii)+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Real(vyarr_in(rx,ry,rz))
               temp_array3(xii,yii,zii) =
     &    temp_array3(xii,yii,zii)+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Real(vyarr_in(rx,ry,rz))
             
        end if              ! for boundary collision test

            end do
          end do
        end do
       
        !$OMP BARRIER
        
       end do !subz
       end do !suby
       end do !subx
      end do !i

       write(*,*) "finish constructing ZA density array"
        time2 =  time()
        write(*,*) "time in sec:", time2-time1

! save za particle cic in temp_array4
!      allocate(temp_array4(grid_num,grid_num,grid_num))  !left: vxarr_in, vyarr_in, vzarr_in, temp_array4
      allocate(temp_array(grid_num,grid_num,grid_num))  !left: vxarr_in, vyarr_in, vzarr_in, temp_array4
C$omp parallel do private(rx,ry,rz)
      do rx = 1, grid_num
        do ry = 1, grid_num
          do rz = 1, grid_num
!         temp_array4(rx,ry,rz) = aimag(vyarr_in(rx,ry,rz)) !NGP
         temp_array(rx,ry,rz) = aimag(vxarr_in(rx,ry,rz))/cell_size**3 !CIC
          end do
        end do
      end do


C$omp parallel do private(rx,ry,rz)
      do rx = 1, grid_num
        do ry = 1, grid_num
           do rz = 1, grid_num
              if(temp_array(rx,ry,rz) .gt. 0.0001)
     &        vyarr_in(rx,ry,rz) = cmplx(real(vyarr_in(rx,ry,rz)),
     &     aimag(vyarr_in(rx,ry,rz))/temp_array(rx,ry,rz)/cell_size**3) !CIC
          end do
        end do
      end do


C$omp parallel do private(rx,ry,rz)
      do rx = 1, grid_num
        do ry = 1, grid_num
          do rz = 1, grid_num
              if(temp_array(rx,ry,rz) .gt. 0.0001)             
     &        temp_array2(rx,ry,rz) =
     &        temp_array2(rx,ry,rz)/temp_array(rx,ry,rz)/cell_size**3 !CIC
          end do
        end do
      end do


C$omp parallel do private(rx,ry,rz)
      do rx = 1, grid_num
        do ry = 1, grid_num
          do rz = 1, grid_num
              if(temp_array(rx,ry,rz) .gt. 0.0001)             
     &        temp_array3(rx,ry,rz) =
     &        temp_array3(rx,ry,rz)/temp_array(rx,ry,rz)/cell_size**3 !CIC
          end do
        end do
      end do

      
!C$omp parallel do private(rx,ry,rz)
!      do rx = 1,grid_pdf
!         do ry = 1, grid_pdf
!            do rz = 1, grid_pdf
!              vxarr_in(rx,ry,rz) = cmplx(real(vxarr_in(rx,ry,rz)),0)
!!              vyarr_in(rx,ry,rz) = cmplx(real(vyarr_in(rx,ry,rz)),0) !image used by NGP density
!!              vzarr_in(rx,ry,rz) = cmplx(real(vzarr_in(rx,ry,rz)),0)
!            end do
!         end do
!      end do
      
      write(*,*) "going to introduce scatter in ZA density field"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

!      density_total = 0
!      density_cut_total=0


C$omp parallel private(i,j,k,temp,ak,bk,
C$omp+ thread_id) reduction(+:density_total,density_cut_total)
      thread_id = omp_get_thread_num()+1

C$omp do
      do i=1, grid_pdf
         do j=1, grid_pdf
            do k=1, grid_pdf
!               temp_array(i,j,k) = temp_array3(i,j,k)
!               density_total = density_total + temp_array(i,j,k)
! cut density
                 if (temp_array(i,j,k) .le. density_cut) then
!               if (Aimag(vxarr_in(i,j,k) .le. density_cut) then
!!!                 density_cut_total = density_cut_total+temp_array(i,j,k)
                   temp_array(i,j,k) = 0
!                   vxarr_in(i,j,k) = cmplx(real(vxarr_in(i,j,k)),0)
                 ak = 0
              else
                 call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)
               if (temp_array(i,j,k) .ge. density_sat) then
!              if (Aimag(vxarr_in(i,j,k)) .ge. density_sat) then
                temp_array(i,j,k) = density_sat
!              vxarr_in(i,j,k) = cmplx(Real(vxarr_in(i,j,k)),density_sat)
              end if

              end if
! scatter model 1: m1
              ak = ak*scatter !*temp_array(i,j,k)**scatter2 ! switch off scatter2
! scatter model 2: m2
!             ak = ak*scatter + bk*scatter2*temp**2
              
              if(ak .lt. 0) ak = exp(ak)-1
!                temp_array(i,j,k)=temp_array(i,j,k)*(1.+ak)
! add scatter2 to boost low density for modifying bik
      temp_array(i,j,k)=10.*(1-exp(-temp_array(i,j,k)/scatter2))*(1.+ak)
!              vxarr_in(i,j,k)=
!     &  cmplx(real(vxarr_in(i,j,k)),Aimag(vxarr_in(i,j,k))*(1.+ak))
            end do
         end do
      end do
!$omp end do
!$omp end parallel

!      write(*,*) "total density=", density_total
!      write(*,*) "total cut-out density=", density_cut_total
      write(*,*) "going to normalize pdf from scattered za field"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

      temp =  maxval(temp_array) 
      norm1 = 1.*(max_num-1)/temp 
      count = 0
!test
      write(*,*) "test 0104", temp, norm1
C$omp parallel private(i,j,k,temp_int,
C$omp+ thread_id) reduction(+:pdf)
      thread_id = omp_get_thread_num()+1
      do i=0, max_num
         pdf(i) = 0
      end do
!test
      write(*,*) "test 0105", thread_id
C$omp do
      do i=1, grid_pdf
         do j=1, grid_pdf
            do k=1, grid_pdf

              temp_array(i,j,k) = temp_array(i,j,k)*norm1
              temp_int = int(temp_array(i,j,k))
              count = count + temp_int
              if(temp_int .gt. max_num) then
               write(*,*)  temp_int,
     &     "is larger than max_num, NOT OK.."
!               goto 500
               temp_int = max_num
              end if
             pdf(temp_int) =    
     &         pdf(temp_int) +1
            end do
         end do
      end do
!$omp end do
!$omp end parallel
      write(*,*) "done normalizing scattered-density form ZA" 
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

!!
!
!!      open(10, file=Trim(output_file))
!      do i = 0,max_num
!!!         write(10,*) i, pdf(i)
!         if(pdf(i) .gt. 0) write(*,*) i, pdf(i)
!      end do
!!!      close(10)

      write(*,*) "read pdf of simulation"
      open(10,file=Trim(sim_pdf_file)//'.new',status='old')
      sim_pdf_max = 0
      do i=0, max_num
        read(10,*,end=201) temp, sim_pdf(i)
!        write(*,*) temp, sim_pdf(i)
        sim_pdf_max =sim_pdf_max+1
      end do
 201  close(10)
!      write(*,*) "sim_pdf_max=",sim_pdf_max-1


      sim_pdf_max = sim_pdf_max - 1 !the first index is 0
      za_pdf_index = max_num
      za_sum = 0       
      do i = sim_pdf_max,1, -1
        if(sim_pdf(i).gt. 0) then
!           write(*,*) "sim pdf(",i,")=", sim_pdf(i)
           do j=za_pdf_index,1,-1
             za_sum = za_sum + pdf(j)
!              write(*,*) "za:",j, za_sum, pdf(j)
             if (za_sum .ge. sim_pdf(i)) then
                za_temp = za_sum - sim_pdf(i)
                pdf_convert(j) = i
                pdf_convert_ratio(j) = 1.0*za_temp/pdf(j)
!            write(*,*) "partial:",j, za_temp, "from", pdf(j),
!     &     "will assinged to", i-1
!            write(*,*) "next sim pdf", sim_pdf(i-1)
!            write(*,*) "pdf_convert_ratio", j, pdf_convert_ratio(j)
                za_sum = za_temp
                goto 301
             else
                pdf_convert(j) = i
             end if
           end do 
 301       za_pdf_index = j - 1
        end if
      end do
      
!      do i=0, max_num
!         write(*,*) pdf_convert(i), pdf_convert_ratio(i)
!      end do
      
      write(*,*) "start producing mock using particles"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1
      open(17,file=Trim(output_file))
      do i=1, total_num_thread
         data_index(i) = 1          
      end do

      write(*,*) "zdist_rate", zdist_rate, zdist_fog

      allocate(dataxyz_omp(6,max_data/total_num_thread,
     &     total_num_thread))
      allocate(dataxyz(6,max_data))
      allocate(temp_int_array(grid_pdf,grid_pdf,grid_pdf))

C$omp parallel do private(rx,ry,rz)
      do rx = 1, grid_num
        do ry = 1, grid_num
          do rz = 1, grid_num
            temp_int_array(rx,ry,rz)=0
          end do
        end do
      end do



      num_to_make = 0
C$omp parallel private(rx,ry,rz,np,np_ratio,temp,temp2,temp3,temp4,
C$omp+ thread_id,ak,bk,xx,yy,zz,x,y,z,ix,iy,iz,xi,yi,zi,xii,yii,zii)
C$omp+ reduction(+:num_to_make)       
      thread_id = omp_get_thread_num()+1
C$omp do
      do rx = 1, grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
               np = pdf_convert(int(temp_array(rx,ry,rz)))
               if(np .gt. 0) then
                 np_ratio =
     &             pdf_convert_ratio(int(temp_array(rx,ry,rz)))
                 if(ran3_omp(iseedxx(thread_id),thread_id).lt.np_ratio)
     &              np = np - 1
               end if
               temp_int_array(rx,ry,rz) = np
               num_to_make = num_to_make + np
!                 vxarr_in(ix,iy,iz)=cmplx(real(vxarr_in(ix,iy,iz)),np)
            end do
         end do
      end do
!$omp end do
!      data_index(thread_id) = data_index(thread_id) -1
!$omp end parallel

      write(*,*) "sum of np=", num_to_make
      

!     
!!assign ZA particles to grid using NGP
C$omp parallel do private(rx,ry,rz,xx,yy,zz,ix,iy,iz,
C$omp+ x,y,z,xi,yi,zi,xii,yii,zii,i,rxx,ryy,rzz,rxi,ryi,rzi,
C$omp+ np,np_ratio,thread_id,ak,bk,temp,temp2,temp3,temp4)
      !!reduction(+:pdf2) 

      
      do i=0, total_num_thread-1
         thread_id = omp_get_thread_num()+1
         rxx = i/(y_div_omp*z_div_omp)
         ryy = mod(i,y_div_omp*z_div_omp)/z_div_omp
         rzz = mod(i,z_div_omp)
!      write(*,*) "thread id", thread_id
        do rxi = 1, grid_num/x_div_omp
          do ryi = 1, grid_num/y_div_omp
            do rzi = 1, grid_num/z_div_omp

              rx = rxx*grid_num/x_div_omp+rxi
              ry = ryy*grid_num/y_div_omp+ryi
              rz = rzz*grid_num/z_div_omp+rzi

               xx = real(vxarr_in(rx,ry,rz)) + gridsize*rx
               if (xx .ge. 0) then
                   xx = xx - int(xx/boxsize)*boxsize
               else
                   xx = xx - (int(xx/boxsize)-1)*boxsize
               end if
               yy = real(vyarr_in(rx,ry,rz)) + gridsize*ry
               if (yy .ge. 0) then
                   yy = yy - int(yy/boxsize)*boxsize
               else
                   yy = yy - (int(yy/boxsize)-1)*boxsize
               end if
               zz = real(vzarr_in(rx,ry,rz)) + gridsize*rz
               if (zz .ge. 0) then
                   zz = zz - int(zz/boxsize)*boxsize
               else
                   zz = zz - (int(zz/boxsize)-1)*boxsize
               end if

!test generate raw za particles
!               if(ran3(iseed) .lt. density) then
!                  write(17,*) xx,yy,zz
!               end if

!     find closest grid point
               ix = int(xx/cell_size+0.5)+1
               if (ix .gt. grid_pdf) ix=1
               iy = int(yy/cell_size+0.5)+1
               if (iy .gt. grid_pdf) iy=1
               iz = int(zz/cell_size+0.5)+1
               if (iz .gt. grid_pdf) iz=1

!               np = pdf_convert(int(temp_array(ix,iy,iz)))
!               if(np .gt. 0) then
!                 np_ratio =
!     &             pdf_convert_ratio(int(temp_array(ix,iy,iz)))
!                 if(ran3_omp(iseedxx(thread_id),thread_id).lt.np_ratio)
!     &              np = np - 1

!                 vxarr_in(ix,iy,iz)=cmplx(real(vxarr_in(ix,iy,iz)),np)
                 
                  if(temp_int_array(ix,iy,iz) .gt. 0) then

                  dataxyz_omp(1,data_index(thread_id),thread_id)=xx
                  dataxyz_omp(2,data_index(thread_id),thread_id)=yy
                  dataxyz_omp(3,data_index(thread_id),thread_id)=zz
                  
!! add velocity with NGP
!                  
!                  temp = Aimag(vyarr_in(ix,iy,iz))
!                  temp2= temp_array2(ix,iy,iz)
!                  temp3= temp_array3(ix,iy,iz)

! add velocity with CIC
                  
               ix = int(xx/cell_size)+1
               if (ix .gt. grid_pdf) ix=grid_pdf
               iy = int(yy/cell_size)+1
               if (iy .gt. grid_pdf) iy=grid_pdf
               iz = int(zz/cell_size)+1
               if (iz .gt. grid_pdf) iz=grid_pdf

            xi = ix
            yi = iy
            zi = iz
             x = xx
             y = yy
             z = zz

               xii = xi + 1
               if(xii.eq.grid_pdf+1) xii=1
               yii = yi + 1
               if(yii.eq.grid_pdf+1) yii=1
               zii = zi + 1
               if(zii.eq.grid_pdf+1) zii=1

!compute density at the location of each halo using cic
               temp4 = Aimag(vxarr_in(xi,yi,zi))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))+
     &          Aimag(vxarr_in(xii,yi,zi))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))+
     &          Aimag(vxarr_in(xi,yii,zi))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))+
     &          Aimag(vxarr_in(xi,yi,zii))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))+
     &          Aimag(vxarr_in(xii,yii,zi))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))+
     &          Aimag(vxarr_in(xi,yii,zii))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))+
     &          Aimag(vxarr_in(xii,yi,zii))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))+
     &          Aimag(vxarr_in(xii,yii,zii))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))


!compute vz, vx, vy using cic
               temp = Aimag(vyarr_in(xi,yi,zi))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yi,zi))+
     &          Aimag(vyarr_in(xii,yi,zi))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yi,zi))+
     &          Aimag(vyarr_in(xi,yii,zi))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yii,zi))+
     &          Aimag(vyarr_in(xi,yi,zii))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yi,zii))+
     &          Aimag(vyarr_in(xii,yii,zi))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yii,zi))+
     &          Aimag(vyarr_in(xi,yii,zii))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yii,zii))+
     &          Aimag(vyarr_in(xii,yi,zii))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yi,zii))+
     &          Aimag(vyarr_in(xii,yii,zii))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yii,zii))

               temp = temp / temp4
 

              temp2 = (temp_array2(xi,yi,zi))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yi,zi))+
     &          (temp_array2(xii,yi,zi))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yi,zi))+
     &          (temp_array2(xi,yii,zi))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yii,zi))+
     &          (temp_array2(xi,yi,zii))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yi,zii))+
     &          (temp_array2(xii,yii,zi))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yii,zi))+
     &          (temp_array2(xi,yii,zii))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yii,zii))+
     &          (temp_array2(xii,yi,zii))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yi,zii))+
     &          (temp_array2(xii,yii,zii))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yii,zii))
               temp2 = temp2 / temp4
 
              temp3 = (temp_array3(xi,yi,zi))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yi,zi))+
     &          (temp_array3(xii,yi,zi))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yi,zi))+
     &          (temp_array3(xi,yii,zi))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yii,zi))+
     &          (temp_array3(xi,yi,zii))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yi,zii))+
     &          (temp_array3(xii,yii,zi))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yii,zi))+
     &          (temp_array3(xi,yii,zii))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yii,zii))+
     &          (temp_array3(xii,yi,zii))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yi,zii))+
     &          (temp_array3(xii,yii,zii))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yii,zii))
               temp3 = temp3 / temp4



                  
                  call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)
                  dataxyz_omp(4,data_index(thread_id),thread_id)=
     &       temp2*zdist_rate + ak*zdist_fog                  
!     &      temp_array2(rx,ry,rz)*zdist_rate + ak*zdist_fog
!     &             real(vxarr_in(rx,ry,rz))*zdist_rate + ak*zdist_fog
                  dataxyz_omp(5,data_index(thread_id),thread_id)=
     &       temp3*zdist_rate + bk*zdist_fog
!     &      temp_array3(rx,ry,rz)*zdist_rate + bk*zdist_fog
!     &             real(vyarr_in(rx,ry,rz))*zdist_rate + bk*zdist_fog                  
                  call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)
                  dataxyz_omp(6,data_index(thread_id),thread_id)=
     &       temp*zdist_rate + ak*zdist_fog
!     &  aimag(vyarr_in(rx,ry,rz))*zdist_rate + ak*zdist_fog
!     &             real(vzarr_in(rx,ry,rz))*zdist_rate + ak*zdist_fog
                  
                  data_index(thread_id) = data_index(thread_id) +1
                  if(data_index(thread_id).eq.max_data/total_num_thread)
     &             then
             write(*,*) "too many particles saved. need to modify code!"
                     stop
                  end if


                  end if
               

         end do
         end do
         end do
      data_index(thread_id) = data_index(thread_id) -1

      end do !i



      mock_data_num = 0
      do i=1, total_num_thread
        mock_data_num = mock_data_num + data_index(i)
      end do 
      write(*,*)
     & "total number of particles kept for assigning galaxies",
     &  mock_data_num

!copy all selected particles to one array (x,y,z,vx,vy,vz) and create a logical array for flag      
      count=0
      do i=1, total_num_thread
         do j=1, data_index(i)
            count=count+1
            dataxyz(1:6,count)=dataxyz_omp(1:6,j,i)
         end do
      end do
!      write(*,*) count, "should equal", mock_data_num
!      time2 =  time()
!      write(*,*) "time in sec:", time2-time1

!     shuffle it      
      do i= mock_data_num, 2, -1

            k=int(ran3(iseed)*i)+1
            dataxyz_temp(1:6) = dataxyz(1:6,i)
            dataxyz(1:6,i) = dataxyz(1:6,k)
            dataxyz(1:6,k) = dataxyz_temp(1:6)

      end do
      write(*,*) "done shuffling array"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

      temp_int = 0
      do rx = 1, grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
               temp_int = temp_int +temp_int_array(rx,ry,rz)
            end do
         end do
      end do
      write(*,*) "temp_int_array sum1= ", temp_int
      
! assign galaxies to particles -> dataxyz_keep(i) = .true.      
       num_data = 0
      do i=1, mock_data_num
         
         xx = dataxyz(1,i)
         yy = dataxyz(2,i)
         zz = dataxyz(3,i)
         ix = int(xx/cell_size+0.5)+1
         if (ix .gt. grid_pdf) ix=1
         iy = int(yy/cell_size+0.5)+1
         if (iy .gt. grid_pdf) iy=1
         iz = int(zz/cell_size+0.5)+1
         if (iz .gt. grid_pdf) iz=1
         
!     if(aimag(vxarr_in(ix,iy,iz)) .gt. 0.1) then
         if(temp_int_array(ix,iy,iz).gt.0) then
         dataxyz_keep(i) = .true.
!            vxarr_in(ix,iy,iz)=
!     & cmplx(real(vxarr_in(ix,iy,iz)),aimag(vxarr_in(ix,iy,iz))-1)
         temp_int_array(ix,iy,iz) = temp_int_array(ix,iy,iz) - 1
            num_data = num_data + 1
         else
            dataxyz_keep(i) = .false.
         end if
      end do
      write(*,*) "total mock produced", num_data


      temp_int = 0
      do rx = 1, grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
               temp_int = temp_int +temp_int_array(rx,ry,rz)
            end do
         end do
      end do
      write(*,*) "temp_int_array sum2= ", temp_int

      

!              if(AIMAG(vyarr_in(i,j,k)) .gt. max_num) then
!               write(*,*)  "increase max_num! ijij" 
!               stop
!              end if
!
!             pdf2(int(AIMAG(vyarr_in(i,j,k)))) =
!     &         pdf2(int(AIMAG(vyarr_in(i,j,k)))) +1
!
!            end do
!         end do
!      end do
!
!      temp = 0
!!!      open(10, file=Trim(output_file)//'.pdf.cic')
!      do i = 0, max_num
!!!         write(10,*) i, pdf2(i)
!         if(pdf2(i) .gt. 0) write(*,*)"converted pdf:",i, pdf2(i)
!         temp = temp + pdf2(i)*i
!      end do
!!!      close(10)
!      write(*,*) "total number from converted pdf", temp

      
!     write mock data
      do i = 1, mock_data_num
!     do j = 1, data_index(i)
         if(dataxyz_keep(i) .eq. .true.) then
           write(17,'(6F10.3)')
     &     real(dataxyz(1,i)),real(dataxyz(2,i)),
     &         real(dataxyz(3,i)), real(dataxyz(4,i)),
     &           real(dataxyz(5,i)),real(dataxyz(6,i))
         end if
!          write(17,*) dataxyz_omp(1:6,j,i)
       ! end do
      end do
      close(17)
      
!      do i = 1, total_num_thread
!        do j = 1, data_index(i)
!         write(17,'(6F10.3)')
!     &  real(dataxyz_omp(1,j,i)),real(dataxyz_omp(2,j,i)),
!     &         real(dataxyz_omp(3,j,i)), real(dataxyz_omp(4,j,i)),
!     &              real(dataxyz_omp(5,j,i)),real(dataxyz_omp(6,j,i))
!!          write(17,*) dataxyz_omp(1:6,j,i)
!        end do
!      end do
!      close(17)
      write(*,*) "finish saving mock from particles,
     &   now working on cic for the rest"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

!!!------------------------------------------ assign the rest of the halos with CIC -----
      write(*,*) "time in sec:", time2-time1
      open(17,file=Trim(output_file),access='append')
      do i=1, total_num_thread
         data_index(i) = 1          
      end do


C$omp parallel private(rx,ry,rz,np,np_ratio,temp,temp2,temp3,temp4,
C$omp+ thread_id,ak,bk,xx,yy,zz,x,y,z,ix,iy,iz,xi,yi,zi,xii,yii,zii) 
      thread_id = omp_get_thread_num()+1
C$omp do
      do rx = 1, grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
!              np = pdf_convert(int(AIMAG(vzarr_in2(rx,ry,rz))))
!              np_ratio = 
!     &             pdf_convert_ratio(int(AIMAG(vzarr_in2(rx,ry,rz))))
!              if(ran3_omp(iseedxx(thread_id),thread_id) .lt. np_ratio)
!     &       np = np - 1
               np = temp_int_array(rx,ry,rz)
!            vyarr_in2(rx,ry,rz) = cmplx(real(vyarr_in2(rx,ry,rz)),1.*np) !??

               do i =1, np
!                  temp = ran3_omp(iseedxx(thread_id),thread_id) - 0.5
                  temp2 = ran3_omp(iseedxx(thread_id),thread_id)
                  temp = 1.*int(100.*temp2)/100-0.5

                  if(temp.ge.0) then
                    dataxyz_omp(1,data_index(thread_id),thread_id)=
     &                  (1.-sqrt(2.*temp)+rx-1)*cell_size
                  else
                    dataxyz_omp(1,data_index(thread_id),thread_id)=
     &                  (-1.+sqrt(-2.*temp)+rx-1)*cell_size
                  end if
                  if(dataxyz_omp(1,data_index(thread_id),thread_id).gt.
     &                boxsize) then
                    dataxyz_omp(1,data_index(thread_id),thread_id)=
     &              dataxyz_omp(1,data_index(thread_id),thread_id)-
     &              boxsize
                  end if
                  if(dataxyz_omp(1,data_index(thread_id),thread_id).lt.
     &                0) then
                    dataxyz_omp(1,data_index(thread_id),thread_id)=
     &              dataxyz_omp(1,data_index(thread_id),thread_id)+
     &              boxsize
                  end if

!                  temp = ran3_omp(iseedxx(thread_id),thread_id) - 0.5
                  temp = 1.*mod(int(temp2*10000),100)/100-0.5
                  if(temp.ge.0) then
                    dataxyz_omp(2,data_index(thread_id),thread_id)=
     &              (1.-sqrt(2.*temp)+ry-1)*cell_size
                  else
                    dataxyz_omp(2,data_index(thread_id),thread_id)=
     &              (-1.+sqrt(-2.*temp)+ry-1)*cell_size
                  end if
                  if(dataxyz_omp(2,data_index(thread_id),thread_id).gt.
     &              boxsize) then
                    dataxyz_omp(2,data_index(thread_id),thread_id)=
     &              dataxyz_omp(2,data_index(thread_id),thread_id)-
     &              boxsize
                  end if
                  if(dataxyz_omp(2,data_index(thread_id),thread_id).lt.
     &              0) then
                    dataxyz_omp(2,data_index(thread_id),thread_id)=
     &              dataxyz_omp(2,data_index(thread_id),thread_id)+
     &              boxsize
                  end if

!                  temp = ran3_omp(iseedxx(thread_id),thread_id) - 0.5
                  temp = 1.*mod(int(temp2*1000000),100)/100-0.5
                  if(temp.ge.0) then
                    dataxyz_omp(3,data_index(thread_id),thread_id)=
     &              (1.-sqrt(2.*temp)+rz-1)*cell_size
                  else
                    dataxyz_omp(3,data_index(thread_id),thread_id)=
     &              (-1.+sqrt(-2.*temp)+rz-1)*cell_size
                  end if
                  if(dataxyz_omp(3,data_index(thread_id),thread_id).gt.
     &               boxsize) then
                    dataxyz_omp(3,data_index(thread_id),thread_id) =
     &              dataxyz_omp(3,data_index(thread_id),thread_id)-
     &              boxsize
                  end if
                  if(dataxyz_omp(3,data_index(thread_id),thread_id).lt.
     &               0) then
                    dataxyz_omp(3,data_index(thread_id),thread_id)=
     &              dataxyz_omp(3,data_index(thread_id),thread_id)+
     &              boxsize
                  end if
                 xx = dataxyz_omp(1,data_index(thread_id),thread_id)
                 yy = dataxyz_omp(2,data_index(thread_id),thread_id)
                 zz = dataxyz_omp(3,data_index(thread_id),thread_id)
! add velocity
               ix = int(xx/cell_size)+1
               if (ix .gt. grid_pdf) ix=grid_pdf
               iy = int(yy/cell_size)+1
               if (iy .gt. grid_pdf) iy=grid_pdf
               iz = int(zz/cell_size)+1
               if (iz .gt. grid_pdf) iz=grid_pdf

            xi = ix
            yi = iy
            zi = iz
             x = xx
             y = yy
             z = zz

               xii = xi + 1
               if(xii.eq.grid_pdf+1) xii=1
               yii = yi + 1
               if(yii.eq.grid_pdf+1) yii=1
               zii = zi + 1
               if(zii.eq.grid_pdf+1) zii=1

!compute density at the location of each halo using cic
               temp4 = Aimag(vxarr_in(xi,yi,zi))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))+
     &          Aimag(vxarr_in(xii,yi,zi))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))+
     &          Aimag(vxarr_in(xi,yii,zi))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))+
     &          Aimag(vxarr_in(xi,yi,zii))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))+
     &          Aimag(vxarr_in(xii,yii,zi))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))+
     &          Aimag(vxarr_in(xi,yii,zii))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))+
     &          Aimag(vxarr_in(xii,yi,zii))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))+
     &          Aimag(vxarr_in(xii,yii,zii))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))


!compute vz, vx, vy using cic
               temp = Aimag(vyarr_in(xi,yi,zi))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yi,zi))+
     &          Aimag(vyarr_in(xii,yi,zi))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yi,zi))+
     &          Aimag(vyarr_in(xi,yii,zi))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yii,zi))+
     &          Aimag(vyarr_in(xi,yi,zii))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yi,zii))+
     &          Aimag(vyarr_in(xii,yii,zi))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yii,zi))+
     &          Aimag(vyarr_in(xi,yii,zii))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yii,zii))+
     &          Aimag(vyarr_in(xii,yi,zii))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yi,zii))+
     &          Aimag(vyarr_in(xii,yii,zii))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yii,zii))

               temp = temp / temp4
 

              temp2 = (temp_array2(xi,yi,zi))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yi,zi))+
     &          (temp_array2(xii,yi,zi))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yi,zi))+
     &          (temp_array2(xi,yii,zi))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yii,zi))+
     &          (temp_array2(xi,yi,zii))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yi,zii))+
     &          (temp_array2(xii,yii,zi))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yii,zi))+
     &          (temp_array2(xi,yii,zii))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yii,zii))+
     &          (temp_array2(xii,yi,zii))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yi,zii))+
     &          (temp_array2(xii,yii,zii))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yii,zii))
               temp2 = temp2 / temp4
 
              temp3 = (temp_array3(xi,yi,zi))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yi,zi))+
     &          (temp_array3(xii,yi,zi))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yi,zi))+
     &          (temp_array3(xi,yii,zi))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xi,yii,zi))+
     &          (temp_array3(xi,yi,zii))*
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yi,zii))+
     &          (temp_array3(xii,yii,zi))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi))*Aimag(vxarr_in(xii,yii,zi))+
     &          (temp_array3(xi,yii,zii))*
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xi,yii,zii))+
     &          (temp_array3(xii,yi,zii))*
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yi,zii))+
     &          (temp_array3(xii,yii,zii))*
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1))*Aimag(vxarr_in(xii,yii,zii))
               temp3 = temp3 / temp4



                  
                  call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)
                  dataxyz_omp(4,data_index(thread_id),thread_id)=
     &       temp2*zdist_rate + ak*zdist_fog                  
!     &      temp_array2(rx,ry,rz)*zdist_rate + ak*zdist_fog
!     &             real(vxarr_in(rx,ry,rz))*zdist_rate + ak*zdist_fog
                  dataxyz_omp(5,data_index(thread_id),thread_id)=
     &       temp3*zdist_rate + bk*zdist_fog
!     &      temp_array3(rx,ry,rz)*zdist_rate + bk*zdist_fog
!     &             real(vyarr_in(rx,ry,rz))*zdist_rate + bk*zdist_fog                  
                  call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)
                  dataxyz_omp(6,data_index(thread_id),thread_id)=
     &       temp*zdist_rate + ak*zdist_fog
!     &  aimag(vyarr_in(rx,ry,rz))*zdist_rate + ak*zdist_fog
!     &             real(vzarr_in(rx,ry,rz))*zdist_rate + ak*zdist_fog
                  
                  data_index(thread_id) = data_index(thread_id) +1
                  if(data_index(thread_id).eq.max_data/total_num_thread)
     &             then
             write(*,*) "too many cic halos saved. need to modify code!"
                     stop
                  end if
               end do !i


            end do
         end do
      end do
!$omp end do
      data_index(thread_id) = data_index(thread_id) -1
!$omp end parallel


      mock_data_num2 = 0
      do i=1, total_num_thread
        mock_data_num2 = mock_data_num2 + data_index(i)
      end do  
      write(*,*) "particle galaxies +  CIC galaxies=",
     & num_data + mock_data_num2


      do i = 1, total_num_thread
        do j = 1, data_index(i)
         write(17,'(6F10.3)')
     &  real(dataxyz_omp(1,j,i)),real(dataxyz_omp(2,j,i)),
     &         real(dataxyz_omp(3,j,i)), real(dataxyz_omp(4,j,i)),
     &              real(dataxyz_omp(5,j,i)),real(dataxyz_omp(6,j,i))
!          write(17,*) dataxyz_omp(1:6,j,i)
        end do
      end do
      close(17)
      write(*,*) "finish saving the rest of mock"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1
      
!!!!------------------- finish writing the rest of the halos -------      

!     compute pdf of mock cata of cic
!$omp parallel do private(rx,ry,rz)
      do rx = 1,grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
!              count_in_cell(rx,ry,rz) = 0
            vzarr_in(rx,ry,rz)=cmplx(real(vzarr_in(rx,ry,rz)),0)
            vxarr_in(rx,ry,rz)=cmplx(real(vxarr_in(rx,ry,rz)),0)
 !           vwarr_in(rx,ry,rz)=cmplx(0,0)
            end do
         end do
      end do

      do i=0, max_num
         pdf2(i) = 0
      end do

      do i=-max_vx_pdf, max_vx_pdf
         za_vx_pdf(i)=0
      end do 

!     might be openmp while generating the data
!BUT, it is for post-process so that the time used is not so important
!      do i = 1, total_num_thread
!     do j = 1, data_index(i)
      do i = 1, mock_data_num
         if(dataxyz_keep(i) .eq. .true.) then
           x=dataxyz(1,i)
           y=dataxyz(2,i)
           z=dataxyz(3,i)
           vz=dataxyz(6,i)

           zz=z+vz*inverse_aH !(1+z)/Hubble
           if(zz.gt.boxsize) zz = zz - boxsize
           if(zz.lt.0.) zz = zz + boxsize
           if(zz.gt.boxsize) zz = zz - boxsize
           if(zz.lt.0.) zz = zz + boxsize
           if(zz.gt.boxsize) zz = zz - boxsize
           if(zz.lt.0.) zz = zz + boxsize
           if(zz.gt.boxsize) zz = zz - boxsize
           if(zz.lt.0.) zz = zz + boxsize
           if(zz.gt.boxsize) zz = zz - boxsize
           if(zz.lt.0.) zz = zz + boxsize
           if(zz.gt.boxsize) then
                  write(*,*) "vz is too large:", vz 
                  stop
           end if
           if(zz.lt.0.) then 
                  write(*,*) "vz is too small:", vz
                  stop
           end if
! these lines seem resulting segmentation fault, so take them out            
!           if(vz.ge.0) then
!             za_vx_pdf(int(vz)) = za_vx_pdf(int(vz))+1
!           else
!             za_vx_pdf(int(vz)-1) = za_vx_pdf(int(vz)-1)+1
!           end if


           if (CIC_assignment .eq. .false.) then !NGP

             xi = int(x/cell_size+0.5)+1
             yi = int(y/cell_size+0.5)+1
             zi = int(z/cell_size+0.5)+1
             zzi = int(zz/cell_size+0.5)+1

             if(xi .eq. grid_pdf+1) xi = 1
             if(yi .eq. grid_pdf+1) yi = 1
             if(zi .eq. grid_pdf+1) zi = 1
             if(zzi .eq. grid_pdf+1) zzi = 1

               vzarr_in(xi,yi,zi) =cmplx(real(vzarr_in(xi,yi,zi)),
     &    AIMAG(vzarr_in(xi,yi,zi))+1)
               vxarr_in(xi,yi,zzi) =cmplx(real(vxarr_in(xi,yi,zzi)),
     &    AIMAG(vxarr_in(xi,yi,zzi))+1) 
             
           else   

             xi = int(x/cell_size)+1
             yi = int(y/cell_size)+1
             zi = int(z/cell_size)+1
             zzi = int(zz/cell_size)+1

              
              
               if(xi .eq. grid_pdf+1) xi = xi-1
               if(yi .eq. grid_pdf+1) yi = yi-1
               if(zi .eq. grid_pdf+1) zi = zi-1
               if(zzi .eq. grid_pdf+1) zzi = zzi-1

               xii = xi + 1
               if(xii.eq.grid_pdf+1) xii=1
               yii = yi + 1
               if(yii.eq.grid_pdf+1) yii=1
               zii = zi + 1
               if(zii.eq.grid_pdf+1) zii=1
               zzii = zzi + 1
               if(zzii.eq.grid_pdf+1) zzii=1


               vzarr_in(xi,yi,zi) =cmplx(real(vzarr_in(xi,yi,zi)),
     &    AIMAG(vzarr_in(xi,yi,zi))+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &              (-z+cell_size*(zi)))
               vzarr_in(xii,yi,zi) =cmplx(real(vzarr_in(xii,yi,zi)),
     &   AIMAG(vzarr_in(xii,yi,zi))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi)))
               vzarr_in(xi,yii,zi) = cmplx(real(vzarr_in(xi,yii,zi)),
     &    AIMAG(vzarr_in(xi,yii,zi))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi)))
               vzarr_in(xi,yi,zii) = cmplx(real(vzarr_in(xi,yi,zii)),
     &    AIMAG(vzarr_in(xi,yi,zii))+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1)))
             vzarr_in(xii,yii,zi) = cmplx(real(vzarr_in(xii,yii,zi)),
     &    AIMAG(vzarr_in(xii,yii,zi))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi)))
             vzarr_in(xi,yii,zii) = cmplx(real(vzarr_in(xi,yii,zii)),
     &    AIMAG(vzarr_in(xi,yii,zii))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1)))
             vzarr_in(xii,yi,zii) = cmplx(real(vzarr_in(xii,yi,zii)),
     &    AIMAG(vzarr_in(xii,yi,zii))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1)))
             vzarr_in(xii,yii,zii) =cmplx(real(vzarr_in(xii,yii,zii)),
     &    AIMAG(vzarr_in(xii,yii,zii))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1)))
             

             vxarr_in(xi,yi,zzi) =cmplx(real(vxarr_in(xi,yi,zzi)),
     &    AIMAG(vxarr_in(xi,yi,zzi))+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-zz+cell_size*(zzi)))
             vxarr_in(xii,yi,zzi) =cmplx(real(vxarr_in(xii,yi,zzi)),
     &   AIMAG(vxarr_in(xii,yi,zzi))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-zz+cell_size*(zzi)))
             vxarr_in(xi,yii,zzi) = cmplx(real(vxarr_in(xi,yii,zzi)),
     &    AIMAG(vxarr_in(xi,yii,zzi))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-zz+cell_size*(zzi)))
             vxarr_in(xi,yi,zzii) = cmplx(real(vxarr_in(xi,yi,zzii)),
     &    AIMAG(vxarr_in(xi,yi,zzii))+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (zz-cell_size*(zzi-1)))
           vxarr_in(xii,yii,zzi) = cmplx(real(vxarr_in(xii,yii,zzi)),
     &    AIMAG(vxarr_in(xii,yii,zzi))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-zz+cell_size*(zzi)))
           vxarr_in(xi,yii,zzii) = cmplx(real(vxarr_in(xi,yii,zzii)),
     &    AIMAG(vxarr_in(xi,yii,zzii))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (zz-cell_size*(zzi-1)))
           vxarr_in(xii,yi,zzii) = cmplx(real(vxarr_in(xii,yi,zzii)),
     &    AIMAG(vxarr_in(xii,yi,zzii))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (zz-cell_size*(zzi-1)))
           vxarr_in(xii,yii,zzii) =cmplx(real(vxarr_in(xii,yii,zzii)),
     &    AIMAG(vxarr_in(xii,yii,zzii))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (zz-cell_size*(zzi-1)))

      end if                    !cic or not

      end if ! keep or not
!        end do
      end do

!---------add the rest of mocks to pdf gird for computing pk ---------- Albert

      do i = 1, total_num_thread
       do j = 1, data_index(i)
!      do i = 1, mock_data_num
!     if(dataxyz_keep(i) .eq. .true.) then
           x=dataxyz_omp(1,j,i)
           y=dataxyz_omp(2,j,i)
           z=dataxyz_omp(3,j,i)
           vz=dataxyz_omp(6,j,i)

           zz=z+vz*inverse_aH !(1+z)/Hubble
           if(zz.gt.boxsize) zz = zz - boxsize
           if(zz.lt.0.) zz = zz + boxsize;
           if(zz.gt.boxsize) zz = zz - boxsize
           if(zz.lt.0.) zz = zz + boxsize;
           if(zz.gt.boxsize) zz = zz - boxsize
           if(zz.lt.0.) zz = zz + boxsize;
           if(zz.gt.boxsize) zz = zz - boxsize
           if(zz.lt.0.) zz = zz + boxsize;
           if(zz.gt.boxsize) zz = zz - boxsize
           if(zz.lt.0.) zz = zz + boxsize;
           if(zz.gt.boxsize) then
                  write(*,*) "vz is too large:", vz 
                  stop
           end if
           if(zz.lt.0.) then 
                  write(*,*) "vz is too small:", vz
                  stop
           end if

! these lines seem resulting segmentation fault, so take them out            
!           if(vz.ge.0) then
!             za_vx_pdf(int(vz)) = za_vx_pdf(int(vz))+1
!           else
!             za_vx_pdf(int(vz)-1) = za_vx_pdf(int(vz)-1)+1
!           end if


           if (CIC_assignment .eq. .false.) then !NGP

             xi = int(x/cell_size+0.5)+1
             yi = int(y/cell_size+0.5)+1
             zi = int(z/cell_size+0.5)+1
             zzi = int(zz/cell_size+0.5)+1

             if(xi .eq. grid_pdf+1) xi = 1
             if(yi .eq. grid_pdf+1) yi = 1
             if(zi .eq. grid_pdf+1) zi = 1
             if(zzi .eq. grid_pdf+1) zzi = 1

               vzarr_in(xi,yi,zi) =cmplx(real(vzarr_in(xi,yi,zi)),
     &    AIMAG(vzarr_in(xi,yi,zi))+1)
               vxarr_in(xi,yi,zzi) =cmplx(real(vxarr_in(xi,yi,zzi)),
     &    AIMAG(vxarr_in(xi,yi,zzi))+1) 
             
           else   

             xi = int(x/cell_size)+1
             yi = int(y/cell_size)+1
             zi = int(z/cell_size)+1
             zzi = int(zz/cell_size)+1

              
              
               if(xi .eq. grid_pdf+1) xi = xi-1
               if(yi .eq. grid_pdf+1) yi = yi-1
               if(zi .eq. grid_pdf+1) zi = zi-1
               if(zzi .eq. grid_pdf+1) zzi = zzi-1

               xii = xi + 1
               if(xii.eq.grid_pdf+1) xii=1
               yii = yi + 1
               if(yii.eq.grid_pdf+1) yii=1
               zii = zi + 1
               if(zii.eq.grid_pdf+1) zii=1
               zzii = zzi + 1
               if(zzii.eq.grid_pdf+1) zzii=1


               vzarr_in(xi,yi,zi) =cmplx(real(vzarr_in(xi,yi,zi)),
     &    AIMAG(vzarr_in(xi,yi,zi))+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &              (-z+cell_size*(zi)))
               vzarr_in(xii,yi,zi) =cmplx(real(vzarr_in(xii,yi,zi)),
     &   AIMAG(vzarr_in(xii,yi,zi))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi)))
               vzarr_in(xi,yii,zi) = cmplx(real(vzarr_in(xi,yii,zi)),
     &    AIMAG(vzarr_in(xi,yii,zi))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi)))
               vzarr_in(xi,yi,zii) = cmplx(real(vzarr_in(xi,yi,zii)),
     &    AIMAG(vzarr_in(xi,yi,zii))+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1)))
             vzarr_in(xii,yii,zi) = cmplx(real(vzarr_in(xii,yii,zi)),
     &    AIMAG(vzarr_in(xii,yii,zi))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi)))
             vzarr_in(xi,yii,zii) = cmplx(real(vzarr_in(xi,yii,zii)),
     &    AIMAG(vzarr_in(xi,yii,zii))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1)))
             vzarr_in(xii,yi,zii) = cmplx(real(vzarr_in(xii,yi,zii)),
     &    AIMAG(vzarr_in(xii,yi,zii))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1)))
             vzarr_in(xii,yii,zii) =cmplx(real(vzarr_in(xii,yii,zii)),
     &    AIMAG(vzarr_in(xii,yii,zii))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1)))
             

             vxarr_in(xi,yi,zzi) =cmplx(real(vxarr_in(xi,yi,zzi)),
     &    AIMAG(vxarr_in(xi,yi,zzi))+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-zz+cell_size*(zzi)))
             vxarr_in(xii,yi,zzi) =cmplx(real(vxarr_in(xii,yi,zzi)),
     &   AIMAG(vxarr_in(xii,yi,zzi))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-zz+cell_size*(zzi)))
             vxarr_in(xi,yii,zzi) = cmplx(real(vxarr_in(xi,yii,zzi)),
     &    AIMAG(vxarr_in(xi,yii,zzi))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-zz+cell_size*(zzi)))
             vxarr_in(xi,yi,zzii) = cmplx(real(vxarr_in(xi,yi,zzii)),
     &    AIMAG(vxarr_in(xi,yi,zzii))+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (zz-cell_size*(zzi-1)))
           vxarr_in(xii,yii,zzi) = cmplx(real(vxarr_in(xii,yii,zzi)),
     &    AIMAG(vxarr_in(xii,yii,zzi))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-zz+cell_size*(zzi)))
           vxarr_in(xi,yii,zzii) = cmplx(real(vxarr_in(xi,yii,zzii)),
     &    AIMAG(vxarr_in(xi,yii,zzii))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (zz-cell_size*(zzi-1)))
           vxarr_in(xii,yi,zzii) = cmplx(real(vxarr_in(xii,yi,zzii)),
     &    AIMAG(vxarr_in(xii,yi,zzii))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (zz-cell_size*(zzi-1)))
           vxarr_in(xii,yii,zzii) =cmplx(real(vxarr_in(xii,yii,zzii)),
     &    AIMAG(vxarr_in(xii,yii,zzii))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (zz-cell_size*(zzi-1)))

      end if                    !cic or not

!      end if ! keep or not
        end do
      end do

      
      write(*,*) "finish assigning mock to pdf grid"

      time2 =  time()
      write(*,*) "time in sec:", time2-time1

      

!---  compute pk ---

      num_data = num_data + mock_data_num2

      write(*,*) "start compute pk"
      do i=1, k_index_max
        pk(i)=0   
        pk_quad(i)=0
        k_num_sep(i)=0
        quad_modes(i)=0
      end do
!$omp parallel do private(i,j,k)
      do i=1, grid_pdf
        do j=1, grid_pdf
           do k=1, grid_pdf
             if(CIC_assignment.eq..true.) then
            vzarr_in(i,j,k)=
     &  cmplx(AIMAG(vzarr_in(i,j,k))/cell_size**3-
     &             (1.*num_data)/(1.*grid_pdf)**3,0)
             else
            vzarr_in(i,j,k)=
     &  cmplx(AIMAG(vzarr_in(i,j,k))-
     &             (1.*num_data)/(1.*grid_pdf)**3,0)
             end if
          end do
        end do
      end do

      write(*,*) "starting FFT for pk"
!      FFTW3
      call sfftw_plan_dft_3d(planz2, grid_num, grid_num, grid_num,
     &  vzarr_in, vzarr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call sfftw_execute_dft(planz2,vzarr_in,vzarr_in)
      call sfftw_destroy_plan(planz2)

      time2 =  time()
      write(*,*) "finish FFT, time in sec:", time2-time1

      box_size = boxsize

!$omp parallel do private(i,j,k,kx,ky,kz,temp,k_index,mu,pk_temp)
!$omp+ reduction(+:k_num_sep,pk_quad,quad_modes,pk)
      do i=1, grid_num
        do j=1, grid_num
          do k=1, grid_num
            kx= 1./box_size*i*2*pi
            ky= 1./box_size*j*2*pi
            kz= 1./box_size*k*2*pi
            if (i .gt. grid_num/2) kx=1./box_size*(i-grid_num)*2*pi
            if (j .gt. grid_num/2) ky=1./box_size*(j-grid_num)*2*pi
            if (k .gt. grid_num/2) kz=1./box_size*(k-grid_num)*2*pi
! CIC binning
           if(CIC_binning .eq. .true.) then
            temp =sqrt(kx**2+ky**2+kz**2)/k_size+0.5
            k_index =int(temp)+1
            mu = 1.*kz/sqrt(1.*kx**2+ky**2+kz**2)
         pk_temp =REAL(vzarr_in(i,j,k))**2+AIMAG(vzarr_in(i,j,k))**2
            if(k_index .eq. 1) then
               pk(k_index) = pk(k_index) + pk_temp
               k_num_sep(k_index) = k_num_sep(k_index) + 1
               if(mu.lt.1 .and. mu.gt.-1) then
            pk_quad(k_index) = pk_quad(k_index) + pk_temp*(3.*mu**2-1)/2
            quad_modes(k_index) = quad_modes(k_index) + (3.*mu**2-1)/2
               end if
            else
               pk(k_index) = pk(k_index) + pk_temp*(temp-k_index+1)
               pk(k_index-1) = pk(k_index-1) + pk_temp*(-temp+k_index)
             k_num_sep(k_index) = k_num_sep(k_index) + (temp-k_index+1)
             k_num_sep(k_index-1)=k_num_sep(k_index-1)+(-temp+k_index)
               if(mu.lt.1 .and. mu.gt.-1) then
            pk_quad(k_index) = pk_quad(k_index) + 
     &                  pk_temp*(3.*mu**2-1)/2*(temp-k_index+1)
            pk_quad(k_index-1) = pk_quad(k_index-1) + 
     &                  pk_temp*(3.*mu**2-1)/2*(-temp+k_index)
            quad_modes(k_index) = quad_modes(k_index) + 
     &                  (3.*mu**2-1)/2*(temp-k_index+1)
            quad_modes(k_index) = quad_modes(k_index) +
     &                  (3.*mu**2-1)/2*(-temp+k_index)
               end if
            end if
           else
            k_index =int(sqrt(1.*kx**2+ky**2+kz**2)/k_size)+1
            mu = 1.*kz/sqrt(1.*kx**2+ky**2+kz**2)
         pk_temp =REAL(vzarr_in(i,j,k))**2+AIMAG(vzarr_in(i,j,k))**2
            pk(k_index) = pk(k_index) + pk_temp
            k_num_sep(k_index) = k_num_sep(k_index) + 1
               if(mu.lt.1 .and. mu.gt.-1) then
            pk_quad(k_index) = pk_quad(k_index) + pk_temp*(3.*mu**2-1)/2
            quad_modes(k_index) = quad_modes(k_index) + (3.*mu**2-1)/2
               end if
           end if
          end do
        end do
      end do


      do i=1, k_index_max
        if(k_num_sep(i) .eq. 0) goto 111
        pk(i) = pk(i)/k_num_sep(i)/num_data/num_data*(1.*box_size)**3
        pk_quad(i) = pk_quad(i)/
     &          k_num_sep(i)/num_data/num_data*(1.*box_size)**3

! shot noise (NGP)
        if(CIC_assignment .eq. .false.) then
          pk(i) = pk(i)-1./num_data*(1.*box_size)**3
          pk_quad(i) = pk_quad(i)-
     &     1./num_data*(1.*box_size)**3*quad_modes(i)/k_num_sep(i)
          pk_quad(i) = 5.*pk_quad(i)
        else
          kny=3.1415927*grid_num/box_size
          pk(i) = pk(i)/
     &     (1.-2./3*SIN(3.1415627*(-0.5+i)*k_size/2/kny)**2)- 
     &         1./num_data*(1.*box_size)**3
          
          pk_quad(i) = pk_quad(i)/
     &     (1.-2./3*SIN(3.1415627*(-0.5+i)*k_size/2/kny)**2)- 
     &     1./num_data*(1.*box_size)**3*quad_modes(i)/k_num_sep(i)
          pk_quad(i) = 5.*pk_quad(i)
        end if
      end do
 111  nozero_k_index_max = i-1

      if(CIC_assignment .eq. .true.) then
       mock_pkfile = Trim(output_file)//".CICassign"
      else
       mock_pkfile = Trim(output_file)//".NGPassign"
      end if

      if(CIC_binning .eq. .true.) then
       mock_pkfile = Trim(mock_pkfile)//".CICbin.pk"
      else
       mock_pkfile = Trim(mock_pkfile)//".NGPbin.pk"
      end if

      open(20,file=Trim(mock_pkfile)//".mono")
      do i=1, 2.*pi/box_size*grid_num/2/k_size ! up to k_nyquist
        write(20,*) (-0.5+i)*k_size, pk(i), k_num_sep(i)
      end do
      close(20)

!write quadrupole
!      open(20,file=Trim(mock_pkfile)//".quad")
!      do i=1, 2.*pi/box_size*grid_num/2/k_size ! up to k_nyquist
!        write(20,*) (-0.5+i)*k_size, pk_quad(i), k_num_sep(i)
!      end do
!      close(20)

!--- compute pk in redshift space---

      write(*,*) "start compute pk in redshift space"
      do i=1, k_index_max
        pk(i)=0   
        pk_quad(i)=0
        k_num_sep(i)=0
        quad_modes(i)=0
      end do
!$omp parallel do private(i,j,k)
      do i=1, grid_pdf
        do j=1, grid_pdf
           do k=1, grid_pdf
              if(CIC_assignment.eq..true.) then
            vxarr_in(i,j,k)=
     &  cmplx(AIMAG(vxarr_in(i,j,k))/cell_size**3-
     &                (1.*num_data)/(1.*grid_pdf)**3,0)
              else
            vxarr_in(i,j,k)=
     &  cmplx(AIMAG(vxarr_in(i,j,k))-
     &                (1.*num_data)/(1.*grid_pdf)**3,0)
              end if
            
          end do
        end do
      end do

      write(*,*) "starting FFT for pk"
!      FFTW3
      call sfftw_plan_dft_3d(planw2, grid_num, grid_num, grid_num,
     &  vxarr_in, vxarr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call sfftw_execute_dft(planw2,vxarr_in,vxarr_in)
      call sfftw_destroy_plan(planw2)

      time2 =  time()
      write(*,*) "finish FFT, time in sec:", time2-time1

      box_size = boxsize
!      num_data = mock_data_num

!$omp parallel do private(i,j,k,kx,ky,kz,temp,k_index,mu,pk_temp)
!$omp+ reduction(+:k_num_sep,pk_quad,quad_modes,pk)
      do i=1, grid_num
        do j=1, grid_num
          do k=1, grid_num
            kx= 1./box_size*i*2*pi
            ky= 1./box_size*j*2*pi
            kz= 1./box_size*k*2*pi
            if (i .gt. grid_num/2) kx=1./box_size*(i-grid_num)*2*pi
            if (j .gt. grid_num/2) ky=1./box_size*(j-grid_num)*2*pi
            if (k .gt. grid_num/2) kz=1./box_size*(k-grid_num)*2*pi
! CIC binning
           if(CIC_binning .eq. .true.) then
            temp =sqrt(kx**2+ky**2+kz**2)/k_size+0.5
            k_index =int(temp)+1
            mu = 1.*kz/sqrt(1.*kx**2+ky**2+kz**2)
         pk_temp =REAL(vxarr_in(i,j,k))**2+AIMAG(vxarr_in(i,j,k))**2
            if(k_index .eq. 1) then
               pk(k_index) = pk(k_index) + pk_temp
               k_num_sep(k_index) = k_num_sep(k_index) + 1
               if(mu.lt.1 .and. mu.gt.-1) then
            pk_quad(k_index) = pk_quad(k_index) + pk_temp*(3.*mu**2-1)/2
            quad_modes(k_index) = quad_modes(k_index) + (3.*mu**2-1)/2
               end if
            else
               pk(k_index) = pk(k_index) + pk_temp*(temp-k_index+1)
               pk(k_index-1) = pk(k_index-1) + pk_temp*(-temp+k_index)
             k_num_sep(k_index) = k_num_sep(k_index) + (temp-k_index+1)
             k_num_sep(k_index-1)=k_num_sep(k_index-1)+(-temp+k_index)
               if(mu.lt.1 .and. mu.gt.-1) then
            pk_quad(k_index) = pk_quad(k_index) + 
     &                  pk_temp*(3.*mu**2-1)/2*(temp-k_index+1)
            pk_quad(k_index-1) = pk_quad(k_index-1) + 
     &                  pk_temp*(3.*mu**2-1)/2*(-temp+k_index)
            quad_modes(k_index) = quad_modes(k_index) + 
     &                  (3.*mu**2-1)/2*(temp-k_index+1)
            quad_modes(k_index) = quad_modes(k_index) +
     &                  (3.*mu**2-1)/2*(-temp+k_index)
               end if
            end if
           else
            k_index =int(sqrt(1.*kx**2+ky**2+kz**2)/k_size)+1
            mu = 1.*kz/sqrt(1.*kx**2+ky**2+kz**2)
         pk_temp =REAL(vxarr_in(i,j,k))**2+AIMAG(vxarr_in(i,j,k))**2
            pk(k_index) = pk(k_index) + pk_temp
            k_num_sep(k_index) = k_num_sep(k_index) + 1
               if(mu.lt.1 .and. mu.gt.-1) then
            pk_quad(k_index) = pk_quad(k_index) + pk_temp*(3.*mu**2-1)/2
            quad_modes(k_index) = quad_modes(k_index) + (3.*mu**2-1)/2
               end if
           end if
          end do
        end do
      end do


      do i=1, k_index_max
        if(k_num_sep(i) .eq. 0) goto 112
        pk(i) = pk(i)/k_num_sep(i)/num_data/num_data*(1.*box_size)**3
        pk_quad(i) = pk_quad(i)/
     &          k_num_sep(i)/num_data/num_data*(1.*box_size)**3

! shot noise (NGP)
        if(CIC_assignment .eq. .false.) then
          pk(i) = pk(i)-1./num_data*(1.*box_size)**3
          pk_quad(i) = pk_quad(i)-
     &     1./num_data*(1.*box_size)**3*quad_modes(i)/k_num_sep(i)
          pk_quad(i) = 5.*pk_quad(i)
        else
          kny=3.1415927*grid_num/box_size
          pk(i) = pk(i)/
     &     (1.-2./3*SIN(3.1415627*(-0.5+i)*k_size/2/kny)**2)- 
     &      1./num_data*(1.*box_size)**3
          pk_quad(i) = pk_quad(i)/
     &     (1.-2./3*SIN(3.1415627*(-0.5+i)*k_size/2/kny)**2)- 
     &     1./num_data*(1.*box_size)**3*quad_modes(i)/k_num_sep(i)
          pk_quad(i) = 5.*pk_quad(i)
        end if
      end do
 112  nozero_k_index_max = i-1

      if(CIC_assignment .eq. .true.) then
       mock_pkfile = Trim(output_file)//".CICassign"
      else
       mock_pkfile = Trim(output_file)//".NGPassign"
      end if

      if(CIC_binning .eq. .true.) then
       mock_pkfile = Trim(mock_pkfile)//".CICbin.pk"
      else
       mock_pkfile = Trim(mock_pkfile)//".NGPbin.pk"
      end if

      open(20,file=Trim(mock_pkfile)//".zdist.mono")
      do i=1, 2.*pi/box_size*grid_num/2/k_size ! up to k_nyquist
        write(20,*) (-0.5+i)*k_size, pk(i), k_num_sep(i)
      end do
      close(20)

!write quadrupole
      open(20,file=Trim(mock_pkfile)//".zdist.quad")
      do i=1, 2.*pi/box_size*grid_num/2/k_size ! up to k_nyquist
        write(20,*) (-0.5+i)*k_size, pk_quad(i), k_num_sep(i)
      end do
      close(20)

      deallocate(dataxyz)
      deallocate(dataxyz_omp)
      deallocate(vxarr_in)
      deallocate(vyarr_in)
      deallocate(vzarr_in)
      deallocate(temp_array)
      deallocate(temp_array2)
      deallocate(temp_array3)
      deallocate(temp_int_array)      
      
      end subroutine zeldovich_ic
! -------------------------------------------------------
      subroutine convert_corr2pk(r_array,corr_array,r_num,
     &   k_array,pk_array,k_num,volume,boxsize)
      implicit none
      integer,intent(in) :: r_num
      real, intent(in) :: r_array(r_num)
      real, intent(in) :: corr_array(r_num)
      integer,intent(in) :: k_num
      real, intent(in) :: k_array(k_num)
      real, intent(out) :: pk_array(k_num)
      real, intent(in) :: volume
      integer,intent(in) :: boxsize
 
      integer :: k
      integer :: r_index

      real :: rmin, rmax
      integer, parameter :: r_num_max = 100000
      integer, parameter :: r_index_max = 1000000
      real, parameter :: pi = 3.1415927
      real    :: corr(r_num_max), r(r_num_max), y2(r_num_max)
      real    :: thisr, nextr, this_corr

      rmin = r_array(1)
      rmax = r_array(r_num)
      do k=1, k_num
          pk_array(k) = 0
      end do

      do r_index = 1, r_num
         r(r_index) = r_array(r_index)
         corr(r_index) = corr_array(r_index)
      end do

      call spline(r,corr,r_num,1.0e30,1.0e30,y2)
      nextr = 10**(1.0*(log10(rmax/rmin)/r_index_max))*rmin     

      do r_index=1,r_index_max
!        if(nextr .lt. boxsize/2) then
          thisr=nextr
	  nextr=10**(1.0*(r_index+1)*(log10(rmax/rmin)/r_index_max))*rmin
          call splint(r,corr,y2,r_num,thisr,this_corr)
          do k=1, k_num
              pk_array(k)= pk_array(k)+this_corr*sin(thisr*k_array(k))*
     &                  thisr/k_array(k)*(nextr-thisr)*4.*pi/volume
          end do
!        end if
      end do

      end subroutine convert_corr2pk

!-----------------------------------------------------------

      subroutine poisson_rand(mean, iseed, output,thread_id)
      implicit none
      real,intent(in)     :: mean
      integer,intent(inout)  :: iseed
      integer,intent(out) :: output

      integer,parameter :: kmax = 300
      integer :: k, thread_id!, omp_get_thread_num
      real    :: el, p, ran3_omp
!!!C$OMP PARALLEL PRIVATE(k,el,p) firstprivate(mean) 
!!!C$OMP+lastprivate(output)      

      k  = 0
      el = exp(-mean)
      p  = 1

      do while(p > el .and. k < kmax)
          p = ran3_omp(iseed, thread_id)*p
          k = k + 1
      end do

      output = k - 1
      if(output .lt. 0) output = 0 
!!!C$OMP END PARALLEL

      end subroutine poisson_rand

!---------------------------------------------

      subroutine gauss_ran(iseed, y1, y2,thread_id)
      implicit none
      integer :: iseed,thread_id
      real    :: y1, y2

      real :: w
      real :: x1, x2, ran3_omp

!!!C$OMP PARALLEL PRIVATE(w,x1,x2)      
      w = 2
      do while ( w .ge. 1.0)
         x1 = 2.0 * ran3_omp(iseed,thread_id) - 1.0
         x2 = 2.0 * ran3_omp(iseed,thread_id) - 1.0
          w = x1 * x1 + x2 * x2
      end do

      w = sqrt( (-2.0 * log( w ) ) / w )
      y1 = x1 * w
      y2 = x2 * w
!!!C$omp end parallel

      end subroutine !gauss_ran
!--------------------------------------------------------

!      subroutine poisson_rand_catalog(boxsize,grid_num,density,
!     &  randfile,iseed)
!      implicit none
!! input
!      real,intent(in)    :: boxsize
!      integer,intent(in) :: grid_num
!      real,intent(in)    :: density
!      character*500,intent(in) :: randfile
!      integer            :: iseed
!
!      integer :: rx, ry, rz, pr, count=0, i, j
!      real    :: mean
!      real,allocatable :: datax(:),datay(:),dataz(:)
!      integer,allocatable :: datanp(:)
!      integer :: data_num_max
!
!      data_num_max = int(1.5*grid_num**3*density)
!      allocate(datax(data_num_max))
!      allocate(datay(data_num_max))
!      allocate(dataz(data_num_max))
!      allocate(datanp(data_num_max))
!
!      mean = density
!
!      do rx = 1, grid_num
!        do ry = 1, grid_num
!          do rz = 1, grid_num
!            call poisson_rand(mean, iseed, pr)
!            if (pr .gt. 0) then
!               count = count + 1
!               datax(count) = 1.*rx*boxsize/grid_num
!               datay(count) = 1.*ry*boxsize/grid_num
!               dataz(count) = 1.*rz*boxsize/grid_num
!               datanp(count)= pr
!            end if
!          end do
!        end do
!      end do
!      write(*,*) "finish generating ", count,
!     &  " random data in box ", boxsize, "cube"
!      open(unit = 10, file= Trim(randfile))
!      do i = 1, count
!         do j = 1, datanp(i)
!           write(10,*) datax(i), datay(i), dataz(i)
!         end do
!      end do
!      close(10)
!
!      deallocate(datax)
!      deallocate(datay)
!      deallocate(dataz)
!      deallocate(datanp)
!   
!      end subroutine
!
!-------------------------------------------------
      subroutine convert_pk2corr(k_array,pk_array,k_num,corr_array,
     &  dr,r_max)
      implicit none
      integer,intent(in) :: k_num
      real, intent(in) :: k_array(k_num)
      real, intent(in) :: pk_array(k_num)
      real, intent(out) :: corr_array(2000)
      real, intent(in) :: dr, r_max
 
      real    :: r
      integer :: r_num
      integer :: k_index, r_index

!      real, parameter    :: kmin = 0.0001, kmax = 100
      real    :: kmin, kmax
      integer, parameter :: k_num_max = 1000000!5000
      integer, parameter :: k_index_max = 1000000!5000
      real    :: pk(k_num_max), k(k_num_max), y2(k_num_max)
      real    :: thisk, nextk, this_pk
      real,parameter :: pi=3.1415927

      if (k_num .gt. k_num_max) then
        write(*,*) "increase k_num_max in convert_pk2corr()"
        stop
      end if

      r_num = int(r_max/dr)

      do r_index=1, r_num
          corr_array(r_index) = 0
      end do

      open(18,file="../fortran_output/testpk.dat")
      do k_index = 1, k_num
         k(k_index) = k_array(k_index)
         pk(k_index) = pk_array(k_index)
!test
         write(18,*) k(k_index), pk(k_index)
      end do

      kmin = k_array(1)
      kmax = 1!k_array(k_num)

      call spline(k,pk,k_num,1.0e30,1.0e30,y2)
      nextk = 10**(1.0*(log10(kmax/kmin)/k_index_max))*kmin     
      do k_index=1,k_index_max-1
          thisk=nextk
	  nextk=10**(1.0*(k_index+1)*(log10(kmax/kmin)/k_index_max))*kmin
          call splint(k,pk,y2,k_num,thisk,this_pk)
          do r_index=1, r_num
              r = dr*(1.0*r_index-0.5)
              corr_array(r_index)= corr_array(r_index)+ 
     &                this_pk*sin(thisk*r)*thisk/r*
     &          (nextk-thisk)
          end do
      end do

      do r_index=1, r_num
         corr_array(r_index) = corr_array(r_index)/2/pi**2
      end do

      end subroutine convert_pk2corr
!-------------------------------------
c-------------- ran3_omp_omp() ---------------------------------- 
      function ran3_omp(idum, thread_id) 
! random number generater for openmp with max_threads = 512
c modified by albert, May 1 2014
      INTEGER idum
      INTEGER MBIG,MSEED,MZ
      REAL    ran3_omp,FAC
      PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1./MBIG)
      INTEGER i,iik
      INTEGER mj,mk,thread_id
      INTEGER,parameter :: max_thread = 512
      INTEGER iff(max_thread),inext(max_thread),inextp(max_thread),
     &   ma(max_thread,55)
      SAVE iff,inext,inextp,ma
      DATA iff /max_thread*0/

!C$OMP PARALLEL PRIVATE(thread_id,mj,mk,i,ii,k)

!      thread_id = omp_get_thread_num()+1

      if(idum.lt.0.or.iff(thread_id).eq.0)then
        iff(thread_id)=1
        mj=MSEED-iabs(idum)
        mj=mod(mj,MBIG)
        ma(thread_id,55)=mj
        mk=1
        do 11 i=1,54
          ii=mod(21*i,55)
          ma(thread_id,ii)=mk
          mk=mj-mk
          if(mk.lt.MZ)mk=mk+MBIG
          mj=ma(thread_id,ii)
11      continue
        do 13 k=1,4
          do 12 i=1,55
            ma(thread_id,i)=ma(thread_id,i)-ma(thread_id,1+mod(i+30,55))
            if(ma(thread_id,i).lt.MZ) ma(thread_id,i)=
     &                                  ma(thread_id,i)+MBIG
12        continue
13      continue
        inext(thread_id)=0
        inextp(thread_id)=31
        idum=1
      endif
      inext(thread_id)=inext(thread_id)+1
      if(inext(thread_id).eq.56)inext(thread_id)=1
      inextp(thread_id)=inextp(thread_id)+1
      if(inextp(thread_id).eq.56)inextp(thread_id)=1
      mj=ma(thread_id,inext(thread_id))-ma(thread_id,inextp(thread_id))
      if(mj.lt.MZ)mj=mj+MBIG
      ma(thread_id,inext(thread_id))=mj
      ran3_omp=mj*FAC

!C$OMP END PARALLEL

      return
      END function ! ran3_omp_omp()
c------------------------------------------
!---------------------------------------------------
      subroutine fast_code_CF_box_linear(x,y,w,wt,data_num,
     &  box_size,num_r,num_mu,ddnofr,rrnofr)
      integer,intent(in) :: data_num
      integer,intent(in) :: num_r
      integer,intent(in) :: num_mu
      real,intent(in) :: x(data_num), y(data_num),
     &                   w(data_num), wt(data_num)
      real,intent(inout) :: DDnofr(num_r,num_mu), RRnofr(num_r,num_mu)
      real          :: box_size

      integer       :: iseed = -4321
      integer       :: count
      real          :: ra, dec, z, r1,r2,r3,r4
      integer       :: i, j, ii, jj,kk
      real          :: totalwt_data, totalwt_rand, corr
      real          :: x_min,x_max,y_min,y_max,z_min,z_max
      real          :: ddbin, d2rbin,rrbin, rrbin_de
      integer       :: nx, ny, nz, xindex,yindex,zindex
      integer       :: max_n_cell
      integer,parameter :: max_cells=300000
      integer :: n_cell_data(max_cells), n_cell_rand(max_cells)
      integer :: ijk2n,nindex,rr,mm
      real,allocatable :: data_cell_array(:,:,:), rand_cell_array(:,:,:)
      logical :: is_neighbor
      real,parameter :: pi=3.1415927

      do i=1, max_cells
        n_cell_data(i)=0
      end do

!find x_min, x_max, y...,z...
      x_min=0
      y_min=0
      z_min=0

!determine numbers of cells nx x ny x nz
      nx = int(box_size/num_r)
      ny = nx
      nz = nx
      write(*,*) "the number of cells are",nx,"x",nx,"x",nx,"=",nx*ny*nz
!!determine max number in cells and construct arrays
      if(max_cells .lt. nx*ny*nz) then
        write(*,*) "should increase max_cells"
        stop
      end if
      do i=1,nx*ny*nz
        n_cell_data(i)=0
      end do

      do i=1, data_num
        xindex=int((x(i))/num_r)
        yindex=int((y(i))/num_r)
        zindex=int((w(i))/num_r)

        n_cell_data(ijk2n(xindex,yindex,zindex,nx,ny,nz))=     
     &    n_cell_data(ijk2n(xindex,yindex,zindex,nx,ny,nz))+1
      end do
      max_n_cell=maxval(n_cell_data)
      write(*,*) "max number in the cells is", max_n_cell, "(data)"
      allocate(data_cell_array(nx*ny*nz,max_n_cell,4)) !4 for x,y,w,wt


!assign data to cells
      do i=1,nx*ny*nz
        n_cell_data(i)=0
      end do

      do i=1,data_num
        xindex=int((x(i)-x_min)/num_r)
        yindex=int((y(i)-y_min)/num_r)
        zindex=int((w(i)-z_min)/num_r)
        nindex=ijk2n(xindex,yindex,zindex,nx,ny,nz) !index = ix+iy*nx+iz*nx*ny +1
        n_cell_data(nindex)=n_cell_data(nindex)+1
        data_cell_array(nindex,n_cell_data(nindex),1)=x(i)     
        data_cell_array(nindex,n_cell_data(nindex),2)=y(i)     
        data_cell_array(nindex,n_cell_data(nindex),3)=w(i)     
        data_cell_array(nindex,n_cell_data(nindex),4)=wt(i)
      end do     


!loop to go over all the pairs of cells to compute counts of galaxy pairs

      do i=1,num_r
        do j=1,num_mu
          ddnofr(i,j)=0
          RRnofr(i,j)=0
        end do
      end do 

!!--------------compute RR--------------------!!
      do i=1,num_r
        do j=1,num_mu
          RRnofr(i,j)=4.*pi/3*(i**3-(i-1)**3)*1./num_mu/box_size**3/2
        end do
      end do

      count = 0
      do i=1, nx*ny*nz
         do j=i, nx*ny*nz

           if(i.eq.j) then

            call twoD_XX_corr_M_box(data_cell_array(i,:,1),
     &                          data_cell_array(i,:,2),
     &                          data_cell_array(i,:,3),
     &                          data_cell_array(i,:,4),
     &                          n_cell_data(i),DDnofr,num_r,num_mu)
           elseif(is_neighbor(i,j,nx,ny,nz).eqv..true.) then
             count = count +1
            call twoD_XY_corr_M_box(data_cell_array(i,:,1),
     &                          data_cell_array(i,:,2),
     &                          data_cell_array(i,:,3),
     &                          data_cell_array(i,:,4),
     &                          n_cell_data(i),
     &                          data_cell_array(j,:,1),
     &                          data_cell_array(j,:,2),
     &                          data_cell_array(j,:,3),
     &                          data_cell_array(j,:,4),
     &                          n_cell_data(j),DDnofr,num_r,num_mu,
     &                          box_size) 
           end if
         end do !j
!         write(*,*) "done with cell",i,"of",nx*ny*nz

       end do !i

       write(*,*) "number of neighbor pairs is", count

       deallocate(data_cell_array)
       end subroutine
c---------------------------------------
c------------------------------------------
      function ijk2n(xindex,yindex,zindex,nx,ny,nz)
      implicit none
      integer :: xindex,yindex,zindex,nx,ny,nz,ijk2n

      ijk2n=xindex+yindex*nx+zindex*ny*nx+1
      return
      end function
C---------------------------------------------
!-----------------------------------------------------------------
      SUBROUTINE twoD_XX_corr_M_box(x,y,w,wt,data_num,nofr,max_r,num_mu)
c calculate DD or RR nofr with matrix input and output
c  -------------  blueprint ------------
c LRGsourcefile_xyz is xyz data file
c nofrfile is nofr file(DD or RR)
c --------------------------------------

      IMPLICIT NONE
      integer, intent(in) :: data_num
      REAL,    intent(in) :: x(data_num),y(data_num)
      REAL,    intent(in) :: w(data_num),wt(data_num)
      integer, intent(in) :: max_r, num_mu
      REAL,intent(inout)  :: nofr(max_r,num_mu)
      integer             :: i, j, mu_index
      real                :: mu,max_r2,dr
      real                :: xc,yc,zc,rc,dx,dy,dz,dr2,dr_r,dr_a
      REAL                :: nofr_temp(max_r,num_mu)
     
      max_r2 = max_r**2
      do i=1, max_r
        do j=1, num_mu
           nofr_temp(i,j) = 0
        end do
      end do

!openmp start
C$omp parallel
C$omp+private(j,xc,yc,zc,dx,dy,dz,rc,dr2,dr_r,dr,mu,mu_index)
C$omp+reduction(+:nofr_temp)
C$omp do
      
      Do i=1,data_num
         Do j=i+1,data_num
            xc = (x(i)+x(j))/2
            yc = (y(i)+y(j))/2
            zc = (w(i)+w(j))/2
            dx = x(i)-x(j)
            dy = y(i)-y(j)
            dz = w(i)-w(j)
            dr2 = dx**2+dy**2+dz**2
            if(dr2 .lt. max_r2 .and. dr2 .gt. 0.000001) then

               dr = sqrt(dr2)
               mu = abs(dz)/dr
               mu_index = int(mu*num_mu)+1
               if(mu_index .gt. num_mu) then
!                 write(*,*) "on the edge:mu"
!                 stop
                  mu_index = num_mu
               end if
               if(int(dr)+1 .gt. max_r) then
!                 write(*,*) "on the edge:r"
!                 stop
                  dr = dr-1
               end if


              nofr_temp(int(dr)+1,mu_index)=
     &          nofr_temp(int(dr)+1,mu_index)+wt(i)*wt(j)

            end if
         end do
      end do
C$omp end do
C$omp end parallel
      do i=1, max_r
        do j=1, num_mu
           nofr(i,j) = nofr(i,j) + nofr_temp(i,j)
        end do
      end do

      END SUBROUTINE ! twoD_XX_corr_M()
c---------------------------------------

      SUBROUTINE twoD_XY_corr_M_box(
     & x,y,w,wt,data_num1,x2,y2,w2,wt2,data_num2,nofr,max_r,num_mu,
     & box_size)
c calculate number density of dr (n(r)) for DR with matrix input and output
c  -------------  blueprint ------------
c LRGsourcefile_xyz,RANsourcefile_xyz are LRG and Random data
c nofrfile is filename of nofr file of DR
c --------------------------------------

      IMPLICIT NONE
      integer, intent(in) :: data_num1, data_num2
      REAL,    intent(in) :: x(data_num1),y(data_num1)
      REAL,    intent(in) :: w(data_num1),wt(data_num1)
      REAL,    intent(in) :: x2(data_num2),y2(data_num2)
      REAL,    intent(in) :: w2(data_num2),wt2(data_num2)
      integer, intent(in) :: max_r, num_mu
      REAL, intent(inout) :: nofr(max_r,num_mu)
      real, intent(in)    :: box_size

      integer             :: i, j, mu_index
      real                :: mu, dr, max_r2, dx2,dy2,dz2
      real                :: xc,yc,zc,dx,rc,dy,dz,dr2,dr_r,dr_a
      real                :: nofr_temp(max_r,num_mu)
     
      max_r2 = max_r**2
      do i=1, max_r
        do j=1, num_mu
           nofr_temp(i,j) = 0
        end do
      end do

!openmp start
C$omp parallel
C$omp+private(j,xc,yc,zc,dx,dy,dz,rc,dr2,dr,dr_r,mu,mu_index)
C$omp+reduction(+:nofr_temp)
C$omp do
      Do i=1,data_num1
         Do j=1,data_num2
            xc = (x(i)+x2(j))/2
            yc = (y(i)+y2(j))/2
            zc = (w(i)+w2(j))/2
            dx = min(abs(x(i)-x2(j)),(box_size-abs(x(i)-x2(j))))
            dy = min(abs(y(i)-y2(j)),(box_size-abs(y(i)-y2(j))))
            dz = min(abs(w(i)-w2(j)),(box_size-abs(w(i)-w2(j))))
            dr2 = dx**2+dy**2+dz**2
            if(dr2 .lt. max_r2 .and. dr2 .gt. 0.000001) then
!               rc = sqrt(xc**2+yc**2+zc**2)
               dr = sqrt(dr2)
!               dr_r = ABS(dx*xc+dy*yc+dz*zc)/rc
               mu = abs(dz)/dr
               mu_index = int(mu*num_mu)+1

               if(mu_index .gt. num_mu) then
!                 write(*,*) "on the edge:mu"
!                 stop
                 mu_index = num_mu
               end if
               if(int(dr)+1 .gt. max_r) then
!                 write(*,*) "on the edge:r"
!                 stop
                 dr = dr -1
               end if

           !!$OMP CRITICAL
              nofr_temp(int(dr)+1,mu_index)=
     &          nofr_temp(int(dr)+1,mu_index)+wt(i)*wt2(j)
           !!$OMP END CRITICAL
            end if
            
         end do
      end do
C$omp end do
C$omp end parallel

      do i=1, max_r
        do j=1, num_mu
           nofr(i,j) = nofr(i,j) + nofr_temp(i,j)
        end do
      end do

      END SUBROUTINE ! twoD_XY_corr_M()  

c -------------------------------------------------------------------
      function is_neighbor(i,j,nx,ny,nz)
!for box
      implicit none
      integer :: i,j,nx,ny,nz
      logical :: is_neighbor
      
      integer :: ix,iy,iz,jx,jy,jz

      call n2ijk(i,ix,iy,iz,nx,nx,nx)
      call n2ijk(j,jx,jy,jz,nx,nx,nx)
      if((abs(ix-jx).le.1 .or. abs(ix-jx).eq.(nx-1)) .and.
     &   (abs(iy-jy).le.1 .or. abs(iy-jy).eq.(nx-1)).and.
     &   (abs(iz-jz).le.1 .or. abs(iz-jz).eq.(nx-1))) then
         is_neighbor = .true.
      else
         is_neighbor = .false.
      end if
      return
      end function
C-----------------------------------------------------
      subroutine n2ijk(i,ix,iy,iz,nx,ny,nz)
      implicit none
      integer :: i,ix,iy,iz,nx,ny,nz

      ix = mod(i-1,nx)
      iy = mod((i-1-ix)/nx,ny)
      iz = ((i-1-ix)/nx-iy)/ny

      end subroutine
      
