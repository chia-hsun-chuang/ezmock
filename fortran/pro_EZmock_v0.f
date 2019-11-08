!note:
      program EZmock
! Compute statistics of sim and generate EZmock catalog,  using white noise to generate IC
! produce ZA mock with scatter
! make clean
! make
      implicit none
!MPI
!      include 'mpif.h'
!      include 'fftw_f77.i'
!output file prefix
      character*200     :: datafile_xyz =
!     & "../"//
!     & "Box_HAM_z0.561800_nbar4.271000e-04_scat0.1901_noincom.dat"
!     & "real_Box_HAM_z0.561800_nbar3.128785e-04.dat"
!     & "Planck_box2500_FOF02_snap29_no_masscut10000.dat"
     & "../reference_cats/"//
!     & "cat_BMD_Planck_nmean_3.5e-4_snap0029.out"
     & "cat_BMD_Planck_nmean_2.e-4_snap0029.out"

      character*200     :: datafile_prefix = 
!     &"FOF02_A-0d09_scat10_0_modpk190_gr0d3_an0d1_den2d5_4_zdist102" 
!     &"FOF_A-0d08_sc10_0_pk170_gr0d3_an0d1_den2d5_3d7_zd101_0_wn960" !test constant fog (no dependence on velocity)
!
!     &"BDM3d5E4_A-0d0_sc10_0_pk170_gr0d3_an0d1_den2d4_3d9_zd101_130"
     &"BDM2d0E4_A-0d04_sc10_0_pk170_gr0d5_an0d1_den2d5_4d6_zd78_100ss"

!      integer,parameter :: start_rank = 0
!      integer,parameter :: total_catalog = 1 !total --> not end_rank

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
      character*200  :: whitenoise_file = 
     & "../BigMD/WhiteNoise_ICS/"//
     & "BigMD_960_wn_ascii.dat"
      logical        :: use_whitenoise_file = .true. 

      character*200     :: pkfile =
     & '../'//
     & 'PlanckDM.linear.pk'
      character*200     :: pknwfile =
     & '../'//
     & 'PlanckDM.nowiggle.pk'

      character*200     :: datafile_path =
     &  "../mock0704/"

      character*200 :: output_suffix = ".dat"

      character*200 :: output_file
      integer       ::  data_num

c... for mpi
!      integer :: my_rank
!      integer,parameter :: server = 0
 !     integer tag
 !     integer size
 !     integer status(MPI_STATUS_SIZE)
 !     integer ierr

      CHARACTER*5  :: rank_string='00000'
      integer      :: i, num_loop,j,k,l, ip,grid_pdf
      real :: temp,ran3_omp
      integer :: iseed0 = -1234
      integer :: iseed
      integer,parameter :: max_data=100000000
      real :: ran3
      integer :: omp_get_max_threads
      real,parameter    :: data_density = 0. !for output raw particle data !per grid cell
!-- 2PCF
!      character*200 :: datafile
      logical :: compute_CF, compute_CF_zdist
      integer       :: skiplines = 6
      character*200 :: twod_corr_suffix = '.bin5.corr'
      real          :: max_r = 250
      integer       :: bin_size = 5
      real          :: redshift = 0.5618
      real          :: om = 0.307115 !planck
      real :: inverse_aH
!!
      character*200 :: twod_nofr_file
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
      integer :: time, time1, time2

!read input file
       namelist /EZmock_v0_input/
     &  datafile_xyz, datafile_prefix, boxsize, grid_num,
     & grow2z0, scatter, scatter2, modify_pk, modify_pdf,
     & antidamping, density_cut, density_sat, zdist_rate,
     & zdist_fog, whitenoise_file, use_whitenoise_file,
     & pkfile, pknwfile, datafile_path, output_suffix,
     & compute_CF, compute_CF_zdist,
     & skiplines, twod_corr_suffix, max_r, bin_size,
     & redshift, om

!input                                                                                                                 
      read(*,EZmock_v0_input)

      time1=time()
      grid_pdf=grid_num
!      call MPI_Init(ierr)
!      call MPI_Comm_rank(MPI_COMM_WORLD, my_rank, ierr)
!      call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

!      write(*,*) "MPI size=",size
!      num_loop = total_catalog/size
!      write(*,*) "num of loop=",num_loop
!      write(*,*) "MPI my_rank=",my_rank

      z = redshift
      Hubble = 100.*(sqrt((1.+z)*(1.+z)*(1.+z)*om+1.-om))! h units
      inverse_aH=(1.+z)/Hubble
      write(*,*) "1/aH =", inverse_aH
!--
!      do ip = 1, num_loop !generate num_loop catalogs per core
!        Write(rank_string,'(I5)') 10000 + my_rank+start_rank+size*(ip-1)
!        write(*,*) "file rank:",rank_string
!        iseed = -iseed0 - (my_rank+start_rank+size*(ip-1))*
!     &      omp_get_max_threads()
        iseed = -iseed0+
     &      omp_get_max_threads()
          write(*,*) "number of threads is",omp_get_max_threads()

        output_file = Trim(datafile_path)//Trim(datafile_prefix)
!     &   //Trim(rank_string)//Trim(output_suffix)
     &   //Trim(output_suffix)

        call zeldovich_ic(datafile_xyz,skiplines,pkfile,pknwfile,
     &  grow2z0, boxsize,grid_num,grid_pdf,
     &  data_density, iseed,data_num, output_file,max_data,scatter,
     &  scatter2,modify_pk,modify_pdf,antidamping,density_cut,
     &  density_sat, zdist_rate, zdist_fog, inverse_aH,
     &  whitenoise_file, use_whitenoise_file)
        write(*,*)
     &      "finish generating catalog in cubic box. Size:",boxsize,
     &  "output file:", Trim(output_file)

!test
!       write(*,*) "do not generate CF for testing"
!       stop

!--- insert 2PCF code here
! initialize
      num_r=int(max_r)
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

      call fast_code_CF_box_linear(x,y,w,wt,data_num,
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
          ddbin=ddbin/data_num/data_num
          write(11,*) 0.5*bin_size+rr-1, ddbin/rrbin-1
      end do   
      close(11)

! redshift space
      open(unit=30,file=Trim(output_file),status='old')
      do i=1,skiplines
        read(30,*)
      end do

      do i=1,data_num
          read(30,*) x(i),y(i),w(i),r1,r2,vz(i)
          wt(i)=1
!          w(i)=w(i)+vz(i)*inverse_aH

          Hubble = 100*(sqrt((1+z)*(1+z)*(1+z)*om+1.-om))! h units
          w(i)=w(i)+vz(i)*(1+z)/Hubble

          if(w(i).gt.boxsize) w(i) = w(i) - boxsize
          if(w(i).lt.0.) w(i) = w(i) + boxsize;
      end do

      close(30)
      write(*,*) "finish reading data to arrays"

      call fast_code_CF_box_linear(x,y,w,wt,data_num,
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
          ddbin=ddbin/data_num/data_num
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
        ddbin=ddbin/data_num/data_num
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
!      end do
      time2 =  time()
      write(*,*) "time in sec:", time2-time1


!      call MPI_Finalize(ierr)
      end program

!-----------------------------------------------------------------

      subroutine zeldovich_ic(datafile_xyz, skip_line,pkfile,
     &     pknwfile,grow2z0,boxsize,grid_num,grid_pdf,density,iseed,
     &     data_num, output_file,max_data,scatter,scatter2,modify_pk,
     &     modify_pdf, antidamping, density_cut,density_sat,zdist_rate,
     &     zdist_fog, inverse_aH,whitenoise_file, use_whitenoise_file)
      implicit none
!      include 'fftw_f77.i'
!      FFTW3
      include 'fftw3.f'
      character*200,intent(in) :: datafile_xyz
      integer,intent(in)          :: skip_line
      character*200,intent(in) :: pkfile
      character*200,intent(in) :: pknwfile
      real,intent(in)          :: grow2z0
      real,intent(in)          :: boxsize
      integer,intent(in)       :: grid_num
      real,intent(in)          :: density
      integer,intent(in)       :: iseed
      integer,intent(out)      :: data_num
      character*200,intent(in) :: output_file
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
      character*200,intent(in) :: whitenoise_file
      logical,intent(in)       :: use_whitenoise_file
!!test
!      real :: scatter2=0.
      logical,parameter :: zeldovich = .true.
!compute pk
      real    :: k_size = 0.005
      logical :: CIC_assignment = .true.
      character*200 :: pkfile_path =
     & "../pk0421/"
      character*200 :: pkfile_prefix =
     & "dk0.005"
!
      logical :: CIC_binning = .false. !I have shown that it's useless
      integer,parameter :: r_num_max = 10000, k_num_max = 10000
      real,parameter    :: kmin = 0.0001, kmax = 100
      real,parameter    :: pi = 3.1415927
!      integer,parameter :: skip_line=2
      character*200     :: sim_pdf_file
      character*200     :: mock_pkfile
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
      real    :: kx, ky, kz, k_mag, pk_out, temp, ak, bk,pk_outo
      real    :: kx1, ky1, kz1, re_in, im_in
!      real,allocatable :: temp_array(:,:,:,:)
      integer :: np, data_index(100)
      integer :: nn(3), ndim=3, isign=-1
      real    :: gridsize, ran3_omp, xx, yy, zz, ran3, temp3
      integer :: sirko
      integer :: totalcount = 0
      real    :: rescaleG, rescaleLN, temp2, sigma2_densityG, tempo
      integer ::  readtime, xi, yi, zi,xii,yii,zii, iseed2=-1234
c... for fftw
      integer*8  :: plan,planx,plany,planz,planm,palno,planw
      double complex, dimension(grid_num,grid_num,grid_num) ::
     & vxarr_in, vyarr_in, vzarr_in, vwarr_in
     & , arr_in
      real :: temp_array(grid_num,grid_num,grid_num)
      real :: temp_array2(grid_num,grid_num,grid_num)
      real :: temp_array3(grid_num,grid_num,grid_num)
      integer,parameter :: max_num = 200
!      real :: f_exp = 1. 
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
      integer :: pdf2(0:max_num), temp_int
      integer :: pdf_convert(0:max_num)
      real :: pdf_convert_ratio(0:max_num)
      real :: np_ratio, exp
      real :: dataxyz(6,max_data)
      real :: dataxyz_omp(6,max_data/100,100)
      integer :: time, time1, time2

      integer,parameter :: k_index_max = 10000
      real              :: pk(k_index_max),pk_quad(k_index_max)
      integer ::  num_data, mock_data_num, i1,i2,i3
      real    :: mu,kny, pk_temp, dr, norm1
      integer :: k_index,nozero_k_index_max, box_size
      real    :: k_num_sep(k_index_max), quad_modes(k_index_max)
!openmp
      integer :: fftw_init_threads, omp_get_max_threads, iret
      integer :: omp_get_thread_num, thread_id, iseedxx(512) !iseed for threads
      integer :: rxx, ryy, rzz, rxi,ryi,rzi, zzi, zzii
      integer :: x_div_omp, y_div_omp, z_div_omp, total_num_thread
!adjust pdf
      real :: A_pdf, B_pdf, sum_pdf,sum_pdf_A

!intitalize for doing cic with openmp (total_num_threads = x_div_omp*y_div_omp*z_div_omp)
      total_num_thread = omp_get_max_threads()
      x_div_omp = 2
      y_div_omp = 2

      if(mod(total_num_thread, x_div_omp*y_div_omp) .ne. 0) then
         write(*,*) "total threads is not integer times 2*2"
         stop  
      end if

      z_div_omp = total_num_thread/x_div_omp/y_div_omp
      if(mod(grid_num,x_div_omp).ne.0 .or.
     & mod(grid_num,y_div_omp).ne.0 .or.
     & mod(grid_num,z_div_omp).ne.0) then
          write(*,*) 
     & "grid size is not proper divided by x_div_omp(or y or z)"
       stop
      end if

      write(*,*)
      time1 = time() 
      write(*,*) "start count time in sec:", 0 

!initial iseed for  openmp
C$OMP PARALLEL PRIVATE(thread_id)
      thread_id =  omp_get_thread_num()+1
      iseedxx(thread_id) =  iseed + thread_id
C$OMP END PARALLEL

! initialize
C$OMP PARALLEL DO PRIVATE(rx,ry,rz)
      do rx=1,grid_num
        do ry=1,grid_num
          do rz=1,grid_num
            arr_in(rx,ry,rz)=cmplx(0,0)
            vxarr_in(rx,ry,rz)=cmplx(0,0)
            vyarr_in(rx,ry,rz)=cmplx(0,0)
            vzarr_in(rx,ry,rz)=cmplx(0,0)
           end do
        end do
      end do

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

!read data from data file and compute pdf in CIC
      open(unit=10,file=Trim(datafile_xyz),status='old')
      do i=1, skip_line
          read(10,*)
      end do
      count = 1
  901 read(10,*,end=902) dataxyz(1,count), dataxyz(2,count),
     &    dataxyz(3,count), dataxyz(4,count)
      if(count .ge. max_data) then
        write(*,*) "need to adjust max_data for simulation"
        stop
      end if
! modify boxsize
      if(dataxyz(1,count) .lt. boxsize .and.
     &    dataxyz(2,count) .lt. boxsize .and.
     &    dataxyz(3,count) .lt. boxsize) then
         count = count + 1
      end if
      goto 901
  902 close(10)
!      num_data=count-1
      write(*,*) "there are", count-1, "data"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

!
      do i=-max_vx_pdf, max_vx_pdf
         sim_vx_pdf(i) = 0
      end do

      do i = 1, count-1
           x=dataxyz(1,i)
           y=dataxyz(2,i)
           z=dataxyz(3,i)
           vx=dataxyz(4,i)

           if (vx .ge. 0) then
             sim_vx_pdf(int(vx)) = sim_vx_pdf(int(vx))+1
           else
             sim_vx_pdf(int(vx)-1) = sim_vx_pdf(int(vx)-1)+1
           end if

           xi = int(x/cell_size)+1
           yi = int(y/cell_size)+1
           zi = int(z/cell_size)+1

               if(xi .eq. grid_pdf+1) xi = xi-1
               if(yi .eq. grid_pdf+1) yi = yi-1
               if(zi .eq. grid_pdf+1) zi = zi-1

               xii = xi + 1
               if(xii.eq.grid_pdf+1) xii=1
               yii = yi + 1
               if(yii.eq.grid_pdf+1) yii=1
               zii = zi + 1
               if(zii.eq.grid_pdf+1) zii=1


! --- original cic assignment
               vyarr_in(xi,yi,zi) = cmplx(0,
     &    AIMAG(vyarr_in(xi,yi,zi))+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi)))
               vyarr_in(xii,yi,zi) =  cmplx(0,
     &    AIMAG(vyarr_in(xii,yi,zi))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-z+cell_size*(zi)))
               vyarr_in(xi,yii,zi) =  cmplx(0,
     &     AIMAG(vyarr_in(xi,yii,zi))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi)))
               vyarr_in(xi,yi,zii) =  cmplx(0,
     &     AIMAG(vyarr_in(xi,yi,zii))+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1)))
               vyarr_in(xii,yii,zi) =  cmplx(0,
     &     AIMAG(vyarr_in(xii,yii,zi))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-z+cell_size*(zi)))
               vyarr_in(xi,yii,zii) =  cmplx(0,
     &     AIMAG(vyarr_in(xi,yii,zii))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1)))
               vyarr_in(xii,yi,zii) =  cmplx(0,
     &     AIMAG(vyarr_in(xii,yi,zii))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (z-cell_size*(zi-1)))
               vyarr_in(xii,yii,zii) = cmplx(0,
     &     AIMAG(vyarr_in(xii,yii,zii))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (z-cell_size*(zi-1)))


        end do

      write(*,*) "finish assigning data to pdf grid"
      write(*,*) "write vx pdf from simulation"
      open(15,file=Trim(sim_pdf_file)//'.vx')
      do i = -max_vx_pdf, max_vx_pdf
        write(15,*) i, sim_vx_pdf(i)
      end do
      close(15)

      time2 =  time()
      write(*,*) "time in sec:", time2-time1
!test
      write(*,*) "going to convert CIC density to integer by poisson"


      count=0
      count2=0
      count_real=0

C$omp parallel private(thread_id,temp,temp_int,i,j,k)
C$omp+reduction(+:sim_pdf2)
      thread_id = omp_get_thread_num()+1
      write(*,*) "thread_id=", thread_id
C$omp do
      do i=1, grid_pdf
         do j=1, grid_pdf
            do k=1, grid_pdf
              temp = AIMAG(vyarr_in(i,j,k))/cell_size**3
!               count_real = count_real + temp
               call  poisson_rand(temp, iseedxx(thread_id), temp_int,
     &                                   thread_id)

              vyarr_in(i,j,k) = cmplx(0,temp_int)
!               count2 = count2+ temp_int
               if(temp_int .gt. max_num) then
               write(*,*)  "increase max_num!", temp_int
               stop
              end if
              sim_pdf2(temp_int) =sim_pdf2(temp_int) +1
            end do
         end do
      end do
C$omp end do
C$omp end parallel

      write(*,*) "finish computing pdf for simulation"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1


      Write(*,*) "PDF from cic assignment (sim)"
      open(10, file=Trim(sim_pdf_file))
      do i = 0, max_num
         write(10,*) i, sim_pdf2(i)
         if(sim_pdf2(i) .gt. 0) write(*,*) i, sim_pdf2(i)
      end do
      close(10)

!!test albert
      write(*,*)" modifying pdf from simulation"
      A_pdf = 1.*sim_pdf2(2)/sim_pdf2(1)
      B_pdf = 1.*sim_pdf2(1)/A_pdf
      write(*,*) "A, B of PDF are", A_pdf, B_pdf
      sum_pdf = 0
      sum_pdf_A = 0
      do i=1, 20
        sum_pdf_A = sum_pdf_A + A_pdf**i*i
      end do
      write(*,*) "sum of iA^i, sum of iA^i*B",sum_pdf_A, sum_pdf_A*B_pdf
! modify A_pdf
      sum_pdf = sum_pdf_A*B_pdf
      A_pdf = A_pdf - modify_pdf
      sum_pdf_A = 0
      do i=1, 20
        sum_pdf_A = sum_pdf_A + A_pdf**i*i
      end do
      B_pdf = sum_pdf/sum_pdf_A
      write(*,*) "new sum of iA^i, A, B", sum_pdf_A, A_pdf, B_pdf

      sum_pdf=0
      temp_int = 2 !initial with any number > 1
      i=1
      do while (temp_int .ge. 1.)
        temp_int = int(B_pdf*A_pdf**i)
        write(*,*) i, temp_int
        sim_pdf2(i) = temp_int
        sum_pdf = sum_pdf + temp_int*i
        i = i+1
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



!-- remove temperary
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
!test -- albert
      write(*,*) "modify input pk here"      
      do i=1,k_num
!          pk_array(i) = pk_arrayo(i)
        if(antidamping .lt. 1) then
         if(k_array(i).lt. 0.5)
     &    pk_arrayo(i) = (pk_arrayo(i)-pk_nw_array(i))*
     &          exp(k_array(i)**2/antidamping) + pk_nw_array(i)
         pk_arrayo(i) = pk_arrayo(i)*(1.+(k_array(i)*modify_pk)**(-1))
        end if
      end do

! density on grid array
! to initial the pk function
      call spline(k_array,pk_arrayo,k_num,1.0e30,1.0e30,y2o)

!   pk(0) = 0
      vxarr_in(1,1,1)=cmplx(0,0)
      vyarr_in(1,1,1)=cmplx(0,0)
      vzarr_in(1,1,1)=cmplx(0,0)
!initial
C$omp parallel do private(rx,ry,rz)
      do rx = 1, grid_num
        do ry = 1, grid_num
          do rz = 1, grid_num
        vyarr_in(rx,ry,rz) = cmplx(0,0)
        vxarr_in(rx,ry,rz) = cmplx(0,0)
          end do
        end do
      end do

!generate white noise on arr_in or read the white noise file
      if(use_whitenoise_file .eq. .false.) then

        thread_id = 1
        do rx = 1, grid_num
          do ry = 1, grid_num
            do rz = 1, grid_num,2  !might use every 2 since gaussion generator provide 2
              call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)
              arr_in(rx,ry,rz) = cmplx(ak,0)
              arr_in(rx,ry,rz+1) = cmplx(bk,0)
            end do
          end do
        end do
        write(*,*) "finish generating white noise"
        time2 =  time()
        write(*,*) "time in sec:", time2-time1

      else
        write(*,*) "reading white noise file"
        open(10,file=Trim(whitenoise_file))
!        read (10) np1,np2,np3, temp_int
!        print*, 'np1,np2,np3,iseed'
!        print*, np1,np2,np3,temp_int
!        if (np1 .ne. grid_num) then
!          write(*,*)"the grid size does not match white noise file:",np1 
!          stop
!        end if
!        allocate(f(np1,np2,np3))

        do i3=1,grid_num
          do i2=1,grid_num
         read(10,*) (temp_array(i1,i2,i3),i1=1,grid_num)
          end do
        enddo
        close(10)
        write(*,*) "finish reading white noise file"
        time2 =  time()
        write(*,*) "time in sec:", time2-time1

!move white noise to arr_in
C$omp parallel do private(rx,ry,rz)
        do rx = 1, grid_num
          do ry = 1, grid_num
            do rz = 1, grid_num
              arr_in(rx,ry,rz) = cmplx(temp_array(rx,ry,rz),0)
            end do
          end do
        end do

      end if

      
!!for openmp
      call dfftw_init_threads(iret)
      if(iret .eq. 0) then
          write(*,*) "something wrong with using fftw and openmp"
          stop
      end if
      write(*,*) "number of threads=",omp_get_max_threads()

      call dfftw_plan_with_nthreads(omp_get_max_threads())

      call dfftw_plan_dft_3d(plan, grid_num, grid_num, grid_num,
     &  arr_in, arr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(plan,arr_in,arr_in)
      call dfftw_destroy_plan(plan)

     
      temp = sqrt(1.0*grid_num**3/2) !to normalize white noise with given grid
!      vxarr_in(1,1,1) = cmplx(Real(arr_in(1,1,1))/temp,0)
C$omp parallel private(nx,ny,nz,kx,ky,kz,k_mag,pk_outo,tempo,ak,bk,
C$omp+nx_conj,ny_conj,nz_conj,thread_id) reduction(+:sigma_G2)
      thread_id = omp_get_thread_num()+1
      sigma_G2 = 0
C$omp do
      do nx = 0, grid_num-1
        do ny = 0, grid_num-1
          do nz = 0, grid_num-1
            if((nx .ne. 0) .and. (ny .ne. 0) .and. (nz .ne. 0)) then
              kx = 2.0*pi/boxsize*nx
              ky = 2.0*pi/boxsize*ny
              kz = 2.0*pi/boxsize*nz
              if (nx .gt. grid_num/2) kx = 2.0*pi/boxsize*(nx-grid_num)
              if (ny .gt. grid_num/2) ky = 2.0*pi/boxsize*(ny-grid_num)
              if (nz .gt. grid_num/2) kz = 2.0*pi/boxsize*(nz-grid_num)
              k_mag = sqrt(kx**2 + ky**2 + kz**2)
              call splint(k_array,pk_arrayo,y2o,k_num,k_mag,pk_outo)
              tempo = sqrt(pk_outo/boxsize**3/2)
              sigma_G2 = sigma_G2 + pk_outo*2 !*2 for nz < 0
!              call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)

 !following is for velocity field, one need to figure out issue of conj if use white noise from sim
              ak = Real(arr_in(nx+1,ny+1,nz+1))/temp
              bk = AIMAG(arr_in(nx+1,ny+1,nz+1))/temp
              vxarr_in(nx+1,ny+1,nz+1)=
     &              cmplx(tempo*bk*kx/k_mag**2,-tempo*ak*kx/k_mag**2)
              vyarr_in(nx+1,ny+1,nz+1)=
     &              cmplx(tempo*bk*ky/k_mag**2,-tempo*ak*ky/k_mag**2)
              vzarr_in(nx+1,ny+1,nz+1)=
     &              cmplx(tempo*bk*kz/k_mag**2,-tempo*ak*kz/k_mag**2)

!              nx_conj = grid_num - nx
!              ny_conj = grid_num - ny
!              nz_conj = grid_num - nz
!              if (nx .eq. 0) nx_conj = 0
!              if (ny .eq. 0) ny_conj = 0
!              if (nz .eq. 0) nz_conj = 0
!              vxarr_in(nx_conj+1,ny_conj+1,nz_conj+1)=
!     &              cmplx(tempo*bk*kx/k_mag**2,tempo*ak*kx/k_mag**2)
!              vyarr_in(nx_conj+1,ny_conj+1,nz_conj+1)=
!     &              cmplx(tempo*bk*ky/k_mag**2,tempo*ak*ky/k_mag**2)
!              vzarr_in(nx_conj+1,ny_conj+1,nz_conj+1)=
!     &              cmplx(tempo*bk*kz/k_mag**2,tempo*ak*kz/k_mag**2)
!
            end if
          end do !kz
        end do !ky
      end do !kx
C$omp end do
C$omp end parallel

      sigma_G2 = sigma_G2/boxsize**3
      write(*,*) "finish assigning pk array"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1


      write(*,*) "start FFT"
!!for openmp
!      call dfftw_init_threads(iret)
!      if(iret .eq. 0) then
!          write(*,*) "something wrong with using fftw and openmp"
!          stop
!      end if
!      write(*,*) "number of threads=",omp_get_max_threads()
!      call dfftw_plan_with_nthreads(omp_get_max_threads())

!      FFTW3
      call dfftw_plan_dft_3d(planx, grid_num, grid_num, grid_num,
     &  vxarr_in, vxarr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(planx,vxarr_in,vxarr_in)
      call dfftw_destroy_plan(planx)

      time2 =  time()
      write(*,*) "1st FFT, time in sec:", time2-time1

!      FFTW3
      call dfftw_plan_dft_3d(plany, grid_num, grid_num, grid_num,
     &  vyarr_in, vyarr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(plany,vyarr_in,vyarr_in)
      call dfftw_destroy_plan(plany)

      time2 =  time()
      write(*,*) "2nd FFT, time in sec:", time2-time1

!      FFTW3
      call dfftw_plan_dft_3d(planz, grid_num, grid_num, grid_num,
     &  vzarr_in, vzarr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(planz,vzarr_in,vzarr_in)
      call dfftw_destroy_plan(planz)

      time2 =  time()
      write(*,*) "3rd FFT, time in sec:", time2-time1

!--- displacement field using CIC ----------------------- raw zeldovich particles
! using image part of array for CIC counting
      write(*,*) "count initial CIC for za particles"
      open(17,file=Trim(output_file)//'.raw')

C$omp parallel do private(rx,ry,rz)
      do rx = 1,grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
              vxarr_in(rx,ry,rz) = cmplx(real(vxarr_in(rx,ry,rz)),0)
              vyarr_in(rx,ry,rz) = cmplx(real(vyarr_in(rx,ry,rz)),0)
              vzarr_in(rx,ry,rz) = cmplx(real(vzarr_in(rx,ry,rz)),0)
            end do
         end do
      end do

!test -- albert
!      do i = -max_vx_pdf ,max_vx_pdf
!         za_vx_pdf(i) = 0
!      end do
!
!      do rx = 1,grid_pdf
!         do ry = 1, grid_pdf
!            do rz = 1, grid_pdf
!              temp = real(vzarr_in(rx,ry,rz))
!              za_vx_pdf(int(temp*100)) = za_vx_pdf(int(temp*100))+1
!            end do
!         end do
!      end do
!               
!      open(11,file="../output/test_displace_z.dat")
!      do i= -max_vx_pdf, max_vx_pdf
!         write(11,*) i, za_vx_pdf(i)
!      end do
!      close(11)
!      stop
!
!assign ZA particles to grid using CIC
C$omp parallel do private(rx,ry,rz,xx,yy,zz,ix,iy,iz,
C$omp+ x,y,z,xi,yi,zi,xii,yii,zii,i,rxx,ryy,rzz,rxi,ryi,rzi)
      do i=0, total_num_thread-1
         rxx = i/(y_div_omp*z_div_omp)
         ryy = mod(i,y_div_omp*z_div_omp)/z_div_omp
         rzz = mod(i,z_div_omp)

        do rxi = 1, grid_num/x_div_omp
          do ryi = 1, grid_num/y_div_omp
            do rzi = 1, grid_num/z_div_omp
                
              rx = rxx*grid_num/x_div_omp + rxi
              ry = ryy*grid_num/y_div_omp + ryi
              rz = rzz*grid_num/z_div_omp + rzi

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
!
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


            end do
          end do
        end do
      end do !i

! compute average vx on grid
!      do i=-max_vx_pdf, max_vx_pdf
!        za_vx_pdf(i) = 0
!      end do

C$omp parallel do private(rx,ry,rz)
      do rx = 1, grid_num
        do ry = 1, grid_num
          do rz = 1, grid_num
        vzarr_in(rx,ry,rz) = cmplx(0,aimag(vzarr_in(i,j,k)))
!        temp_array2(rx,ry,rz) = 0
          end do
        end do
      end do

C$omp parallel private(thread_id,temp,temp_int,i,j,k)
C$omp+reduction(+:za_vx_pdf)
      thread_id = omp_get_thread_num()+1
      write(*,*) "thread_id=", thread_id
C$omp do
      do i=1, grid_pdf
         do j=1, grid_pdf
            do k=1, grid_pdf
              if(AIMAG(vxarr_in(i,j,k)).gt.0.00001) then
                temp =
     &    AIMAG(vyarr_in(i,j,k))/AIMAG(vxarr_in(i,j,k))
                temp_int = int(temp*100)
                vzarr_in(i,j,k)=cmplx(temp,aimag(vzarr_in(i,j,k)))
          temp_array2(i,j,k) = temp_array2(i,j,k)/AIMAG(vxarr_in(i,j,k))
          temp_array3(i,j,k) = temp_array3(i,j,k)/AIMAG(vxarr_in(i,j,k))

!                temp_array2(i,j,k) = temp_int
              else
!          write(*,*) "density is too small (za)", aimag(vxarr_in(i,j,k))
                temp_int=0
                vzarr_in(i,j,k)=cmplx(0,aimag(vzarr_in(i,j,k)))
                temp_array2(i,j,k) = 0
                temp_array3(i,j,k) = 0
              end if
              if(abs(temp_int) .lt. max_vx_pdf) then
               za_vx_pdf(temp_int) =za_vx_pdf(temp_int)+1
              end if
            end do
         end do
      end do
C$omp end do
C$omp end parallel

!      open(10, file="../output/test_za_vx_pdf.dat")
!      do i=-max_vx_pdf, max_vx_pdf
!        write(10,*) i, za_vx_pdf(i)
!      end do
!      close(10)
!      write(*,*) "finish vx pdf from za particles"
!      stop
! albert albert -- no deviation here
!      temp = 0
!      temp2 =0
!      do i=-max_vx_pdf, max_vx_pdf
!        temp = temp + i*za_vx_pdf(i)
!        temp2 = temp2 + za_vx_pdf(i)
!      end do
!      temp = temp/temp2
!      do i=1, grid_pdf
!         do j=1, grid_pdf
!            do k=1, grid_pdf
!               temp = temp+ temp_array2(i,j,k)/cell_size**3
!            end do
!         end do
!      end do
!      write(*,*) "test vzarr", temp
!      stop
!
      close(17)
      write(*,*) "finish assigning data to pdf grid & generate raw data"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

      write(*,*) "going to convert CIC density to integer by poisson"
! 500  count=0
!      f_exp = f_exp*0.5
!      write(*,*) "reduce normalization factor of za density:",
!     &           " new f_exp=",f_exp 

!! testing normalization for za density
!! cannot openmp because might jump out the loop (goto 500)
!      thread_id = 1
!      do i=1, grid_pdf
!         do j=1, grid_pdf
!            do k=1, grid_pdf
!              temp = AIMAG(vxarr_in(i,j,k))/cell_size**3 
!              call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)
!! scatter model 1: m1
!              ak = ak*scatter*temp**scatter2
!              if(ak .lt. 0) ak = exp(ak)-1
!               temp = temp*(1.+ak)*f_exp
!! scatter model 2: m2
!!              ak=ak*scatter*temp**scatter2
!!              if(ak .lt. 0) ak = exp(ak)-1
!!          temp = temp*(1.+ak)*f_exp*(1.-exp(-(temp/scatter3)**scatter4))
!
!
!              call  poisson_rand(temp, iseedxx(thread_id), 
!     &                   temp_int,thread_id)
!!              vzarr_in(i,j,k) = cmplx(real(vzarr_in(i,j,k)),temp_int)
!!              count = count+ temp_int
!              if(temp_int .gt. max_num) then
!               write(*,*)  temp_int,"is larger than max_num"
!               write(*,*) ak, temp, AIMAG(vxarr_in(i,j,k))/cell_size**3
!               goto 500
!              end if
!
!!             pdf(temp_int) =    
!!     &         pdf(temp_int) +1
!
!            end do
!         end do
!      end do

C$omp parallel private(i,j,k,temp,ak,bk,
C$omp+ thread_id)
      thread_id = omp_get_thread_num()+1
C$omp do
      do i=1, grid_pdf
         do j=1, grid_pdf
            do k=1, grid_pdf
              temp_array(i,j,k) = AIMAG(vxarr_in(i,j,k))/cell_size**3 
! cut density
              if (temp_array(i,j,k) .le. density_cut) then
                 temp_array(i,j,k) = 0
                 ak = 0
              else
                 call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)
!test !need to remove after testing !!!!!!!!!!!! albert
           if (temp_array(i,j,k) .ge. density_sat) then
        temp_array(i,j,k) = density_sat
           end if

              end if
! scatter model 1: m1
              ak = ak*scatter*temp_array(i,j,k)**scatter2
! scatter model 2: m2
!              ak = ak*scatter + bk*scatter2*temp**2
              if(ak .lt. 0) ak = exp(ak)-1
              temp_array(i,j,k)=temp_array(i,j,k)*(1.+ak)
            end do
         end do
      end do
!$omp end do
!$omp end parallel

      write(*,*) "going to compute pdf from normalized za field"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

      temp =  maxval(temp_array) 
      norm1 = max_num/(temp + 3.*sqrt(temp))
      count = 0
C$omp parallel private(i,j,k,temp_int,
C$omp+ thread_id) reduction(+:pdf)
      thread_id = omp_get_thread_num()+1
      do i=0, max_num
         pdf(i) = 0
      end do
C$omp do
      do i=1, grid_pdf
         do j=1, grid_pdf
            do k=1, grid_pdf

              temp_array(i,j,k) = temp_array(i,j,k)*norm1
              call  poisson_rand(temp_array(i,j,k), iseedxx(thread_id), 
     &                                temp_int, thread_id)
              vzarr_in(i,j,k) = cmplx(real(vzarr_in(i,j,k)),temp_int)
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
      write(*,*) "total raw za mock from CIC", count
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

!!
!
!      open(10, file=Trim(output_file))
      do i = 0,max_num
!         write(10,*) i, pdf(i)
         if(pdf(i) .gt. 0) write(*,*) i, pdf(i)
      end do
      close(10)

      write(*,*) "read pdf of simulation"
      open(10,file=Trim(sim_pdf_file)//'.new',status='old')
      sim_pdf_max = 0
      do i=0, max_num
        read(10,*,end=201) temp, sim_pdf(i)
!        write(*,*) temp, sim_pdf(i)
        sim_pdf_max =sim_pdf_max+1
      end do
 201  close(10)
      write(*,*) "sim_pdf_max=",sim_pdf_max-1


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

      write(*,*) "start producing mock using cic"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1
      open(17,file=Trim(output_file))
      do i=1, total_num_thread
         data_index(i) = 1          
      end do


C$omp parallel private(rx,ry,rz,np,np_ratio,temp,
C$omp+ thread_id,ak,bk) 
      thread_id = omp_get_thread_num()+1
C$omp do
      do rx = 1, grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
              np = pdf_convert(int(AIMAG(vzarr_in(rx,ry,rz))))
              np_ratio = 
     &             pdf_convert_ratio(int(AIMAG(vzarr_in(rx,ry,rz))))
              if(ran3_omp(iseedxx(thread_id),thread_id) .lt. np_ratio)
     &       np = np - 1
              vyarr_in(rx,ry,rz) = cmplx(real(vyarr_in(rx,ry,rz)),1.*np)

               do i =1, np
                  temp = ran3_omp(iseedxx(thread_id),thread_id) - 0.5
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

                  temp = ran3_omp(iseedxx(thread_id),thread_id) - 0.5
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

                  temp = ran3_omp(iseedxx(thread_id),thread_id) - 0.5
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
!add velocity - only z direction is correct so far -- testing
                 call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)

                  temp = real(vzarr_in(rx,ry,rz))  *zdist_rate

                  dataxyz_omp(6,data_index(thread_id),thread_id) = 
     &            temp + zdist_fog*ak!*temp/abs(temp)*abs(temp)**(0.5)

                  temp =  temp_array2(rx,ry,rz)  *zdist_rate
                  dataxyz_omp(4,data_index(thread_id),thread_id) = 
     &            temp + zdist_fog*bk!*temp/abs(temp)*abs(temp)**(0.5)

                 call gauss_ran(iseedxx(thread_id),ak,bk,thread_id)

                  temp = temp_array3(rx,ry,rz)  *zdist_rate
                  dataxyz_omp(5,data_index(thread_id),thread_id) = 
     &            temp + zdist_fog*ak!*temp/abs(temp)*abs(temp)**(0.5)

                  data_index(thread_id) = data_index(thread_id) + 1
               end do

            end do
         end do
      end do
!$omp end do
      data_index(thread_id) = data_index(thread_id) -1
!$omp end parallel


      mock_data_num = 0
      do i=1, total_num_thread
        mock_data_num = mock_data_num + data_index(i)
      end do  

!$omp parallel do private(i,j,k) reduction(+:pdf2)
      do i=1, grid_pdf
         do j=1, grid_pdf
            do k=1, grid_pdf
              if(AIMAG(vyarr_in(i,j,k)) .gt. max_num) then
               write(*,*)  "increase max_num! ijij" 
               stop
              end if

             pdf2(int(AIMAG(vyarr_in(i,j,k)))) =
     &         pdf2(int(AIMAG(vyarr_in(i,j,k)))) +1

            end do
         end do
      end do

      temp = 0
!      open(10, file=Trim(output_file)//'.pdf.cic')
      do i = 0, max_num
!         write(10,*) i, pdf2(i)
         if(pdf2(i) .gt. 0) write(*,*)"converted pdf:",i, pdf2(i)
         temp = temp + pdf2(i)*i
      end do
!      close(10)
      write(*,*) "total number from converted pdf", temp
      write(*,*) "total mock produced", mock_data_num
      
! write mock data
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
      write(*,*) "finish saving mock, now working on cic pdf"
      time2 =  time()
      write(*,*) "time in sec:", time2-time1

! compute pdf of mock cata of cic
!$omp parallel do private(rx,ry,rz)
      do rx = 1,grid_pdf
         do ry = 1, grid_pdf
            do rz = 1, grid_pdf
!              count_in_cell(rx,ry,rz) = 0
            vzarr_in(rx,ry,rz)=cmplx(real(vzarr_in(rx,ry,rz)),0)
            vwarr_in(rx,ry,rz)=cmplx(0,0)
            end do
         end do
      end do

      do i=0, max_num
         pdf2(i) = 0
      end do

      do i=-max_vx_pdf, max_vx_pdf
         za_vx_pdf(i)=0
      end do 
!might be openmp while generating the data
!BUT, it is for post-process so that the time used is not so important
      do i = 1, total_num_thread
        do j = 1, data_index(i)
           x=dataxyz_omp(1,j,i)
           y=dataxyz_omp(2,j,i)
           z=dataxyz_omp(3,j,i)
           vz=dataxyz_omp(6,j,i)

           zz=z+vz*inverse_aH !(1+z)/Hubble
           if(zz.gt.boxsize) zz = zz - boxsize
           if(zz.lt.0.) zz = zz + boxsize;

           if(vz.ge.0) then
             za_vx_pdf(int(vz)) = za_vx_pdf(int(vz))+1
           else
             za_vx_pdf(int(vz)-1) = za_vx_pdf(int(vz)-1)+1
           end if

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
     &   (-z+cell_size*(zi)))
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

             vwarr_in(xi,yi,zzi) =cmplx(real(vwarr_in(xi,yi,zzi)),
     &    AIMAG(vwarr_in(xi,yi,zzi))+ 
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (-zz+cell_size*(zzi)))
             vwarr_in(xii,yi,zzi) =cmplx(real(vwarr_in(xii,yi,zzi)),
     &   AIMAG(vwarr_in(xii,yi,zzi))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (-zz+cell_size*(zzi)))
             vwarr_in(xi,yii,zzi) = cmplx(real(vwarr_in(xi,yii,zzi)),
     &    AIMAG(vwarr_in(xi,yii,zzi))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (-zz+cell_size*(zzi)))
             vwarr_in(xi,yi,zzii) = cmplx(real(vwarr_in(xi,yi,zzii)),
     &    AIMAG(vwarr_in(xi,yi,zzii))+
     &   (-x+cell_size*(xi))*(-y+cell_size*(yi))*
     &   (zz-cell_size*(zzi-1)))
             vwarr_in(xii,yii,zzi) = cmplx(real(vwarr_in(xii,yii,zzi)),
     &    AIMAG(vwarr_in(xii,yii,zzi))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (-zz+cell_size*(zzi)))
             vwarr_in(xi,yii,zzii) = cmplx(real(vwarr_in(xi,yii,zzii)),
     &    AIMAG(vwarr_in(xi,yii,zzii))+
     &   (-x+cell_size*(xi))*(y-cell_size*(yi-1))*
     &   (zz-cell_size*(zzi-1)))
             vwarr_in(xii,yi,zzii) = cmplx(real(vwarr_in(xii,yi,zzii)),
     &    AIMAG(vwarr_in(xii,yi,zzii))+
     &   (x-cell_size*(xi-1))*(-y+cell_size*(yi))*
     &   (zz-cell_size*(zzi-1)))
             vwarr_in(xii,yii,zzii) =cmplx(real(vwarr_in(xii,yii,zzii)),
     &    AIMAG(vwarr_in(xii,yii,zzii))+
     &   (x-cell_size*(xi-1))*(y-cell_size*(yi-1))*
     &   (zz-cell_size*(zzi-1)))

        end do
      end do

      write(*,*) "finish assigning mock to pdf grid"
      write(*,*) "write vx pdf from data"
      temp=0
      temp2=0
      open(15,file=Trim(output_file)//'.pdf.vx')
      do i = -max_vx_pdf, max_vx_pdf
        write(15,*) i, za_vx_pdf(i)
      end do
      close(15)

      time2 =  time()
      write(*,*) "time in sec:", time2-time1

      write(*,*) "going to convert CIC density to integer by poisson"

!test
      temp = 0
      temp2 = 0
      do i=1, grid_pdf
         do j=1, grid_pdf
            do k=1, grid_pdf
              temp = temp+AIMAG(vzarr_in(i,j,k))/cell_size**3
              temp2 = temp2+AIMAG(vwarr_in(i,j,k))/cell_size**3
            end do
         end do
      end do
      write(*,*) "test: total number from cic of mock", temp,temp2
!      stop

C$omp parallel private(i,j,k,temp,temp_int,
C$omp+ thread_id) reduction(+:pdf2,temp2)
      thread_id = omp_get_thread_num()+1
      temp2 = 0
C$omp do
      do i=1, grid_pdf
         do j=1, grid_pdf
            do k=1, grid_pdf
              temp = AIMAG(vzarr_in(i,j,k))/cell_size**3
              temp2 = temp2 + temp
              call  poisson_rand(temp, iseedxx(thread_id),
     &           temp_int,thread_id)
              if(temp_int .gt. max_num) then
               write(*,*)  "increase max_num!", temp_int
               stop
              end if

             pdf2(temp_int) =
     &         pdf2(temp_int) +1

            end do
         end do
      end do
!$omp end do
!$omp end parallel

      temp = 0
      open(10, file=Trim(output_file)//'.pdf.cic')
      do i = 0, max_num
         write(10,*) i, pdf2(i)
        if(pdf2(i) .gt. 0) write(*,*) i, pdf2(i)
        temp = temp + pdf2(i)*i
      end do
      close(10)
      write(*,*) "total number from mock pdf", temp
      write(*,*) "total temp",temp2

      time2 =  time()
      write(*,*) "Done computing pdf of mock, time in sec:", time2-time1

!--- compute pk ---

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
            vzarr_in(i,j,k)=
     &  cmplx(AIMAG(vzarr_in(i,j,k))/cell_size**3-
     &                  (1.*mock_data_num)/(1.*grid_pdf)**3,0)
          end do
        end do
      end do

      write(*,*) "starting FFT for pk"
!      FFTW3
      call dfftw_plan_dft_3d(planz, grid_num, grid_num, grid_num,
     &  vzarr_in, vzarr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(planz,vzarr_in,vzarr_in)
      call dfftw_destroy_plan(planz)

      time2 =  time()
      write(*,*) "finish FFT, time in sec:", time2-time1

      box_size = boxsize
      num_data = mock_data_num

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
        pk(i) = pk(i)/k_num_sep(i)/num_data/num_data*(1.*box_size)**3*2
        pk_quad(i) = pk_quad(i)/
     &          k_num_sep(i)/num_data/num_data*(1.*box_size)**3*2

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
            vwarr_in(i,j,k)=
     &  cmplx(AIMAG(vwarr_in(i,j,k))/cell_size**3-
     &                  (1.*mock_data_num)/(1.*grid_pdf)**3,0)
          end do
        end do
      end do

      write(*,*) "starting FFT for pk"
!      FFTW3
      call dfftw_plan_dft_3d(planw, grid_num, grid_num, grid_num,
     &  vwarr_in, vwarr_in, FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute_dft(planw,vwarr_in,vwarr_in)
      call dfftw_destroy_plan(planw)

      time2 =  time()
      write(*,*) "finish FFT, time in sec:", time2-time1

      box_size = boxsize
      num_data = mock_data_num

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
         pk_temp =REAL(vwarr_in(i,j,k))**2+AIMAG(vwarr_in(i,j,k))**2
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
         pk_temp =REAL(vwarr_in(i,j,k))**2+AIMAG(vwarr_in(i,j,k))**2
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
        pk(i) = pk(i)/k_num_sep(i)/num_data/num_data*(1.*box_size)**3*2
        pk_quad(i) = pk_quad(i)/
     &          k_num_sep(i)/num_data/num_data*(1.*box_size)**3*2

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
!      character*200,intent(in) :: randfile
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
      
