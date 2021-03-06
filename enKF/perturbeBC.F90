!
! Copyright (C) 2020, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
!
!------------------------------------------------------------------------------
! Creates an ensemble of files with perturbed BC from an initial BC file.
! Can handle both lateral and surface BC.
!------------------------------------------------------------------------------
program perturbeBC

  use m_sample2D
  implicit none

  character(len=6) :: arg1
  character(len=3) :: arg2
  character(len=3) :: arg3
  character(len=3) :: arg4
  character(len=12) :: arg5
  character(len=12) :: arg6
  character(len=12) :: arg7
  character(len=12) :: arg8
  character(len=12) :: arg9
  character(len=80) :: arg10

  integer nrens, var_dim
  character(len=80) :: filein
  character(len=3) :: filety
  real var_std
  double precision mem_time
  double precision time_std	
  integer pert_type

  character(len=80) :: inbname

  ! FEM files
  integer iunit,iformat
  integer nvers           !version of file format
  integer np              !size of data (horizontal, nodes or elements)
  integer nvar            !number of variables to write
  integer ntype           !type of information contained
  integer datetime(2)     !date and time information
  integer ierr            !return error code
  real*4 regpar(7)        !regular array params
  integer lmax            !vertical values
  real*4,allocatable :: hlv(:)     !vertical structure
  character(len=50),allocatable :: vstring(:)
  integer,allocatable :: ilhkv(:)
  real*4,allocatable :: hd(:)
  real x0,y0,dx,dy,flag
  integer nx,ny
  integer nlvddi
  !---

  integer           :: nrec_file

  real              :: var
  double precision  :: tnew,told,dtime
  double precision  :: tramp,tfact,enstime
  logical	    :: blast
  real, parameter   :: dt = 60.   		!time step for rounding time in secs
  real, parameter   :: time_ramp = 43200.  	!extension of time ramp for temporal shift
  double precision, allocatable  :: tens(:)	!in seconds
  character(len=80) :: dstring
  real              :: var_min,var_max

  real, allocatable :: pvec(:),pvec1(:,:)
  real, allocatable :: amat(:,:,:)
  real, allocatable :: pmat(:,:,:,:)

  ! 2D perturbation parameters
  real rx,ry,theta
  integer fmult
  logical samp_fix,verbose

  ! geostrophic vars
  real               :: flat
  logical, parameter :: bpress = .true.

  real, allocatable :: var3d(:,:,:,:),var3d_ens(:,:,:,:)

  ! SHYFEM single precision variables
  real*4, allocatable :: femdata(:,:,:)

  integer n,i,ne,l,ix,iy
  integer fid

!********************************
! Get input from stdin
!********************************
  ! read input
  call get_command_argument(1, arg1)
  call get_command_argument(2, arg2)
  call get_command_argument(3, arg3)
  call get_command_argument(4, arg4)
  call get_command_argument(5, arg5)
  call get_command_argument(6, arg6)
  call get_command_argument(7, arg7)
  call get_command_argument(8, arg8)
  call get_command_argument(9, arg9)
  call get_command_argument(10, arg10)

  if (trim(arg10) .eq. '') then
      write(*,*) ''
      write(*,*) 'Usage:'
      write(*,*) ''
      write(*,*) 'perturbeBC [nrens] [file_type] [pert_type] [var_dim]' // &
                           ' [var_std] [var_min] [var_max] [mem_time] [time_std] [input]'
      write(*,*) ''
      write(*,*) '[nrens] is the n. of ens members, control included.'
      write(*,*) '[file_type] can be fem or ts (timeseries).'
      write(*,*) '[pert_type] type of perturbation, see below.'
      write(*,*) '[var_dim] spatial dimension of the variable (0->3).'
      write(*,*) '[var_std] standard deviation of the ensemble distribution.'
      write(*,*) '[var_min] minimum value for the variable (-999 to disable, not used for wind).'
      write(*,*) '[var_max] maximum value for the variable (-999 to disable, if wind this is wind speed).'
      write(*,*) '[mem_time] time correlation of the perturbations (red noise) in seconds.'
      write(*,*) '[time_std] standard deviation for perturbation in time [in seconds].'
      write(*,*) '[input] is the name of the input unperturbed BC file.'
      write(*,*) ''
      write(*,*) 'pert_type can be:'
      write(*,*) ''
      write(*,*) '1- Spatially constant perturbation of each variable (1D, 2D, 3D)'
      write(*,*) '2- 2D wind pseudo Gaussian perturbation of u and v, no press.'
      write(*,*) '3- 2D wind pseudo Gaussian perturbation of press, geostrophic u and v (press var_std).'
      write(*,*) '4- 2D wind-speed pseudo Gaussian perturbation, no press.'
      write(*,*) ''
      write(*,*) 'Use var_std = 0 to disable this type of perturbations'
      write(*,*) 'Use time_std = 0 to disable the time-shift perturbations'
      write(*,*) ''
      write(*,*) 'The time-shift perturbations shift the ens members in time, according to a Gaussian'
      write(*,*) 'distribution (time_std). A time-ramp of half a day is used in order to have all the'
      write(*,*) 'members at the initial times.'
      write(*,*) ''
      stop
  end if

  read(arg1,*) nrens
  filety = arg2
  read(arg3,*) pert_type
  read(arg4,*) var_dim
  read(arg5,*) var_std
  read(arg6,*) var_min
  read(arg7,*) var_max
  read(arg8,*) mem_time
  read(arg9,*) time_std
  filein = arg10

  flag = -999.

!********************************
! some checks on inputs
!********************************
  if (( nrens < 3 ) .or. (mod(nrens,2) == 0) .or. (nrens > 1000)) error stop 'perturbeBC: bad nrens'
  if (( var_dim > 3 ) .or. ( var_dim < 0 )) error stop 'perturbeBC: bad dimension'
  if (( trim(filety) /= 'ts' ) .and. ( trim(filety) /= 'fem' )) error stop 'perturbeBC: bad file_type'
  if (( trim(filety) == 'ts' ) .and. ( var_dim /= 0 )) error stop 'perturbeBC: file_type and dimension not compatible'
  if (( trim(filety) == 'fem' ) .and. ( var_dim < 1 )) error stop 'perturbeBC: file_type and dimension not compatible'
  if (mem_time <= 0) write(*,*) 'Warning: zero or negative mem_time. Setting to 0 (white noise).'
  if ( pert_type > 6 ) error stop 'perturbeBC: bad pert_type'
  if (( var_std < 0 ) .or. ( time_std < 0 )) error stop 'perturbeBC: bad var_std and time_std'
  if (( var_std == 0 ) .and. ( time_std == 0 )) error stop 'perturbeBC: bad var_std and time_std'

!********************************
! Compute time perturbations. tens(1) = 0.
!********************************
  tramp = 0.
  if (.not.allocated(tens)) allocate(tens(nrens))
  call perturbe_time(nrens,time_std,tens)

!********************************
! open and close to read informations
!********************************
  select case(trim(filety))
  case('ts')
     call read_ts(.true.,filein,nrec_file,var,tnew,dstring)
  case('fem')
     np = 0
     y0 = 0.
     x0 = 0.
     dx = 0.
     dy = 0.
     nx = 0
     ny = 0
     call fem_file_read_open(filein,np,iformat,iunit)
  end select

  iunit = 20
  call open_files(iunit,nrens,filein,filety,inbname)

  lmax = 1
  told = -999.

!********************************
! time loop start
!********************************
  n = 0
  do

    n = n + 1

!********************************
! read record
!********************************
    select case(trim(filety))
!----------------------------------------------
    case('ts')
!----------------------------------------------

	if (n > nrec_file) exit
        call read_ts(.false.,filein,nrec_file,var,tnew,dstring)

!----------------------------------------------
    case('fem')
!----------------------------------------------

	! read 1st header
	call fem_file_read_params(iformat,iunit,dtime &
               ,nvers,np,lmax,nvar,ntype,datetime,ierr)
        if( ierr .lt. 0 ) exit

	call dts_convert_to_atime(datetime,dtime,tnew)

	! read 2nd header
	if (.not.allocated(hlv)) allocate(hlv(lmax))
        if (.not.allocated(ilhkv)) allocate(ilhkv(np))
	if (.not.allocated(hd)) allocate(hd(np))
	if (.not.allocated(vstring)) allocate(vstring(nvar))

        nlvddi = lmax
        call fem_file_read_2header(iformat,iunit,ntype,lmax &
             ,hlv,regpar,ierr)
        nx = nint(regpar(1))
        ny = nint(regpar(2))
	x0 = regpar(3)
	y0 = regpar(4)
        dx = regpar(5)
        dy = regpar(6)
        flag = regpar(7)

	! read variables
	if (.not.allocated(femdata)) allocate(femdata(lmax,nx,ny))
	if (.not. allocated(var3d)) allocate(var3d(nvar,nx,ny,lmax))
	do i = 1,nvar

	   femdata = flag
           call fem_file_read_data(iformat,iunit,nvers,np,lmax, &
		   vstring(i),ilhkv,hd,nlvddi,femdata,ierr)
	   do l = 1,lmax
	    do iy = 1,ny
	    do ix = 1,nx
	       var3d(i,ix,iy,l) = femdata(l,ix,iy) 
	    end do
	    end do
	   end do

	   if (n==1) write(*,*) '  Reading: ',vstring(i)
        end do
	
	! check if it is the last record
        call fem_check_last(iformat,iunit,blast)

    end select

!----------------------------------------------
    !smooth ramp [0-1] for time perturbation
!----------------------------------------------
    if (n > 1 ) tramp = (tramp + (tnew-told))
    tfact = min(tramp/time_ramp,1.)
    

!********************************
! generate perturbations and write the files
!********************************

    if (.not.allocated(pvec)) allocate(pvec(nrens-1))
    select case(var_dim)
!----------------------------------------------
    case(0)	! 0D variable
!----------------------------------------------

        call perturbe_0d(nrens-1,pvec)
        call red_noise_0d(told,tnew,pvec,nrens-1,mem_time,1,1)
        call write_record_0d(iunit,told,nrens,dstring,var,var_std,var_min,var_max,pvec,flag)

!----------------------------------------------
    case(1)	! 1D variable
!----------------------------------------------

        print*, 'todo'
	stop
        !call perturbe_1d

!----------------------------------------------
    case(2)	! 2D variable
!----------------------------------------------

	if (trim(filety) /= 'fem') error stop 'Bad file format'
	if ( lmax /= 1 ) error stop 'Bad vertical dimensions'

        write(*,*) "**********************"
        write(*,*) "  Field n. ",n
        write(*,*) "**********************"

	select case(pert_type)
	!-----------------------
	case(1)
	!-----------------------

            if (var_std > 0.) write(*,*) '  Case 1: spatially constant perturbations'

            if (.not.allocated(pvec1)) allocate(pvec1(nvar,nrens-1))
	    do i=1,nvar
               call perturbe_0d(nrens-1,pvec)
               call red_noise_0d(told,tnew,pvec,nrens-1,mem_time,i,nvar)
               pvec1(i,:) = pvec
	    end do

	    do ne=1,nrens
	      ! create member
	      if (.not.allocated(var3d_ens)) allocate(var3d_ens(nvar,nx,ny,lmax))
	      var3d_ens = flag
              call make_member_2D1(nvar,nrens,nx,ny,lmax,ne,pvec1,var3d,var3d_ens,var_std,flag)

	      ! write file
              fid = iunit + 10 + ne
	      enstime = dtime + nint(tens(ne)*tfact/dt)*dt
              if (blast) enstime = max(enstime,dtime)

              call fem_file_write_header(iformat,fid,enstime,nvers,np,lmax &
                     ,nvar,ntype,nlvddi,hlv,datetime,regpar)
              do i = 1,nvar
	         femdata = flag
	         do iy = 1,ny
	          do ix = 1,nx
		     var = var3d_ens(i,ix,iy,1)
		     if (nint(var_min) /= nint(flag)) call var_limit_min(var,var_min,flag)
		     if (nint(var_max) /= nint(flag)) call var_limit_max(var,var_max,flag)
	             femdata(1,ix,iy) = var
	          end do
                 end do
                 call fem_file_write_data(iformat,fid,nvers,np,lmax &
                       ,vstring(i),ilhkv,hd,nlvddi,femdata)
              end do
	    end do

	!-----------------------
	case(2)
	!-----------------------

	    if (var_std > 0.) write(*,*) '  Case 2: 2D wind pseudo-Gaussian perturbations'

	    if (nvar /= 3) error stop 'Dimension error'
	    if (trim(vstring(1)) /= 'wind velocity - x') error stop 'Invalid variables'

	    ! set parameters
	    verbose = .false.
	    samp_fix = .true.
	    theta = 0.
	    fmult = 5
	    call set_decorrelation(nx,dx,rx)
	    call set_decorrelation(ny,dy,ry)
	    flat = y0 + (ny/2 * dy)

	    write(*,*) 'nx,ny,dx,dy,rx,ry,theta,samp_fix: ',nx,ny,dx,dy,rx,ry,theta,samp_fix

            ! Make the random fields
            allocate(amat(nx,ny,nrens-1),pmat(nvar-1,nx,ny,nrens-1))
	    do i=1,nvar-1
               call sample2D(amat,nx,ny,nrens-1,fmult,dx,dy,rx,ry,theta &
                    ,samp_fix,verbose)
               pmat(i,:,:,:) = amat
	    end do
            deallocate(amat)

            call red_noise_2d(nx,ny,nvar-1,told,tnew,pmat,nrens-1,mem_time)

	    do ne = 1,nrens

	      if (.not.allocated(var3d_ens)) allocate(var3d_ens(nvar,nx,ny,lmax))

              call make_member_2D2(nvar,nrens,nx,ny,lmax,ne,pmat,var3d,var3d_ens,var_std,flag)

	      call limit_wind(nvar,nx,ny,lmax,var3d_ens,var_max,flag)

	      ! write file
              fid = iunit + 10 + ne
	      enstime = dtime + nint(tens(ne)*tfact/dt)*dt	
              if (blast) enstime = max(enstime,dtime)

              call fem_file_write_header(iformat,fid,enstime,nvers,np,lmax &
                     ,nvar,ntype,nlvddi,hlv,datetime,regpar)
              do i = 1,nvar
	         femdata = flag
	         do iy = 1,ny
	          do ix = 1,nx
		     var = var3d_ens(i,ix,iy,1)
	             femdata(1,ix,iy) = var
	          end do
                 end do
                 call fem_file_write_data(iformat,fid,nvers,np,lmax &
                       ,vstring(i),ilhkv,hd,nlvddi,femdata)
              end do

	    end do
	    deallocate(pmat)

	!-----------------------
	case(3)
	!-----------------------

	    if (var_std > 0.) write(*,*) '  Case 2: 2D wind-press pseudo-Gaussian geostrophic perturbations'

	    if (nvar /= 3) error stop 'Dimension error'
	    if (trim(vstring(1)) /= 'wind velocity - x') error stop 'Invalid variables'

	    ! set parameters
	    verbose = .false.
	    samp_fix = .true.
	    theta = 0.
	    fmult = 5
	    call set_decorrelation(nx,dx,rx)
	    call set_decorrelation(ny,dy,ry)
	    flat = y0 + (ny/2 * dy)

	    write(*,*) 'nx,ny,dx,dy,rx,ry,theta,samp_fix: ',nx,ny,dx,dy,rx,ry,theta,samp_fix

            ! Make the random fields
            allocate(amat(nx,ny,nrens-1),pmat(1,nx,ny,nrens-1))
            call sample2D(amat,nx,ny,nrens-1,fmult,dx,dy,rx,ry,theta &
                    ,samp_fix,verbose)
            pmat(1,:,:,:) = amat
            deallocate(amat)

            call red_noise_2d(nx,ny,1,told,tnew,pmat,nrens-1,mem_time)

	    do ne = 1,nrens

	      if (.not.allocated(var3d_ens)) allocate(var3d_ens(nvar,nx,ny,lmax))

	      call make_geo_field(nvar,ne,nrens,nx,ny,lmax,dx,dy,pmat,var3d,var3d_ens,flag,var_std, &
                                flat,bpress)

	      call limit_wind(nvar,nx,ny,lmax,var3d_ens,var_max,flag)

	      ! write file
              fid = iunit + 10 + ne
	      enstime = dtime + nint(tens(ne)*tfact/dt)*dt	
              if (blast) enstime = max(enstime,dtime)

              call fem_file_write_header(iformat,fid,enstime,nvers,np,lmax &
                     ,nvar,ntype,nlvddi,hlv,datetime,regpar)
              do i = 1,nvar
	         femdata = flag
	         do iy = 1,ny
	          do ix = 1,nx
		     var = var3d_ens(i,ix,iy,1)
	             femdata(1,ix,iy) = var
	          end do
                 end do
                 call fem_file_write_data(iformat,fid,nvers,np,lmax &
                       ,vstring(i),ilhkv,hd,nlvddi,femdata)
              end do

	    end do
	    deallocate(pmat)

	!-----------------------
	case(4)
	!-----------------------

	    if (var_std > 0.) write(*,*) '  Case 2: 2D wind speed pseudo-Gaussian perturbations'

	    if (nvar /= 3) error stop 'Dimension error'
	    if (trim(vstring(1)) /= 'wind velocity - x') error stop 'Invalid variables'

	    ! set parameters
	    verbose = .false.
	    samp_fix = .true.
	    theta = 0.
	    fmult = 5
	    call set_decorrelation(nx,dx,rx)
	    call set_decorrelation(ny,dy,ry)
	    flat = y0 + (ny/2 * dy)

	    write(*,*) 'nx,ny,dx,dy,rx,ry,theta,samp_fix: ',nx,ny,dx,dy,rx,ry,theta,samp_fix

            ! Make the random fields
            allocate(amat(nx,ny,nrens-1),pmat(1,nx,ny,nrens-1))
            call sample2D(amat,nx,ny,nrens-1,fmult,dx,dy,rx,ry,theta &
                    ,samp_fix,verbose)
            pmat(1,:,:,:) = amat
            deallocate(amat)

            call red_noise_2d(nx,ny,1,told,tnew,pmat,nrens-1,mem_time)

	    do ne = 1,nrens

	      if (.not.allocated(var3d_ens)) allocate(var3d_ens(nvar,nx,ny,lmax))

	      call make_ws_pert(nvar,ne,nrens,nx,ny,lmax,pmat,var3d,var3d_ens,flag,var_std)

	      call limit_wind(nvar,nx,ny,lmax,var3d_ens,var_max,flag)

	      ! write file
              fid = iunit + 10 + ne
	      enstime = dtime + nint(tens(ne)*tfact/dt)*dt	
              if (blast) enstime = max(enstime,dtime)

              call fem_file_write_header(iformat,fid,enstime,nvers,np,lmax &
                     ,nvar,ntype,nlvddi,hlv,datetime,regpar)
              do i = 1,nvar
	         femdata = flag
	         do iy = 1,ny
	          do ix = 1,nx
		     var = var3d_ens(i,ix,iy,1)
	             femdata(1,ix,iy) = var
	          end do
                 end do
                 call fem_file_write_data(iformat,fid,nvers,np,lmax &
                       ,vstring(i),ilhkv,hd,nlvddi,femdata)
              end do

	    end do
	    deallocate(pmat)

	end select
!----------------------------------------------
    case(3)	! 3D variable
!----------------------------------------------

	if (trim(filety) /= 'fem') error stop 'Bad file format'
	if ( lmax < 2 ) error stop 'Bad vertical dimensions'

	select case(pert_type)
	case(1)

            if (n==1) write(*,*) '  Case 1: spatially constant perturbations'

            if (.not.allocated(pvec1)) allocate(pvec1(nvar,nrens-1))
	    do i=1,nvar
               call perturbe_0d(nrens-1,pvec)
               call red_noise_0d(told,tnew,pvec,nrens-1,mem_time,i,nvar)
               pvec1(i,:) = pvec
	    end do

	    do ne = 1,nrens
	      ! create member
	      if (.not.allocated(var3d_ens)) allocate(var3d_ens(nvar,nx,ny,lmax))
              call make_member_3D1(nvar,nrens,nx,ny,lmax,ne,pvec1,var3d,var3d_ens,var_std,flag)

	      ! write file
              fid = iunit + 10 + ne
	      enstime = dtime + nint(tens(ne)*tfact/dt)*dt	
              if (blast) enstime = max(enstime,dtime)

              call fem_file_write_header(iformat,fid,enstime,nvers,np,lmax &
                     ,nvar,ntype,nlvddi,hlv,datetime,regpar)
              do i = 1,nvar
	         femdata = flag
	         do l = 1,lmax
	            do iy = 1,ny
	              do ix = 1,nx
		         var = var3d_ens(i,ix,iy,l)
		         if (nint(var_min) /= nint(flag)) call var_limit_min(var,var_min,flag)
		         if (nint(var_max) /= nint(flag)) call var_limit_max(var,var_max,flag)
		         femdata(l,ix,iy) = var
	              end do
                    end do
		 end do
                 call fem_file_write_data(iformat,fid,nvers,np,lmax &
                       ,vstring(i),ilhkv,hd,nlvddi,femdata)
              end do
	      
	    end do
        end select

    end select

    told = tnew

!********************************
! time loop end
!********************************
  end do

!********************************
! close and rename files
!********************************
  call close_and_rename(filety,iunit,nrens,inbname)

end program perturbeBC

!********************************************************************************************
!********************************************************************************************
!********************************************************************************************
!********************************************************************************************

!-----------------------------------------------
  subroutine read_ts(linit,filein,kdim,v,vtime,dstring)
!-----------------------------------------------

  use iso8601
  implicit none

  logical, intent(in)          :: linit
  character(len=*),intent(in)  :: filein
  integer, intent(inout)       :: kdim
  real, intent(out)            :: v
  double precision, intent(out):: vtime
  character*80, intent(out) :: dstring
  integer ios
  integer ierr
  integer date, time
  integer k

  v = -999.
  vtime = -999.

  select case(linit)
  
     ! just check values and see the lenght
     !
     case (.true.)
          open(26,file=trim(filein), status = 'old', form = 'formatted', iostat = ios)
          if (ios /= 0) error stop 'read_ts: error opening file'

          k = 0
 90       read(26,*,end=100) dstring,v
          if (isnan(v)) error stop 'perturbeBC: input file contains nans'
          call string2date(trim(dstring),date,time,ierr)
          if (ierr /= 0) error stop "read_ts: error reading string"
          k = k + 1
          goto 90

 100      continue
          kdim = k
          rewind(26, iostat = ios)
          if (ios /= 0) error stop 'read_ts: error in file'
          return

     ! read the values
     !
     case (.false.)

          read(26,*,end=101) dstring,v
          call string2date(trim(dstring),date,time,ierr)
          if (ierr /= 0) error stop "read_ts: error reading string"
          call dts_to_abs_time(date,time,vtime)

  end select

  return

  101  close(26)

  end subroutine read_ts


!-----------------------------------------------
  subroutine open_files(iunit,nrens,fin,ftype,bname)
!-----------------------------------------------
  implicit none

  integer, intent(in) :: nrens,iunit
  character(len=80), intent(in) :: fin
  character(len=3), intent(in) :: ftype
  character(len=80),intent(out) :: bname

  integer n
  integer fid
  integer :: ppos
  character(len=80) :: fname
  character(len=3) :: nlab

  ! find basename
  ppos = scan(trim(fin),".", BACK= .true.)
  if ( ppos > 0 ) bname = fin(1:ppos-1)
  
  if (trim(ftype) == 'ts') then
     do n = 1,nrens
        fid = iunit + 10 + n
        write(nlab,'(i3.3)') n-1
        fname = trim(bname)//'_'//nlab//'.dat'
        open(fid,file=fname,status='unknown')
     end do
  end if

  end subroutine open_files

!-----------------------------------------------
  subroutine close_and_rename(filety,iunit,nrens,bname)
!-----------------------------------------------
  implicit none
  character(len=3), intent(in) :: filety
  integer, intent(in) :: iunit
  integer, intent(in) :: nrens
  character(len=80), intent(in) :: bname

  character(len=90) :: filein,fileout
  character(len=3) :: nlab,lfid
  integer n

  ! close files
  do n = 1,nrens
     close(iunit + 10 + n)
  end do

  ! rename fem files
  if (trim(filety) == 'fem') then
	  do n = 1,nrens
	     write(nlab,'(i3.3)') n-1
	     write(lfid,'(i3)') (iunit + 10 + n)

	     filein = 'fort.'//adjustl(trim(lfid))
	     fileout = trim(bname)//'_'//nlab//'.fem'
	     call rename(trim(filein),trim(fileout))
	  end do
  end if

  end subroutine close_and_rename

!-----------------------------------------------
  subroutine perturbe_0d(nrensp,pvec)
!-----------------------------------------------
  use m_random
  implicit none
  integer, intent(in) :: nrensp
  real, intent(out) :: pvec(nrensp)

  integer n
  real aaux,ave

  call random(pvec,nrensp)

  ! remove outlayers
  do n = 1,nrensp
     aaux = pvec(n)
     if( abs(aaux).ge.3. ) then
        aaux = sign(1.,aaux) * (abs(aaux)-floor(abs(aaux)) + 1.)
     end if
     pvec(n) = aaux
  end do
  ! set mean eq to zero
  ave = sum(pvec)/float(nrensp)
  pvec = pvec - ave

  end subroutine perturbe_0d

!-----------------------------------------------
  subroutine red_noise_0d(told,tnew,pvec,nrensp,tau,n,nvar)
!-----------------------------------------------
  implicit none
  double precision, intent(in) :: tnew,tau
  double precision, intent(inout) :: told
  integer, intent(in) :: nrensp,nvar,n
  real, intent(inout) :: pvec(nrensp)
  real, allocatable, save :: pveco(:,:)

  double precision alpha

  if (n > nvar) error stop 'Dimension error.'

  if (told < 0) then
     if (.not. allocated(pveco)) allocate(pveco(nrensp,nvar))
     pveco(:,n) = pvec
  else
     if (tau > 0) then
        alpha = 1. -  (tnew - told)/tau
     else
        alpha = 0.
     end if
     if (alpha < 0.) alpha = 0.

     pveco(:,n) = alpha * pveco(:,n) + sqrt(1 - alpha**2) * pvec
  end if

  pvec = pveco(:,n)

  if (n == nvar) told = tnew

  end subroutine red_noise_0d

!-----------------------------------------------
  subroutine red_noise_2d(nx,ny,nvar,told,tnew,pmat,nrensp,tau)
!-----------------------------------------------
  implicit none
  integer, intent(in) :: nx,ny,nvar
  double precision, intent(in) :: tnew,tau
  double precision, intent(inout) :: told
  integer, intent(in) :: nrensp
  real, intent(inout) :: pmat(nvar,nx,ny,nrensp)
  real, allocatable, save :: pmato(:,:,:,:)

  double precision alpha

  if (told < 0) then
     if (.not. allocated(pmato)) allocate(pmato(nvar,nx,ny,nrensp))
     pmato = pmat
  else
    if (tau > 0) then
       alpha = 1. -  (tnew - told)/tau
    else
       alpha = 0.
    end if
    if (alpha < 0.) alpha = 0.

    pmato = alpha * pmato + sqrt(1 - alpha**2) * pmat
  end if

  pmat = pmato

  told = tnew
  
  end subroutine red_noise_2d


!-----------------------------------------------
  subroutine perturbe_time(nrens,time_std,tens)
!-----------------------------------------------
  use m_random
  implicit none
  integer, intent(in) :: nrens
  double precision, intent(in) :: time_std
  double precision, intent(out) :: tens(nrens)

  real tmpv(nrens)
  integer n
  real aaux,ave

  tens = 0.d0
  if ( nrens == 1) return
  if ( time_std <= 0. ) return
  write(*,*)'  Perturbing forcing in time'

  call random(tmpv,nrens)

  ! remove outlayers
  do n = 1,nrens
     aaux = tmpv(n)
     if( abs(aaux).ge.3. ) then
        aaux = sign(1.,aaux) * (abs(aaux)-floor(abs(aaux)) + 1.)
     end if
     tmpv(n) = aaux
  end do
  ! set mean eq to zero
  ave = sum(tmpv)/float(nrens)
  tmpv = tmpv - ave
  tmpv(1) = 0.

  tens = tmpv * time_std

  end subroutine perturbe_time

!-----------------------------------------------
  subroutine write_record_0d(iunit,told,nrens,dstring,var0,var_std,var_min,var_max,pvec,flag)
!-----------------------------------------------
  implicit none
  double precision, intent(in) :: told
  integer, intent(in) :: nrens,iunit
  character(len=80), intent(in) :: dstring
  real, intent(in) :: var0,var_std,var_min,var_max
  real, intent(in) :: pvec(nrens-1)
  real, intent(in) :: flag

  integer n,fid
  real var


  do n = 1,nrens
     fid = iunit + 10 + n

     if (n > 1) then
        var = (pvec(n-1) * var_std) + var0
     else
        var = var0
     end if

     if (told < 0) var = var0

     if (nint(var_min) /= nint(flag)) call var_limit_min(var,var_min,flag)
     if (nint(var_max) /= nint(flag)) call var_limit_max(var,var_max,flag)

     write(fid,*) trim(dstring),var
  end do

  end subroutine write_record_0d

!-----------------------------------------------
  subroutine make_member_2D1(nvar,nrens,nx,ny,lmax,ne,vec1,var3d,var3d_ens,var_std,flag)
!-----------------------------------------------
  implicit none
  integer, intent(in) :: nvar,nrens,nx,ny,lmax,ne
  real, intent(in) :: vec1(nvar,nrens-1),var_std
  real, intent(in) :: var3d(nvar,nx,ny,lmax)
  real, intent(out) :: var3d_ens(nvar,nx,ny,lmax)
  real, intent(in) :: flag

  integer i
  real vec11(nvar,nrens)

  if ( lmax /= 1 ) error stop 'Bad vertical dimension'

  vec11(:,1) = 0.
  vec11(:,2:nrens) = vec1
  var3d_ens = flag

  do i = 1,nvar
     where (nint(var3d(i,:,:,1)).ne.nint(flag))
	  var3d_ens(i,:,:,1) = var3d(i,:,:,1) + var_std * vec11(i,ne)
     end where
  end do

  end subroutine make_member_2D1

!-----------------------------------------------
  subroutine make_member_2D2(nvar,nrens,nx,ny,lmax,ne,pmat,var3d,var3d_ens,var_std,flag)
!-----------------------------------------------
  implicit none
  integer, intent(in) :: nvar,nrens,nx,ny,lmax,ne
  real, intent(in) :: pmat(nvar-1,nx,ny,nrens-1),var_std
  real, intent(in) :: var3d(nvar,nx,ny,lmax)
  real, intent(out) :: var3d_ens(nvar,nx,ny,lmax)
  real, intent(in) :: flag

  integer i,ix,iy
  real pmat1(nvar-1,nx,ny,nrens),v

  if ( lmax /= 1 ) error stop 'Bad vertical levels'

  pmat1(:,:,:,1) = 0.
  pmat1(:,:,:,2:nrens) = pmat
  var3d_ens = flag

  ! u,v wind perturbed
  do i=1,nvar-1
    where (nint(var3d(i,:,:,1)).ne.nint(flag))
	    var3d_ens(i,:,:,1) = var3d(i,:,:,1) + var_std * pmat1(i,:,:,ne)
    end where
  end do

  ! pressure unperturbed
  var3d_ens(nvar,:,:,1) = var3d(nvar,:,:,1)

  end subroutine make_member_2D2

!-----------------------------------------------
  subroutine make_member_3D1(nvar,nrens,nx,ny,lmax,ne,vec1,var3d,var3d_ens,var_std,flag)
!-----------------------------------------------
  implicit none
  integer, intent(in) :: nvar,nrens,nx,ny,lmax,ne
  real, intent(in) :: vec1(nvar,nrens-1),var_std
  real, intent(in) :: var3d(nvar,nx,ny,lmax)
  real, intent(out) :: var3d_ens(nvar,nx,ny,lmax)
  real, intent(in) :: flag

  integer i,ix,iy,l
  real vec11(nvar,nrens),v

  vec11(:,1) = 0.
  vec11(:,2:nrens) = vec1
  var3d_ens = flag

  do i=1,nvar
    do l = 1,lmax
       where (nint(var3d(i,:,:,l)).ne.nint(flag))
	       var3d_ens(i,:,:,l) = var3d(i,:,:,l) + var_std * vec11(i,ne)
       end where
    end do
  end do

  end subroutine make_member_3D1

!-----------------------------------------------
  subroutine set_decorrelation(ntot,delta,lrange)
!-----------------------------------------------
! For pseudo 2D fields, set rx,ry, the horizontal decorrelation length. Suppose to use geo degree
  implicit none
  real delta,lrange,lrange_low
  integer ntot
  integer ires
  integer, parameter :: nmin = 12	!minimum number of deltas to resolve the range

  if ( ntot < nmin ) error stop 'Field too small'

  ! minimum range allowed
  lrange_low = delta * nmin

  ! range of 3 degrees
  lrange = 4.

  end subroutine set_decorrelation


!--------------------------------------------------
	subroutine make_geo_field(nvar,iens,nrens,nx,ny,lmax,dx,dy, &
     		pmat,datain,dataout,flag,sigmaP,flat,bpress)
!--------------------------------------------------
	implicit none

	integer,intent(in) :: nvar,iens
	integer,intent(in) :: nrens,nx,ny,lmax
	real,intent(in) :: dx,dy
	real,intent(in) :: sigmaP,flag
	real,intent(in) :: flat
        logical,intent(in) :: bpress
	real,intent(in) :: pmat(1,nx,ny,nrens-1)
	real,intent(in) :: datain(nvar,nx,ny,lmax)
	real,intent(out) :: dataout(nvar,nx,ny,lmax)

	real :: mat(1,nx,ny,nrens)

	integer ix,iy,ivar
	real Up,Vp
	real dxm,dym

	real, parameter :: pi = acos(-1.)
	real, parameter :: rhoa = 1.2041
	real, parameter :: er1 = 6378137. !max earth radius
	real, parameter :: er2 = 6356752. !min earth radius
	real fcor,er,theta,DPY,DPX

	if ( lmax /= 1 ) error stop 'Bad vertical levels'

        mat(:,:,:,1) = 0.
        mat(:,:,:,2:nrens) = pmat
  
	theta = flat * pi/180.
	fcor = 2. * sin(theta) * (2.* pi / 86164.)
	! earth radius with latitude
	er = sqrt( ( (er1**2 * cos(theta))**2 + &
     			(er2**2 * sin(theta))**2 ) / &
     			( (er1 * cos(theta))**2 + &
     			(er2 * sin(theta))**2 ) ) 

	dxm = dx * pi/180. * er
	dym = dy * pi/180. * er

	do ivar = 1,nvar

          select case(ivar)
	      case(1)	!u-wind

		do ix = 1,nx
		do iy = 2,ny

		  DPY = sigmaP * (mat(1,ix,iy,iens) - mat(1,ix,iy-1,iens))

		  Up = - (DPY/dym) * 1./(rhoa * fcor)

		  !dataout(ivar,ix,iy,1) = Up
		  dataout(ivar,ix,iy,1) = datain(ivar,ix,iy,1) + Up

		  if (datain(ivar,ix,iy,1) == flag) dataout(ivar,ix,iy,1) = flag

		end do
		end do

		! First row unchanged
		dataout(ivar,:,1,1) = datain(ivar,:,1,1)
		
	    case(2)	!v-wind

		do iy = 1,ny
		do ix = 2,nx

		  DPX = sigmaP * (mat(1,ix,iy,iens) - mat(1,ix-1,iy,iens))

		  Vp = (DPX/dxm) * 1./(rhoa * fcor)

		  !dataout(ivar,ix,iy,1) =  Vp
		  dataout(ivar,ix,iy,1) = datain(ivar,ix,iy,1) + Vp

		  if (datain(ivar,ix,iy,1) == flag) dataout(ivar,ix,iy,1) = flag

		end do
		end do

		! First row unchanged
		dataout(ivar,1,:,1) = datain(ivar,1,:,1)

	    case(3)	!pressure

	      if (bpress) then
		do iy = 1,ny
		do ix = 1,nx
		     dataout(ivar,ix,iy,1) = datain(ivar,ix,iy,1) + sigmaP &
     					* mat(1,ix,iy,iens)
		  if (datain(ivar,ix,iy,1) == flag) dataout(ivar,ix,iy,1) = flag
		end do
		end do
              else
		!write(*,*) 'pressure not perturbed'
		dataout = datain ! no perturbation
              end if

	  end select
	end do
	
	end subroutine make_geo_field


!--------------------------------------------------
	subroutine make_ws_pert(nvar,iens,nrens,nx,ny,lmax,pmat, &
		datain,dataout,flag,err)
!--------------------------------------------------
	implicit none
        integer,intent(in) :: nvar,iens
        integer,intent(in) :: nrens,nx,ny,lmax
        real,intent(in) :: pmat(1,nx,ny,nrens-1)
        real,intent(in) :: datain(nvar,nx,ny,lmax)
        real,intent(out) :: dataout(nvar,nx,ny,lmax)
        real,intent(in) :: err,flag

	real :: mat(1,nx,ny,nrens)

	real wso(nx,ny),wse(nx,ny)
	integer ivar,ix,iy

	if ( lmax /= 1 ) error stop 'Bad vertical dimensions'

        mat(:,:,:,1) = 0.
        mat(:,:,:,2:nrens) = pmat

	wso = sqrt(datain(1,:,:,1)**2 + datain(2,:,:,1)**2)
	wse = wso + err * mat(1,:,:,iens)
	where (wse < 0)
		wse = 0
	end where

	do ivar = 1,nvar
	   select case (ivar)
	     case default	!wind
		do ix = 1,nx
		do iy = 1,ny
	     	  dataout(ivar,ix,iy,1) = datain(ivar,ix,iy,1) * wse(ix,iy)/wso(ix,iy)
		  if (datain(ivar,ix,iy,1) == flag) dataout(ivar,ix,iy,1) = flag
		end do
		end do
	     case (3)	!pressure
	       dataout(ivar,:,:,1) = datain(ivar,:,:,1)
	   end select
	end do


	end subroutine make_ws_pert

!--------------------------------------------------
  subroutine var_limit_min(var,var_min,flag)
!--------------------------------------------------
  implicit none
  real,intent(inout) :: var
  real,intent(in) :: var_min
  real,intent(in) :: flag

  if ((var < var_min).and.(nint(var).ne.nint(flag))) then
	  var = var_min
  end if
    
  end subroutine var_limit_min


!--------------------------------------------------
  subroutine var_limit_max(var,var_max,flag)
!--------------------------------------------------
  implicit none
  real,intent(inout) :: var
  real,intent(in) :: var_max
  real,intent(in) :: flag

  if ((var > var_max).and.(nint(var).ne.nint(flag))) var = var_max
    
  end subroutine var_limit_max

!--------------------------------------------------
  subroutine limit_wind(nv,nx,ny,lmax,wdata,wsmax,flag)
!--------------------------------------------------
  implicit none
  integer, intent(in) :: nv,nx,ny,lmax
  real, intent(inout) ::  wdata(nv,nx,ny,lmax)
  real, intent(in)    :: wsmax
  real, intent(in)    :: flag
  real u,v,ws,k
  integer ix,iy

  if (nint(wsmax) == nint(flag)) return

  if ( lmax /= 1 ) error stop 'Bad vertical dimension'
  if ((nv /= 3).or.(wsmax < 1.)) error stop 'Wrong wind speed limit'

  do ix = 1,nx
  do iy = 1,ny
     u = wdata(1,ix,iy,1)
     v = wdata(2,ix,iy,1)
     ws = sqrt(u**2 + v**2)
     if (ws > wsmax) then
        k = wsmax/ws
	u = u * k
	v = v * k
     end if
     wdata(1,ix,iy,1) = u
     wdata(2,ix,iy,1) = v
  end do
  end do
  
  end subroutine limit_wind

!--------------------------------------------------
  subroutine fem_check_last(iformat,iunit,blast)
!--------------------------------------------------
  implicit none
  integer, intent(in) :: iformat,iunit
  logical, intent(inout) ::  blast

  double precision dtime
  integer nvers,np,lmax,nvar,ntype,ierr
  integer datetime(2)     !date and time information

  blast = .false.

  call fem_file_peek_params(iformat,iunit,dtime &
         ,nvers,np,lmax,nvar,ntype,datetime,ierr)

  if( ierr /= 0 ) blast = .true.

  end subroutine fem_check_last

