!
! Copyright (C) 2017, Marco Bajo, CNR-ISMAR Venice, All rights reserved.
! 
!---------------------------------------------------------------------------------
! Convert an EnKF simulation saved in RST format into an EnKS analysis.
! This program reads the ensemble RST files from the KF run and the file
! "X5_tot.uf" (the ensemble transform matrix X5 at each timestep) and produces
! EnKS mean/STD (and optionally the full smoothed ensemble).
!---------------------------------------------------------------------------------
program enKF2enKS
  use iso_fortran_env, only: dp => real64
  use mod_restart
  use levels,     only: nlvdi, nlv, hlv, ilhv, ilhkv
  use shympi
  use basin
  implicit none

  !---------------------------- CLI arguments ----------------------------
  character(len=80) :: basinf           ! basin filename (with extension)
  character(len=3)  :: lnnlv            ! number of vertical levels (as text)
  character(len=6)  :: lnrens           ! number of ensemble members (as text)
  character(len=4)  :: lout             ! output mode: 'full' or 'norm'
  character(len=6)  :: lnlag            ! smoother lag steps, -1 => all

  integer :: nnlv                       ! parsed nlv
  integer :: nrens, nre                 ! parsed ensemble size, loop index
  character(len=5) :: nrel              ! member index formatted (00000..)
  integer :: ierr, iflag
  real(dp) :: atime                     ! record time
  integer :: fid                        ! file unit id

  !---------------------------- Ensemble storage ----------------------------
  real(dp), allocatable :: Astate(:,:), AmeanKF(:), AstdKF(:), AmeanKS(:), AstdKS(:)
  integer :: rrec, sdim                 ! record counter; state vector size
  integer :: nlag                       ! (re)parsed smoother lag

  !---------------------------- Read & check inputs ----------------------------
  call get_command_argument(1, basinf)
  call get_command_argument(2, lnnlv)
  call get_command_argument(3, lnrens)
  call get_command_argument(4, lout)
  call get_command_argument(5, lnlag)

  if (lnlag == '') then
     write(*,*) ''
     write(*,*) 'Usage: enKF2enKS [basinf] [nlv] [nrens] [output] [nlag]'
     write(*,*) ''
     write(*,*) '[basinf] bas file with the extension.'
     write(*,*) '[nlv]     number of vertical levels.'
     write(*,*) '[nrens]   number of ensemble members, control included.'
     write(*,*) "[output]  full (full) or just mean and std (norm)."
     write(*,*) '[nlag]    number of forward analysis steps to consider (-1 = all)'
     write(*,*) ''
     stop
  end if

  read(lnnlv, *) nnlv
  if (nnlv < 1) error stop 'enKF2enKS: bad vertical levels'

  read(lnrens, *) nrens
  if (nrens < 3) error stop 'enKF2enKS: not a valid nrens number'

  if ((lout /= 'full') .and. (lout /= 'norm')) error stop 'enKF2enKS: bad output specification'

  read(lnlag, *) nlag
  if ((nlag < 2) .and. (nlag /= -1)) error stop 'enKF2enKS: bad nlag'

  write(*,*) 'Limited Kalman smoother. Number of time steps: ', nlag

  call init_shyfem(basinf, nnlv)

  !---------------------------- Open input RST files ----------------------------
  do nre = 1, nrens
     fid = 20 + nre
     call num2str(nre-1, nrel)
     open(fid, file='analKF_en'//nrel//'.rst', status='old', form='unformatted', iostat=ierr)
     if (ierr /= 0) error stop 'enKF2enKS: error opening file'
  end do

  !---------------------------- Open outputs ----------------------------
  if (lout == 'full') then
     write(*,*) 'Writing the full output'
     do nre = 1, nrens
        fid = 20 + nrens + nre
        call num2str(nre-1, nrel)
        open(fid, file='analKS_en'//nrel//'.rst', status='unknown', form='unformatted')
     end do
  end if

  open(18, file='analKS_mean.rst', status='unknown', form='unformatted')
  open(19, file='analKS_std.rst',  status='unknown', form='unformatted')

  !---------------------------- Loop on time records ----------------------------
  rrec = 0
89 continue

  nkn_global = nkn; nel_global = nel; nlv_global = nlv

  ! Read ensemble record across all members, same time index
  do nre = 1, nrens
     fid = 20 + nre
     call rst_read_record(fid, atime, iflag, ierr)
     if (ierr == -1) goto 90
     if (ierr > 0) error stop 'Error in reading the restart files'

     if (.not. allocated(Astate)) then
        sdim = nkn + 2*nlv*nel + 2*nlv*nkn
        allocate(Astate(sdim, nrens))
        allocate(AmeanKF(sdim), AmeanKS(sdim))
        allocate(AstdKF(sdim),  AstdKS(sdim))
        Astate = 0.0_dp
     end if

     call push_matrix(sdim, nrens, nre, Astate)
  end do

  ! Add RST parameters (first record only)
  if (rrec == 0) then
     call addpar('ibarcl', ibarcl_rst)
     call addpar('iconz',  iconz_rst)
     call addpar('iwvert', iwvert_rst)  ! maybe not
     call addpar('ieco',   ieco_rst)    ! maybe not
     call addpar('ibio',   0)
     call addpar('ibfm',   0)
     call addpar('imerc',  imerc_rst)
     call addpar('iturb',  iturb_rst)
     call daddpar('date',  0.0_dp)
     call daddpar('time',  0.0_dp)
  end if

  ! Perform analysis (method 2 is faster; method 1 kept for reference)
  ! call make_analysis1(atime, sdim, nrens, Astate)
  call make_analysis2(atime, sdim, nrens, Astate, nlag)

  ! Compute EnKS mean and std of the ensemble
  call make_mn_std(sdim, nrens, Astate, AmeanKS, AstdKS)

  call pull_state(sdim, AmeanKS)
  call rst_write_record(atime, 18)

  call pull_state(sdim, AstdKS)
  call rst_write_record(atime, 19)

  ! Optional: write full smoothed ensemble
  if (lout == 'full') then
     do nre = 1, nrens
        fid = 20 + nrens + nre
        call pull_matrix(sdim, nrens, nre, Astate)
        call rst_write_record(atime, fid)
     end do
  end if

  rrec = rrec + 1
  write(*,*) 'N. of record: ', rrec
  goto 89

90 continue

  !---------------------------- Close files ----------------------------
  do nre = 1, nrens
     fid = 20 + nre
     close(fid)
  end do

  close(18)
  close(19)

  if (lout == 'full') then
     do nre = 1, nrens
        fid = 20 + nrens + nre
        close(fid)
     end do
  end if

end program enKF2enKS

!*******************************************************************************
! Utility: format an integer with leading zeros into a 5-char string
!*******************************************************************************
subroutine num2str(num, str)
  implicit none
  integer, intent(in)        :: num
  character(len=5), intent(out) :: str

  if ((num >= 0).and.(num < 10)) then
     write(str,'(a4,i1)') '0000', num
  elseif ((num >= 10).and.(num < 100)) then
     write(str,'(a3,i2)') '000',  num
  elseif ((num >= 100).and.(num < 1000)) then
     write(str,'(a2,i3)') '00',   num
  elseif ((num >= 1000).and.(num < 10000)) then
     write(str,'(a1,i4)') '0',    num
  elseif ((num >= 10000).and.(num < 100000)) then
     write(str,'(i5)')     num
  else
     error stop 'num2str: num out of range'
  end if
end subroutine num2str

!*******************************************************************************
! Push model state from SHYFEM fields into the ensemble matrix column nre
! Layout (vectorized): u, v, z, T, S
!*******************************************************************************
subroutine push_matrix(sdim, nrens, nre, Amat)
  use iso_fortran_env, only : dp => real64
  use basin
  use levels, only : nlv
  use mod_hydro
  use mod_hydro_vel
  use mod_ts
  use mod_conz
  implicit none
  integer,  intent(in)    :: sdim, nrens, nre
  real(dp), intent(inout) :: Amat(sdim,nrens)

  integer :: dimz, dimuv, dimts

  dimz  = nkn
  dimuv = nlv*nel
  dimts = nlv*nkn

  Amat(1:dimuv, nre)                             = reshape(real(utlnv, dp), (/dimuv/))
  Amat(dimuv+1:2*dimuv, nre)                     = reshape(real(vtlnv, dp), (/dimuv/))
  Amat(2*dimuv+1:2*dimuv+dimz, nre)              = real(znv, dp)
  Amat(2*dimuv+dimz+1:2*dimuv+dimz+dimts, nre)   = reshape(real(tempv, dp), (/dimts/))
  Amat(2*dimuv+dimz+dimts+1:2*dimuv+dimz+2*dimts, nre) = reshape(real(saltv, dp), (/dimts/))
end subroutine push_matrix

!*******************************************************************************
! Pull one ensemble member back into SHYFEM fields from the ensemble matrix
!   PATCH: added 'use iso_fortran_env, only: dp => real64' (dp visible here)
!*******************************************************************************
subroutine pull_matrix(sdim, nrens, nre, Amat)
  use iso_fortran_env, only : dp => real64
  use basin
  use levels, only : nlv
  use mod_hydro
  use mod_hydro_vel
  use mod_ts
  use mod_conz
  implicit none
  integer,  intent(in) :: sdim, nrens, nre
  real(dp), intent(in) :: Amat(sdim,nrens)

  integer :: dimz, dimuv, dimts

  dimz  = nkn
  dimuv = nlv*nel
  dimts = nlv*nkn

  utlnv = reshape(real(Amat(1:dimuv, nre)),                    (/nlv, nel/))
  vtlnv = reshape(real(Amat(dimuv+1:2*dimuv, nre)),            (/nlv, nel/))
  znv   = real(Amat(2*dimuv+1:2*dimuv+dimz, nre))
  tempv = reshape(real(Amat(2*dimuv+dimz+1:2*dimuv+dimz+dimts, nre)), (/nlv, nkn/))
  saltv = reshape(real(Amat(2*dimuv+dimz+dimts+1:2*dimuv+dimz+2*dimts, nre)), (/nlv, nkn/))
end subroutine pull_matrix

!*******************************************************************************
! Pull a full state vector back into SHYFEM fields (mean or std)
!   PATCH: added 'use iso_fortran_env, only: dp => real64' (dp visible here)
!*******************************************************************************
subroutine pull_state(nvec, Avec)
  use iso_fortran_env, only : dp => real64
  use basin
  use levels, only : nlv
  use mod_hydro
  use mod_hydro_vel
  use mod_ts
  use mod_conz
  implicit none
  integer,  intent(in) :: nvec
  real(dp), intent(in) :: Avec(nvec)

  integer :: dimz, dimuv, dimts

  dimz  = nkn
  dimuv = nlv*nel
  dimts = nlv*nkn

  utlnv = reshape(real(Avec(1:dimuv)),                         (/nlv, nel/))
  vtlnv = reshape(real(Avec(dimuv+1:2*dimuv)),                 (/nlv, nel/))
  znv   = real(Avec(2*dimuv+1:2*dimuv+dimz))
  tempv = reshape(real(Avec(2*dimuv+dimz+1:2*dimuv+dimz+dimts)), (/nlv, nkn/))
  saltv = reshape(real(Avec(2*dimuv+dimz+dimts+1:2*dimuv+dimz+2*dimts)), (/nlv, nkn/))
end subroutine pull_state

!*******************************************************************************
! Initialize SHYFEM environment from basin file and set levels
!*******************************************************************************
subroutine init_shyfem(basinf, nnlv)
  use basin
  use levels
  use mod_geom_dynamic
  use mod_hydro
  use mod_hydro_vel
  use mod_ts
  use mod_conz
  use mod_restart
  implicit none
  character(len=80), intent(in) :: basinf
  integer,          intent(in) :: nnlv

  ! Read basin and check dimensions
  open(21, file=basinf, status='old', form='unformatted')
  call basin_read_by_unit(21)
  close(21)

  ! Set vertical levels
  nlv   = nnlv
  nlvdi = nnlv

  ! Initialize SHYFEM modules (sizes from basin)
  call mod_geom_dynamic_init(nkn, nel)
  call mod_hydro_init(      nkn, nel, nlv)
  call mod_hydro_vel_init(  nkn, nel, nlv)
  call mod_ts_init(         nkn, nlv)
  call levels_init(         nkn, nel, nlv)

  ! Concentrations initialization (left commented as in original)
  ! call mod_conz_init(1, nkn, nlvdi)
end subroutine init_shyfem

!*******************************************************************************
! Apply X5 transform read from file (method 1: per-time multiplication)
!*******************************************************************************
subroutine make_analysis1(atime, sdim, nrens, Amat)
  use iso_fortran_env, only: dp => real64
  implicit none
  real(dp), intent(in)    :: atime
  integer,  intent(in)    :: sdim, nrens
  real(dp), intent(inout) :: Amat(sdim, nrens)

  real(dp) :: tt
  character(len=6) :: alabel
  character(len=2) :: tag
  integer :: nren
  real(dp), allocatable :: Aaux(:,:), X5(:,:)

  external :: dgemm

  allocate(Aaux(sdim, nrens))
  allocate(X5(nrens, nrens))

  open(15, status='old', form='unformatted', file='X5_tot.uf')
45 read(15, end=44) tt, alabel, tag
  read(15) nren
  if ((tag == 'X3') .or. (nren /= nrens)) then
     error stop 'Error reading X5_tot.uf.'
  end if
  read(15) X5

  Aaux = 0.0_dp
  if (tt >= atime) then
     call dgemm('N','N', sdim, nrens, nrens, 1.0_dp, Amat, sdim, X5, nrens, 0.0_dp, Aaux, sdim)
     Amat = Aaux
  end if
  goto 45

44 close(15)
  deallocate(Aaux, X5)
end subroutine make_analysis1

!*******************************************************************************
! Accumulate X5 transforms forward up to nlag steps (or all if nlag=-1)
! and apply the resulting product to the ensemble (faster method).
!*******************************************************************************
subroutine make_analysis2(atime, sdim, nrens, Amat, nlag)
  use iso_fortran_env, only: dp => real64
  implicit none
  real(dp), intent(in)    :: atime
  integer,  intent(in)    :: sdim, nrens
  real(dp), intent(inout) :: Amat(sdim, nrens)
  integer,  intent(inout) :: nlag

  real(dp) :: tt
  character(len=6) :: alabel
  character(len=2) :: tag
  integer :: nren
  real(dp), allocatable :: Aaux(:,:), X5(:,:), X5old(:,:), X5aux(:,:)
  integer :: n, k

  external :: dgemm

  allocate(Aaux(sdim, nrens))
  allocate(X5(nrens,nrens), X5old(nrens,nrens), X5aux(nrens,nrens))

  if (nlag == -1) nlag = 1000000000

  ! Initialize X5old = I
  X5old = 0.0_dp
  do n = 1, nrens
     X5old(n,n) = 1.0_dp
  end do

  open(15, status='old', form='unformatted', file='X5_tot.uf')
  k = 0
45 read(15, end=44) tt, alabel, tag
  read(15) nren
  read(15) X5

  if ((tag == 'X3') .or. (nren /= nrens)) error stop 'X3 not implemented or bad n. of ens members'

  if ((tt > atime) .and. (k < nlag)) then
     X5aux = 0.0_dp
     call dgemm('N','N', nrens, nrens, nrens, 1.0_dp, X5old, nrens, X5, nrens, 0.0_dp, X5aux, nrens)
     X5old = X5aux
     k = k + 1
  end if
  goto 45

44 continue
  close(15)

  Aaux = 0.0_dp
  call dgemm('N','N', sdim, nrens, nrens, 1.0_dp, Amat, sdim, X5old, nrens, 0.0_dp, Aaux, sdim)
  Amat = Aaux

  deallocate(Aaux, X5, X5old, X5aux)
end subroutine make_analysis2

!*******************************************************************************
! Ensemble mean and standard deviation (per state component)
!*******************************************************************************
subroutine make_mn_std(ndim, nens, Amat, Am, Astd)
  use iso_fortran_env, only: dp => real64
  implicit none
  integer,  intent(in)  :: ndim, nens
  real(dp), intent(in)  :: Amat(ndim, nens)
  real(dp), intent(out) :: Am(ndim), Astd(ndim)

  real(dp) :: esum, esumsq
  integer  :: n

  if (nens < 2) error stop 'Ensemble too small'

  do n = 1, ndim
     esum   = sum(Amat(n,:))
     esumsq = sum(Amat(n,:)*Amat(n,:))
     Am(n)   = esum / real(nens, dp)
     Astd(n) = sqrt( (esumsq - esum*esum/real(nens,dp)) / (real(nens,dp) - 1.0_dp) )
  end do
end subroutine make_mn_std
