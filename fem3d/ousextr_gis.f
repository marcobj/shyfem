c
c $Id: ousextr.f,v 1.3 2009-09-14 08:20:58 georg Exp $
c
c extracts records from OUS files
c
c revision log :
c
c 02.09.2003	ggu	adapted to new OUS format
c 24.01.2005	ggu	computes maximum velocities for 3D (only first level)
c 23.03.2010	ggu	extracts reocrds
c 26.03.2010	ggu	bug fix: set nkn and nel
c 23.11.2010	ggu	new for 3D gis output
c
c***************************************************************

	program ousextr_gis

c reads ous file and writes extracted records in ascii to new file

	implicit none

        include 'param.h'

        integer nrdim
        parameter ( nrdim = 2000 )

        integer irec(nrdim)

	character*80 title
	integer nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw
	common /nkonst/ nkn,nel,nrz,nrq,nrb,nbc,ngr,mbw

	real xgv(nkndim), ygv(nkndim)
	real hm3v(3,neldim)
	integer nen3v(3,neldim)
	integer ipev(neldim), ipv(nkndim)
	integer iarv(neldim)
	common /xgv/xgv, /ygv/ygv
	common /hm3v/hm3v
	common /nen3v/nen3v
	common /ipev/ipev, /ipv/ipv
	common /iarv/iarv

	integer ilhv(neldim)
	integer ilhkv(nkndim)
	real hlv(nlvdim)
        real utlnv(nlvdim,neldim)
        real vtlnv(nlvdim,neldim)
	common /ilhv/ilhv
	common /ilhkv/ilhkv
	common /hlv/hlv
        common /utlnv/utlnv
        common /vtlnv/vtlnv

	real hev(neldim)

	real znv(nkndim)
	real zenv(3,neldim)

	real uprv(nlvdim,nkndim)
	real vprv(nlvdim,nkndim)
	real ut2v(neldim)
	real vt2v(neldim)
	real u2v(neldim)
	real v2v(neldim)
	real weight(nlvdim,nkndim)

	integer ilnv(nlvdim,nkndim)

	character*80 name,file
        integer nvers,nin,nlv
        integer itanf,itend,idt,idtous
	integer it,ie,i
        integer ierr,nread,nextr,nb
        integer nknous,nelous,nlvous
        real href,hzoff,hlvmin
	real zmin,zmax
	real umin,umax
	real vmin,vmax
	real xe,ye
	integer k,ke,ivar
	integer lmax,l

	integer iapini,ideffi,ifileo
	logical berror,ball,bwrite

c-------------------------------------------------------------------
c initialize params
c-------------------------------------------------------------------

	ball = .true.		!write all records
	ball = .false.

	nread=0
	nextr=0

c-------------------------------------------------------------------
c get simulation
c-------------------------------------------------------------------

	if(iapini(3,nkndim,neldim,0).eq.0) then
		stop 'error stop : iapini'
	end if

        nin=ideffi('datdir','runnam','.ous','unform','old')
        if(nin.le.0) goto 100

c--------------------------------------------------------------------
c open OUS file and read header
c--------------------------------------------------------------------

	nvers=1
        call rfous(nin
     +			,nvers
     +			,nknous,nelous,nlvous
     +			,href,hzoff
     +			,title
     +			,ierr)
        if(ierr.ne.0) goto 100

        write(6,*)
        write(6,*)   title
        write(6,*)
        write(6,*) ' nvers        : ',nvers
        write(6,*) ' href,hzoff   : ',href,hzoff
        write(6,*) ' nkn,nel      : ',nknous,nelous
        write(6,*) ' nlv          : ',nlvous
        write(6,*)

	nkn=nknous
	nel=nelous
	nlv=nlvous
	call dimous(nin,nkndim,neldim,nlvdim)

	call rsous(nin,ilhv,hlv,hev,ierr)
        if(ierr.ne.0) goto 100

	call level_e2k(nkn,nel,nen3v,ilhv,ilhkv)

        write(6,*) 'Available levels: ',nlv
        write(6,*) (hlv(l),l=1,nlv)

c-------------------------------------------------------------------
c get records to extract from STDIN
c-------------------------------------------------------------------

        if( .not. ball ) then
          call get_records_from_stdin(nrdim,irec)
        end if


c-------------------------------------------------------------------
c open OUS output file
c-------------------------------------------------------------------

        call mkname(' ','extract_ous','.gis',file)
        write(6,*) 'writing file ',file(1:50)
        nb = ifileo(55,file,'form','new')
        if( nb .le. 0 ) goto 98

c-------------------------------------------------------------------
c loop on input records
c-------------------------------------------------------------------

  300   continue

        call rdous(nin,it,nlvdim,ilhv,znv,zenv,utlnv,vtlnv,ierr)

        if(ierr.gt.0) write(6,*) 'error in reading file : ',ierr
        if(ierr.ne.0) goto 100

	nread=nread+1
	if( nread .gt. nrdim ) goto 100
	write(6,*) 'time : ',nread,it,irec(nread)

        bwrite = ball .or. irec(nread) .ne. 0

        if( bwrite ) then
	  call transp2vel(nel,nkn,nlv,nlvdim,hev,zenv,nen3v
     +                          ,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv,weight)
          call wrgis_3d(nb,it,nkn,ilhkv,znv,uprv,vprv)
          nextr = nextr + 1
        end if

	goto 300

  100	continue

c-------------------------------------------------------------------
c end of loop
c-------------------------------------------------------------------

	write(6,*)
	write(6,*) nread,' records read'
        write(6,*) nextr,' records written to file extract.ous'
	write(6,*)

        if( nextr .le. 0 ) stop 'no file written'

c-------------------------------------------------------------------
c end of routine
c-------------------------------------------------------------------

	stop
   98   continue
        write(6,*) 'error opening file'
        stop 'error stop ousextr_records'
   99   continue
        write(6,*) 'error writing file'
        stop 'error stop ousextr_records'
	end

c******************************************************************

        subroutine get_records_from_stdin(ndim,irec)

c gets records to extract from stdin

        implicit none

        integer ndim
        integer irec(ndim)

        integer i,ir

        do i=1,ndim
          irec(i) = 0
        end do

        write(6,*) 'Please enter the record numbers to be extracted.'
        write(6,*) 'Enter every record on a single line.'
        write(6,*) 'Finish with 0 on the last line.'
        write(6,*) 'example:'
        write(6,*) '   5'
        write(6,*) '  10'
        write(6,*) '  15'
        write(6,*) '  0'
        write(6,*) ' '

        do while(.true.)
          write(6,*) 'Enter record to extract (0 to end): '
          ir = 0
          read(5,'(i10)') ir

          if( ir .le. 0 ) then
            return
          else if( ir .gt. ndim ) then
            write(6,*) 'Cannot extract records higher than ',ndim
            write(6,*) 'Please change ndim and recompile.'
          else
            irec(ir) = 1
          end if
        end do

        end

c******************************************************************

        subroutine wrgis_3d(nb,it,nkn,ilhkv,znv,uprv,vprv)

c writes one record to file nb (3D)

        implicit none

        include 'param.h'

        integer nb,it,nkn
        integer ilhkv(nkndim)
        real znv(nkndim)
        real uprv(nlvdim,nkndim)
        real vprv(nlvdim,nkndim)

        double precision x0,y0
        !parameter ( x0 = 2330000.-50000., y0 = 5000000. )
        parameter ( x0 = 0., y0 = 0. )

        real xgv(nkndim), ygv(nkndim)
        common /xgv/xgv, /ygv/ygv

        integer k,l,lmax
        real x,y

        write(nb,*) it,nkn

        do k=1,nkn
          lmax = ilhkv(k)
          x = xgv(k) + x0
          y = ygv(k) + y0

          write(nb,*) x,y,lmax,znv(k)
          write(nb,*) (uprv(l,k),l=1,lmax)
          write(nb,*) (vprv(l,k),l=1,lmax)
        end do

        end

c******************************************************************

        subroutine transp2vel(nel,nkn,nlv,nlvdim,hev,zenv,nen3v
     +				,ilhv,hlv,utlnv,vtlnv
     +                          ,uprv,vprv,weight)

c transforms transports at elements to velocities at nodes

        implicit none

        integer nel
        integer nkn
	integer nlv
        integer nlvdim
        real hev(1)
        real zenv(3,1)
	integer nen3v(3,1)
	integer ilhv(1)
	real hlv(1)
        real utlnv(nlvdim,1)
        real vtlnv(nlvdim,1)
        real uprv(nlvdim,1)
        real vprv(nlvdim,1)
        real weight(nlvdim,1)

        integer ie,ii,k,l,lmax
        real zmed,hmed,u,v,w
	real hbot,htop

	do k=1,nkn
	  do l=1,nlv
	    weight(l,k) = 0.
	    uprv(l,k) = 0.
	    vprv(l,k) = 0.
	  end do
	end do
	      
        do ie=1,nel
          zmed = 0.
          do ii=1,3
            zmed = zmed + zenv(ii,ie)
          end do
          zmed = zmed / 3.

	  htop = -zmed
	  lmax = ilhv(ie)
	  do l=1,lmax
	    hbot = hlv(l)
	    if( l .eq. lmax ) hbot = hev(ie)
	    hmed = hbot - htop
	    u = utlnv(l,ie) / hmed
	    v = vtlnv(l,ie) / hmed
	    do ii=1,3
	      k = nen3v(ii,ie)
	      uprv(l,k) = uprv(l,k) + u
	      vprv(l,k) = vprv(l,k) + v
	      weight(l,k) = weight(l,k) + 1.
	    end do
	    htop = hbot
	  end do
	end do

	do k=1,nkn
	  do l=1,nlv
	    w = weight(l,k)
	    if( w .gt. 0. ) then
	      uprv(l,k) = uprv(l,k) / w
	      vprv(l,k) = vprv(l,k) / w
	    end if
	  end do
	end do
	      
	end

c******************************************************************

	subroutine level_e2k(nkn,nel,nen3v,ilhv,ilhkv)

c computes ilhkv from ilhv

	implicit none

	integer nkn,nel
	integer nen3v(3,1)
	integer ilhv(1)
	integer ilhkv(1)

	integer k,ie,ii,lmax

	do k=1,nkn
	  ilhkv(k) = 1
	end do

	do ie=1,nel
	  lmax = ilhv(ie)
	  do ii=1,3
	    k = nen3v(ii,ie)
	    if( ilhkv(k) .lt. lmax ) ilhkv(k) = lmax
	  end do
	end do

	end

c******************************************************************
