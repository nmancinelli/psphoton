program post
 !An example program that postprocesses the output
 ! of psphonon and writes is as an ascii array to
 ! be read by plotting program.
 ! 
 !This program corrects for changes in surface area in each range bin and takes the square
 ! root, which converts power/energy to amplitude.
 !
	implicit none
	real, dimension(:,:), allocatable   :: rbin
	real, dimension(10)                 :: dummy
	character (len=200)                 :: arrayfile, outfile
	integer                             :: it,ix

	real :: x1,x2,t1,t2,stnmin,stnmax,zsource
	integer :: ntdim,nxdim,iwstart,nray,nsurf

	real :: deldeg,delrad,degrad,sang,xbinwidth
	
    degrad=180./3.1415927

	print *, 'Enter binary array file name: '
	read *, arrayfile

	print *, 'Enter output file name: '
	read *, outfile


	open (11, file=trim(arrayfile),form='unformatted')
	read (11) x1,x2,nxdim,t1,t2,ntdim,iwstart,stnmin,stnmax,zsource,zsource,nray,nsurf
	!print *, 'ntdim, ndxim', ntdim, nxdim
	allocate ( rbin ( ntdim , nxdim ) )
	rbin=0. !init to zero
	read (11) dummy
	read (11) ((rbin(it,ix),it=1,ntdim),ix=1,nxdim)
	close (11)
	
	open (12, file=trim(outfile))

	write (12,*) x1,x2,nxdim,t1,t2,ntdim,iwstart,stnmin,stnmax,zsource,zsource,nray,nsurf

    xbinwidth = (x2-x1)/real(nxdim)
	do ix=1,nxdim
			deldeg=(real(ix)-0.5)*xbinwidth + x1
			delrad=deldeg/degrad
			sang=sin(delrad)
			rbin(:,ix)=rbin(:,ix)/sang
			rbin(:,ix)=sqrt(rbin(:,ix))
			write(12,*) rbin(:, ix)
	end do

	close (12)

end program post
