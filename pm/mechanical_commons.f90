!####################################################################
!####################################################################
!####################################################################
module mechanical_commons
   use amr_commons, ONLY:dp,ndim,ncpu
   integer(kind=4),parameter::nvarSN=12

   ! Important Note: SN stands for SN cell, not SN particle (SNp)

   ! Array to define neighbors for SN
   integer, parameter::nSNnei=48  ! number of neighboring cells to deposit mass/momentum/energy
   real(dp),parameter::nSNcen=4   ! number of cells corresponding to the central cell to deposit mass
   real(dp),dimension(1:3,1:nSNnei)::xSNnei
   real(dp),dimension(1:3,1:nSNnei)::vSNnei
   real(dp)::f_LOAD,f_LEFT,f_ESN
   ! SN cells that are needed to communicate across different CPUs
   ! note that SNs are processed on a cell-by-cell basis
   ! hence the central position is taken as the central leaf cell
   integer ::nSN_comm   ! the number of cells to be communicated (specific to myid)
   integer ::uidSN_comm ! the unique id of SN comm

  ! momentum input
  ! p_sn = A_SN*nH**(alpha)*ESN**(beta)*ZpSN**(gamma)
  ! ex) Thornton et al.
  !     A_SN = 3e5, alphaN = -2/17, beta = 16/17, gamma = -0.14
  ! ex) Kim & Ostriker (2015) uniform case
  !     A_SN = 2.17e5, alpha = -0.13, beta = 0.93
   real(dp),parameter:: expE_SN=+16d0/17d0
   real(dp),parameter:: expZ_SN=-0.14

#ifndef WITHOUTMPI
   integer, parameter ::ncomm_max = 50000
   integer ,dimension(1:ncomm_max)::iSN_comm  ! cpu list
   integer ,dimension(1:ncomm_max)::idSN_comm  ! id of SNe in each cpu
   !integer ,dimension(1:ncomm_max)::lSN_comm  ! level
   real(dp),dimension(1:ncomm_max)::mSN_comm          ! gas mass of SNe
   real(dp),dimension(1:ncomm_max)::mZSN_comm         ! metal mass of SNe
   real(dp),dimension(1:ncomm_max)::mZdSN_comm        ! dust mass of SNe 
   real(dp),dimension(1:ncomm_max)::mloadSN_comm      ! ejecta mass + gas entrained
   real(dp),dimension(1:ncomm_max)::eloadSN_comm      ! kinetic energy of ejecta + gas entrained
   real(dp),dimension(1:ncomm_max)::mZloadSN_comm     ! metals ejected + entrained
   real(dp),dimension(1:ncomm_max)::mZdloadSN_comm    ! dust ejected + entrained
   real(dp),dimension(1:3,1:ncomm_max)::xSN_comm      ! pos of SNe host cell (leaf)
   real(dp),dimension(1:3,1:ncomm_max)::pSN_comm      ! total momentum of total SNe in each leaf cell
   real(dp),dimension(1:3,1:ncomm_max)::ploadSN_comm  ! momentum from original star + gas entrained
   real(dp),dimension(1:ncomm_max)::floadSN_comm      ! fraction of gas to be loaded from the central cell
   integer,dimension(:,:),allocatable::icpuSN_comm,icpuSN_comm_mpi
   integer,dimension(:)  ,allocatable::nSN_comm_cpu,nSN_comm_mpi
#endif

   ! refinement for resolved feedback (used 'ring' instead of 'shell' to be more catchy)
   integer::ncshell3                            ! the maximum number of cells within a (1+2*nshell_re
   integer,dimension(:,:),allocatable::xrefnei  ! relative position of the neighboring cells
   integer,dimension(:),allocatable::irefnei    ! index of the nei cells
   integer,dimension(:),allocatable::lrefnei    ! level of the neighboring cells
   integer,dimension(:),allocatable::icellnei   ! cell index of the neighbors
   real(dp),dimension(:),allocatable::mrefnei   ! gas mass within each cell in Msun
   real(dp),dimension(:),allocatable::mzrefnei  ! metal mass within each cell in Msun
   integer,dimension(:),allocatable::nrefnei_ring  ! cumulative number of nei cells per ring - useful
   real(dp),dimension(:),allocatable::mrefnei_ring ! gas mass within each shell in Msun
   real(dp),dimension(:),allocatable::mzrefnei_ring! metal mass within each shell in Msun

   integer,dimension(:),allocatable::icommr     ! for communication

   ! Added by Joki (or rather, moved from amr_parameters.f90, since not used)
   integer ::nshell_resolve=3  ! r_shell will be resolved on 3 cells in one direction

end module
!####################################################################
!####################################################################
!####################################################################
subroutine init_mechanical
   use amr_commons
   use mechanical_commons
   implicit none
   integer::i,j,k,ind,indall
   real(kind=dp)::x,y,z,r
   logical::ok
   integer,allocatable::ltmp(:)
   integer::ncshell,iring,i2,j2,k2


   !------------------------------------------------
   ! Warning messages
   !------------------------------------------------
   ok=.false.
   if(.not.metal.and.mechanical_feedback>0)then
      print *, '>>> mechanical Err: Please turn on metal'
      ok=.true.
   endif
   if(ok) call clean_stop

#ifndef WITHOUTMPI
   allocate(nSN_comm_cpu(1:ncpu))
   allocate(nSN_comm_mpi(1:ncpu))
#endif

   ! some parameters
   f_LOAD = nSNnei / dble(nSNcen + nSNnei)
   f_LEFT = nSNcen / dble(nSNcen + nSNnei)
   f_ESN  = 0.676   ! Blondin+(98) at t=trad

   ! Arrays to define neighbors (center=[0,0,0])
   ! normalized to dx = 1 = size of the central leaf cell in which a SN particle sits
   ! from -0.75 to 0.75
   ind=0
   do k=1,4
   do j=1,4
   do i=1,4
      ok=.true.
      if((i==1.or.i==4).and.&
         (j==1.or.j==4).and.&
         (k==1.or.k==4)) ok=.false. ! edge
      if((i==2.or.i==3).and.&
         (j==2.or.j==3).and.&
         (k==2.or.k==3)) ok=.false. ! centre
      if(ok)then
         ind=ind+1
         x = (i-1)+0.5d0 - 2
         y = (j-1)+0.5d0 - 2
         z = (k-1)+0.5d0 - 2
         r = dsqrt(dble(x*x+y*y+z*z))
         xSNnei(1,ind) = x/2d0
         xSNnei(2,ind) = y/2d0
         xSNnei(3,ind) = z/2d0
         vSNnei(1,ind) = x/r
         vSNnei(2,ind) = y/r
         vSNnei(3,ind) = z/r
         !indall(i+(j-1)*4+(k-1)*4*4) = ind
      endif
   enddo
   enddo
   enddo

   !=======================================================
   ! For careful refinement (for resolved feedback)
   !=======================================================
   ncshell  = (1+2*nshell_resolve)
   ncshell3 = ncshell**3

   allocate(xrefnei(1:3,1:ncshell3))
   allocate(irefnei(1:ncshell3))
   ! notice that xrefnei,irefnei have a different ordering than the rest of variables
   ! xrefnei <- irefnei <- mrefnei
   allocate(mrefnei(1:ncshell3))
   allocate(mzrefnei(1:ncshell3))
   allocate(icellnei(1:ncshell3))
   allocate(lrefnei(1:ncshell3))
   allocate(ltmp   (1:ncshell3))

   allocate(nrefnei_ring (0:nshell_resolve))
   allocate(mrefnei_ring (0:nshell_resolve))
   allocate(mzrefnei_ring(0:nshell_resolve))

   allocate(icommr (1:ncshell3))

   xrefnei=0;irefnei=0;nrefnei_ring=1;ltmp=0
   mrefnei=0d0;mrefnei_ring=0d0;icellnei=0;lrefnei=0

   ind=1
   do iring=0,nshell_resolve
      do k=-iring,iring
      do j=-iring,iring
      do i=-iring,iring

         i2=i+nshell_resolve  ! [0,nshell_resolve*2]
         j2=j+nshell_resolve
         k2=k+nshell_resolve

         indall=i2+(j2+k2*ncshell)*ncshell+1

         if(ltmp(indall)==0)then
           ltmp(indall)=1
           irefnei(ind)=indall
           nrefnei_ring(iring)=ind   ! just keep track of the last index
           ind=ind+1
         endif

         xrefnei(1,indall)=i ! [-3,3]
         xrefnei(2,indall)=j ! [-3,3]
         xrefnei(3,indall)=k ! [-3,3]

      end do
      end do
      end do

   end do

   deallocate(ltmp)

end subroutine init_mechanical
!####################################################################
!####################################################################
!####################################################################
