!####################################################################
!####################################################################
!####################################################################
module mechanical_commons
   use amr_commons, ONLY:dp,ndim,ncpu,nchem
   integer(kind=4),parameter::nvarSN=13+nchem

   ! Important Note: SN stands for SN cell, not SN particle (SNp)

   ! Array to define neighbors for SN
   integer, parameter::nSNnei=48  ! number of neighboring cells to deposit mass/momentum/energy
   real(dp),parameter::nSNcen=4   ! number of cells corresponding to the central cell to deposit mass
   real(dp),dimension(1:3,1:nSNnei)::xSNnei
   real(dp),dimension(1:3,1:nSNnei)::vSNnei
   real(dp)::f_LOAD,f_LEFT,f_ESN,f_PCAN
   ! SN cells that are needed to communicate across different CPUs
   ! note that SNs are processed on a cell-by-cell basis
   ! hence the central position is taken as the central leaf cell
   integer ::ncomm_SN   ! the number of cells to be communicated (specific to myid)
   integer ::uidSN_comm ! the unique id of SN comm

  ! momentum input
  ! p_sn = A_SN*nH**(alpha)*ESN**(beta)*ZpSN**(gamma)
  ! ex) Thornton et al.
  !     A_SN = 3e5, alphaN = -2/17, beta = 16/17, gamma = -0.14
  ! ex) Kim & Ostriker (2015) uniform case
  !     A_SN = 2.17e5, alpha = -0.13, beta = 0.93
   real(dp),parameter:: expE_SN=+16d0/17d0
   real(dp),parameter:: expZ_SN=-0.14
   real(dp),parameter:: expN_SN_boost=-0.15
   real(dp),parameter:: expE_SN_boost=16d0/17d0

#ifndef WITHOUTMPI
  ! common lists for SNe across different cpus
   integer, parameter ::ncomm_max = 50000
   integer ,dimension(1:ncomm_max)::iSN_comm  ! cpu list
   integer ,dimension(1:ncomm_max)::idSN_comm  ! id of SNe in each cpu
   real(dp),dimension(1:ncomm_max)::nSN_comm          ! number of SNe
   real(dp),dimension(1:ncomm_max)::mSN_comm          ! gas mass of SNe
   real(dp),dimension(1:ncomm_max)::mZSN_comm         ! metal mass of SNe
   real(dp),dimension(1:ncomm_max)::mloadSN_comm      ! ejecta mass + gas entrained
   real(dp),dimension(1:ncomm_max)::eloadSN_comm      ! kinetic energy of ejecta + gas entrained
   real(dp),dimension(1:ncomm_max)::mZloadSN_comm     ! metals ejected + entrained
   real(dp),dimension(1:3,1:ncomm_max)::xSN_comm      ! pos of SNe host cell (leaf)
   real(dp),dimension(1:3,1:ncomm_max)::pSN_comm      ! total momentum of total SNe in each leaf cell
   real(dp),dimension(1:3,1:ncomm_max)::ploadSN_comm  ! momentum from original star + gas entrained
   real(dp),dimension(1:ncomm_max)::floadSN_comm      ! fraction of gas to be loaded from the central cell
   integer,dimension(:,:),allocatable::icpuSN_comm,icpuSN_comm_mpi
   integer,dimension(:)  ,allocatable::ncomm_SN_cpu,ncomm_SN_mpi
   real(dp),dimension(1:ncomm_max)::rSt_comm          ! Stromgren radius in pc
   real(dp),dimension(1:ncomm_max,1:nchem)::mchloadSN_comm     ! chemical species ejected + entrained
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

   ! For chemical abundance due to SN II
   real(dp)::Zejecta_chem_II(1:nchem)

   ! For chemical abundance due to SN Ia
   real(dp)::mejecta_Ia
   real(dp)::Zejecta_chem_Ia(1:nchem)

end module
!####################################################################
!####################################################################
!####################################################################
subroutine init_mechanical
   use amr_commons
   use mechanical_commons
   use hydro_parameters, ONLY:ichem
   implicit none
   integer::i,j,k,ind,indall
   real(kind=dp)::x,y,z,r
   logical::ok
   integer,allocatable::ltmp(:)
   integer::ncshell,iring,i2,j2,k2,nrad,irad,ich
   character(len=2)::element_name


   !------------------------------------------------
   ! Warning messages
   !------------------------------------------------
   ok=.false.
   if(.not.metal.and.mechanical_feedback==2)then
      print *, '>>> mechanical Err: Please turn on metal'
      ok=.true.
   endif
   if(ok) call clean_stop

#ifndef WITHOUTMPI
   allocate(ncomm_SN_cpu(1:ncpu))
   allocate(ncomm_SN_mpi(1:ncpu))
#endif

   ! some parameters
   f_LOAD = nSNnei / dble(nSNcen + nSNnei)
   f_LEFT = nSNcen / dble(nSNcen + nSNnei)
   f_ESN  = 0.676   ! Blondin+(98) at t=trad
   f_PCAN = 0.9387  ! correction due to direct momentum cancellation
                    ! due to the assumption of 48 neighboring cells
                    ! even in the uniform case where there are 18 immediate neighbors

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

   if(.not.variable_yield_SNII)then ! assume solar case
      do ich=1,nchem
         element_name=chem_list(ich)
         select case (element_name)
            case ('H ')
               Zejecta_chem_II(ich) = 10.**(-0.30967822)
            case ('He')
               Zejecta_chem_II(ich) = 10.**(-0.40330181)
            case ('C ')
               Zejecta_chem_II(ich) = 10.**(-1.9626259)
            case ('N ')
               Zejecta_chem_II(ich) = 10.**(-2.4260355)
            case ('O ')
               Zejecta_chem_II(ich) = 10.**(-1.1213435)
            case ('Mg')
               Zejecta_chem_II(ich) = 10.**(-2.3706062)
            case ('Si')
               Zejecta_chem_II(ich) = 10.**(-2.0431845)
            case ('S ')
               Zejecta_chem_II(ich) = 10.**(-2.2964020)
            case ('Fe')
               Zejecta_chem_II(ich) = 10.**(-2.2126987)
            case default
               Zejecta_chem_II(ich)=0
         end select

      end do
   endif

  if (snIa) call init_snIa_yield

  ! This is not to double count SN feedback
  if (mechanical_feedback==2) f_w = -1


end subroutine init_mechanical
!####################################################################
!####################################################################
!####################################################################
subroutine init_snIa_yield
   use amr_commons
   use mechanical_commons, ONLY: mejecta_Ia, Zejecta_chem_Ia
   implicit none
   integer::ich
   real(kind=8)::yield_snIa(1:66)
   character(len=2)::element_name
!----------------------------------------------------------------------------
!  Nomoto et al. (1997,1984) W7 (carbon-deflagration model)
!            (better fit with Tycho observation)
!----------------------------------------------------------------------------
!       12C ,     13C ,   14N ,    15N ,    16O ,    17O ,    18O ,    19F ,
!       20Ne,     21Ne,   22Ne,    23Na,    24Mg,    25Mg,    26Mg,    27Al,
!       28Si,     29Si,   30Si,    31P ,    32S ,    33S ,    34S ,    36S ,
!       35Cl,     37Cl,   36Ar,    38Ar,    40Ar,    39K ,    41K ,    40Ca,
!       42Ca,     43Ca,   44Ca,    46Ca,    48Ca,    45Sc,    46Ti,    47Ti,
!       48Ti,     49Ti,   50Ti,    50V ,    51V ,    50Cr,    52Cr,    53Cr,
!       54Cr,     55Mn,   54Fe,    56Fe,    57Fe,    58Fe,    59Co,    58Ni,
!       60Ni,     61Ni,   62Ni,    64Ni,    63Cu,    65Cu,    64Zn,    66Zn,
!       67Zn,     68Zn
!----------------------------------------------------------------------------

   yield_snIa = (/ &  ! Msun per SN
    &4.83E-02,1.40E-06,1.16E-06,1.32E-09,1.43E-01,3.54E-08,8.25E-10,5.67E-10,&
    &2.02E-03,8.46E-06,2.49E-03,6.32E-05,8.50E-03,4.05E-05,3.18E-05,9.86E-04,&
    &1.50E-01,8.61E-04,1.74E-03,4.18E-04,8.41E-02,4.50E-04,1.90E-03,3.15E-07,&
    &1.34E-04,3.98E-05,1.49E-02,1.06E-03,1.26E-08,8.52E-05,7.44E-06,1.23E-02,&
    &3.52E-05,1.03E-07,8.86E-06,1.99E-09,7.10E-12,2.47E-07,1.71E-05,6.04E-07,&
    &2.03E-04,1.69E-05,1.26E-05,8.28E-09,5.15E-05,2.71E-04,5.15E-03,7.85E-04,&
    &1.90E-04,8.23E-03,1.04E-01,6.13E-01,2.55E-02,9.63E-04,1.02E-03,1.28E-01,&
    &1.05E-02,2.51E-04,2.66E-03,1.31E-06,1.79E-06,6.83E-07,1.22E-05,2.12E-05,&
    &1.34E-08,1.02E-08/)

   mejecta_Ia=sum(yield_snIa)
   do ich=1,nchem
      element_name=chem_list(ich)
      select case (element_name)
         case ('C ')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(1:2))/mejecta_Ia
         case ('N ')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(3:4))/mejecta_Ia
         case ('O ')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(5:7))/mejecta_Ia
         case ('Mg')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(13:15))/mejecta_Ia
         case ('Si')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(17:19))/mejecta_Ia
         case ('S ')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(21:24))/mejecta_Ia
         case ('Fe')
            Zejecta_chem_Ia(ich)=sum(yield_snIa(51:54))/mejecta_Ia
         case default
            Zejecta_chem_Ia(ich)=0d0
      end select
   enddo

end subroutine init_snIa_yield
