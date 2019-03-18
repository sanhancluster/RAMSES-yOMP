subroutine init_sink
  use amr_commons
  use pm_commons
  use clfind_commons
#ifdef RT
  use rt_parameters,only: rt_AGN, nGroups
#endif
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer,parameter::tag=1112,tag2=1113
  integer::dummy_io,info2
#endif
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer::idim,ilevel
  integer::isink,itype
  integer::ilun
  integer::nsinkold,nsinknew
  real(dp)::xx1,xx2,xx3,vv1,vv2,vv3,mm1,ll1,ll2,ll3
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::isp
  logical::ic_sink=.false.
  character(LEN=80)::filename
  character(LEN=80)::fileloc
  character(LEN=5)::nchar,ncharcpu

  allocate(total_volume(1:nsinkmax))
  allocate(wdens(1:nsinkmax))
  allocate(wvol(1:nsinkmax))
  allocate(wmom(1:nsinkmax,1:ndim))
  allocate(wc2(1:nsinkmax))
  allocate(wdens_new(1:nsinkmax))
  allocate(wvol_new(1:nsinkmax))
  allocate(wmom_new(1:nsinkmax,1:ndim))
  allocate(wc2_new(1:nsinkmax))
  allocate(msink(1:nsinkmax))
  allocate(msink_new(1:nsinkmax))
  allocate(msink_all(1:nsinkmax))
  allocate(idsink(1:nsinkmax))
  ! Important to set nindsink
  idsink=0
  allocate(idsink_new(1:nsinkmax))
  allocate(idsink_all(1:nsinkmax))
  allocate(tsink(1:nsinkmax))
  allocate(tsink_new(1:nsinkmax))
  allocate(tsink_all(1:nsinkmax))
  allocate(vsink(1:nsinkmax,1:ndim))
  allocate(xsink(1:nsinkmax,1:ndim))
  allocate(vsink_new(1:nsinkmax,1:ndim))
  allocate(vsink_all(1:nsinkmax,1:ndim))
  allocate(xsink_new(1:nsinkmax,1:ndim))
  allocate(xsink_all(1:nsinkmax,1:ndim))
  allocate(dMBHoverdt(1:nsinkmax))
  allocate(dMEdoverdt(1:nsinkmax))
  allocate(r2sink(1:nsinkmax))
  allocate(r2k(1:nsinkmax))
  allocate(v2sink(1:nsinkmax))
  allocate(c2sink(1:nsinkmax))
  allocate(v2sink_new(1:nsinkmax))
  allocate(c2sink_new(1:nsinkmax))
  allocate(v2sink_all(1:nsinkmax))
  allocate(c2sink_all(1:nsinkmax))
  allocate(weighted_density(1:nsinkmax,1:nlevelmax))
  allocate(weighted_volume (1:nsinkmax,1:nlevelmax))
  allocate(weighted_momentum(1:nsinkmax,1:nlevelmax,1:ndim))
  allocate(weighted_c2 (1:nsinkmax,1:nlevelmax))
  allocate(oksink_new(1:nsinkmax))
  allocate(oksink_all(1:nsinkmax))
  allocate(jsink(1:nsinkmax,1:ndim))
  allocate(jsink_new(1:nsinkmax,1:ndim))
  allocate(jsink_all(1:nsinkmax,1:ndim))
  allocate(dMBH_coarse    (1:nsinkmax))
  allocate(dMEd_coarse    (1:nsinkmax))
  allocate(dMsmbh         (1:nsinkmax))
  allocate(Esave          (1:nsinkmax))
  allocate(Efeed          (1:nsinkmax))  !Feedback energy during a given feedback event, collected to be written to file
  allocate(dMBH_coarse_new(1:nsinkmax))
  allocate(dMEd_coarse_new(1:nsinkmax))
  allocate(dMsmbh_new     (1:nsinkmax))
  allocate(Esave_new      (1:nsinkmax))
  allocate(Efeed_new      (1:nsinkmax))
  allocate(dMBH_coarse_all(1:nsinkmax))
  allocate(dMEd_coarse_all(1:nsinkmax))
  allocate(dMsmbh_all     (1:nsinkmax))
  allocate(Esave_all      (1:nsinkmax))
  allocate(Efeed_all      (1:nsinkmax))
  allocate(sink_stat      (1:nsinkmax,levelmin:nlevelmax,1:ndim*2+1))
  allocate(sink_stat_all  (1:nsinkmax,levelmin:nlevelmax,1:ndim*2+1))
  allocate(v_avgptr(1:nsinkmax))
  allocate(c_avgptr(1:nsinkmax))
  allocate(d_avgptr(1:nsinkmax))
  allocate(spinmag(1:nsinkmax),bhspin(1:nsinkmax,1:ndim))
  allocate(spinmag_new(1:nsinkmax),bhspin_new(1:nsinkmax,1:ndim))
  allocate(spinmag_all(1:nsinkmax),bhspin_all(1:nsinkmax,1:ndim))
  allocate(eps_sink(1:nsinkmax))

  ! Initialize all to 0
  total_volume=0; wdens=0; wvol=0; wmom=0; wc2=0; wdens_new=0; wvol_new=0; wmom_new=0
  wc2_new=0; msink=0; msink_new=0; msink_all=0; idsink=0; idsink_new=0; idsink_all=0
  tsink=0; tsink_new=0; tsink_all=0; vsink=0; xsink=0; vsink_new=0; vsink_all=0; xsink_new=0
  xsink_all=0; dMBHoverdt=0; dMEdoverdt=0; r2sink=0; r2k=0; v2sink=0; c2sink=0; v2sink_new=0;
  c2sink_new=0; v2sink_all=0; c2sink_all=0; weighted_density=0; weighted_volume =0; weighted_momentum=0
  weighted_c2 =0; oksink_new=0; oksink_all=0; jsink=0; jsink_new=0; jsink_all=0; dMBH_coarse=0
  dMEd_coarse=0; dMsmbh=0; Esave=0; Efeed=0; dMBH_coarse_new=0; dMEd_coarse_new=0; dMsmbh_new =0
  Esave_new=0; Efeed_new=0; dMBH_coarse_all=0; dMEd_coarse_all=0; dMsmbh_all =0; Esave_all=0
  Efeed_all=0; sink_stat=0; sink_stat_all=0; v_avgptr=0; c_avgptr=0; d_avgptr=0; spinmag=0
  spinmag_new=0; spinmag_all=0; eps_sink=0

#ifdef RT
  ! Only allocate some variables if AGNRT is activated
  if (rt_AGN) then
     allocate(LAGN_coarse     (1:nsinkmax))  ! AGNRT
     allocate(dMeff_coarse    (1:nsinkmax))  ! AGNRT
     allocate(dMeff_coarse_new(1:nsinkmax))  ! AGNRT
     allocate(dMeff_coarse_all(1:nsinkmax))  ! AGNRT
     allocate(lumfrac_AGN(1:nsinkmax,1:nGroups))
     LAGN_coarse=0; dMeff_coarse=0; dMeff_coarse_new=0; dMeff_coarse_all=0; lumfrac_AGN=0
  end if
#endif
  allocate(rg_scale(1:nsinkmax))
  rg_scale=0
  ! Dynamical friction from particles (HP)
  if (drag_part) then
     ! quantities needed to compute DF at
     ! measured around each BHs (1:nsinkmax)
     ! for each level (levelmin:nlevelmax)
     ! for stars and DM (1:2).
     allocate(v_background(1:nsinkmax, 1:ndim, 1:2))
     allocate(n_background(1:nsinkmax, 1:2))
     allocate(m_background(1:nsinkmax, 1:2))
     allocate(mass_lowspeed_background(1:nsinkmax, 1:2))
     allocate(fact_fast_background(1:nsinkmax, 1:2))
     allocate(vrel_sink(1:nsinkmax, 1:ndim, 1:2))
     allocate(vrel_sink_norm(1:nsinkmax, 1:2))
     allocate(v_DFnew(1:nsinkmax, 1:ndim, 1:2))
     allocate(v_DFall(1:nsinkmax, 1:ndim, 1:2))
     allocate(v_DF(1:nsinkmax, levelmin:nlevelmax, 1:ndim, 1:2))
     allocate(v_DFnew_all(1:nsinkmax, levelmin:nlevelmax, 1:ndim, 1:2))
     allocate(mass_DFnew(1:nsinkmax, 1:2))
     allocate(mass_DFall(1:nsinkmax, 1:2))
     allocate(mass_DF(1:nsinkmax, levelmin:nlevelmax, 1:2))
     allocate(mass_DFnew_all(1:nsinkmax, levelmin:nlevelmax, 1:2))
     allocate(fact_fastnew(1:nsinkmax, 1:2))
     allocate(fact_fastall(1:nsinkmax, 1:2))
     allocate(fact_fast(1:nsinkmax, levelmin:nlevelmax, 1:2))
     allocate(fact_fastnew_all(1:nsinkmax, levelmin:nlevelmax, 1:2))
     allocate(mass_lowspeednew(1:nsinkmax, 1:2))
     allocate(mass_lowspeedall(1:nsinkmax, 1:2))
     allocate(mass_lowspeed(1:nsinkmax, levelmin:nlevelmax, 1:2))
     allocate(mass_lowspeednew_all(1:nsinkmax, levelmin:nlevelmax, 1:2))
     allocate(n_partnew(1:nsinkmax, 1:2))
     allocate(n_partall(1:nsinkmax, 1:2))
     allocate(n_part(1:nsinkmax, levelmin:nlevelmax, 1:2))
     allocate(n_partnew_all(1:nsinkmax, levelmin:nlevelmax, 1:2))

     ! Set initial value to 0
     v_background=0; n_background=0; m_background=0; mass_lowspeed_background=0; fact_fast_background=0
     vrel_sink=0; vrel_sink_norm=0; v_DFnew=0; v_DFall=0; v_DF=0; v_DFnew_all=0; mass_DFnew=0
     mass_DFall=0; mass_DF=0; mass_DFnew_all=0; fact_fastnew=0; fact_fastall=0; fact_fast=0;
     fact_fastnew_all=0; mass_lowspeednew=0; mass_lowspeedall=0; mass_lowspeed=0; mass_lowspeednew_all=0
     n_partnew=0; n_partall=0; n_part=0; n_partnew_all=0

     ! This table is needed to properly fill the above one
     ! after a merger
     allocate(most_massive_sink(1:nsinkmax))
     most_massive_sink=0

     ! This tables is here to speed up the computation
     allocate(sink_cell(1:nsinkmax))
     sink_cell=0
  end if

  eps_sink=0.057190958d0

  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  call compute_ncloud_sink

  if(nrestart>0)then
     ilun=4*ncpu+myid+10
     call title(nrestart,nchar)

     if(IOGROUPSIZEREP>0)then
        call title(((myid-1)/IOGROUPSIZEREP)+1,ncharcpu)
        fileloc='output_'//TRIM(nchar)//'/group_'//TRIM(ncharcpu)//'/sink_'//TRIM(nchar)//'.out'
     else
        fileloc='output_'//TRIM(nchar)//'/sink_'//TRIM(nchar)//'.out'
     endif


     call title(myid,nchar)
     fileloc=TRIM(fileloc)//TRIM(nchar)

     inquire(file=fileloc, exist=ic_sink)

     nsink=0
     nindsink=0
     if (ic_sink) then

        ! Wait for the token
#ifndef WITHOUTMPI
        if(IOGROUPSIZE>0) then
           if (mod(myid-1,IOGROUPSIZE)/=0) then
              call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
                   & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
           end if
        endif
#endif

        open(unit=ilun,file=fileloc,form='unformatted')
        rewind(ilun)
        read(ilun)nsink
        read(ilun)nindsink

        if(nsink>0)then
           allocate(xdp(1:nsink))
           allocate(isp(1:nsink))
           read(ilun)isp
           idsink(1:nsink)=isp
           ! Important for the indexation of sinks
           nindsink=MAXVAL(idsink)
           deallocate(isp)
           read(ilun)xdp
           msink(1:nsink)=xdp
           do idim=1,ndim
              read(ilun)xdp
              xsink(1:nsink,idim)=xdp
           end do
           do idim=1,ndim
              read(ilun)xdp
              vsink(1:nsink,idim)=xdp
           end do
           read(ilun)xdp
           tsink(1:nsink)=xdp
           read(ilun)xdp
           dMsmbh(1:nsink)=xdp
           read(ilun)xdp
           dMBH_coarse(1:nsink)=xdp
           read(ilun)xdp
           dMEd_coarse(1:nsink)=xdp
           read(ilun)xdp
           Esave(1:nsink)=xdp
           do idim=1,ndim
              read(ilun)xdp
              jsink(1:nsink,idim)=xdp
           end do
           do idim=1,ndim
              read(ilun)xdp
              bhspin(1:nsink,idim)=xdp
           end do
           read(ilun)xdp
           spinmag(1:nsink)=xdp
           read(ilun)xdp
           eps_sink(1:nsink)=xdp
           do idim=1,ndim*2+1
              do ilevel=levelmin,nlevelmax
                 read(ilun)xdp
                 sink_stat(1:nsink,ilevel,idim)=xdp
              enddo
           enddo
           ! AGNRT
#ifdef RT
           if(rt_AGN) then
              ! Read AGN radiation to be released
              read(ilun)xdp
              LAGN_coarse(1:nsink)=xdp
              !! Uncomment this to restart from a non-AGNRT simulation
              !! LAGN_coarse(1:nsink) = 0.d0
           endif
#endif
           !/AGNRT
           deallocate(xdp)
        end if
        close(ilun)
        ! Send the token
#ifndef WITHOUTMPI
        if(IOGROUPSIZE>0) then
           if(mod(myid,IOGROUPSIZE)/=0 .and.(myid.lt.ncpu))then
              dummy_io=1
              call MPI_SEND(dummy_io,1,MPI_INTEGER,myid-1+1,tag, &
                   & MPI_COMM_WORLD,info2)
           end if
        endif
#endif

     end if
  end if

  if (nrestart>0)then
     nsinkold=nsink
     if(TRIM(initfile(levelmin)).NE.' ')then
        filename=TRIM(initfile(levelmin))//'/ic_sink_restart'
     else
        filename='ic_sink_restart'
     end if
     INQUIRE(FILE=filename, EXIST=ic_sink)
     if (myid==1)write(*,*)'Looking for file ic_sink_restart: ',filename
  else
     nsink=0
     nindsink=0
     nsinkold=0
     if(TRIM(initfile(levelmin)).NE.' ')then
        filename=TRIM(initfile(levelmin))//'/ic_sink'
     else
        filename='ic_sink'
     end if
     INQUIRE(FILE=filename, EXIST=ic_sink)
     if (myid==1)write(*,*)'Looking for file ic_sink: ',filename
  end if

  if (ic_sink)then
#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if (mod(myid-1,IOGROUPSIZE)/=0) then
           dummy_io=1
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag2,&
                &MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif
     open(10, file=filename, form='formatted')
     read(10, *) nsinknew
     nsink = nsink + nsinknew

     ! FORMAT FOR THE SINKFILE
     ! 1st line:  Number of sink particles to be read from file
     ! Every other line: sink mass [Msun], x , y, z, vx [km/s], vy [km/s], vz [km/s], lx, ly, lz
     ! The positions x,y and z are given in units of code_length, so boxlen/2 is the middle of the box

     !EXAMPLE for a case with 1 10Msun black hole at rest in the centre of the box
     ! boxlen=10
     ! ------- Beginning of file ic_sink
     ! 1
     ! 10,5,5,5,0,0,0,0,0,0
     ! ------ End of file ic_sink

     do isink=1, nsinknew
        read(10,*,end=102) mm1, xx1, xx2, xx3, vv1, vv2, vv3, ll1, ll2, ll3
        nindsink=nindsink+1
        idsink(nindsink)=nindsink
        msink(nindsink)=mm1*2.d33/(scale_d*scale_l**3)
        xsink(nindsink,1)=xx1
        xsink(nindsink,2)=xx2
        xsink(nindsink,3)=xx3
        vsink(nindsink,1)=vv1*1.d5/scale_v
        vsink(nindsink,2)=vv2*1d5/scale_v
        vsink(nindsink,3)=vv3*1d5/scale_v
        jsink(nindsink,1)=ll1
        jsink(nindsink,2)=ll2
        jsink(nindsink,3)=ll3
        tsink(nindsink)=t
        dMsmbh(nindsink)=0d0
        Esave(nindsink)=0d0
        dMBH_coarse(nindsink)=0d0
        dMEd_coarse(nindsink)=0d0
        spinmag(nindsink)=0d0
        eps_sink(nindsink)=0d0
        bhspin(nindsink,1:3)=jsink(nindsink,1:3)
        ! AGNRT
#ifdef RT
        if (rt_AGN) LAGN_coarse(nindsink) = 0.d0
#endif
        !/AGNRT
        ! Particles dynamical friction (HP)
        ! itype=1 for stars, itype=2 for DM
        if (drag_part) then
           do itype = 1, 2
              v_DF(nindsink, levelmin:nlevelmax, 1, itype) = vsink(nindsink, 1)
              v_DF(nindsink, levelmin:nlevelmax, 2, itype) = vsink(nindsink, 2)
              v_DF(nindsink, levelmin:nlevelmax, 3, itype) = vsink(nindsink, 3)
              mass_DF(nindsink, levelmin:nlevelmax, itype) = 0d0
              mass_lowspeed(nindsink, levelmin:nlevelmax, itype) = 0d0
              fact_fast(nindsink, levelmin:nlevelmax, itype) = 0d0
              n_part(nindsink, levelmin:nlevelmax, itype) = 0
           end do
        end if
        !/Particles dynamical friction (HP)
     end do
102  continue
     close(10)

#ifndef WITHOUTMPI
     if(IOGROUPSIZE>0) then
        if ((mod(myid-1,IOGROUPSIZE)/=0) .and. myid .lt. ncpu) then
           dummy_io=1
           call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag2,&
                &MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
        end if
     endif
#endif

     if (myid==1.and.nsink-nsinkold>0)then
        write(*,*)'sinks read from file '//filename
        write(*,"(999(A))")'   Id           M [code mass]          x             y         &
             &z            vx            vy            vz            lx            &
             &ly            lz       '
        write(*,"(999(A))") ' ==========================================================&
             &==========================================================================================='
        do isink=nsinkold+1,nsink
           write(*,'(I8,2X,10(2X,E12.5))')idsink(isink),msink(isink),xsink(isink,1:ndim),&
                vsink(isink,1:ndim), jsink(isink,1:ndim)
        end do
     end if
  end if

end subroutine init_sink
!################################################################
!################################################################
!################################################################
!################################################################
subroutine compute_ncloud_sink
  use amr_commons, only:dp,myid
  use pm_commons, only:ir_cloud,ncloud_sink
  real(dp)::xx,yy,zz,rr
  integer::ii,jj,kk

  ! Compute number of cloud particles
  ncloud_sink=0
  do kk=-2*ir_cloud,2*ir_cloud
     zz=dble(kk)/2.0
     do jj=-2*ir_cloud,2*ir_cloud
        yy=dble(jj)/2.0
        do ii=-2*ir_cloud,2*ir_cloud
              xx=dble(ii)/2.0
           rr=sqrt(xx*xx+yy*yy+zz*zz)
           if(rr<=dble(ir_cloud))ncloud_sink=ncloud_sink+1
        end do
     end do
  end do

  if(myid==1)write(*,*)"Number of cloud particles per sink:",ncloud_sink
end subroutine compute_ncloud_sink
