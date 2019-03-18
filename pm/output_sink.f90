subroutine backup_sink(filename, filename_desc)
  use amr_commons
  use pm_commons
  use dump_utils, only : dump_header_info, generic_dump, dim_keys
#ifdef RT
  use rt_parameters,only: rt_AGN
#endif
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer,parameter::tag=1135
  integer::dummy_io,info2
#endif

  character(LEN=*), intent(in)::filename, filename_desc

  integer::idim,i,ilevel, unit_out, unit_info
  character(LEN=80)::fileloc
  character(LEN=5)::nchar
  character(len=100) :: field_name
  real(dp),allocatable,dimension(:)::xdp
  integer,allocatable,dimension(:)::ii

  logical :: dump_info
  integer :: ivar

  if(.not. sink) return

  if(verbose) write(*,*)'Entering backup_sink'

  call title(myid,nchar)
  fileloc=TRIM(filename)//TRIM(nchar)

  ! Set ivar to 1 for first variable
  ivar = 1

  ! Wait for the token
#ifndef WITHOUTMPI
  if(IOGROUPSIZE>0) then
     if (mod(myid-1,IOGROUPSIZE)/=0) then
        call MPI_RECV(dummy_io,1,MPI_INTEGER,myid-1-1,tag,&
             & MPI_COMM_WORLD,MPI_STATUS_IGNORE,info2)
     end if
  endif
#endif

  fileloc = TRIM(filename) // TRIM(nchar)
  if (myid == 1) then
     open(newunit=unit_info, file=trim(filename_desc), form='formatted')
     call dump_header_info(unit_info)
     dump_info = .true.
  else
     dump_info = .false.
  end if

  open(newunit=unit_out,file=TRIM(fileloc),form='unformatted')
  rewind(unit_out)

  write(unit_out)nsink
  write(unit_out)nindsink
  if(nsink>0)then
     allocate(ii(1:nsink))
     ! Write identity sink
     do i=1,nsink
        ii(i)=idsink(i)
     end do
     call generic_dump("identity", ivar, ii, unit_out, dump_info, unit_info)
     deallocate(ii)
     allocate(xdp(1:nsink))
     ! Write mass
     do i=1,nsink
        xdp(i)=msink(i)
     end do
     call generic_dump("mass", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write position
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=xsink(i,idim)
        end do
        call generic_dump("position_"//dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
     enddo
     ! Write velocity
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=vsink(i,idim)
        end do
        call generic_dump("velocity_"//dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
     enddo
     ! Write time
     do i=1,nsink
        xdp(i)=tsink(i)
     end do
     call generic_dump("birth_time", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write real accretion
     do i=1,nsink
        xdp(i)=dMsmbh(i)
     enddo
     call generic_dump("dMsmbh", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write Bondi accretion
     do i=1,nsink
        xdp(i)=dMBH_coarse(i)
     end do
     call generic_dump("dMBH_coarse", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write Eddington accretion
     do i=1,nsink
        xdp(i)=dMEd_coarse(i)
     end do
     call generic_dump("dMEd_coarse", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write Esave
     do i=1,nsink
        xdp(i)=Esave(i)
     end do
     call generic_dump("Esave", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write gas spin axis
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=jsink(i,idim)
        end do
        call generic_dump("jsink_" // dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
     enddo
     ! Write BH spin axis
     do idim=1,ndim
        do i=1,nsink
           xdp(i)=bhspin(i,idim)
        end do
        call generic_dump("spin_" // dim_keys(idim), ivar, xdp, unit_out, dump_info, unit_info)
     enddo
     ! Write BH spin amplitude (signed)
     do i=1,nsink
        xdp(i)=spinmag(i)
     end do
     call generic_dump("spin_magnitude", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write BH efficiency
     do i=1,nsink
        xdp(i)=eps_sink(i)
     end do
     call generic_dump("eps_sink", ivar, xdp, unit_out, dump_info, unit_info)
     ! Write sink_stat
     do idim=1,ndim*2+1
        do ilevel=levelmin,nlevelmax
           do i=1,nsink
              xdp(i)=sink_stat(i,ilevel,idim)
           end do
           write(field_name, "('sink_stat_', i0.2, '_', i0.2)") ilevel, idim
           call generic_dump(field_name, ivar, xdp, unit_out, dump_info, unit_info)
        enddo
     enddo
     ! AGNRT
#ifdef RT
     if(rt_AGN) then
        ! Write AGN radiation to be released
        do i=1,nsink
           xdp(i)=LAGN_coarse(i)
        end do
        call generic_dump("LAGN_coarse", ivar, xdp, unit_out, dump_info, unit_info)
     endif
#endif
     !/AGNRT

     deallocate(xdp)
  endif
  close(unit_out)
  if (myid == 1) close(unit_info)

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

end subroutine backup_sink



subroutine output_sink(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename

  integer::isink
  integer::nx_loc,ilun
  real(dp)::scale,dx_min
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  character(LEN=80)::fileloc

  if(verbose)write(*,*)'Entering output_sink'

  ilun=myid+10

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_m=scale_d*scale_l**3d0
  nx_loc=(icoarse_max-icoarse_min+1)
  scale=boxlen/dble(nx_loc)
  dx_min=scale*0.5D0**(nlevelmax-nlevelsheld)/aexp

  if(verbose)write(*,*)'Entering output_sink'

  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',status='replace')
  !======================
  ! Write sink properties
  !======================
  write(ilun,*)'Number of sink = ',nsink

  write(ilun,'(" ================================================================================================================================== ")')
  write(ilun,'("        Id       Mass(Msol)             x                y                z               vx               vy               vz      ")')
  write(ilun,'(" ================================================================================================================================== ")')

  do isink=1,nsink
     write(ilun,'(I10,7(2X,E15.7))')idsink(isink),msink(isink)*scale_m/2d33,xsink(isink,1:ndim),vsink(isink,1:ndim)
  end do
  write(ilun,'(" ================================================================================================================================== ")')
  close(ilun)

end subroutine output_sink

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine output_sink_csv(filename)
  use amr_commons
  use pm_commons
  implicit none
  character(LEN=80)::filename,fileloc

  integer::ilun,isink

  if(verbose)write(*,*)'Entering output_sink_csv'

  ilun=2*ncpu+myid+10

  fileloc=TRIM(filename)
  open(unit=ilun,file=TRIM(fileloc),form='formatted',status='replace', recl=500)
  !======================
  ! Write sink properties
  !======================
  do isink=1,nsink
     write(ilun,'(I10,9(A1,ES20.10))')idsink(isink),',',msink(isink),&
          ',',xsink(isink,1),',',xsink(isink,2),',',xsink(isink,3),&
          ',',vsink(isink,1),',',vsink(isink,2),',',vsink(isink,3),&
          ',',t-tsink(isink),',',dMBHoverdt(isink)
  end do

  close(ilun)

end subroutine output_sink_csv
