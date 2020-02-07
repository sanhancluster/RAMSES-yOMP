subroutine parallel_link_all(head,ngrid,head_thr,ngrid_thr)
  use amr_commons
  use pm_commons
  implicit none
  !----------------------------------------------------------------------
  ! This routine interates over grid and particle linked list and divide
  ! them by number of threads weighted by specific type of particles and
  ! put head information for each thread.
  !----------------------------------------------------------------------
  integer :: head,ngrid
  integer, dimension(1:nthr) :: head_thr,ngrid_thr
  integer, dimension(1:ngrid) :: npart_grid

  integer :: igrid,jgrid,ipart,jpart,next_part,npart1,npart_tot,npart_now,ngrid_now,ngrid_res,npart_thr,ithr

  if(ngrid==0) then
     do ithr=1,nthr
        ngrid_thr(ithr)=0
     end do
     return
  end if

  npart_grid=0
  ! Count the total number of particles
  igrid=head
  npart_tot=0
  do jgrid=1,ngrid
     npart1=numbp(igrid)
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle  <---- Very important !!!
          next_part=nextp(ipart)
          ! Skip tracers (except "classic" tracers)
          npart_tot=npart_tot+1
          npart_grid(jgrid)=npart_grid(jgrid)+1
          ! End MC Tracer patch
          ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)
  end do

  ! Calculate rough number of particles per grids
  npart_thr=npart_tot/nthr

  igrid=head
  head_thr(1)=igrid
  npart_now=0
  ngrid_now=0
  ngrid_res=ngrid
  do ithr=1,nthr
     ngrid_thr(ithr)=0
  end do
  ithr=1

  do jgrid=1,ngrid
     if(ithr == nthr) then
        ngrid_thr(ithr)=ngrid_res
        return
     end if
     ngrid_now=ngrid_now+1
     npart_now=npart_now+npart_grid(jgrid)

     igrid=next(igrid)
     ! Save state and move to next thread
     if(npart_now>=npart_thr) then

        ngrid_thr(ithr)=ngrid_now
        ngrid_res=ngrid_res-ngrid_now
        ngrid_now=0

        ithr=ithr+1
        head_thr(ithr)=igrid

        ! Correct the weight
        npart_tot=npart_tot-npart_now
        npart_thr=npart_tot/(nthr-ithr+1)
        npart_now=0
     end if
  end do
end subroutine parallel_link_all
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine parallel_link_notracer(head,ngrid,head_thr,ngrid_thr)
  use amr_commons
  use pm_commons
  implicit none
  !----------------------------------------------------------------------
  ! This routine interates over grid and particle linked list and divide
  ! them by number of threads weighted by specific type of particles and
  ! put head information for each thread.
  !----------------------------------------------------------------------
  integer :: head,ngrid
  integer, dimension(1:nthr) :: head_thr,ngrid_thr
  integer, dimension(1:ngrid) :: npart_grid

  integer :: igrid,jgrid,ipart,jpart,next_part,npart1,npart_tot,npart_now,ngrid_now,ngrid_res,npart_thr,ithr

  if(ngrid==0) then
     do ithr=1,nthr
        ngrid_thr(ithr)=0
     end do
     return
  end if

  npart_grid=0
  ! Count the total number of particles
  igrid=head
  npart_tot=0
  do jgrid=1,ngrid
     npart1=numbp(igrid)
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle  <---- Very important !!!
          next_part=nextp(ipart)
          ! Skip tracers (except "classic" tracers)
          if (.not. (MC_tracer .and. is_tracer(typep(ipart)))) then
             npart_tot=npart_tot+1
             npart_grid(jgrid)=npart_grid(jgrid)+1
          end if
          ! End MC Tracer patch
          ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)
  end do

  ! Calculate rough number of particles per grids
  npart_thr=npart_tot/nthr

  igrid=head
  head_thr(1)=igrid
  npart_now=0
  ngrid_now=0
  ngrid_res=ngrid
  do ithr=1,nthr
     ngrid_thr(ithr)=0
  end do
  ithr=1

  do jgrid=1,ngrid
     if(ithr == nthr) then
        ngrid_thr(ithr)=ngrid_res
        return
     end if
     ngrid_now=ngrid_now+1
     npart_now=npart_now+npart_grid(jgrid)

     igrid=next(igrid)
     ! Save state and move to next thread
     if(npart_now>=npart_thr) then

        ngrid_thr(ithr)=ngrid_now
        ngrid_res=ngrid_res-ngrid_now
        ngrid_now=0

        ithr=ithr+1
        head_thr(ithr)=igrid

        ! Correct the weight
        npart_tot=npart_tot-npart_now
        npart_thr=npart_tot/(nthr-ithr+1)
        npart_now=0
     end if
  end do
end subroutine parallel_link_notracer
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine parallel_link_tracer(head,ngrid,head_thr,ngrid_thr)
  use amr_commons
  use pm_commons
  implicit none
  !----------------------------------------------------------------------
  ! This routine interates over grid and particle linked list and divide
  ! them by number of threads weighted by specific type of particles and
  ! put head information for each thread.
  !----------------------------------------------------------------------
  integer :: head,ngrid
  integer, dimension(1:nthr) :: head_thr,ngrid_thr
  integer, dimension(1:ngrid) :: npart_grid

  integer :: igrid,jgrid,ipart,jpart,next_part,npart1,npart_tot,npart_now,ngrid_now,ngrid_res,npart_thr,ithr

  if(ngrid==0) then
     do ithr=1,nthr
        ngrid_thr(ithr)=0
     end do
     return
  end if

  npart_grid=0
  ! Count the total number of particles
  igrid=head
  npart_tot=0
  do jgrid=1,ngrid
     npart1=numbp(igrid)
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle  <---- Very important !!!
          next_part=nextp(ipart)
          ! Skip tracers (except "classic" tracers)
          if (is_tracer(typep(ipart))) then
             npart_tot=npart_tot+1
             npart_grid(jgrid)=npart_grid(jgrid)+1
          end if
          ! End MC Tracer patch
          ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)
  end do

  ! Calculate rough number of particles per grids
  npart_thr=npart_tot/nthr

  igrid=head
  head_thr(1)=igrid
  npart_now=0
  ngrid_now=0
  ngrid_res=ngrid
  do ithr=1,nthr
     ngrid_thr(ithr)=0
  end do
  ithr=1

  do jgrid=1,ngrid
     if(ithr == nthr) then
        ngrid_thr(ithr)=ngrid_res
        return
     end if
     ngrid_now=ngrid_now+1
     npart_now=npart_now+npart_grid(jgrid)

     igrid=next(igrid)
     ! Save state and move to next thread
     if(npart_now>=npart_thr) then

        ngrid_thr(ithr)=ngrid_now
        ngrid_res=ngrid_res-ngrid_now
        ngrid_now=0

        ithr=ithr+1
        head_thr(ithr)=igrid

        ! Correct the weight
        npart_tot=npart_tot-npart_now
        npart_thr=npart_tot/(nthr-ithr+1)
        npart_now=0
     end if
  end do
end subroutine parallel_link_tracer
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine parallel_link_startracer(head,ngrid,head_thr,ngrid_thr)
  use amr_commons
  use pm_commons
  implicit none
  !----------------------------------------------------------------------
  ! This routine interates over grid and particle linked list and divide
  ! them by number of threads weighted by specific type of particles and
  ! put head information for each thread.
  !----------------------------------------------------------------------
  integer :: head,ngrid
  integer, dimension(1:nthr) :: head_thr,ngrid_thr
  integer, dimension(1:ngrid) :: npart_grid

  integer :: igrid,jgrid,ipart,jpart,next_part,npart1,npart_tot,npart_now,ngrid_now,ngrid_res,npart_thr,ithr

  if(ngrid==0) then
     do ithr=1,nthr
        ngrid_thr(ithr)=0
     end do
     return
  end if

  npart_grid=0
  ! Count the total number of particles
  igrid=head
  npart_tot=0
  do jgrid=1,ngrid
     npart1=numbp(igrid)
     if(npart1>0)then
        ipart=headp(igrid)
        ! Loop over particles
        do jpart=1,npart1
          ! Save next particle  <---- Very important !!!
          next_part=nextp(ipart)
          ! Skip tracers (except "classic" tracers)
          if (is_star_tracer(typep(ipart))) then
             npart_tot=npart_tot+1
             npart_grid(jgrid)=npart_grid(jgrid)+1
          end if
          ! End MC Tracer patch
          ipart=next_part  ! Go to next particle
        end do
        ! End loop over particles
     end if
     igrid=next(igrid)
  end do

  ! Calculate rough number of particles per grids
  npart_thr=npart_tot/nthr

  igrid=head
  head_thr(1)=igrid
  npart_now=0
  ngrid_now=0
  ngrid_res=ngrid
  do ithr=1,nthr
     ngrid_thr(ithr)=0
  end do
  ithr=1

  do jgrid=1,ngrid
     if(ithr == nthr) then
        ngrid_thr(ithr)=ngrid_res
        return
     end if
     ngrid_now=ngrid_now+1
     npart_now=npart_now+npart_grid(jgrid)

     igrid=next(igrid)
     ! Save state and move to next thread
     if(npart_now>=npart_thr) then

        ngrid_thr(ithr)=ngrid_now
        ngrid_res=ngrid_res-ngrid_now
        ngrid_now=0

        ithr=ithr+1
        head_thr(ithr)=igrid

        ! Correct the weight
        npart_tot=npart_tot-npart_now
        npart_thr=npart_tot/(nthr-ithr+1)
        npart_now=0
     end if
  end do
end subroutine parallel_link_startracer
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cpu_parallel_link(ilevel,head_thr,ngrid_thr,icpu_thr,jgrid_thr,npart_thr,target)
  use amr_commons
  use pm_commons
  implicit none
  !----------------------------------------------------------------------
  ! This routine interates over grid and particle linked list and divide
  ! them by number of threads weighted by specific type of particles and
  ! put head information for each thread.
  ! Additionally, this routine is designed to run on multiple cpu loops,
  ! so it stores current status of the loop (icpu, jgrid)
  !----------------------------------------------------------------------
  integer :: ilevel,target
  integer, dimension(1:nthr) :: head_thr,ngrid_thr,icpu_thr,jgrid_thr,npart_thr
  integer, dimension(:), allocatable :: npart_grid

  integer :: icpu,ngrid_tot,igrid,jgrid,kgrid,ipart,jpart,next_part,npart1,npart_tot,npart_now,ngrid_now,weight_eql,ithr

  ngrid_tot=0
  do icpu=1,ncpu
     ngrid_tot=ngrid_tot+numbl(icpu,ilevel)
  end do

  allocate(npart_grid(ngrid_tot))

  if(ngrid_tot==0) then
     ngrid_thr=0
     return
  end if

!$omp parallel do
  do kgrid=1,ngrid_tot
     npart_grid(kgrid)=0 ! base weight 1
  end do

  if(target==1) then
     call get_partmap_notracer(ilevel,npart_grid,ngrid_tot,npart_tot)
  else if(target==2) then
     call get_partmap_cloud(ilevel,npart_grid,ngrid_tot,npart_tot)
  end if

  ! Calculate the rough number of particles per grids
  weight_eql=(npart_tot+ngrid_tot)/nthr

  head_thr(1)=0

  jgrid_thr(1)=1
  kgrid=1
  npart_now=0
  ngrid_now=0

  ngrid_thr=0
  npart_thr=0

  ithr=1
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     do jgrid=1,numbl(icpu,ilevel)
        if(head_thr(1)==0) then
           head_thr(1)=igrid
           icpu_thr(1)=icpu
        end if
        if(ithr == nthr) then
           ngrid_thr(ithr)=ngrid_tot
           npart_thr(ithr)=npart_tot
           return
        end if
        ngrid_now=ngrid_now+1
        npart_now=npart_now+npart_grid(kgrid)
        kgrid=kgrid+1
        igrid=next(igrid)
        if(npart_now+ngrid_now>=weight_eql) then
           ! Save state and move to the next thread
           ngrid_thr(ithr)=ngrid_now
           npart_thr(ithr)=npart_now
           ngrid_tot=ngrid_tot-ngrid_now
           ngrid_now=0

           ithr=ithr+1
           head_thr(ithr)=igrid
           icpu_thr(ithr)=icpu
           jgrid_thr(ithr)=jgrid+1

           ! Correct the weight
           npart_tot=npart_tot-npart_now
           weight_eql=(npart_tot+ngrid_tot)/(nthr-ithr+1)
           npart_now=0
        end if
     end do
  end do
end subroutine cpu_parallel_link
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cpu_parallel_link_all(ilevel,head_thr,ngrid_thr,icpu_thr,jgrid_thr,npart_thr)
  use amr_commons
  use pm_commons
  implicit none
  !----------------------------------------------------------------------
  ! This routine interates over grid and particle linked list and divide
  ! them by number of threads weighted by specific type of particles and
  ! put head information for each thread.
  ! Additionally, this routine is designed to run on multiple cpu loops,
  ! so it stores current status of the loop (icpu, jgrid)
  !----------------------------------------------------------------------
  integer :: ilevel
  integer, dimension(1:nthr) :: head_thr,ngrid_thr,icpu_thr,jgrid_thr,npart_thr

  integer :: icpu,ngrid_tot,igrid,jgrid,kgrid,ipart,jpart,next_part,npart1,npart_tot,npart_now,ngrid_now,weight_eql,ithr

  ngrid_tot=0
  do icpu=1,ncpu
     ngrid_tot=ngrid_tot+numbl(icpu,ilevel)
  end do

  if(ngrid_tot==0) then
     ngrid_thr=0
     return
  end if

  npart_tot=0
  ! Count the total number of particles
!$omp parallel private(igrid) reduction(+:npart_tot)
  do icpu=1,ncpu
!$omp do
     do jgrid=1,numbl(icpu,ilevel)
        if(icpu==myid)then
           igrid=active(ilevel)%igrid(jgrid)
        else
           igrid=reception(icpu,ilevel)%igrid(jgrid)
        end if
        npart_tot=npart_tot+numbp(igrid)
     end do
!$omp end do nowait
  end do
!$omp end parallel

  ! Calculate rough number of particles per grids
  weight_eql=(npart_tot+ngrid_tot)/nthr

  head_thr(1)=0; jgrid_thr(1)=1
  npart_now=0; ngrid_now=0; ngrid_thr=0; npart_thr=0

  ithr=1
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     do jgrid=1,numbl(icpu,ilevel)
        if(head_thr(1)==0) then
           head_thr(1)=igrid
           icpu_thr(1)=icpu
        end if
        if(ithr == nthr) then
           ngrid_thr(ithr)=ngrid_tot
           npart_thr(ithr)=npart_tot
           return
        end if
        ngrid_now=ngrid_now+1
        npart_now=npart_now+numbp(igrid)
        igrid=next(igrid)
        ! Save state and move to next thread
        if(npart_now+ngrid_now>=weight_eql) then

           ngrid_thr(ithr)=ngrid_now
           npart_thr(ithr)=npart_now

           ngrid_tot=ngrid_tot-ngrid_now
           ngrid_now=0

           ithr=ithr+1

           head_thr(ithr)=igrid
           icpu_thr(ithr)=icpu
           jgrid_thr(ithr)=jgrid+1

           ! Correct the weight
           npart_tot=npart_tot-npart_now
           weight_eql=(npart_tot+ngrid_tot)/(nthr-ithr+1)
           npart_now=0
        end if
     end do
  end do
end subroutine cpu_parallel_link_all
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cpu_parallel_reception(ilevel,icpu_thr,igrid_thr,ngrid_thr)
  use amr_commons
  use pm_commons
  implicit none
  !----------------------------------------------------------------------
  ! This routine is designed to run on multiple cpu loops,
  ! so it stores current status of the loop (icpu, igrid)
  !----------------------------------------------------------------------
  integer :: ilevel
  integer, dimension(1:nthr) :: icpu_thr,igrid_thr,ngrid_thr

  integer :: icpu,ithr,ngrid_tot,igrid,ngrid_now,weight_eql,ncache,ngrid


  ngrid_tot=0
!$omp parallel do private(icpu) reduction(+:ngrid_tot)
  do icpu=1,ncpu
     ngrid_tot=ngrid_tot+reception(icpu,ilevel)%ngrid
  end do

  if(ngrid_tot==0) then
     icpu_thr=1
     igrid_thr=1
     ngrid_thr=0
     return
  end if

  weight_eql=ngrid_tot/nthr

  ngrid_now=0

  icpu_thr(1)=1
  igrid_thr(1)=1
  ngrid_thr=0

  ithr=2
  cpu_loop: do icpu=1,ncpu
     ncache=reception(icpu,ilevel)%ngrid
     do igrid=1,ncache,nvector
        ! Gather nvector grids
        ngrid=MIN(nvector,ncache-igrid+1)
        if(ngrid_now>=weight_eql) then
           icpu_thr(ithr)=icpu
           igrid_thr(ithr)=igrid

           ngrid_thr(ithr-1)=ngrid_now
           ngrid_tot=ngrid_tot-ngrid_now
           ngrid_now=0
           weight_eql=ngrid_tot/(nthr-ithr+1)

           ithr=ithr+1
        end if
        if(ithr>nthr) then
           exit cpu_loop
        end if
        ngrid_now=ngrid_now+ngrid
     end do
  end do cpu_loop
  ngrid_thr(ithr-1)=ngrid_tot

end subroutine cpu_parallel_reception
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine cpu_parallel_reception_part(ilevel,icpu_thr,igrid_thr,ngrid_thr,ipcom_thr)
  use amr_commons
  use pm_commons
  implicit none
  !----------------------------------------------------------------------
  ! This routine is designed to run on multiple cpu loops,
  ! so it stores current status of the loop (icpu, igrid)
  !----------------------------------------------------------------------
  integer :: ilevel
  integer, dimension(1:nthr) :: icpu_thr,igrid_thr,ngrid_thr,ipcom_thr

  integer :: icpu,ithr,ngrid_tot,npart_tot,igrid,npart_now,ngrid_now,weight_eql,ncache,ngridm,ipcom,npart1

  ngrid_tot=0
!$omp parallel do private(icpu) reduction(+:ngrid_tot)
  do icpu=1,ncpu
     ngrid_tot=ngrid_tot+reception(icpu,ilevel)%ngrid
  end do

  if(ngrid_tot==0) then
     icpu_thr=1
     igrid_thr=1
     ngrid_thr=0
     return
  end if

  npart_tot=0
  ! Count the total number of particles
!$omp parallel reduction(+:npart_tot)
  do icpu=1,ncpu
!$omp do
     do igrid=1,reception(icpu,ilevel)%ngrid
        npart_tot=npart_tot+numbp(reception(icpu,ilevel)%igrid(igrid))
     end do
!$omp end do nowait
  end do
!$omp end parallel

  if(npart_tot==0) then
     icpu_thr=1
     igrid_thr=1
     ngrid_thr=0
     return
  end if

  weight_eql=(ngrid_tot+npart_tot)/nthr

  npart_now=0; ngrid_now=0

  icpu_thr(1)=1; igrid_thr(1)=1; ngrid_thr=0; ipcom_thr(1)=0

  ithr=2
  cpu_loop: do icpu=1,ncpu
     ipcom=0
     do igrid=1,reception(icpu,ilevel)%ngrid
        ! Gather nvector grids
        if(npart_now+ngrid_now>=weight_eql) then
           icpu_thr(ithr)=icpu
           igrid_thr(ithr)=igrid
           ipcom_thr(ithr)=ipcom

           ngrid_thr(ithr-1)=ngrid_now

           npart_tot=npart_tot-ngrid_now
           npart_now=0

           ngrid_tot=ngrid_tot-ngrid_now
           ngrid_now=0

           weight_eql=(npart_tot+ngrid_tot)/(nthr-ithr+1)

           ithr=ithr+1
        end if
        if(ithr>nthr) then
           exit cpu_loop
        end if

        ngrid_now=ngrid_now+1

        npart1=numbp(reception(icpu,ilevel)%igrid(igrid))
        npart_now=npart_now+npart1
        ipcom=ipcom+npart1

     end do
  end do cpu_loop
  ngrid_thr(ithr-1)=ngrid_tot

end subroutine cpu_parallel_reception_part
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_partmap_notracer(ilevel,npart_grid,ngrid_tot,npart_tot)
  use amr_commons
  use pm_commons
  implicit none
  integer :: ilevel,icpu,igrid,jgrid,kgrid,ipart,jpart
  integer :: npart1,npart_now
  integer,intent(in) :: ngrid_tot
  integer,dimension(1:ngrid_tot) :: npart_grid
  integer,intent(out) :: npart_tot

  npart_tot=0
  ! Count the total number of particles
!$omp parallel private(igrid,npart1,npart_now,ipart,kgrid) reduction(+:npart_tot)
  kgrid=0
  do icpu=1,ncpu
!$omp do
     do jgrid=1,numbl(icpu,ilevel)
        if(icpu==myid)then
           igrid=active(ilevel)%igrid(jgrid)
        else
           igrid=reception(icpu,ilevel)%igrid(jgrid)
        end if
        npart1=numbp(igrid)
        npart_now=0
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              if (is_not_tracer(typep(ipart))) then
                 npart_now=npart_now+1
              end if
              ipart=nextp(ipart) ! Go to next particle
           end do
           ! End loop over particles
        end if
        npart_grid(kgrid+jgrid)=npart_now
        npart_tot=npart_tot+npart_now
     end do
!$omp end do nowait
     kgrid=kgrid+numbl(icpu,ilevel)
  end do
!$omp end parallel

end subroutine get_partmap_notracer
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################
subroutine get_partmap_cloud(ilevel,npart_grid,ngrid_tot,npart_tot)
  use amr_commons
  use pm_commons
  implicit none
  integer :: ilevel,icpu,igrid,jgrid,kgrid,ipart,jpart
  integer :: npart1,npart_now
  integer,intent(in) :: ngrid_tot
  integer,dimension(1:ngrid_tot) :: npart_grid
  integer,intent(out) :: npart_tot

  npart_tot=0
  ! Count the total number of particles
!$omp parallel private(igrid,npart1,npart_now,ipart,kgrid) reduction(+:npart_tot)
  kgrid=0
  do icpu=1,ncpu
!$omp do
     do jgrid=1,numbl(icpu,ilevel)
        if(icpu==myid)then
           igrid=active(ilevel)%igrid(jgrid)
        else
           igrid=reception(icpu,ilevel)%igrid(jgrid)
        end if
        npart1=numbp(igrid)
        npart_now=0
        if(npart1>0)then
           ipart=headp(igrid)
           ! Loop over particles
           do jpart=1,npart1
              if (is_cloud(typep(ipart))) then
                 npart_now=npart_now+1
              end if
              ipart=nextp(ipart) ! Go to next particle
           end do
           ! End loop over particles
        end if
        npart_grid(kgrid+jgrid)=npart_now
        npart_tot=npart_tot+npart_now
     end do
!$omp end do nowait
     kgrid=kgrid+numbl(icpu,ilevel)
  end do
!$omp end parallel

end subroutine get_partmap_cloud
!#########################################################################
!#########################################################################
!#########################################################################
!#########################################################################


subroutine pthreadLinkedList(ipart,nworks,nthreads,nparticles,ptrhead,nextlink)
  implicit none
  integer, intent(in):: ipart, nworks
  integer:: eqlwrks, nthreads,icount,ithread,i,jpart
  integer, dimension(1:*), intent(in):: nextlink
  integer, dimension(0:nthreads-1), intent(out):: nparticles, ptrhead
  integer, dimension(0:nthreads-1):: beforenwrks
  integer:: runningwrks,remainwrks,allocworks,k
  eqlwrks = (nworks+nthreads-1)/ nthreads
  runningwrks = 0
  remainwrks = nworks
  jpart = ipart

  !   k = 1000
  !   do while( k.ge.0)
  !      k = 1000
  !   enddo

  do i = 0, nthreads-1
     if(remainwrks .ge. eqlwrks) then
        allocworks = eqlwrks
     else
        allocworks = remainwrks
     endif
     nparticles(i) = allocworks
     runningwrks = runningwrks + allocworks
     remainwrks = remainwrks - allocworks
  enddo
  beforenwrks(0) = 0
  do i = 1, nthreads-1
     beforenwrks(i) = nparticles(i-1) + beforenwrks(i-1)
  enddo
  !  ptrhead(0) = jpart
  !  ithread = 1;
  ptrhead = 0

  ithread = 0

  do i=1,nworks
     if(i .eq. beforenwrks(ithread)+1) then
        ptrhead(ithread) = jpart
        ithread = ithread + 1
        if(ithread .ge. nthreads) then
          return
        endif
     endif
     jpart=nextlink(jpart)
  enddo
  return
end subroutine pthreadLinkedList


