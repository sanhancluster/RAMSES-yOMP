module sink_particle_tracer
  use amr_parameters, only : dp, ndim, nvector, verbose, aexp, nlevelmax, boxlen, icoarse_max, icoarse_min, nlevelsheld
  use amr_parameters, only : X_floor
  use amr_commons, only: ncpu, MC_tracer, metal, myid, star, sink, twondim, use_initial_mass, &
       write_stellar_densities, communicator
  ! Particle list and stuff
  use pm_commons, only: numbp_free, numbp_free_tot, headp, itmpp, nextp, numbp, headp, nextp, levelp, &
       typep, tracer_seed
  ! Particle data
  use pm_commons, only: idp, itmpp, typep, xp, vp, mp, tp, zp, st_n_tp, st_n_SN, st_e_SN, mp0, partp, move_flag
  ! Particle functions
  use pm_commons, only: is_gas_tracer, int2part, part2int

  ! Util function
  use utils, only : cgs, cell2grid, grid_level, normalize_position, pi, draw_normal, find_grid_containing_max
  use random, only : ranf, gaussdev
  use mpi_mod
  use geom_utils, only : basis_from_axis


  implicit none

  private
  real(dp), allocatable :: xAGN(:, :), jAGN(:, :), mAGN(:), dAGNcell(:)
  integer, allocatable :: ind_blast(:)
  logical, allocatable :: jet_mode_AGN(:)

  integer :: nAGN

  public :: MC_tracer_to_jet, prepare_MC_tracer_to_jet

contains
  subroutine prepare_MC_tracer_to_jet(xAGN_, jAGN_, mAGN_, dAGNcell_, &
       X_radio_, ind_blast_, nAGN_)
    ! This routines stores all the quantities required for jet/tracer interaction.
    ! It should be called each time the jet quantities are computed.
    real(dp), dimension(1:nAGN_, 1:3), intent(in) :: xAGN_, jAGN_
    real(dp), dimension(1:nAGN_), intent(in) :: mAGN_, X_radio_, dAGNcell_
    integer, dimension(1:nAGN_), intent(in) :: ind_blast_

    integer, intent(in) :: nAGN_

    nAGN = nAGN_
    if (allocated(xAGN)) then
       deallocate(xAGN, jAGN, mAGN, dAGNcell, ind_blast, jet_mode_AGN)
    end if

    allocate(xAGN(1:nAGN, 1:ndim), jAGN(1:nAGN, 1:ndim), mAGN(1:nAGN), dAGNcell(1:nAGN), &
         ind_blast(1:nAGN), jet_mode_AGN(1:nAGN))
    xAGN(:, :)      = xAGN_(:, :)           ! Position of the AGN
    jAGN(:, :)      = jAGN_(:, :)           ! Spin of the AGN
    mAGN(:)         = mAGN_(:)              ! Mass send in feedback
    dAGNcell(:)     = dAGNcell_(:)          ! Mass of the central cell
    ind_blast(:)    = ind_blast_(:)         ! Index of the cell
    jet_mode_AGN(:) = X_radio_(:) < X_floor ! True when in jet mode
  end subroutine prepare_MC_tracer_to_jet

  subroutine MC_tracer_to_jet(ilevel)
    ! This routines treats the MC tracers to move them in the jet.
    ! It is called at each level and only moves particles at *this
    ! level*

    use amr_parameters, only: rAGN
    integer, intent(in) :: ilevel

#ifdef WITHOUTMPI
    integer :: MPI_STATUS_SIZE = 1
#endif
    integer :: iAGN, ipart, i,  ii, jj, ncache
    integer, dimension(1:nvector), save :: ind_AGN, ind_part, ind_grid, ind_cell
    real(dp) :: proba, rand
    real(dp), dimension(1:nvector, 1:ndim), save ::AGN_pos, AGN_j, buffer_j, buffer_pos
    real(dp), dimension(1:nvector), save ::AGN_mass, AGN_cell_mass
    real(dp) :: rmax, scale, scale_l, scale_t, scale_d, scale_v, scale_nH, scale_T2, scale_m
    real(dp) :: dx_min
    integer :: nx_loc

    logical, dimension(1:nvector), save :: AGN_ok, grid_ok
    logical :: ok

    integer :: j
    integer, dimension(:), allocatable, save :: sendbuf, recvbuf
    integer, dimension(:), allocatable, save :: reqsend, reqrecv
    integer, dimension(:, :), allocatable, save :: statuses
    if (.not. allocated(sendbuf)) then
       allocate(sendbuf(1:ncpu), recvbuf(1:ncpu))
       allocate(reqsend(1:6*ncpu), reqrecv(1:6*ncpu))
       allocate(statuses(1:MPI_STATUS_SIZE, 1:4*ncpu))
    end if

    if (verbose) write(*, *) ' Entering treat_MC_tracers at level', ilevel
    call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
    scale_m=scale_d*scale_l**3d0
    nx_loc = (icoarse_max - icoarse_min + 1)
    scale = boxlen / dble(nx_loc)
    dx_min=scale*0.5d0**(nlevelmax-nlevelsheld) / aexp
    rmax = max(dx_min * scale_l, rAGN * cgs%kpc)
    rmax = rmax/scale_l

    !----------------------------------------
    ! Select the AGNs in jet mode, and get their particles
    !----------------------------------------
    j = 0

    ! Loop over AGNs in jet mode, gather tracer particles and move
    ! them in the direction of the jet

    ! NOTE
    ! ----
    ! We only select the grids at the current level, so that the
    ! particles are detached from the correct grid.
    ii = 0
    ! Look on AGN (vectorized)
    do jj = 1, nAGN, nvector
       ncache = min(nvector, nAGN - jj + 1)
       do ii = 1, ncache
          iAGN = ii + jj - 1
          AGN_pos(ii, :) = xAGN(iAGN, :)
          AGN_mass(ii)   = mAGN(iAGN)
          ! Which one ?
          AGN_ok(ii) = jet_mode_AGN(iAGN) .and. ind_blast(iAGN) > 0
          ind_cell(ii) = ind_blast(iAGN)
          ind_grid(ii) = cell2grid(ind_cell(ii))

          ! This is the mass the cell has before removing the AGN mass
          AGN_cell_mass(ii) = dAGNcell(iAGN)
          AGN_j(ii, 1:ndim) = jAGN(iAGN, 1:ndim)

          ! Only keep grids that are at the current level
          if (AGN_ok(ii)) then
             grid_ok(ii) = (grid_level(ind_grid(ii)) == ilevel)
          else
             grid_ok(ii) = .false.
          end if
       end do

       ! This array is used to store the target CPU for tracer particles
       do ii = 1, ncache
          iAGN = ii + jj - 1
          if (AGN_ok(ii) .and. grid_ok(ii)) then
             ! Compute amount of gas moved w.r.t. cell mass
             ! TODO: check me here (should be uold before mass removal,
             ! no? if so, move *before* mass removal
             if (AGN_cell_mass(ii) == 0d0) then
                proba = 0
             else
                proba = AGN_mass(ii) / AGN_cell_mass(ii)
             end if

             ! Loop over particles in grid
             ipart = headp(ind_grid(ii))
             do i = 1, numbp(ind_grid(ii))
                ok = .false.
                if (is_gas_tracer(typep(ipart))) then
                   ! Reset "CPU" flag
                   itmpp(ipart) = 0
                   ok = (partp(ipart) == ind_cell(ii)) .and. move_flag(ipart) == 0
                end if

                if (ok) then
                   ! Draw random number
                   call ranf(tracer_seed, rand)

                   if (rand < proba) then
                      j = j + 1
                      ind_AGN(j) = iAGN
                      ind_part(j) = ipart
                      buffer_pos(j, :) = AGN_pos(ii, :)
                      buffer_j(j, :)   = AGN_j(ii, :)

                      if (j == nvector) then
                         call tracer2jet(ind_AGN, ind_part, &
                              buffer_pos, buffer_j, &
                              nvector, rmax)
                         j = 0
                      end if
                   end if
                end if
                ipart = nextp(ipart)
             end do ! end loop on particles
          end if
       end do ! end cache loop
    end do ! end loop on AGN

    if (j > 0) then
       call tracer2jet(ind_AGN, ind_part, &
            buffer_pos, buffer_j,&
            j, rmax)
    end if

    call tracer_particle_send_fine(ilevel)

    if (verbose) write (*, *) 'Exiting treat_MC_tracers'

  end subroutine MC_tracer_to_jet

  subroutine tracer2jet(ind_AGN, ind_part, AGN_pos, AGN_j, npart, rmax)
    ! Compute the position of the particle within the jet
    ! On exit, the array itmpp contains the target CPU
    integer, intent(in), dimension(1:nvector) :: ind_AGN, ind_part
    real(dp), intent(in), dimension(1:nvector, 1:3) :: AGN_pos, AGN_j
    real(dp), intent(in) :: rmax
    integer, intent(in) :: npart

    real(dp) :: radius2, hh, ux(3), uy(3), uz(3), rmax2, xx, yy
    real(dp), dimension(1:nvector, 1:3), save :: newPos
    integer, dimension(1:nvector), save :: cpus
    integer :: i
    logical :: ok

    rmax2 = rmax**2
    !----------------------------------------
    ! Generate new positions
    !----------------------------------------
    do i = 1, npart
       ok = .false.
       do while (.not. ok)
          call ranf(tracer_seed, hh)
          ! Draw a position following a normal law between -1 and 1 (in rmax unit)
          xx = 2
          yy = 2
          do while (norm2([xx, yy]) > 1)
             call gaussdev(tracer_seed, xx)
             call gaussdev(tracer_seed, yy)
          end do

          hh = (hh - 0.5_dp) * 4 * rmax ! *4 to include the spherical caps

          ! rescale to rmax units
          xx = xx * rmax
          yy = yy * rmax
          radius2 = xx**2 + yy**2

          ok = .true.
          if (abs(hh) > rmax) then ! Spherical cap case
             ! Check the point is in the sphere
             ok = ((abs(hh) - rmax)**2 + radius2 < rmax2)
          end if
       end do

       ! Compute new position
       uz = AGN_j(i, :) / norm2(AGN_j(i, :))
       ! Get fresh basis from AGN spin
       call basis_from_axis(uz, ux, uy)
       newpos(i, :) = AGN_pos(i, :) + &
            xx * ux + &
            yy * uy + &
            hh * uz
    end do

    ! Take care of boundary conditions
    call normalize_position(newpos, npart)

    ! Compute location of particles in CPUs
    call cmp_cpumap(newpos, cpus(1:npart), npart)

    do i = 1, npart
       itmpp(ind_part(i)) = cpus(i)
       xp(ind_part(i), :) = newpos(i, :)
    end do

  end subroutine tracer2jet

  subroutine tracer_particle_send_fine(ilevel)
    ! Communicate tracer particles sent through a jet across CPU boundaries.
    use amr_commons, only: emission, reception
    integer, intent(in) :: ilevel
    !-----------------------------------------------------------------------
    ! This subroutine moves tracer particles across processors boundaries.
    !-----------------------------------------------------------------------
    integer :: igrid, icell, ipart, ncache_tot
    integer :: ip, ipcom, npart1, icpu, ncache
    integer :: info, buf_count, tagf=102, tagu=102
    integer :: countsend, countrecv
#ifndef WITHOUTMPI
    integer, dimension(MPI_STATUS_SIZE, 2*ncpu) :: statuses
    integer, dimension(2*ncpu) :: reqsend, reqrecv
    integer, dimension(ncpu) :: sendbuf, recvbuf
#endif
    integer, dimension(1:nvector), save :: ind_part, ind_list, ind_com
    logical :: ok_free, ok
    integer :: particle_data_width
    integer :: iAGN

    integer :: i, next_ipart

    if(verbose) write(*, 111) ilevel

#ifdef WITHOUTMPI
    return
#else
    ! Count tracer particle to be sent
    do icpu = 1, ncpu
       reception(icpu, ilevel)%npart = 0
    end do

    !------------------------------------------------------------
    ! Count the number of particles to send
    !------------------------------------------------------------
    do iAGN = 1, nAGN
       ! Select AGN in jet mode, with their cell at the current level
       icell = ind_blast(iAGN)
       ok = (icell > 0) .and. jet_mode_AGN(iAGN)
       if (ok) then
          igrid = cell2grid(icell)
          ok = (grid_level(igrid) == ilevel)
       end if

       if (ok) then
          ! Loop on the particles in the grid
          ipart = headp(igrid)
          do i = 1, numbp(igrid)
             icpu = itmpp(ipart)
             if (is_gas_tracer(typep(ipart)) .and. icpu > 0) then
                reception(icpu, ilevel)%npart = reception(icpu, ilevel)%npart + 1
             end if

             ipart = nextp(ipart)
          end do
       end if
    end do

    do icpu = 1, ncpu
       sendbuf(icpu) = reception(icpu, ilevel)%npart
    end do

    !------------------------------------------------------------
    ! Calculate how many particle properties are being transferred
    !------------------------------------------------------------
    particle_data_width = twondim+1
    if(star.or.sink) then
       if(metal) then
          particle_data_width=twondim+6 !3
       else
          particle_data_width=twondim+5 !2
       endif
       if(write_stellar_densities) &
            particle_data_width = particle_data_width + 3
       if(use_initial_mass) particle_data_width = particle_data_width + 1
    endif

#ifdef OUTPUT_PARTICLE_POTENTIAL
    particle_data_width=particle_data_width+1
#endif

    ! Allocate communication buffer in emission
    do icpu=1,ncpu
       ncache=reception(icpu,ilevel)%npart
       if(ncache>0)then
          ! Allocate reception buffer
          allocate(reception(icpu,ilevel)%fp(1:ncache,1:5))
          allocate(reception(icpu,ilevel)%up(1:ncache,1:particle_data_width))
       end if
    end do

    !------------------------------------------------------------
    ! Fill communication buffer
    !------------------------------------------------------------
    do icpu=1, ncpu
       ipcom = 0
       ip = 0
       do iAGN = 1, nAGN
          icell = ind_blast(iAGN)
          ok = (icell > 0) .and. jet_mode_AGN(iAGN)
          if (ok) then
             igrid = cell2grid(icell)
             ok = (grid_level(igrid) == ilevel)
          end if

          if (ok) then
             ! Loop on the particles in the grid
             ipart = headp(igrid)
             do i = 1, numbp(igrid)
                ! Store next particle
                next_ipart = nextp(ipart)
                if (is_gas_tracer(typep(ipart)) .and. itmpp(ipart) == icpu) then
                   ! Very important: reset flag to prevent the particle from being sent later
                   itmpp(ipart) = 0

                   ! Fill the buffer
                   ipcom = ipcom + 1
                   ip = ip + 1
                   ind_com (ip) = ipcom
                   ind_part(ip) = ipart
                   ind_list(ip) = igrid

                   if (ip == nvector) then
                      call fill_tracer_comm(ind_part, ind_com, ind_list, nvector, ilevel, icpu, reception)
                      ip = 0
                   end if
                end if
                ipart = next_ipart
             end do
          end if
       end do
       if (ip > 0) call fill_tracer_comm(ind_part, ind_com, ind_list, ip, ilevel, icpu, reception)
    end do

    ! IMPORTANT
    ! Now reset ind_blast to prevent tracers to be sent at next step
    do iAGN = 1, nAGN
       icell = ind_blast(iAGN)
       if (icell > 0) then
          igrid = cell2grid(icell)
          if (grid_level(igrid) == ilevel) ind_blast(iAGN) = 0
       end if
    end do

    !------------------------------------------------------------
    ! Communicate virtual particle number to parent cpu
    !------------------------------------------------------------
    call MPI_ALLTOALL(sendbuf, 1, MPI_INTEGER, recvbuf, 1, MPI_INTEGER, MPI_COMM_WORLD, info)

    ! Allocate communication buffer in reception
    do icpu = 1, ncpu
       emission(icpu, ilevel)%npart = recvbuf(icpu)
       ncache = emission(icpu, ilevel)%npart
       if (ncache>0) then
          ! Allocate reception buffer
          allocate(emission(icpu,ilevel)%fp(1:ncache,1:5))
          allocate(emission(icpu,ilevel)%up(1:ncache,1:particle_data_width))
       end if
    end do

    ! Receive particles
    countrecv = 0
    do icpu = 1, ncpu
       ncache = emission(icpu,ilevel)%npart
       if(ncache>0)then
          buf_count = ncache*5
          countrecv = countrecv+1
#ifndef LONGINT
          call MPI_IRECV(emission(icpu, ilevel)%fp, buf_count,  &
               & MPI_INTEGER, icpu-1, &
               & tagf, MPI_COMM_WORLD, reqrecv(countrecv), info)
#else
          call MPI_IRECV(emission(icpu, ilevel)%fp, buf_count,  &
               & MPI_INTEGER8, icpu-1, &
               & tagf, MPI_COMM_WORLD, reqrecv(countrecv), info)
#endif
          buf_count = ncache*particle_data_width
          countrecv = countrecv+1
          call MPI_IRECV(emission(icpu, ilevel)%up, buf_count,  &
               & MPI_DOUBLE_PRECISION, icpu-1, &
               & tagu, MPI_COMM_WORLD, reqrecv(countrecv), info)
       end if
    end do

    ! Send particles
    countsend = 0
    do icpu = 1, ncpu
       ncache = reception(icpu, ilevel)%npart
       if(ncache>0)then
          buf_count = ncache*5
          countsend = countsend+1
#ifndef LONGINT
          call MPI_ISEND(reception(icpu, ilevel)%fp, buf_count,  &
               & MPI_INTEGER, icpu-1, &
               & tagf, MPI_COMM_WORLD, reqsend(countsend), info)
#else
          call MPI_ISEND(reception(icpu, ilevel)%fp, buf_count,  &
               & MPI_INTEGER8, icpu-1, &
               & tagf, MPI_COMM_WORLD, reqsend(countsend), info)
#endif
          buf_count = ncache*particle_data_width
          countsend = countsend+1
          call MPI_ISEND(reception(icpu, ilevel)%up, buf_count,  &
               & MPI_DOUBLE_PRECISION, icpu-1, &
               & tagu, MPI_COMM_WORLD, reqsend(countsend), info)
       end if
    end do

    ! Wait for full completion of receives
    call MPI_WAITALL(countrecv, reqrecv, statuses, info)

    ! Compute total number of newly created particles
    ncache_tot = 0
    do icpu = 1, ncpu
       ncache_tot = ncache_tot + emission(icpu, ilevel)%npart
    end do

    ! Wait for full completion of sends
    call MPI_WAITALL(countsend, reqsend, statuses, info)

    call MPI_ALLREDUCE(numbp_free, numbp_free_tot, 1, MPI_INTEGER, MPI_MIN, &
         & MPI_COMM_WORLD, info)
    ok_free = (numbp_free - ncache_tot) >= 0
    if (.not. ok_free) then
       write(*, *) 'No more free memory for particles'
       write(*, *) 'Increase npartmax'
       write(*, *) numbp_free, ncache_tot,  ilevel
       write(*, *) myid
       write(*, *) emission(1:ncpu, ilevel)%npart
       write(*, *) '============================'
       write(*, *) reception(1:ncpu, ilevel)%npart
       call MPI_ABORT(MPI_COMM_WORLD, 1, info)
    end if
    !------------------------------------------------------------
    ! Empty communication buffer
    !------------------------------------------------------------
    do icpu = 1, ncpu
       ! Loop over particles by vector sweeps
       ncache = emission(icpu, ilevel)%npart
       do ipart = 1, ncache, nvector
          npart1 = min(nvector, ncache-ipart+1)
          do ip = 1, npart1
             ind_com(ip) = ipart + ip - 1
          end do
          call empty_tracer_comm(ind_com, npart1, ilevel, icpu, emission)
       end do
    end do

    ! Deallocate temporary communication buffers
    do icpu = 1, ncpu
       ncache = emission(icpu, ilevel)%npart
       if (ncache>0) then
          deallocate(emission(icpu, ilevel)%fp)
          deallocate(emission(icpu, ilevel)%up)
       end if
       ncache = reception(icpu, ilevel)%npart
       if (ncache>0) then
          deallocate(reception(icpu, ilevel)%fp)
          deallocate(reception(icpu, ilevel)%up)
       end if
    end do
#endif

111 format('   Entering tracer_particle_send_fine for level ', I2)
  end subroutine tracer_particle_send_fine

  subroutine fill_tracer_comm(ind_part, ind_com, ind_list, np, ilevel, icpu, &
       reception)
    ! Fill the communicator array and remove particles from the linked list
    integer,  intent(in) :: np, ilevel, icpu
    integer, intent(in), dimension(1:nvector) :: ind_part, ind_com, ind_list
    integer :: current_property
    integer :: i, idim
    logical, dimension(1:nvector), save :: ok=.true.

    type(communicator), intent(inout) :: reception(:, :)

    ! Gather particle level and identity
    do i=1,np
       reception(icpu,ilevel)%fp(ind_com(i),2)=levelp(ind_part(i))
       reception(icpu,ilevel)%fp(ind_com(i),3)=idp   (ind_part(i))
    end do

    ! Gather particle position and velocity
    do idim=1,ndim
       do i=1,np
          reception(icpu,ilevel)%up(ind_com(i),idim     )=xp(ind_part(i),idim)
          reception(icpu,ilevel)%up(ind_com(i),idim+ndim)=vp(ind_part(i),idim)
       end do
    end do

    current_property = twondim+1
    ! Gather particle mass
    do i=1,np
       reception(icpu,ilevel)%up(ind_com(i),current_property)=mp(ind_part(i))
    end do
    current_property = current_property+1

#ifdef OUTPUT_PARTICLE_POTENTIAL
    ! Gather particle potential
    do i=1,np
       reception(icpu,ilevel)%up(ind_com(i),current_property)=ptcl_phi(ind_part(i))
    end do
    current_property = current_property+1
#endif

    ! Gather particle birth epoch
    if(star.or.sink)then
       do i=1,np
          reception(icpu,ilevel)%up(ind_com(i),current_property)=tp(ind_part(i))
       end do
       current_property = current_property+1
       if(metal)then
          do i=1,np
             reception(icpu,ilevel)%up(ind_com(i),current_property)=zp(ind_part(i))
          end do
          current_property = current_property+1
       end if
       if(write_stellar_densities) then
          do i=1,np
             reception(icpu,ilevel)%up(ind_com(i),current_property)  =st_n_tp(ind_part(i))
             reception(icpu,ilevel)%up(ind_com(i),current_property+1)=st_n_SN(ind_part(i))
             reception(icpu,ilevel)%up(ind_com(i),current_property+2)=st_e_SN(ind_part(i))
          end do
          current_property = current_property+3
       endif
       if(use_initial_mass)then
          do i=1,np
             reception(icpu,ilevel)%up(ind_com(i),current_property)=mp0(ind_part(i))
          end do
          current_property = current_property+1
       end if
#ifdef NTRACEGROUPS
       do i=1,np
          reception(icpu,ilevel)%up(ind_com(i),current_property)=ptracegroup(ind_part(i))
       end do
       current_property = current_property + 1
#endif
    end if
    ! Family
    do i=1,np
       reception(icpu,ilevel)%fp(ind_com(i), 4) = part2int(typep(ind_part(i)))
    end do

    if (MC_tracer) then
       do i=1,np
          reception(icpu, ilevel)%fp(ind_com(i), 5) = &
               partp(ind_part(i))
       end do
       ! Only useful if you add extra properties
       current_property = current_property + 2
    end if


    ! Remove particles
    call remove_list(ind_part, ind_list, ok, np)
    call add_free(ind_part, np)

  end subroutine fill_tracer_comm

  subroutine empty_tracer_comm(ind_com, np, ilevel, icpu, emission)
    ! Empty communicators for tracer particles and attach particles to the grid
#ifdef NTRACEGROUPS
    use pm_commons, only: ptracegroup
#endif
    integer, intent(in) ::np, icpu, ilevel
    integer, intent(in), dimension(1:nvector)::ind_com
    type(communicator), intent(in) :: emission(:, :)

    integer :: i, idim
    integer, dimension(1:nvector), save :: ind_list, ind_part, ind_cell, ind_level
    logical, dimension(1:nvector), save :: ok=.true.
    integer :: current_property
    real(dp), dimension(1:nvector, 1:ndim), save :: pos

    ! Compute parent grid index
    do idim = 1, 3
       do i = 1, np
          pos(i, idim) = emission(icpu, ilevel)%up(ind_com(i), idim)
       end do
    end do

    ! Particles will be attached at most as deep as the current level
    call find_grid_containing_max(pos, ind_list, ind_cell, ind_level, np, ilevel)

    ! Add particle to parent linked list
    call remove_free(ind_part, np)
    call add_list(ind_part, ind_list, ok, np)

    ! Scatter particle level and identity
    do i=1,np
       levelp(ind_part(i)) = emission(icpu, ilevel)%fp(ind_com(i), 2)
       idp   (ind_part(i)) = emission(icpu, ilevel)%fp(ind_com(i), 3)
    end do

    ! Scatter particle position and velocity
    do idim = 1, ndim
       do i = 1, np
          xp(ind_part(i), idim) = emission(icpu, ilevel)%up(ind_com(i), idim     )
          vp(ind_part(i), idim) = emission(icpu, ilevel)%up(ind_com(i), idim+ndim)
       end do
    end do

    current_property = twondim + 1

    ! Scatter particle mass
    do i = 1, np
       mp(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
    end do
    current_property = current_property + 1

#ifdef OUTPUT_PARTICLE_POTENTIAL
    ! Scatter particle phi
    do i = 1, np
       ptcl_phi(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
    end do
    current_property = current_property+1
#endif

    ! Scatter particle birth epoch
    if(star.or.sink)then
       do i = 1, np
          tp(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
       end do
       current_property = current_property+1
       if(metal)then
          do i = 1, np
             zp(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
          end do
          current_property = current_property+1
       end if
       if(write_stellar_densities) then
          do i = 1, np
             st_n_tp(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)   !SD
             st_n_SN(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property+1) !SD
             st_e_SN(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property+2) !SD
          end do
          current_property = current_property+3
       endif
       if(use_initial_mass)then
          do i = 1, np
             mp0(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
          end do
          current_property = current_property+1
       end if
#ifdef NTRACEGROUPS
       do i = 1, np
          ptracegroup(ind_part(i)) = emission(icpu, ilevel)%up(ind_com(i), current_property)
       end do
       current_property = current_property+1
#endif
    end if
    ! Add family
    do i = 1, np
       typep(ind_part(i)) = int2part(emission(icpu, ilevel)%fp(ind_com(i), 4))
       typep(ind_part(i))%tag = typep(ind_part(i))%tag + 10
    end do
    ! MC Tracer
    if (MC_tracer) then
       do i = 1, np
          partp(ind_part(i)) = emission(icpu, ilevel)%fp(ind_com(i), 5)
          move_flag(ind_part(i)) = 1
       end do
    end if
  end subroutine empty_tracer_comm

end module sink_particle_tracer

subroutine tracer2othersink(ind_tracer, isink_new_part, xsink_loc, np)
  ! Transfer tracers from one sink to another one
  use pm_commons
  implicit none

  integer, dimension(1:nvector), intent(in) :: ind_tracer, isink_new_part
  integer, intent(in) :: np
  real(dp), dimension(1:nvector, 1:ndim), intent(in) :: xsink_loc

  integer :: j,  isink_new
  do j = 1, np
     ! Retrieve the index of the new sink
     isink_new = isink_new_part(j)

     ! If the sink has changed, reattach to the new one
     if (isink_new /= 0) then
        partp(ind_tracer(j)) = isink_new

        ! Do not move tracers there, they will be moved in move_fine.f90.
        ! xp(ind_tracer(j), :) = xsink(isink_new, :)
     end if
  end do

end subroutine tracer2othersink

subroutine tracer2sink(ind_tracer, proba, xsink_loc, isink, nattach, dx_loc)
  ! Attach tracers to a sink (accretion onto the BH)
  use amr_commons
  use random, only : ranf
  use pm_commons, only : tracer_seed, typep, partp, move_flag, FAM_TRACER_CLOUD
  implicit none

  integer, intent(in) :: nattach
  integer, dimension(1:nvector), intent(in) :: ind_tracer, isink
  real(dp), dimension(1:nvector), intent(in) :: proba
  real(dp), dimension(1:nvector, 1:3), intent(in) :: xsink_loc
  real(dp), intent(in) :: dx_loc

  logical, dimension(1:nvector), save :: attach = .false.
  integer :: i
  real(dp) :: r

  do i = 1, nattach
     call ranf(tracer_seed, r)
     attach(i) = r < proba(i)
  end do

  ! Change particle pointer and kind
  do i = 1, nattach
     if (attach(i)) then
        partp(ind_tracer(i)) = isink(i)
        typep(ind_tracer(i))%family = FAM_TRACER_CLOUD
     end if
  end do

  ! Do not move tracers there, they will be moved in move_fine.f90 specifically.

  ! do idim = 1, ndim
  !    do i = 1, nattach
  !       if (attach(i)) then
  !          ! Move at most 2*dx in the direction
  !          xtmp = xp(ind_tracer(i), idim)
  !          if (abs(xsink_loc(i, idim) - xtmp) > 2*dx_loc) then
  !             xp(ind_tracer(i), idim) = xtmp + (xsink_loc(i, idim)-xtmp) / dx_loc * 2
  !          else
  !             xp(ind_tracer(i), idim) = xsink_loc(i, idim)
  !          end if
  !       end if
  !    end do
  ! end do

  ! Save state of the particle
  do i = 1, nattach
     if (attach(i)) then
        ! Set to 1 to prevent further moves
        move_flag(ind_tracer(i)) = 1
     end if
  end do

end subroutine tracer2sink
