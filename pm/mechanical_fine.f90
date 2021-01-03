!####################################################################
!####################################################################
!####################################################################
subroutine mechanical_feedback_fine(ilevel,icount)
  use pm_commons
  use amr_commons
  use mechanical_commons
  use hydro_commons,only:uold
#ifdef _OPENMP
  use omp_lib
#endif
#ifdef RT
  use rt_parameters,only:group_egy
  use SED_module,only:nSEDgroups,inp_SED_table
#endif
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
#endif
  !------------------------------------------------------------------------
  ! This routine computes the energy liberated from supernova II ,
  ! and inject momentum and energy to the surroundings of the young stars.
  ! This routine is called every fine time step.
  ! ind_pos_cell: position of the cell within an oct
  ! m8,mz8,nph8: temporary variable necessary for an oct to add up 
  !              the mass, metal, etc. on a cell by cell basis
  ! mejecta: mass of the ejecta from SNe
  ! Zejecta: metallicity of the ejecta (not yield)
  ! mZSNe : total metal mass from SNe in each cell
  ! mchSNe: total mass of each chemical element in each cell
  ! nphSNe: total production rate of ionising radiation in each cell.
  !         This is necessary to estimate the Stromgren sphere
  !------------------------------------------------------------------------
  integer::igrid,jgrid,ipart,jpart,next_part
  integer::npart1,npart2,icpu,icount,idim,ip,ich
  integer::ind,ind_son,ind_cell,ilevel,iskip,info
  integer::nSNc,nSNc_mpi
  integer,dimension(1:nvector),save::ind_grid,ind_pos_cell
  real(dp)::nsnII_star,mass0,mass_t,nsnII_tot,nsnII_mpi
  real(dp)::tyoung,current_time,dteff
  real(dp)::skip_loc(1:3),scale,dx,dx_loc,vol_loc,x0(1:3)
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_msun
  real(dp),dimension(1:twotondim),save::m8, mz8, n8, nph8 ! SNe
  real(dp),dimension(1:twotondim,1:3),save:: p8  ! SNe
  real(dp),dimension(1:nvector),save::nSNe, mSNe, mZSNe, nphSNe
  real(dp),dimension(1:nvector,1:3),save::pSNe
  real(dp),dimension(1:nvector,1:nchem),save::mchSNe
  real(dp),dimension(1:twotondim,1:nchem),save::mch8 ! SNe
  real(dp)::mejecta,mejecta_ch(1:nchem),Zejecta,mfrac_snII
  real(dp)::snII_freq_noboost, M_SNII_var=0.0
  real(dp),parameter::msun2g=1.989d33
  real(dp),parameter::myr2s=3.1536000d+13
  real(dp)::ttsta,ttend
  logical::ok,done_star
#ifdef RT
  real(dp),allocatable,dimension(:)::L_star
  real(dp)::Z_star,age,L_star_ion
  integer::iph,igroup
  Z_star=z_ave
  allocate(L_star(1:nSEDgroups))
#endif
  ! MC Tracer =================================================
  real(dp)::rand
  real(dp)::new_xp(1:ndim)
  integer :: irand
  ! End MC Tracer =============================================

  integer,dimension(1:IRandNumSize),save :: ompseed
!$omp threadprivate(ompseed)

  if(icount==2) return
  if(.not.hydro) return
  if(ndim.ne.3)  return
  if(numbtot(1,ilevel)==0)return
  if(nstar_tot==0)return

#ifndef WITHOUTMPI
  if(myid.eq.1) ttsta=MPI_WTIME(info)
#endif 
  nSNc=0; nsnII_tot=0d0

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g

  ! Mesh spacing in that level
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc) 

  ! To filter out old particles and compute individual time steps
  ! NB: the time step is always the coarser level time step, since no feedback for icount=2
  if (ilevel==levelmin)then
     dteff = dtnew(ilevel)
  else
     dteff = dtnew(ilevel-1)
  endif

  if (use_proper_time)then
     tyoung = t_sne*myr2s/(scale_t/aexp**2) 
     current_time=texp
     dteff = dteff*aexp**2
  else
     tyoung = t_sne*myr2s/scale_t 
     current_time=t
  endif
  tyoung = current_time - tyoung

  ncomm_SN=0  ! important to initialize; number of communications (not SNe)
#ifndef WITHOUTMPI
  xSN_comm=0d0;ploadSN_comm=0d0;mSN_comm=0d0
  mloadSN_comm=0d0;mZloadSN_comm=0d0;iSN_comm=0
  floadSN_comm=0d0;eloadSN_comm=0d0
#endif

#ifdef _OPENMP
!$omp parallel
    ! Give slight offsets for each OMP threads
    ompseed=MOD(tracer_seed+omp_get_thread_num()+1,4096)
!$omp end parallel
#else
    ompseed=MOD(tracer_seed+1,4096)
#endif
    call ranf(tracer_seed,rand)

  ! MC Tracer =================================================
  ! Reset tmpp array that contains the probability to be detached from the particle
  if (MC_tracer) then
     do icpu=1,ncpu
        ! Loop over grids
!$omp do private(igrid,npart1,ipart) schedule(dynamic,nchunk)
        do jgrid=1,numbl(icpu,ilevel)
           if(icpu==myid)then
              igrid=active(ilevel)%igrid(jgrid)
           else
              igrid=reception(icpu,ilevel)%igrid(jgrid)
           end if
           npart1=numbp(igrid)  ! Number of particles in the grid
           ! Count star particles
           if(npart1>0)then
              ipart=headp(igrid)
              do jpart=1,npart1
                 if(is_star(typep(ipart))) tmpp(ipart) = 0d0
                 ipart=nextp(ipart)
              end do
           end if
        end do
     end do
  end if
  ! End MC Tracer

  ! Type II Supernova frequency per Msun
  snII_freq_noboost = eta_sn / M_SNII

  if (snII_freq<0) then
     snII_freq = snII_freq_noboost
  endif

  if(.not.SNII_zdep_yield) then
      M_SNII_var = M_SNII * (snII_freq_noboost / snII_freq)  
      ! boostx5 -> M_SNII_var~20/5=4 Msun
  endif

  ! Loop over cpus
  do icpu=1,ncpu
     igrid=headl(icpu,ilevel)
     ip=0

     ! Loop over grids
     do jgrid=1,numbl(icpu,ilevel)
        npart1=numbp(igrid)  ! Number of particles in the grid
        npart2=0

        ! Count star particles
        if(npart1>0)then
           do idim=1,ndim
              x0(idim) = xg(igrid,idim) -dx -skip_loc(idim) 
           end do
 
           ipart=headp(igrid)

           m8=0d0;mz8=0d0;p8=0d0;n8=0d0;nph8=0d0;mch8=0d0
    
           ! Loop over particles
           do jpart=1,npart1
              ! Save next particle   <--- Very important !!!
              next_part=nextp(ipart)
#ifdef POP3
              if(pop3 .and. (zp(ipart).lt.Zcrit_pop3) ) then
                 ok=.false. ! this will be done elsewhere
                 done_star=.false.
              else
#endif
                 nsnII_star=0d0

                 ! active star particles?
                 ok = is_star_active(typep(ipart))    ! SanHan: check this out

                 if(ok)then
                    ! initial mass
                    if(use_initial_mass) then
                       mass0 = mp0(ipart)*scale_msun
                    else
                       mass0 = mp (ipart)*scale_msun
                    endif
                    mass_t = mp (ipart)*scale_msun
 

                    ok=.false.
                    if(sn2_real_delay)then
                       if(tp(ipart).ge.tyoung)then  ! if younger than t_sne

                          if(SNII_zdep_yield)then
                             call SNII_yield (zp(ipart), mfrac_snII, Zejecta, Zejecta_chem_II)
                             ! Adjust M_SNII mass not to double-count the mass loss from massive stars
                             M_SNII_var = mfrac_snII / snII_freq ! ex) 0.1 / 0.01 = 10 Msun
                          else
                             Zejecta = zp(ipart)+(1d0-zp(ipart))*yield
                          endif

                          call get_number_of_sn2  (tp(ipart), dteff, zp(ipart), idp(ipart),&
                                & mass0, mass_t, M_SNII_var, nsnII_star, done_star)
                          if(nsnII_star>0)ok=.true.

                       endif

                    else ! single SN event
                       if(tp(ipart).le.tyoung)then   ! if older than t_sne
                          ok=.true.

                          if(SNII_zdep_yield)then
                             call SNII_yield (zp(ipart), mfrac_snII, Zejecta, Zejecta_chem_II)
                             ! Adjust M_SNII mass not to double-count the mass loss from massive stars
                             M_SNII_var = mfrac_snII / snII_freq ! ex) 0.1 / 0.05 = 2 Msun
                          else
                             Zejecta = zp(ipart)+(1d0-zp(ipart))*yield
                          endif

                          ! number of sn doesn't have to be an integer
                          nsnII_star = mass0*snII_freq
                       endif
                    endif
                 endif
#ifdef POP3
              endif
#endif

              if(ok)then
                 ! Find the cell index to get the position of it
                 ind_son=1
                 do idim=1,ndim
                    ind = int((xp(ipart,idim)/scale - x0(idim))/dx)
                    ind_son=ind_son+ind*2**(idim-1)
                 end do 
                 iskip=ncoarse+(ind_son-1)*ngridmax
                 ind_cell=iskip+igrid
                 if(son(ind_cell)==0)then  ! leaf cell

                    !----------------------------------
                    ! For Type II explosions
                    !----------------------------------
                    ! total ejecta mass in code units
                    mejecta = M_SNII_var/scale_msun*nsnII_star

                    ! number of SNII
                    n8(ind_son) = n8(ind_son) + nsnII_star
                    ! mass return from SNe
                    m8 (ind_son)  = m8(ind_son) + mejecta
                    ! momentum from the original star, not the one generated by SNe
                    p8 (ind_son,1) = p8(ind_son,1) + mejecta*vp(ipart,1)
                    p8 (ind_son,2) = p8(ind_son,2) + mejecta*vp(ipart,2)
                    p8 (ind_son,3) = p8(ind_son,3) + mejecta*vp(ipart,3)
                    ! metal mass return from SNe including the newly synthesised one
                    if(metal)then
                       mz8(ind_son) = mz8(ind_son) + mejecta*Zejecta
                    endif
                    do ich=1,nchem
                       mch8(ind_son,ich) = mch8(ind_son,ich) + mejecta*Zejecta_chem_II(ich)
                    end do 

                    ! subtract the mass return
                    mp(ipart)=mp(ipart)-mejecta
                    if (MC_tracer) tmpp(ipart) = mejecta / mp(ipart)

                    ! mark if we are done with this particle
                    if(sn2_real_delay) then
                       if(done_star)then ! only if all SNe exploded
                          !idp(ipart)=abs(idp(ipart))
                          typep(ipart)%tag = 0 
                       endif 
                    else
                       !idp(ipart)=abs(idp(ipart))
                       typep(ipart)%tag = 0 
                    endif
                  
#ifdef RT
                    ! Enhanced momentum due to pre-processing of the ISM due to radiation
                    if(rt.and.mechanical_geen) then
                       ! Let's count the total number of ionising photons per sec 
                       call getAgeGyr(tp(ipart), age)
                       if(metal) Z_star=zp(ipart)
                       Z_star=max(Z_star,10.d-5)

                       ! compute the number of ionising photons from SED
                       call inp_SED_table(age, Z_star, 1, .false., L_star) ! L_star = [# s-1 Msun-1]
                       L_star_ion = 0d0
                       do igroup=1,nSEDgroups
                          if(group_egy(igroup).ge.13.6) L_star_ion = L_star_ion + L_star(igroup)
                       end do
                       nph8 (ind_son)=nph8(ind_son) + mass0*L_star_ion ! [# s-1]
                    endif
#endif

                 endif
              endif
              ipart=next_part  ! Go to next particle
           end do
           ! End loop over particles
    
           do ind=1,twotondim
              if (abs(n8(ind))>0d0)then
                 ip=ip+1
                 ind_grid(ip)=igrid
                 ind_pos_cell(ip)=ind
     
                 ! collect information 
                 nSNe(ip)=n8(ind)
                 mSNe(ip)=m8(ind)
                 mZSNe(ip)=mz8(ind)
                 pSNe(ip,1)=p8(ind,1)
                 pSNe(ip,2)=p8(ind,2)
                 pSNe(ip,3)=p8(ind,3)
                 nphSNe(ip)=nph8(ind)  ! mechanical_geen
                 do ich=1,nchem
                    mchSNe(ip,ich)=mch8(ind,ich)
                 end do
   
                 ! statistics
                 nSNc=nSNc+1
                 nsnII_tot = nsnII_tot + nsnII_star
 
                 if(ip==nvector)then
                    call mech_fine(ind_grid,ind_pos_cell,ip,ilevel,dteff,nSNe,mSNe,pSNe,mZSNe,nphSNe,mchSNe)
                    ip=0
                 endif 
              endif
           enddo
    
    
        end if
        igrid=next(igrid)   ! Go to next grid
     end do ! End loop over grids
    
     if (ip>0) then
        call mech_fine(ind_grid,ind_pos_cell,ip,ilevel,dteff,nSNe,mSNe,pSNe,mZSNe,nphSNe,mchSNe)
        ip=0
     endif
  end do ! End loop over cpus

  do icpu=1,ncpu
     if (MC_tracer) then
        ! MC Tracer =================================================
        ! Loop over grids
!$omp do private(igrid,npart1,ipart,next_part,rand,irand,new_xp)
        do jgrid = 1, numbl(icpu, ilevel)
           if(icpu==myid)then
              igrid=active(ilevel)%igrid(jgrid)
           else
              igrid=reception(icpu,ilevel)%igrid(jgrid)
           end if
           npart1 = numbp(igrid)  ! Number of particles in the grid
           ipart = headp(igrid)

           ! Loop over tracer particles
           do jpart = 1, npart1
              next_part=nextp(ipart)

              if (is_star_tracer(typep(ipart))) then
                 ! Detach particle if required
                 call ranf(ompseed, rand)

                 ! Detach particles
                 if (rand < tmpp(partp(ipart))) then
                    typep(ipart)%family = FAM_TRACER_GAS

                    ! Change here to tag the particle with the star id
                    typep(ipart)%tag = typep(ipart)%tag + 1
                    move_flag(ipart) = 1

                    ! Detached, now decide where to move it
                    call ranf(ompseed, rand)

                    do idim = 1, ndim
                       ! Draw number between 1 and nSNnei
                       irand = floor(rand * nSNnei) + 1
                       new_xp(idim) = xSNnei(idim, irand)*dx
                    end do

                    do idim = 1, ndim
                       xp(ipart, idim) = xp(ipart, idim) + new_xp(idim)
                    end do
                 end if
              end if
              ipart = next_part ! Go to next particle
           end do ! End loop over particles
        end do ! End loop over grids
!$omp end do nowait
        ! End MC Tracer =============================================
     end if
  end do ! End loop over cpus


#ifndef WITHOUTMPI
  nSNc_mpi=0; nsnII_mpi=0d0
  ! Deal with the stars around the bounary of each cpu (need MPI)
  call mech_fine_mpi(ilevel)
  call MPI_ALLREDUCE(nSNc,nSNc_mpi,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  call MPI_ALLREDUCE(nsnII_tot,nsnII_mpi,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
  nSNc = nSNc_mpi
  nsnII_tot = nSNII_mpi
!  if(myid.eq.1.and.nSNc>0.and.log_mfb) then
!     ttend=MPI_WTIME(info)
!     write(*,*) '--------------------------------------'
!     write(*,*) 'Time elapsed in mechanical_fine [sec]:', sngl(ttend-ttsta), nSNc, sngl(nsnII_tot) 
!     write(*,*) '--------------------------------------'
!  endif
#endif

#ifdef RT
   deallocate(L_star)
#endif

end subroutine mechanical_feedback_fine
!################################################################
!################################################################
!################################################################
subroutine mech_fine(ind_grid,ind_pos_cell,np,ilevel,dteff,nSN,mSN,pSN,mZSN,nphSN,mchSN)
  use amr_commons
  use pm_commons
  use hydro_commons
  use mechanical_commons
  implicit none
  integer::np,ilevel ! actually the number of cells
  integer,dimension(1:nvector)::ind_grid,ind_pos_cell,icellvec
  real(dp),dimension(1:nvector)::nSN,mSN,mZSN,floadSN,nphSN
  real(dp),dimension(1:nvector)::mloadSN,mZloadSN,eloadSN
  real(dp),dimension(1:nvector,1:3)::pSN,ploadSN
  real(dp),dimension(1:nvector,1:nchem)::mchSN,mchloadSN
  !-----------------------------------------------------------------------
  ! This routine is called by subroutine mechanical_feedback_fine 
  !-----------------------------------------------------------------------
  integer::i,j,nwco,nwco_here,idim,icell,igrid,ista,iend,ilevel2
  integer::ind_cell,ncell,irad,ii,ich
  real(dp)::d,u,v,w,e,z,eth,ekk,Tk,d0,u0,v0,w0,dteff
  real(dp)::dx,dx_loc,scale,vol_loc,nH_cen,fleftSN
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_kms
  real(dp)::scale_msun,msun2g=1.989d33
  real(dp)::skip_loc(1:3),Tk0,ekk0,eth0,etot0,T2min
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  ! Grid based arrays
  real(dp),dimension(1:ndim,1:nvector),save::xc2
  real(dp),dimension(1:nvector,1:nSNnei), save::p_solid,ek_solid
  real(dp)::d_nei,Z_nei,Z_neisol,dm_ejecta,vol_nei
  real(dp)::mload,vload,Zload=0d0,f_esn2
  real(dp)::num_sn,nH_nei,f_w_cell,f_w_crit
  real(dp)::t_rad,r_rad,r_shell,m_cen,ekk_ej
  real(dp)::uavg,vavg,wavg,ul,vl,wl,ur,vr,wr
  real(dp)::d1,d2,d3,d4,d5,d6,dtot,pvar(1:nvarMHD)
  real(dp)::vturb,vth,Mach,sig_s2,dratio,mload_cen
  ! For stars affecting across the boundary of a cpu
  integer, dimension(1:nSNnei),save::icpuSNnei
  integer ,dimension(1:nvector,0:twondim):: ind_nbor
  logical,dimension(1:nvector,1:nSNnei),save ::snowplough
  real(dp),dimension(1:nvector)::rStrom ! in pc
  real(dp)::dx_loc_pc,psn_tr,chi_tr,psn_thor98,psn_geen15,fthor
  real(dp)::km2cm=1d5,M_SNII_var,boost_geen_ad=0d0,p_hydro,vload_rad,f_wrt_snow
  ! chemical abundance
  real(dp),dimension(1:nchem)::chload,z_ch
  ! fractional abundances ; for ionisation fraction and ref, etc
  real(dp),dimension(1:NVAR),save::fractions ! not compatible with delayed cooling
  integer::i_fractions
  real(dp)::emag

  ! starting index for passive variables except for imetal and chem
  i_fractions = imetal+nchem+1

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g
  scale_kms  = scale_v/1d5

  ! Mesh variables
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)
  dx_loc_pc = dx_loc*scale_l/3.08d18

  ! Record position of each cell [0.0-1.0] regardless of boxlen
  xc2=0d0
  do i=1,np
     do idim=1,ndim
        xc2(idim,i)=xg(ind_grid(i),idim)-skip_loc(idim)+xc(ind_pos_cell(i),idim)
     end do 
  end do

  !======================================================================
  ! Determine p_solid before redistributing mass 
  !   (momentum along some solid angle or cell) 
  ! - This way is desirable when two adjacent SNe explode simulataenously.
  ! - if the neighboring cell does not belong to myid, this will be done 
  !      in mech_fine_mpi
  !======================================================================
  p_solid=0d0;ek_solid=0d0;snowplough=.false.

  do i=1,np
     ind_cell = ncoarse+ind_grid(i)+(ind_pos_cell(i)-1)*ngridmax

     ! redistribute the mass/metals to the central cell
     call get_icell_from_pos (xc2(1:3,i), ilevel+1, igrid, icell,ilevel2)

     ! Sanity Check
     if((cpu_map(father(igrid)).ne.myid).or.&
        (ilevel.ne.ilevel2).or.&
        (ind_cell.ne.icell))     then 
        print *,'>>> fatal error in mech_fine'
        print *, cpu_map(father(igrid)),myid
        print *, ilevel, ilevel2
        print *, ind_cell, icell
        stop 
     endif
 
     num_sn    = nSN(i) ! doesn't have to be an integer
     M_SNII_var = mSN(i)*scale_msun / num_sn
     nH_cen    = uold(icell,1)*scale_nH
     m_cen     = uold(icell,1)*vol_loc*scale_msun

     emag=0.0d0
#ifdef SOLVERmhd
     do idim=1,ndim
        emag=emag+0.125d0*(uold(icell,idim+ndim+2)+uold(icell,idim+nvar))**2
     enddo
#endif
     d   = uold(icell,1)
     u   = uold(icell,2)/d
     v   = uold(icell,3)/d
     w   = uold(icell,4)/d
     e   = uold(icell,5)
     e   = e-0.5d0*d*(u**2+v**2+w**2) - emag
#if NENER>0
     do irad=1,nener
        e = e - uold(icell,ndim+2+irad) 
     end do
#endif
     Tk  = e/d*scale_T2*(gamma-1.0)
     if(Tk<0)then
        print *,'TKERR : mech fbk (pre-call): TK<0', TK,icell
        print *,'nH [H/cc]= ',d*scale_nH
        print *,'u  [km/s]= ',u*scale_v/1d5
        print *,'v  [km/s]= ',v*scale_v/1d5
        print *,'w  [km/s]= ',w*scale_v/1d5
!        stop
     endif

     ! For stability - T<0 seems to happen if strong shock occurs 
     ! due to mechanical feedback
     T2min=T2_star*(d*scale_nH/n_star)**(g_star-1.0)
     if(Tk<T2min)then
        e=T2min*d/scale_T2/(gamma-1.0)*1.2195
        uold(icell,5) = 0.5d0*d*(u**2+v**2+w**2) + e + emag
#if NENER>0
        do irad=1,nener
           uold(icell,5) = uold(icell,5) + uold(icell,ndim+2+irad)
        end do
#endif
     endif
    
     if(metal)then
        z   = uold(icell,imetal)/d
     else
        z   = z_ave*0.02
     endif

     !==========================================
     ! estimate floadSN(i)
     !==========================================
     ! notice that f_LOAD / f_LEFT is a global parameter for SN ejecta themselves!!!
     if(loading_type.eq.1)then
        ! based on Federrath & Klessen (2012)
        ! find the mach number of this cell
        vth   = sqrt(gamma*1.38d-16*Tk/1.673d-24) ! cm/s
        ncell = 1
        icellvec(1) = icell
        call getnbor(icellvec,ind_nbor,ncell,ilevel)
        u0 = uold(icell        ,2) 
        v0 = uold(icell        ,3) 
        w0 = uold(icell        ,4) 
        ul = uold(ind_nbor(1,1),2)-u0 
        ur = uold(ind_nbor(1,2),2)-u0
        vl = uold(ind_nbor(1,3),3)-v0
        vr = uold(ind_nbor(1,4),3)-v0
        wl = uold(ind_nbor(1,5),4)-w0
        wr = uold(ind_nbor(1,6),4)-w0
        vturb = sqrt(ul**2 + ur**2 + vl**2 + vr**2 + wl**2 + wr**2)
        vturb = vturb*scale_v ! cm/s
    
        Mach  = vturb/vth
 
        ! get volume-filling density for beta=infty (non-MHD gas), b=0.4 (a stochastic mixture of forcing modes)
        sig_s2 = log(1d0 + (0.4*Mach)**2d0)
        ! s = -0.5*sig_s2;   s = ln(rho/rho0)
        dratio = max(exp(-0.5*sig_s2),0.01) ! physicall M~100 would be hard to achieve
        floadSN(i) = min(dratio,f_LOAD)

     else
        dratio     = 1d0
        floadSN(i) = f_LOAD
     endif

     !==========================================
     ! estimate Stromgren sphere (relevant to RHD simulations only)
     ! (mechanical_geen=.true)
     !==========================================
     if(mechanical_geen.and.rt) rStrom(i) = (3d0*nphSN(i)/4./3.141592/2.6d-13/nH_cen**2d0)**(1d0/3d0)/3.08d18 ! [pc] 

     if(log_mfb)then
398     format('MFB z= ',f9.5,' N= ',f7.1,' nH,T,Z= ',3(f7.3,1x),' dx= ',f9.5)
        write(*,398) 1./aexp-1,num_sn,log10(d*scale_nH),log10(Tk),log10(z/0.02),log10(dx_loc*scale_l/3.08d18)
     endif

     dm_ejecta = f_LOAD*mSN(i)/dble(nSNnei)  ! per solid angle
     mload     = f_LOAD*mSN(i) + uold(icell,1)*vol_loc*floadSN(i)  ! total SN ejecta + host cell
     if(metal) Zload = (f_LOAD*mZSN(i) + uold(icell,imetal)*vol_loc*floadSN(i))/mload
     do ich=1,nchem
        chload(ich) = (f_LOAD*mchSN(i,ich) + uold(icell,ichem+ich-1)*vol_loc*floadSN(i))/mload
     end do 
 
     do j=1,nSNnei
        call get_icell_from_pos (xc2(1:3,i)+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell,ilevel2)
        if(cpu_map(father(igrid)).eq.myid) then ! if belong to myid

           Z_nei = z_ave*0.02 ! For metal=.false. 
           if(ilevel>ilevel2)then ! touching level-1 cells
              d_nei     = unew(icell,1)
              if(metal) Z_nei = unew(icell,imetal)/d_nei
           else
              d_nei     = uold(icell,1)
              if(metal) Z_nei = uold(icell,imetal)/d_nei
           endif

           f_w_cell  = (mload/dble(nSNnei) + d_nei*vol_loc/8d0)/dm_ejecta - 1d0
           nH_nei    = d_nei*scale_nH*dratio
           Z_neisol  = max(0.01, Z_nei/0.02)

           ! transition mass loading factor (momentum conserving phase)
           ! psn_tr = sqrt(2*chi_tr*Nsn*Esn*Msn*fe)
           ! chi_tr = (1+f_w_crit)
           psn_thor98   = A_SN * num_sn**(expE_SN) * nH_nei**(expN_SN) * Z_neisol**(expZ_SN)  !km/s Msun
           psn_tr       = psn_thor98 
           if(mechanical_geen)then 
              ! For snowplough phase, psn_tr will do the job
              psn_geen15 = A_SN_Geen * num_sn**(expE_SN)* Z_neisol**(expZ_SN)  !km/s Msun

              if(rt)then
                 fthor   = exp(-dx_loc_pc/rStrom(i))
                 psn_tr  = psn_thor98*fthor + psn_geen15*(1d0-fthor)
              else
                 psn_tr  = psn_geen15
              endif
              psn_tr = max(psn_tr, psn_thor98)

              ! For adiabatic phase
              ! psn_tr =  A_SN * (E51 * boost_geen)**expE_SN_boost * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
              !        =  p_hydro * boost_geen**expE_SN_boost
              p_hydro = A_SN * num_sn**(expE_SN_boost) * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
              boost_geen_ad = (psn_tr / p_hydro)**(1d0/expE_SN_boost)
              boost_geen_ad = max(boost_geen_ad-1d0,0.0)
           endif

           chi_tr   = psn_tr**2d0 / (2d0 * num_sn**2d0 * (E_SNII/msun2g/km2cm**2d0) * M_SNII_var * f_ESN)
           !          (Msun*km/s)^2                 (Msun *km2/s2)              (Msun)
           f_w_crit = max(chi_tr-1d0, 0d0)

           !f_w_crit = (A_SN/1d4)**2d0/(f_ESN*M_SNII)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
           !f_w_crit = max(0d0,f_w_crit)
           vload_rad = dsqrt(2d0*f_ESN*E_SNII*(1d0+f_w_crit)/(M_SNII_var*msun2g))/scale_v/(1d0+f_w_cell)/f_LOAD/f_PCAN
           if(f_w_cell.ge.f_w_crit)then ! radiative phase
              ! ptot = sqrt(2*chi_tr*Mejtot*(fe*Esntot))
              ! vload = ptot/(chi*Mejtot) = sqrt(2*chi_tr*fe*Esntot/Mejtot)/chi = sqrt(2*chi_tr*fe*Esn/Mej)/chi
              vload = vload_rad
              snowplough(i,j)=.true.
           else ! adiabatic phase
              ! ptot = sqrt(2*chi*Mejtot*(fe*Esntot))
              ! vload = ptot/(chi*Mejtot) = sqrt(2*fe*Esntot/chi/Mejtot) = sqrt(2*fe*Esn/chi/Mej)
              f_esn2 = 1d0-(1d0-f_ESN)*f_w_cell/f_w_crit ! to smoothly correct the adibatic to the radiative phase
              vload = dsqrt(2d0*f_esn2*E_SNII/(1d0+f_w_cell)/(M_SNII_var*msun2g))/scale_v/f_LOAD
              if(mechanical_geen) then
                 !f_wrt_snow = (f_esn2 - f_ESN)/(1d0-f_ESN)
                 f_wrt_snow = 2d0-2d0/(1d0+exp(-f_w_cell/f_w_crit/0.3)) ! 0.3 is obtained by calibrating
                 vload = vload * dsqrt(1d0 + boost_geen_ad*f_wrt_snow)
                 ! NB. this sometimes give too much momentum because expE_SN_boost != expE_SN. A limiter is needed
              endif
              ! safety device: limit the maximum velocity so that it does not exceed p_{SN,final}
              if(vload>vload_rad) vload = vload_rad
              snowplough(i,j)=.false.
           endif
           p_solid(i,j)=(1d0+f_w_cell)*dm_ejecta*vload
           ek_solid(i,j)=ek_solid(i,j)+p_solid(i,j)*(vload*f_LOAD)/2d0 !ek=(m*v)*v/2, not (d*v)*v/2

           if(log_mfb_mega)then
             write(*,'(" MFBN nHcen=", f6.2," nHnei=", f6.2, " mej=", f6.2, " mcen=", f6.2, " vload=", f6.2, " lv2=",I3," fwcrit=", f6.2, " fwcell=", f6.2, " psol=",f6.2, " mload/48=",f6.2," mnei/8=",f6.2," mej/48=",f6.2)') &
            & log10(nH_cen),log10(nH_nei),log10(mSN(i)*scale_msun),log10(m_cen),log10(vload*scale_v/1d5),ilevel2-ilevel,&
            & log10(f_w_crit),log10(f_w_cell),log10(p_solid(i,j)*scale_msun*scale_v/1d5),log10(mload*scale_msun/48),&
            & log10(d_nei*vol_loc/8d0*scale_msun),log10(dm_ejecta*scale_msun)
            endif

        endif
        
     enddo ! loop over neighboring cells
  enddo ! loop over SN cells

  !-----------------------------------------
  ! Redistribute mass from the SN cell
  !-----------------------------------------
  do i=1,np
     icell = ncoarse+ind_grid(i)+(ind_pos_cell(i)-1)*ngridmax
     d     = uold(icell,1)
     u     = uold(icell,2)/d
     v     = uold(icell,3)/d
     w     = uold(icell,4)/d
     e     = uold(icell,5)
#if NENER>0
     do irad=1,nener
        e = e - uold(icell,ndim+2+irad)
     enddo
#endif
     emag=0.0d0
#ifdef SOLVERmhd
     do idim=1,ndim
        emag=emag+0.125d0*(uold(icell,idim+ndim+2)+uold(icell,idim+nvar))**2
     enddo
#endif

     ekk   = 0.5*d*(u**2 + v**2 + w**2)
     eth   = e - ekk - emag  ! thermal pressure 

     ! ionisation fractions, ref, etc.
     do ii=i_fractions,nvar
        fractions(ii) = uold(icell,ii)/d
     end do

     mloadSN (i) = mSN (i)*f_LOAD + d*vol_loc*floadSN(i)
     if(metal)then
        z = uold(icell,imetal)/d
        mZloadSN(i) = mZSN(i)*f_LOAD + d*z*vol_loc*floadSN(i)
     endif

     do ich=1,nchem
        z_ch(ich) = uold(icell,ichem+ich-1)/d
        mchloadSN(i,ich) = mchSN(i,ich)*f_LOAD + d*z_ch(ich)*vol_loc*floadSN(i)
     end do 

     ! original momentum by star + gas entrained from the SN cell
     ploadSN(i,1) = pSN(i,1)*f_LOAD + vol_loc*d*u*floadSN(i)
     ploadSN(i,2) = pSN(i,2)*f_LOAD + vol_loc*d*v*floadSN(i)
     ploadSN(i,3) = pSN(i,3)*f_LOAD + vol_loc*d*w*floadSN(i)

     ! update the hydro variable
     fleftSN = 1d0 - floadSN(i)
     uold(icell,1) = uold(icell,1)*fleftSN + mSN(i)  /vol_loc*f_LEFT 
     uold(icell,2) = uold(icell,2)*fleftSN + pSN(i,1)/vol_loc*f_LEFT  ! rho*v, not v
     uold(icell,3) = uold(icell,3)*fleftSN + pSN(i,2)/vol_loc*f_LEFT 
     uold(icell,4) = uold(icell,4)*fleftSN + pSN(i,3)/vol_loc*f_LEFT

     if(metal)then
        uold(icell,imetal) = mZSN(i)/vol_loc*f_LEFT + d*z*fleftSN
     endif
     do ich=1,nchem
        uold(icell,ichem+ich-1) = mchSN(i,ich)/vol_loc*f_LEFT + d*z_ch(ich)*fleftSN
     end do
     do ii=i_fractions,nvar
        uold(icell,ii) = fractions(ii) * uold(icell,1)
     end do

     ! original kinetic energy of the gas entrained
     eloadSN(i) = ekk*vol_loc*floadSN(i) 

     ! original thermal energy of the gas entrained (including the non-thermal part)
     eloadSN(i) = eloadSN(i) + eth*vol_loc*floadSN(i)

     ! reduce total energy as we are distributing it to the neighbours
     !uold(icell,5) = uold(icell,5)*fleftSN 
     uold(icell,5) = uold(icell,5) - (ekk+eth)*floadSN(i)

     ! add the contribution from the original kinetic energy of SN particle
     d = mSN(i)/vol_loc
     u = pSN(i,1)/mSN(i)
     v = pSN(i,2)/mSN(i)
     w = pSN(i,3)/mSN(i)
     uold(icell,5) = uold(icell,5) + 0.5d0*d*(u**2 + v**2 + w**2)*f_LEFT
    
     ! add the contribution from the original kinetic energy of SN to outflow
     eloadSN(i) = eloadSN(i) + 0.5d0*mSN(i)*(u**2 + v**2 + w**2)*f_LOAD

     ! update ek_solid     
     ek_solid(i,:) = ek_solid(i,:) + eloadSN(i)/dble(nSNnei)

  enddo  ! loop over SN cell


  !-------------------------------------------------------------
  ! Find and save stars affecting across the boundary of a cpu
  !-------------------------------------------------------------
  do i=1,np

     ind_cell = ncoarse+ind_grid(i)+(ind_pos_cell(i)-1)*ngridmax

     nwco=0; icpuSNnei=0
     do j=1,nSNnei
        call get_icell_from_pos (xc2(1:3,i)+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell,ilevel2)
 
        if(cpu_map(father(igrid)).ne.myid) then ! need mpi
           nwco=nwco+1
           icpuSNnei(nwco)=cpu_map(father(igrid))
        else  ! can be handled locally
           vol_nei = vol_loc*(2d0**ndim)**(ilevel-ilevel2)
           pvar(:) = 0d0 ! temporary primitive variable
           if(ilevel>ilevel2)then ! touching level-1 cells
              pvar(1:nvarMHD) = unew(icell,1:nvarMHD)
           else
              pvar(1:nvarMHD) = uold(icell,1:nvarMHD)
           endif 
           do ii=i_fractions,nvar ! fractional quantities that we don't want to change
              fractions(ii) = pvar(ii)/pvar(1)
           end do

           emag=0.0d0
#ifdef SOLVERmhd
           do idim=1,ndim
              emag=emag+0.125d0*(pvar(idim+ndim+2)+pvar(idim+nvar))**2
           enddo
#endif
           d0=pvar(1)
           u0=pvar(2)/d0
           v0=pvar(3)/d0
           w0=pvar(4)/d0
           ekk0=0.5d0*d0*(u0**2+v0**2+w0**2)
           eth0=pvar(5)-ekk0-emag
#if NENER>0
           do irad=1,nener
              eth0=eth0-pvar(ndim+2+irad)
           end do
#endif
           ! For stability
           Tk0 =eth0/d0*scale_T2*(gamma-1.0)
           T2min=T2_star*(d0*scale_nH/n_star)**(g_star-1.0)
           if(Tk0<T2min)then
              eth0=T2min*d0/scale_T2/(gamma-1.0)
           endif 

           d= mloadSN(i  )/dble(nSNnei)/vol_nei
           u=(ploadSN(i,1)/dble(nSNnei)+p_solid(i,j)*vSNnei(1,j))/vol_nei/d
           v=(ploadSN(i,2)/dble(nSNnei)+p_solid(i,j)*vSNnei(2,j))/vol_nei/d
           w=(ploadSN(i,3)/dble(nSNnei)+p_solid(i,j)*vSNnei(3,j))/vol_nei/d
           pvar(1)=pvar(1)+d
           pvar(2)=pvar(2)+d*u
           pvar(3)=pvar(3)+d*v
           pvar(4)=pvar(4)+d*w
           ekk_ej = 0.5*d*(u**2 + v**2 + w**2)   
           etot0  = eth0+ekk0+emag+ek_solid(i,j)/vol_nei ! additional energy from SNe+entrained gas

           ! the minimum thermal energy input floor
           d   = pvar(1)
           u   = pvar(2)/d
           v   = pvar(3)/d
           w   = pvar(4)/d
           ekk = 0.5*d*(u**2 + v**2 + w**2)
           pvar(5) = max(etot0, ekk+eth0+emag)

           ! sanity check
           Tk = (pvar(5)-ekk-emag)/d*scale_T2*(gamma-1)
           if(Tk<0)then
              print *,'TKERR: mech (post-call): Tk<0 =',Tk
              print *,'nH [H/cc]= ',d*scale_nH
              print *,'u  [km/s]= ',u*scale_v/1d5
              print *,'v  [km/s]= ',v*scale_v/1d5
              print *,'w  [km/s]= ',w*scale_v/1d5
              print *,'T0 [K]   = ',Tk0
              stop
           endif 

#if NENER>0
           do irad=1,nener
              pvar(5) = pvar(5) + pvar(ndim+2+irad)
           end do
#endif
           if(metal)then
               pvar(imetal)=pvar(imetal)+mzloadSN(i)/dble(nSNnei)/vol_nei
           end if
           do ich=1,nchem
               pvar(ichem+ich-1)=pvar(ichem+ich-1)+mchloadSN(i,ich)/dble(nSNnei)/vol_nei
           end do
           do ii=i_fractions,nvar
               pvar(ii)=fractions(ii)*pvar(1)
           end do

           ! update the hydro variable
           if(ilevel>ilevel2)then ! touching level-1 cells
              unew(icell,1:nvarMHD) = pvar(1:nvarMHD)
           else
              uold(icell,1:nvarMHD) = pvar(1:nvarMHD)
           endif 

        end if
     end do ! loop over 48 neighbors


#ifndef WITHOUTMPI 
     if(nwco>0)then  ! for SNs across different cpu
        if(nwco>1)then
           nwco_here=nwco
           ! remove redundant cpu list for this SN cell
           call redundant_non_1d(icpuSNnei(1:nwco), nwco_here, nwco)
        endif
        ista=ncomm_SN+1
        iend=ista+nwco-1

        if(iend>ncomm_max)then
           write(*,*) 'Error: increase ncomm_max in mechanical_fine.f90', ncomm_max, iend
           call clean_stop
        endif
        iSN_comm (ista:iend)=icpuSNnei(1:nwco)
        nSN_comm (ista:iend )=nSN(i)
        mSN_comm (ista:iend )=mSN(i)
        mloadSN_comm (ista:iend  )=mloadSN(i)
        xSN_comm (1,ista:iend)=xc2(1,i)
        xSN_comm (2,ista:iend)=xc2(2,i)
        xSN_comm (3,ista:iend)=xc2(3,i)
        ploadSN_comm (1,ista:iend)=ploadSN(i,1)
        ploadSN_comm (2,ista:iend)=ploadSN(i,2)
        ploadSN_comm (3,ista:iend)=ploadSN(i,3)
        floadSN_comm (ista:iend  )=floadSN(i)
        eloadSN_comm (ista:iend  )=eloadSN(i)
        if(metal) mZloadSN_comm(ista:iend  )=mZloadSN(i)
        do ich=1,nchem
            mchloadSN_comm(ista:iend,ich)=mchloadSN(i,ich)
        end do
        if(mechanical_geen.and.rt) rSt_comm (ista:iend)=rStrom(i)
        ncomm_SN=ncomm_SN+nwco
     endif
#endif

  end do ! loop over SN cell

end subroutine mech_fine
!################################################################
!################################################################
!################################################################
subroutine mech_fine_mpi(ilevel)
  use amr_commons
  use mechanical_commons
  use pm_commons
  use hydro_commons
  implicit none
#ifndef WITHOUTMPI
  include 'mpif.h'
  integer::i,j,info,nSN_tot,icpu,ncpu_send,ncpu_recv,ncc
  integer::ncell_recv,ncell_send,cpu2send,cpu2recv,tag,np
  integer::isend_sta,irecv_sta,irecv_end
  real(dp),dimension(:,:),allocatable::SNsend,SNrecv,p_solid,ek_solid
  integer ,dimension(:),allocatable::list2recv,list2send
  integer, dimension(:),allocatable::reqrecv,reqsend
  integer, dimension(:,:),allocatable::statrecv,statsend
  ! SN variables
  real(dp)::dx,dx_loc,scale,vol_loc
  real(dp)::mloadSN_i,zloadSN_i,ploadSN_i(1:3),mSN_i,xSN_i(1:3),fload_i
  real(dp)::f_esn2,d_nei,Z_nei,Z_neisol,f_w_cell,f_w_crit,nH_nei,dratio
  real(dp)::num_sn,vload,Tk,vol_nei,dm_ejecta
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_kms
  real(dp)::scale_msun,msun2g=1.989d33, pvar(1:nvarMHD),etot0
  real(dp)::skip_loc(1:3),d,u,v,w,ekk,eth,d0,u0,v0,w0,eth0,ekk0,Tk0,ekk_ej,T2min
  integer::igrid,icell,ilevel,ilevel2,irad,ii,ich
  real(dp),dimension(1:twotondim,1:ndim),save::xc
  logical,allocatable,dimension(:,:)::snowplough
  real(dp),dimension(1:nchem),save::chloadSN_i
  real(dp)::rSt_i,dx_loc_pc,psn_tr,chi_tr,psn_thor98,psn_geen15,fthor
  real(dp)::km2cm=1d5,M_SNII_var,boost_geen_ad=0d0,p_hydro,vload_rad,f_wrt_snow
  ! fractional abundances ; for ionisation fraction and ref, etc
  real(dp),dimension(1:NVAR),save::fractions ! not compatible with delayed cooling
  integer::i_fractions,idim
  real(dp)::emag

  if(ndim.ne.3) return

  ! starting index for passive variables except for imetal and chem
  i_fractions = imetal+nchem+1

  !============================================================
  ! For MPI communication
  !============================================================
  ncpu_send=0;ncpu_recv=0

  ncomm_SN_cpu=0 
  ncomm_SN_mpi=0
  ncomm_SN_mpi(myid)=ncomm_SN
  ! compute the total number of communications needed
  call MPI_ALLREDUCE(ncomm_SN_mpi,ncomm_SN_cpu,ncpu,&
                   & MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)
  nSN_tot = sum(ncomm_SN_cpu)
  if(nSN_tot==0) return

  allocate(icpuSN_comm    (1:nSN_tot,1:2))
  allocate(icpuSN_comm_mpi(1:nSN_tot,1:2))

  ! index for mpi variable
  if(myid==1)then
     isend_sta = 0 
  else
     isend_sta = sum(ncomm_SN_cpu(1:myid-1)) 
  endif

  icpuSN_comm=0
  do i=1,ncomm_SN_cpu(myid)
     icpuSN_comm(isend_sta+i,1)=myid
     icpuSN_comm(isend_sta+i,2)=iSN_comm(i)
     ! iSN_comm:   local variable
     ! icpuSN_comm:  local (but extended) variable to be passed to a mpi variable
     ! icpuSN_comm_mpi: mpi variable
  end do

  ! share the list of communications
  icpuSN_comm_mpi=0
  call MPI_ALLREDUCE(icpuSN_comm,icpuSN_comm_mpi,nSN_tot*2,&
                   & MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,info)

  ncell_send = ncomm_SN_cpu(myid)
  ncell_recv = count(icpuSN_comm_mpi(:,2).eq.myid, 1)

  ! check if myid needs to send anything
  if(ncell_send>0)then
     allocate( SNsend (1:nvarSN,1:ncell_send)) ! x(3),m,mz,p,pr
     allocate( list2send(1:ncell_send)    )
     list2send=0;  SNsend=0d0
     list2send=icpuSN_comm_mpi(isend_sta+1:isend_sta+ncell_send,2)
     ncpu_send=1
     if(ncell_send>1) call redundant_non_1d (list2send,ncell_send,ncpu_send)
     ! ncpu_send = No. of cpus to which myid should send info 
     allocate( reqsend  (1:ncpu_send) )
     allocate( statsend (1:MPI_STATUS_SIZE,1:ncpu_send) )
     reqsend=0; statsend=0
  endif

  ! check if myid needs to receive anything
  if(ncell_recv>0)then
     allocate( SNrecv (1:nvarSN,1:ncell_recv)) ! x(3),m,mz,p,pr
     allocate( list2recv(1:ncell_recv)    )
     list2recv=0;  SNrecv=0d0
     j=0
     do i=1,nSN_tot
        if(icpuSN_comm_mpi(i,2).eq.myid)then
           j=j+1 
           list2recv(j) = icpuSN_comm_mpi(i,1)
        endif
     end do
 
     ncc = j
     if(ncc.ne.ncell_recv)then ! sanity check
        write(*,*) 'Error in mech_fine_mpi: ncc != ncell_recv',ncc,ncell_recv,myid
        call clean_stop
     endif

     ncpu_recv=1 ! No. of cpus from which myid should receive info 
     if(j>1) call redundant_non_1d(list2recv,ncc,ncpu_recv)

     allocate( reqrecv  (1:ncpu_recv) )
     allocate( statrecv (1:MPI_STATUS_SIZE,1:ncpu_recv) )
     reqrecv=0; statrecv=0
  endif

  ! prepare one variable and send
  if(ncell_send>0)then
     do icpu=1,ncpu_send
        cpu2send = list2send(icpu)
        ncc=0 ! number of SN host cells that need communications with myid=cpu2send
        do i=1,ncell_send
           j=i+isend_sta
           if(icpuSN_comm_mpi(j,2).eq.cpu2send)then
              ncc=ncc+1
              SNsend(1:3,ncc)=xSN_comm (1:3,i)
              SNsend(4  ,ncc)=mSN_comm (    i)
              SNsend(5  ,ncc)=mloadSN_comm (i)
              SNsend(6:8,ncc)=ploadSN_comm (1:3,i)
              SNsend(9  ,ncc)=floadSN_comm (i)
              SNsend(10 ,ncc)=eloadSN_comm (i)
              SNsend(11 ,ncc)=nSN_comm     (i)
              if(metal)SNsend(12,ncc)=mZloadSN_comm(i)
              if(mechanical_geen.and.rt)SNsend(13,ncc)=rSt_comm(i)
              do ich=1,nchem
                 SNsend(13+ich,ncc)=mchloadSN_comm(i,ich)
              end do
           endif
        end do ! i

        tag = myid + cpu2send + ncc
        call MPI_ISEND (SNsend(1:nvarSN,1:ncc),ncc*nvarSN,MPI_DOUBLE_PRECISION, &
                      & cpu2send-1,tag,MPI_COMM_WORLD,reqsend(icpu),info) 
     end do ! icpu

  endif ! ncell_send>0


  ! receive one large variable
  if(ncell_recv>0)then
     irecv_sta=1
     do icpu=1,ncpu_recv
        cpu2recv = list2recv(icpu)
        ncc=0 ! number of SN host cells that need communications with cpu2recv
        do i=1,nSN_tot
           if((icpuSN_comm_mpi(i,1)==cpu2recv).and.&
             &(icpuSN_comm_mpi(i,2)==myid)      )then
              ncc=ncc+1
           endif
        end do
        irecv_end=irecv_sta+ncc-1
        tag = myid + cpu2recv + ncc
        
        call MPI_IRECV (SNrecv(1:nvarSN,irecv_sta:irecv_end),ncc*nvarSN,MPI_DOUBLE_PRECISION,&
                     & cpu2recv-1,tag,MPI_COMM_WORLD,reqrecv(icpu),info)

        irecv_sta=irecv_end+1
     end do ! icpu 

  endif ! ncell_recv >0

  if(ncpu_send>0)call MPI_WAITALL(ncpu_send,reqsend,statsend,info)
  if(ncpu_recv>0)call MPI_WAITALL(ncpu_recv,reqrecv,statrecv,info)


  !============================================================
  ! inject mass/metal/momentum
  !============================================================

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  scale_msun = scale_l**3*scale_d/msun2g
  scale_kms  = scale_v/1d5
 
  ! Mesh variables
  call mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)
  dx_loc_pc = dx_loc*scale_l/3.08d18

  np = ncell_recv
  if(ncell_recv>0) then
     allocate(p_solid(1:np,1:nSNnei))
     allocate(ek_solid(1:np,1:nSNnei))
     allocate(snowplough(1:np,1:nSNnei))
     p_solid=0d0;ek_solid=0d0;snowplough=.false.
  endif

  ! Compute the momentum first before redistributing mass
  do i=1,np

     xSN_i(1:3) = SNrecv(1:3,i)
     mSN_i      = SNrecv(4,i)
     mloadSN_i  = SNrecv(5,i)
     fload_i    = SNrecv(9,i)
     dm_ejecta  = f_LOAD*mSN_i/dble(nSNnei)
     num_sn     = SNrecv(11,i)
     ek_solid(i,:) = SNrecv(10,i)/dble(nSNnei) ! kinetic energy of the gas mass entrained from the host cell + SN
     if(metal) ZloadSN_i = SNrecv(12,i)/SNrecv(5,i)
     if(mechanical_geen.and.rt) rSt_i = SNrecv(13,i) ! Stromgren sphere in pc 

     M_SNII_var = mSN_i*scale_msun/num_sn

     do j=1,nSNnei
        call get_icell_from_pos (xSN_i+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell,ilevel2)
        if(cpu_map(father(igrid)).eq.myid) then ! if belong to myid
           Z_nei = z_ave*0.02 ! for metal=.false.
           if(ilevel>ilevel2)then ! touching level-1 cells
              d_nei     = unew(icell,1)
              if(metal) Z_nei = unew(icell,imetal)/d_nei
           else
              d_nei     = uold(icell,1)
              if(metal) Z_nei = uold(icell,imetal)/d_nei
           endif
           if(loading_type.eq.1.and.fload_i<f_LOAD)then
              dratio = fload_i
           else
              dratio = 1d0
           endif
           f_w_cell = (mloadSN_i/dble(nSNnei) + d_nei*vol_loc/8d0)/dm_ejecta - 1d0
           nH_nei   = d_nei*scale_nH*dratio
           Z_neisol = max(0.01,Z_nei/0.02)

           ! transition mass loading factor (momentum conserving phase)
           ! psn_tr = sqrt(2*chi_tr*Nsn*Esn*Msn*fe)
           ! chi_tr = (1+f_w_crit)
           psn_thor98   = A_SN * num_sn**(expE_SN) * nH_nei**(expN_SN) * Z_neisol**(expZ_SN)  !km/s Msun
           psn_tr       = psn_thor98 
           if(mechanical_geen)then
              ! For snowplough phase, psn_tr will do the job
              psn_geen15= A_SN_Geen * num_sn**(expE_SN) * Z_neisol**(expZ_SN)  !km/s Msun
              if(rt)then
                 fthor   = exp(-dx_loc_pc/rSt_i)
                 psn_tr  = psn_thor98*fthor + psn_geen15*(1d0-fthor)
              else
                 psn_tr  = psn_geen15
              endif
              psn_tr = max(psn_tr, psn_thor98)

              ! For adiabatic phase
              ! psn_tr =  A_SN * (E51 * boost_geen)**expE_SN_boost * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
              !        =  p_hydro * boost_geen**expE_SN_boost
              p_hydro = A_SN * num_sn**(expE_SN_boost) * nH_nei**(expN_SN_boost) * Z_neisol**(expZ_SN)
              boost_geen_ad = (psn_tr / p_hydro)**(1d0/expE_SN_boost)
              boost_geen_ad = max(boost_geen_ad-1d0,0.0)
           endif
           chi_tr   = psn_tr**2d0 / (2d0 * num_sn**2d0 * (E_SNII/msun2g/km2cm**2d0) * M_SNII_var * f_ESN)
           !          (Msun*km/s)^2                 (Msun *km2/s2)              (Msun)
           f_w_crit = max(chi_tr-1d0, 0d0)
 
           !f_w_crit = (A_SN/1d4)**2d0/(f_ESN*M_SNII)*num_sn**((expE_SN-1d0)*2d0)*nH_nei**(expN_SN*2d0)*Zdepen - 1d0
           !f_w_crit = max(0d0,f_w_crit)

           vload_rad = dsqrt(2d0*f_ESN*E_SNII*(1d0+f_w_crit)/(M_SNII_var*msun2g))/(1d0+f_w_cell)/scale_v/f_LOAD/f_PCAN
           if(f_w_cell.ge.f_w_crit)then ! radiative phase
              vload = vload_rad
              snowplough(i,j)=.true.
           else ! adiabatic phase
              f_esn2 = 1d0-(1d0-f_ESN)*f_w_cell/f_w_crit
              vload = dsqrt(2d0*f_esn2*E_SNII/(1d0+f_w_cell)/(M_SNII_var*msun2g))/scale_v/f_LOAD
              if(mechanical_geen) then
                 f_wrt_snow = 2d0-2d0/(1d0+exp(-f_w_cell/f_w_crit/0.30))
                 boost_geen_ad = boost_geen_ad * f_wrt_snow
                 vload = vload * dsqrt(1d0+boost_geen_ad) ! f_boost
                 ! NB. this sometimes give too much momentum because expE_SN_boost != expE_SN. A limiter is needed
              endif
              ! safety device: limit the maximum velocity so that it does not exceed p_{SN,final}
              if(vload>vload_rad) vload = vload_rad
              snowplough(i,j)=.false.
           endif
           p_solid(i,j)=(1d0+f_w_cell)*dm_ejecta*vload
           ek_solid(i,j)=ek_solid(i,j)+p_solid(i,j)*(vload*f_LOAD)/2d0 !ek
        endif
       
     enddo ! loop over neighboring cells

  end do



  ! Apply the SNe to this cpu domain
  do i=1,np

     xSN_i(1:3) = SNrecv(1:3,i)
     mSN_i      = SNrecv(4,i)
     mloadSN_i  = SNrecv(5,i)
     ploadSN_i(1:3)=SNrecv(6:8,i)
     if(metal) ZloadSN_i = SNrecv(12,i)/mloadSN_i
     do ich=1,nchem
        chloadSN_i(ich) = SNrecv(13+ich,i)/mloadSN_i
     end do

     do j=1,nSNnei
 
        call get_icell_from_pos (xSN_i(1:3)+xSNnei(1:3,j)*dx, ilevel+1, igrid, icell, ilevel2)
        if(cpu_map(father(igrid))==myid)then

           vol_nei = vol_loc*(2d0**ndim)**(ilevel-ilevel2)
           if(ilevel>ilevel2)then ! touching level-1 cells
              pvar(1:nvarMHD) = unew(icell,1:nvarMHD)
           else
              pvar(1:nvarMHD) = uold(icell,1:nvarMHD)
           endif
           do ii=i_fractions,nvar
              fractions(ii)=pvar(ii)/pvar(1)
           enddo
     
           emag=0.0d0
#ifdef SOLVERmhd
           do idim=1,ndim
              emag=emag+0.125d0*(pvar(idim+ndim+2)+pvar(idim+nvar))**2
           enddo
#endif
           d0=pvar(1)
           u0=pvar(2)/d0
           v0=pvar(3)/d0
           w0=pvar(4)/d0
           ekk0=0.5d0*d0*(u0**2d0 + v0**2d0 + w0**2d0)
           eth0=pvar(5)-ekk0-emag
#if NENER>0
           do irad=1,nener
              eth0 = eth0 - pvar(ndim+2+irad)
           end do
#endif
           Tk0 = eth0/d0*scale_T2*(gamma-1)
           T2min = T2_star*(d0*scale_nH/n_star)**(g_star-1.0)
           if(Tk0<T2min)then
              eth0 = T2min*d0/scale_T2/(gamma-1)
           endif 

           d= mloadSN_i   /dble(nSNnei)/vol_nei
           u=(ploadSN_i(1)/dble(nSNnei)+p_solid(i,j)*vSNnei(1,j))/vol_nei/d
           v=(ploadSN_i(2)/dble(nSNnei)+p_solid(i,j)*vSNnei(2,j))/vol_nei/d
           w=(ploadSN_i(3)/dble(nSNnei)+p_solid(i,j)*vSNnei(3,j))/vol_nei/d

           pvar(1)=pvar(1)+d
           pvar(2)=pvar(2)+d*u
           pvar(3)=pvar(3)+d*v
           pvar(4)=pvar(4)+d*w
           ekk_ej = 0.5*d*(u**2 + v**2 + w**2)
           etot0  = eth0+ekk0+emag+ek_solid(i,j)/vol_nei  ! additional energy from SNe+entrained gas

           ! the minimum thermal energy input floor
           d   = pvar(1)
           u   = pvar(2)/d
           v   = pvar(3)/d
           w   = pvar(4)/d
           ekk = 0.5*d*(u**2 + v**2 + w**2)

           pvar(5) = max(etot0, ekk+eth0+emag)

#if NENER>0
           do irad=1,nener
              pvar(5) = pvar(5) + pvar(ndim+2+irad)
           end do
#endif

           if(metal)then
              pvar(imetal)=pvar(imetal)+mloadSN_i/dble(nSNnei)*ZloadSN_i/vol_nei
           end if
           do ich=1,nchem
              pvar(ichem+ich-1)=pvar(ichem+ich-1)+mloadSN_i/dble(nSNnei)*chloadSN_i(ich)/vol_nei
           end do
           do ii=i_fractions,nvar
              pvar(ii)=fractions(ii)*pvar(1)
           end do

           ! update the hydro variable
           if(ilevel>ilevel2)then ! touching level-1 cells
              unew(icell,1:nvarMHD) = pvar(1:nvarMHD)
           else
              uold(icell,1:nvarMHD) = pvar(1:nvarMHD)
           endif
 
        endif ! if this belongs to me


     end do ! loop over neighbors
  end do ! loop over SN cells


  deallocate(icpuSN_comm_mpi, icpuSN_comm)
  if(ncell_send>0) deallocate(list2send,SNsend,reqsend,statsend)
  if(ncell_recv>0) deallocate(list2recv,SNrecv,reqrecv,statrecv)
  if(ncell_recv>0) deallocate(p_solid,ek_solid,snowplough)

  ncomm_SN=nSN_tot
#endif

end subroutine mech_fine_mpi
!################################################################
!################################################################
!################################################################
subroutine get_number_of_sn2(birth_time,dteff,zp_star,id_star,mass0,mass_t,M_SNII_var,nsn,done_star)
  use amr_commons, ONLY:dp,snII_freq,eta_sn,sn2_real_delay
  use random
  implicit none
  real(kind=dp)::birth_time,zp_star,mass0,mass_t,dteff ! birth_time in code, mass in Msun
  real(kind=dp)::nsn
  integer::nsn_tot,nsn_sofar,nsn_age2
  integer::i,localseed,id_star   ! necessary for the random number
  real(kind=dp)::age1,age2,tMyr,xdum,ydum,logzpsun
  real(kind=dp)::co0,co1,co2,M_SNII_var
  ! fit to the cumulative number fraction for Kroupa IMF
  ! ~kimm/soft/starburst99/output_hires_new/snrate.pro 
  ! ~kimm/soft/starburst99/output_hires_new/fit.pro 
  ! ~kimm/soft/starburst99/output_hires_new/sc.pro 
  real(kind=dp),dimension(1:3)::coa=(/-2.677292E-01,1.392208E-01,-5.747332E-01/)
  real(kind=dp),dimension(1:3)::cob=(/4.208666E-02, 2.152643E-02, 7.893866E-02/)
  real(kind=dp),dimension(1:3)::coc=(/-8.612668E-02,-1.698731E-01,1.867337E-01/)
  real(kind=dp),external::ran1 
  logical:: done_star


  ! determine the SNrateCum curve for Zp
  logzpsun=max(min(zp_star,0.05),0.008) ! no extrapolation
  logzpsun=log10(logzpsun/0.02)

  co0 = coa(1)+coa(2)*logzpsun+coa(3)*logzpsun**2
  co1 = cob(1)+cob(2)*logzpsun+cob(3)*logzpsun**2
  co2 = coc(1)+coc(2)*logzpsun+coc(3)*logzpsun**2

  ! RateCum = co0+sqrt(co1*Myr+co2)

  ! get stellar age
  call getStarAgeGyr(birth_time+dteff, age1)
  call getStarAgeGyr(birth_time      , age2)

  ! convert Gyr -> Myr
  age1 = age1*1d3
  age2 = age2*1d3

  nsn=0d0; done_star=.false.
  if(age2.le.(-co2/co1))then
     return
  endif

  ! total number of SNe
  nsn_tot   = NINT(mass0*snII_freq,kind=4)
  if(nsn_tot.eq.0)then
      write(*,*) 'Fatal error: please increase the mass of your star particle'
      stop
  endif

  ! number of SNe up to this point
  nsn_sofar = NINT((mass0-mass_t)/M_SNII_var,kind=4)
  if(nsn_sofar.ge.nsn_tot)then
     done_star=.true.
     return
  endif

  localseed = -abs(id_star) !make sure that the number is negative

  nsn_age2 = 0
  do i=1,nsn_tot
     xdum =  ran1(localseed)
     ! inverse function for y=co0+sqrt(co1*x+co2)
     ydum = ((xdum-co0)**2.-co2)/co1
     if(ydum.le.age2)then
        nsn_age2=nsn_age2+1
     endif
  end do

  nsn_age2 = min(nsn_age2,nsn_tot)  ! just in case...
  nsn = max(nsn_age2 - nsn_sofar, 0) 

  if(nsn_age2.ge.nsn_tot) done_star=.true.

end subroutine get_number_of_sn2
!################################################################
!################################################################
!################################################################
function ran1(idum)
   implicit none
   integer:: idum,IA,IM,IQ,IR,NTAB,NDIV
   real(kind=8):: ran1,AM,EPS,RNMX
   parameter(IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,&
            &NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
   integer::j,k,iv(NTAB),iy
   save iv,iy
   data iv /NTAB*0/, iy /0/ 
   if (idum.le.0.or.iy.eq.0) then! initialize
      idum=max(-idum,1)
      do j=NTAB+8,1,-1
         k=idum/IQ
         idum=IA*(idum-k*IQ)-IR*k
         if (idum.lt.0) idum=idum+IM
         if (j.le.NTAB) iv(j)=idum
      end do
      iy=iv(1)
   end if
   k=idum/IQ
   idum=IA*(idum-k*IQ)-IR*k
   if (idum.lt.0) idum=idum+IM
   j=1+iy/NDIV
   iy=iv(j)
   iv(j)=idum
   ran1=min(AM*iy,RNMX)
   return
end
!################################################################
!################################################################
!################################################################
!################################################################
subroutine getStarAgeGyr(birth_time,age_star)
   use amr_commons
   implicit none
   real(dp)::age_star,age_star_pt,birth_time

   if (use_proper_time)then
      call getAgeGyr    (birth_time, age_star)
   else
      call getProperTime(birth_time, age_star_pt)
      call getAgeGyr    (age_star_pt,age_star)
   endif

end subroutine getStarAgeGyr
!################################################################
!################################################################
!################################################################
subroutine redundant_non_1d(list,ndata,ndata2)
   implicit none
   integer,dimension(1:ndata)::list,list2
   integer,dimension(1:ndata)::ind
   integer::ndata,i,y1,ndata2

   ! sort first
   call heapsort(list,ind,ndata) 

   list(:)=list(ind)

   ! check the redundancy
   list2(:) = list(:)

   y1=list(1)
   ndata2=1
   do i=2,ndata
       if(list(i).ne.y1)then
          ndata2=ndata2+1
          list2(ndata2)=list(i)
          y1=list(i)
       endif
   end do

   list =list2

end subroutine redundant_non_1d
!################################################################
!################################################################
!################################################################
subroutine heapsort(ain,ind,n)
   implicit none
   integer::n
   integer,dimension(1:n)::ain,aout,ind
   integer::i,j,l,ir,idum,rra

   l=n/2+1
   ir=n
   do i=1,n
      aout(i)=ain(i)                        ! Copy input array to output array
      ind(i)=i                                   ! Generate initial idum array
   end do
   if(n.eq.1) return                            ! Special for only one record
10 continue
   if(l.gt.1)then
      l=l-1
      rra=aout(l)
      idum=ind(l)
   else
      rra=aout(ir)
      idum=ind(ir)
      aout(ir)=aout(1)
      ind(ir)=ind(1)
      ir=ir-1
      if(ir.eq.1)then
        aout(1)=rra
        ind(1)=idum
        return
      endif
    endif
    i=l
    j=l+l
20  if(j.le.ir)then
       if(j.lt.ir)then
          if(aout(j).lt.aout(j+1))j=j+1
       endif
       if(rra.lt.aout(j))then
          aout(i)=aout(j)
          ind(i)=ind(j)
          i=j
          j=j+j
       else
          j=ir+1
       endif
       go to 20
    endif
    aout(i)=rra
    ind(i)=idum
    go to 10
end subroutine heapsort
!################################################################
!################################################################
!################################################################
subroutine mesh_info (ilevel,skip_loc,scale,dx,dx_loc,vol_loc,xc)
  use amr_commons
  implicit none
  integer::ilevel,ind,ix,iy,iz,nx_loc
  real(dp)::skip_loc(1:3),scale,dx,dx_loc,vol_loc
  real(dp),dimension(1:twotondim,1:ndim):: xc

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Cells center position relative to grid center position
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     xc(ind,1)=(dble(ix)-0.5D0)*dx
     xc(ind,2)=(dble(iy)-0.5D0)*dx
     xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

end subroutine mesh_info
!################################################################
!################################################################
!################################################################
subroutine get_icell_from_pos (fpos,ilevel_max,ind_grid,ind_cell,ilevel_out)
  use amr_commons
  implicit none
  real(dp)::fpos(1:3)
  integer ::ind_grid,ind_cell
  !-------------------------------------------------------------------
  ! This routnies find the index of the leaf cell for a given position
  ! fpos: positional info, from [0.0-1.0] (scaled by scale)
  ! ilevel_max: maximum level you want to search
  ! ind_cell: index of the cell
  ! ind_grid: index of the grid that contains the cell
  ! ilevel_out: level of this cell
  ! You can check whether this grid belongs to this cpu
  !     by asking cpu_map(father(ind_grid))
  !-------------------------------------------------------------------
  integer::ilevel_max
  integer::ilevel_out,nx_loc,i,ind,idim
  real(dp)::scale,dx,fpos2(1:3)
  real(dp)::skip_loc(1:3),x0(1:3)
  logical ::not_found

  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  skip_loc(1)=dble(icoarse_min)
  skip_loc(2)=dble(jcoarse_min)
  skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
 
  fpos2=fpos
  if (fpos2(1).gt.1d0) fpos2(1)=fpos2(1)-1d0
  if (fpos2(2).gt.1d0) fpos2(2)=fpos2(2)-1d0
  if (fpos2(3).gt.1d0) fpos2(3)=fpos2(3)-1d0
  if (fpos2(1).lt.0d0) fpos2(1)=fpos2(1)+1d0
  if (fpos2(2).lt.0d0) fpos2(2)=fpos2(2)+1d0
  if (fpos2(3).lt.0d0) fpos2(3)=fpos2(3)+1d0

  not_found=.true.
  ind_grid =1  ! this is level=1 grid
  ilevel_out=1
  do while (not_found)
     dx = 0.5D0**ilevel_out
     x0 = xg(ind_grid,1:ndim)-dx-skip_loc  !left corner of *this* grid [0-1], not cell (in ramses, grid is basically a struture containing 8 cells)

     ind  = 1 
     do idim = 1,ndim
        i = int( (fpos2(idim) - x0(idim))/dx)
        ind = ind + i*2**(idim-1)
     end do

     ind_cell = ind_grid+ncoarse+(ind-1)*ngridmax
!     write(*,'(2(I2,1x),2(I10,1x),3(f10.8,1x))') ilevel_out,ilevel_max,ind_grid,ind_cell,fpos2
     if(son(ind_cell)==0.or.ilevel_out==ilevel_max) return

     ind_grid=son(ind_cell)
     ilevel_out=ilevel_out+1
  end do

end subroutine get_icell_from_pos
!################################################################
!################################################################
!################################################################
subroutine SNII_yield (zp_star, ej_m, ej_Z, ej_chem)
  use amr_commons, ONLY:dp,nchem,chem_list
  use hydro_parameters, ONLY:ichem 
  implicit none
  real(dp)::zp_star,ej_m,ej_Z,ej_chem(1:nchem)
!-----------------------------------------------------------------
! Notice that even if the name is 'yield', 
! the return value is actually a metallicity fraction for simplicity
! These numbers are based on the outputs from Starburst99 
!                   (i.e. essentially Woosley & Weaver 95)
!-----------------------------------------------------------------
  real(dp),dimension(1:5)::log_SNII_m, log_Zgrid, log_SNII_Z
  real(dp),dimension(1:5)::log_SNII_H,log_SNII_He,log_SNII_C,log_SNII_N,log_SNII_O
  real(dp),dimension(1:5)::log_SNII_Mg,log_SNII_Si,log_SNII_S,log_SNII_Fe, dum1d
  real(dp)::log_Zstar,fz
  integer::nz_SN=5, izg, ich
  character(len=2)::element_name

  ! These are the numbers calculated from Starburst99 (Kroupa with 50Msun cut-off)
  ! (check library/make_stellar_winds.pro)
  log_SNII_m = (/-0.85591807,-0.93501857,-0.96138483,-1.0083450,-1.0544419/)
  log_Zgrid = (/-3.3979400,-2.3979400,-2.0969100,-1.6989700,-1.3010300/)
  log_SNII_Z=(/-0.99530662,-0.98262223,-0.96673581,-0.94018599,-0.93181853/)
  log_SNII_H=(/-0.28525316, -0.28988675, -0.29588988, -0.30967822, -0.31088065/)
  log_SNII_He=(/-0.41974152, -0.41688929, -0.41331511, -0.40330181, -0.40426739/)
  log_SNII_C=(/-1.9739731, -1.9726015, -1.9692265, -1.9626259, -1.9635311/)
  log_SNII_N=(/-3.7647616, -3.1325173, -2.8259748, -2.4260355, -2.1260417/)
  log_SNII_O=(/-1.1596291, -1.1491227, -1.1344617, -1.1213435, -1.1202089/)
  log_SNII_Mg=(/-2.4897201, -2.4368979, -2.4043654, -2.3706062, -2.3682933/)
  log_SNII_Si=(/-2.2157073, -2.1895132, -2.1518758, -2.0431845, -2.0444829/)
  log_SNII_S=(/-2.5492508, -2.5045016, -2.4482936, -2.2964020, -2.2988656/)
  log_SNII_Fe=(/-2.0502141, -2.0702598, -2.1074876, -2.2126987, -2.3480877/)

  ! search for the metallicity index
  log_Zstar = log10(zp_star)
  call binary_search(log_Zgrid, log_Zstar, nz_SN, izg )

  fz  = (log_Zgrid(izg+1) - log_Zstar )/( log_Zgrid(izg+1) - log_Zgrid(izg) )
  ! no extraploation
  if (fz  < 0.0) fz  = 0.0
  if (fz  > 1.0) fz  = 1.0


  ej_m = log_SNII_m(izg)*fz + log_SNII_m(izg+1)*(1d0-fz)
  ej_m = 10d0**ej_m

  ej_Z = log_SNII_Z(izg)*fz + log_SNII_Z(izg+1)*(1d0-fz)
  ej_Z = 10d0**ej_Z
 
  do ich=1,nchem
      element_name=chem_list(ich)
      select case (element_name)
         case ('H ')
            dum1d = log_SNII_H
         case ('He')
            dum1d = log_SNII_He
         case ('C ')
            dum1d = log_SNII_C
         case ('N ')
            dum1d = log_SNII_N
         case ('O ')
            dum1d = log_SNII_O
         case ('Mg')
            dum1d = log_SNII_Mg
         case ('Si')
            dum1d = log_SNII_Si
         case ('S ')
            dum1d = log_SNII_S
         case ('Fe')
            dum1d = log_SNII_Fe
         case default
            dum1d = 0d0
      end select

     ej_chem(ich) = dum1d(izg)*fz + dum1d(izg+1)*(1d0-fz)
     ej_chem(ich) = 10d0**ej_chem(ich)
  end do

end subroutine SNII_yield
!################################################################
!################################################################
!################################################################
!################################################################
subroutine binary_search(database,xtarget,ndata,i)
   use amr_commons,ONLY:dp
   implicit none 
   integer::i,j,k
   integer,intent(in)::ndata
   real(dp),intent(in)::database(1:ndata),xtarget

   i=1  
   j=ndata
   do   
     k=(i+j)/2
     if (xtarget<database(k)) then 
         j=k  
     else 
         i=k  
     end if
     if (i+1>=j) exit 
   end do

end subroutine binary_search
