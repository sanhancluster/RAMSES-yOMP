subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
#ifdef grackle
  use grackle_parameters
#endif
  use mpi_mod
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid,ii
  integer,dimension(1:nvector),save::ind_grid
  real(dp),dimension(1:ndust,1:4)::dM_dust_add

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  dM_dust_add=0d0
!$omp parallel do private(ngrid,ind_grid) reduction(+:dM_dust_add) schedule(static)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid,ngrid,ilevel,dM_dust_add)
  end do

  do ii=1,ndust
     !!sum on every cell
     dM_acc(ii)=dM_acc(ii)+dM_dust_add(ii,1)
     dM_spu(ii)=dM_spu(ii)+dM_dust_add(ii,2)
     dM_coa(ii)=dM_coa(ii)+dM_dust_add(ii,3)
     dM_sha(ii)=dM_sha(ii)+dM_dust_add(ii,4)
  enddo

  if((cooling.and..not.neq_chem).and.ilevel==levelmin.and.cosmo)then
#ifdef grackle
     if(use_grackle==0)then
        if(myid==1)write(*,*)'Computing new cooling table'
        call set_table(dble(aexp))
     endif
#else
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
#endif
  endif

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel,dM_dust_add)
  use amr_commons
  use hydro_commons
  use cooling_module
#ifdef grackle
  use grackle_parameters
#endif
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
#ifdef RT
  use rt_parameters, only: nGroups, iGroups
  use rt_hydro_commons
  use rt_cooling_module, only: rt_solve_cooling,iIR,rt_isIRtrap &
       ,rt_pressBoost,iIRtrapVar,kappaSc,kappaAbs,a_r,is_kIR_T,rt_vc
#endif
  use mpi_mod
  implicit none
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ii,ind,iskip,idim,nleaf,nx_loc,ich
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM,damp_factor,cooling_switch,t_blast
  real(dp)::polytropic_constant
  integer,dimension(1:nvector)::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector)::nH,T2,delta_T2,ekk,err,emag
#if defined(RT) || defined(grackle)
  real(kind=8),dimension(1:nvector)::T2_new
#endif
  real(kind=8),dimension(1:nvector)::T2min,Zsolar,boost
  real(kind=8),dimension(1:nvector,1:nchem)::Zchem
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
  real(dp)::mdust,mdustC,mdustSil,sigma2                      ! Dust (YD)
  real(dp)::year,d,T6,myT2,rhoDT0,rhoDT00,rhoDT               ! Dust (YD)
  real(dp)::error_rel,error_rel1,error_rel2,den0,den          ! Dust (YD)
  real(dp)::rhoG0,rhogZ0,rhoZ0                                ! Dust (YD)
  real(dp)::rhoGZ00,dtremain,dtloc,halfdtloc                  ! Dust (YD)
  real(dp),dimension(1:ndust)::dtloc_bin                      !!$dust_dev
  integer ::countmax=10000                                    ! Dust (YD)
  logical ::okdust                                            ! Dust (YD)
  logical ,dimension(1:ndust)::okdt_bin                       !!$dust_dev
  real(dp),dimension(1:ndust)::drhoD,rhoD,rhoD0               !!$dust_dev
  real(dp),dimension(1:ndust)::k1,k2,k3,k4                    !!$dust_dev
  real(dp),dimension(1:ndust)::rhoD0k1,rhoD0k2,rhoD0k3        !!$dust_dev
  real(dp),dimension(1:ndust)::rhoGZ0k1,rhoGZ0k2,rhoGZ0k3     !!$dust_dev
  real(kind=8),dimension(1:nvector,1:ndust)::fdust       ! Dust (YD) !! Dust-to-gas ratio
  real(kind=8),dimension(1:nvector)::sigma               ! Dust (YD) !! gas vel. dispersion
  real(dp),dimension(1:ndust)::t_des,t_acc,t0_des,t0_acc ! Dust (YD)
  real(dp),dimension(1:ndust)::oneovertdes,oneovertacc   ! Dust (YD)
  real(dp),dimension(1:ndust)::rhoD00                    ! Dust (YD)
  real(dp),dimension(1:ndust,1:4)::t0                         !!$dust_dev
  integer,dimension(1:ndust)::icount                          ! Dust (YD) !!$dust_dev
  real(dp)::t0_coa,t0_sha,t_coa,t_sha                         ! Dust (YD) !!$dust_dev
  real(dp)::oneovertsha,oneovertcoa                           ! Dust (YD) !!$dust_dev
  integer::ilow,ihigh                                         ! Dust (YD)
  real(dp),dimension(1:ndust,1:4)::dM_dust_add

#ifdef RT
  integer::ii,ig,iNp,il
  real(kind=8),dimension(1:nvector):: ekk_new
  logical,dimension(1:nvector)::cooling_on=.true.
  real(dp)::scale_Np,scale_Fp,work,Npc,Npnew, kIR, E_rad, TR
  real(dp),dimension(1:ndim)::Fpnew
  real(dp),dimension(nIons, 1:nvector):: xion
  real(dp),dimension(nGroups, 1:nvector):: Np, Np_boost=0d0, dNpdt=0d0
  real(dp),dimension(ndim, nGroups, 1:nvector):: Fp, Fp_boost=0d0, dFpdt
  real(dp),dimension(ndim, 1:nvector):: p_gas, u_gas
  real(kind=8)::f_trap, NIRtot, EIR_trapped, unit_tau, tau, Np2Ep
  real(kind=8)::aexp_loc, f_dust, xHII
  real(dp),dimension(nDim, nDim):: tEdd ! Eddington tensor
  real(dp),dimension(nDim):: flux
#endif
#if NENER>0
  integer::irad
#endif

  fdust=0d0
  !countmax=100000
  year=3600_dp*24_dp*365_dp
  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim


  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
#ifdef RT
  call rt_units(scale_Np, scale_Fp)
#endif

  ! Typical ISM density in H/cc
  nISM = n_star; nCOM=0d0
  if(cosmo)then
     nCOM = del_star*omega_b*rhoc*(h0/100.)**2/aexp**3*X/mH
  endif
  nISM = MAX(nCOM,nISM)

  ! Polytropic constant for Jeans length related polytropic EOS
  if(jeans_ncells>0)then
     polytropic_constant=2d0*(boxlen*jeans_ncells*0.5d0**dble(nlevelmax-nlevelsheld)*scale_l/aexp)**2/ &
          & (twopi)*6.67e-8*scale_d*(scale_t/scale_l)**2
  else
       polytropic_constant=1.
  endif

#ifdef RT
#if NGROUPS>0
  if(rt_isIRtrap) then
     ! For conversion from photon number density to photon energy density:
     Np2Ep = scale_Np * group_egy(iIR) * ev_to_erg                       &
           * rt_pressBoost / scale_d / scale_v**2
  endif
#endif
  aexp_loc=aexp
  ! Allow for high-z UV background in noncosmo sims:
  if(.not. cosmo .and. haardt_madau .and. aexp_ini .le. 1.)              &
       aexp_loc = aexp_ini
#endif

  ! Loop over cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,ngrid
        ind_cell(i)=iskip+ind_grid(i)
     end do

     ! Gather leaf cells
     nleaf=0
     do i=1,ngrid
        if(son(ind_cell(i))==0)then
           nleaf=nleaf+1
           ind_leaf(nleaf)=ind_cell(i)
        end if
     end do
     if(nleaf.eq.0)cycle

     ! Compute rho
     do i=1,nleaf
        nH(i)=MAX(uold(ind_leaf(i),1),smallr)
     end do

     ! Compute metallicity in solar units
     if(metal)then
        if(dust .and. metal_gasonly)then
           do i=1,nleaf
              Zsolar(i)=uold(ind_leaf(i),imetal)/nH(i)/0.02 !! total met (gas+dust)
              do ich=1,nchem
                 Zchem(i,ich)=uold(ind_leaf(i),ichem+ich-1)/nH(i)/0.02
                 if(Zchem(i,ich)<0)write(*,'(A,3I7,es13.5)')'WARNING:',ind_leaf(i),myid,ich,Zchem(i,ich)
              enddo
              if(dust_chem)then
                 mdustC=0.0d0;mdustSil=0.0d0
                 ilow=idust;ihigh=ilow+dndsize
                 mdustC  =SUM(uold(ind_leaf(i),ilow:ihigh))
                 ilow=ihigh+1;ihigh=ilow+dndsize
                 mdustSil=SUM(uold(ind_leaf(i),ilow:ihigh))/SioverSil

!!$                 ndchemtype=ndust/2 ! not bullet proof (only works with ndust=2 or 4): more dust types needs some revision
!!$                 do ii=1,ndchemtype
!!$                    mdustC  =mdustC  +uold(ind_leaf(i),idust-1+ii)
!!$                 enddo
!!$                 do ii=ndchemtype+1,ndust
!!$                    mdustSil=mdustSil+uold(ind_leaf(i),idust-1+ii)/SioverSil
!!$                 enddo
                 mdust=mdustC+mdustSil
                 Zsolar(i)=Zsolar(i)-mdust/nH(i)/0.02 !! gas met
!!$                 if(dustdebug)write(*,'(A,5es13.5)')'Zsolar',Zsolar(i),Zsolar(i)+mdust/nH(i)/0.02,mdust/nH(i)/0.02,mdustC/nH(i)/0.02,mdustSil/nH(i)/0.02

                 Zchem(i,ichC)=Zchem(i,ichC)-mdustC/nH(i)/0.02
                 Zchem(i,ichMg)=Zchem(i,ichMg)-mdustSil*MgoverSil/nH(i)/0.02
                 Zchem(i,ichFe)=Zchem(i,ichFe)-mdustSil*FeoverSil/nH(i)/0.02
                 Zchem(i,ichSi)=Zchem(i,ichSi)-mdustSil*SioverSil/nH(i)/0.02
                 Zchem(i,ichO)=Zchem(i,ichO)-mdustSil* OoverSil/nH(i)/0.02

                 do ich=1,nchem
                    if(dustdebug)then
!!$                    if(TRIM(chem_list(ich))=='C' )then
!!$                       write(*,'(A,3es13.5)')'ZC (tot,gas,dust)=',Zchem(i,ich),Zchem(i,ich)-mdustC/nH(i)/0.02,mdustC/nH(i)/0.02!YD DEBUG
!!$                    endif
                    endif
                    if(Zchem(i,ich).lt.0.0d0)then
                       write(*,*)'Problem with Zchem<0 in cooling_fine'
                       if(TRIM(chem_list(ich))=='C' )write(*,'(A,I2,A,3es13.5,3I7)')'Zchem(',ich,'), C  gas,dust,tot',Zchem(i,ich),mdustC/nH(i)/0.02,Zchem(i,ich)+mdustC/nH(i)/0.02,i,ind_leaf(i),myid
                       if(TRIM(chem_list(ich))=='Mg')write(*,'(A,I2,A,3es13.5,3I7)')'Zchem(',ich,'), Mg gas,dust,tot',Zchem(i,ich),mdustSil*MgoverSil/nH(i)/0.02,Zchem(i,ich)+mdustSil*MgoverSil/nH(i)/0.02,i,ind_leaf(i),myid
                       if(TRIM(chem_list(ich))=='Fe')write(*,'(A,I2,A,3es13.5,3I7)')'Zchem(',ich,'), Fe gas,dust,tot',Zchem(i,ich),mdustSil*FeoverSil/nH(i)/0.02,Zchem(i,ich)+mdustSil*FeoverSil/nH(i)/0.02,i,ind_leaf(i),myid
                       if(TRIM(chem_list(ich))=='Si')write(*,'(A,I2,A,3es13.5,3I7)')'Zchem(',ich,'), Si gas,dust,tot',Zchem(i,ich),mdustSil*SioverSil/nH(i)/0.02,Zchem(i,ich)+mdustSil*SioverSil/nH(i)/0.02,i,ind_leaf(i),myid
                       if(TRIM(chem_list(ich))=='O' )write(*,'(A,I2,A,3es13.5,3I7)')'Zchem(',ich,'), O  gas,dust,tot',Zchem(i,ich),mdustSil*OoverSil /nH(i)/0.02,Zchem(i,ich)+mdustSil*OoverSil /nH(i)/0.02,i,ind_leaf(i),myid
                       write(*,*) dM_dust_add
!!$                       stop
                    endif
                 enddo
              else
                 do ii=1,ndust
                    Zsolar(i)=Zsolar(i)-uold(ind_leaf(i),idust-1+ii)/nH(i)/0.02 !! gas met
                 enddo
              endif
           end do
        else
           do i=1,nleaf
              Zsolar(i)=uold(ind_leaf(i),imetal)/nH(i)/0.02
              do ich=1,nchem
                 Zchem(i,ich)=uold(ind_leaf(i),ichem+ich-1)/nH(i)/0.02
              enddo
           end do
        endif
        if(dust)then
           do i=1,nleaf
              do ii=1,ndust!!$dust_dev
                 fdust(i,ii)=uold(ind_leaf(i),idust-1+ii)/nH(i) ! Dust-to-gas ratio
                 if(fdust(i,ii)<0.0 .and. log_dc)then
                    write(*,*)'for bin :', ii, 'fdust<0 (beg cooling)',ind_leaf(i),uold(ind_leaf(i),idust-1+ii),nH(i)
                    fdust(i,ii)=0.0d0
                 endif
              end do
           end do
        endif
     else
        do i=1,nleaf
           Zsolar(i)=z_ave
        end do
     endif

#ifdef RT
     ! Floor density (prone to go negative with strong rad. pressure):
     do i=1,nleaf
        uold(ind_leaf(i),1) = max(uold(ind_leaf(i),1),smallr)
     end do
     ! Initialise gas momentum and velocity for photon momentum abs.:
     do i=1,nleaf
        p_gas(:,i) = uold(ind_leaf(i),2:ndim+1) * scale_d * scale_v
        u_gas(:,i) = uold(ind_leaf(i),2:ndim+1) &
                     /uold(ind_leaf(i),1) * scale_v
     end do

#if NGROUPS>0
     if(rt_isIRtrap) then  ! Gather also trapped photons for solve_cooling
        iNp=iGroups(iIR)
        do i=1,nleaf
           il=ind_leaf(i)
           rtuold(il,iNp) = rtuold(il,iNp) + uold(il,iIRtrapVar)/Np2Ep
           if(rt_smooth) &
                rtunew(il,iNp)= rtunew(il,iNp) + uold(il,iIRtrapVar)/Np2Ep
        end do
     endif

     if(rt_vc) then      ! Add/remove radiation work on gas. Eq A6 in RT15
        iNp=iGroups(iIR)
        do i=1,nleaf
           il=ind_leaf(i)
           NIRtot = rtuold(il,iNp)
           kIR  = kappaSc(iIR)
           if(is_kIR_T) then                        ! kIR depends on T_rad
              ! For rad. temperature,  weigh the energy in each group by
              ! its opacity over IR opacity (derived from IR temperature)
              E_rad = group_egy(iIR) * ev_to_erg * NIRtot * scale_Np
              TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)     ! IR temp.
              kIR = kappaAbs(iIR) * (TR/10d0)**2
              do ig=1,nGroups
                 if(i .ne. iIR)                                          &
                      E_rad = E_rad + kappaAbs(ig) / kIR                 &
                            * max(rtuold(il,iGroups(ig)),smallNp)        &
                            * ev_to_erg * scale_Np
              end do
              TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)   ! Rad. temp.
              ! Set the IR opacity according to the rad. temperature:
              kIR  = kappaSc(iIR)  * (TR/10d0)**2 * exp(-TR/1d3)
           endif
           kIR = kIR*scale_d*scale_l           !  Convert to code units
           flux = rtuold(il,iNp+1:iNp+ndim)
           xHII = uold(il,iIons-1+ixHII)/uold(il,1)
           f_dust = (1.-xHII)                     ! No dust in ionised gas
           work = scale_v/c_cgs * kIR * sum(uold(il,2:ndim+1)*flux) &
                * Zsolar(i) * f_dust * dtnew(ilevel) !               Eq A6
           ! Add work to gas energy
           uold(il,ndim+2) = uold(il,ndim+2) &
                + work * group_egy(iIR) &
                * ev_to_erg / scale_d / scale_v**2 / scale_l**3

           rtuold(il,iNp) = rtuold(il,iNp) - work !Remove from rad density
           rtuold(il,iNp) = max(rtuold(il,iNp),smallnp)
           call reduce_flux(rtuold(il,iNp+1:iNp+ndim),rtuold(il,iNp)*rt_c)
        enddo
     endif
#endif
#endif

     ! Compute thermal pressure
     do i=1,nleaf
        T2(i)=uold(ind_leaf(i),ndim+2)
     end do
     do i=1,nleaf
        ekk(i)=0.0d0
     end do
     do idim=1,ndim
        do i=1,nleaf
           ekk(i)=ekk(i)+0.5*uold(ind_leaf(i),idim+1)**2/nH(i)
        end do
     end do
     do i=1,nleaf
        err(i)=0.0d0
     end do
#if NENER>0
     do irad=0,nener-1
        do i=1,nleaf
           err(i)=err(i)+uold(ind_leaf(i),inener+irad)
        end do
     end do
#endif
     do i=1,nleaf
        emag(i)=0.0d0
     end do
#ifdef SOLVERmhd
     do idim=1,ndim
        do i=1,nleaf
           emag(i)=emag(i)+0.125d0*(uold(ind_leaf(i),idim+ndim+2)+uold(ind_leaf(i),idim+nvar))**2
        end do
     end do
#endif
     do i=1,nleaf
        T2(i)=(gamma-1.0)*(T2(i)-ekk(i)-err(i)-emag(i))
     end do

     ! Compute T2=T/mu in Kelvin
     do i=1,nleaf
        T2(i)=T2(i)/nH(i)*scale_T2
     end do

     ! Compute nH in H/cc
     do i=1,nleaf
        nH(i)=nH(i)*scale_nH
     end do

!!$     if(dust_chem)then
!!$        do i=1,nleaf
!!$           if(myid==25.and.ind_leaf(i)==504376)write(*,'(A)')'               Zdust/0.02   nH(H/cc)     T(K)         myid Ztot/0.02   Zdustkey/0.02'
!!$           do ii=1,ndchemtype
!!$              if(myid==25.and.ind_leaf(i)==504376)write(*,'(A,I2,3es13.5,I4)')'(0) C  dust',ii,uold(ind_leaf(i),idust-1+ii)/(nH(i)/scale_nH)/0.02,nh(i),T2(i),myid
!!$           enddo
!!$           do ii=ndchemtype+1,ndust
!!$              if(myid==25.and.ind_leaf(i)==504376)write(*,'(A,I2,3es13.5,I4)')'(0) Si dust',ii,uold(ind_leaf(i),idust-1+ii)/SioverSil/(nH(i)/scale_nH)/0.02,nh(i),T2(i),myid
!!$           enddo
!!$        enddo
!!$     endif

     ! Compute radiation boost factor
     if(self_shielding)then
        do i=1,nleaf
           boost(i)=MAX(EXP(-nH(i)/0.01),1.0D-20)
        end do
#ifdef ATON
     else if (aton) then
        do i=1,nleaf
           boost(i)=MAX(Erad(ind_leaf(i))/J0simple(aexp), &
                &                   J0min/J0simple(aexp) )
        end do
#endif
     else
        do i=1,nleaf
           boost(i)=1.0
        end do
     endif

     !==========================================
     ! Compute temperature from polytrope EOS
     !==========================================
     if(jeans_ncells>0)then
        do i=1,nleaf
           T2min(i) = nH(i)*polytropic_constant*scale_T2
        end do
     else
        do i=1,nleaf
           T2min(i) = T2_star*(nH(i)/nISM)**(g_star-1.0)
        end do
     endif
     !==========================================
     ! You can put your own polytrope EOS here
     !==========================================
     if(cooling)then
        ! Compute thermal temperature by subtracting polytrope
        do i=1,nleaf
           T2(i) = min(max(T2(i)-T2min(i),T2_min_fix),T2max)
        end do
     endif

     sigma(1:nvector)=0.0d0
     if(sticking_coef=='subgrid')then
        do i=1,nleaf
           if(nH(i).ge.0.1d0)then
              call cmp_sigma_turb(ind_leaf(i),sigma2,ilevel)
              sigma(i)=sqrt(sigma2)*scale_v
           endif
           if(nvar>idust+ndust-1)uold(ind_leaf(i),nvar)=uold(ind_leaf(i),1)*sigma(i)/scale_v
        enddo
     endif


     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t

#ifdef RT
     if(neq_chem) then
        ! Get the ionization fractions
        do ii=0,nIons-1
           do i=1,nleaf
              xion(1+ii,i) = uold(ind_leaf(i),iIons+ii)/uold(ind_leaf(i),1)
           end do
        end do

        ! Get photon densities and flux magnitudes
        do ig=1,nGroups
           iNp=iGroups(ig)
           do i=1,nleaf
              il=ind_leaf(i)
              Np(ig,i)        = scale_Np * rtuold(il,iNp)
              Fp(1:ndim, ig, i) = scale_Fp * rtuold(il,iNp+1:iNp+ndim)
           enddo
           if(rt_smooth) then                           ! Smooth RT update
              do i=1,nleaf !Calc addition per sec to Np, Fp for current dt
                 il=ind_leaf(i)
                 Npnew = scale_Np * rtunew(il,iNp)
                 Fpnew = scale_Fp * rtunew(il,iNp+1:iNp+ndim)
                 dNpdt(ig,i)   = (Npnew - Np(ig,i)) / dtcool
                 dFpdt(:,ig,i) = (Fpnew - Fp(:,ig,i)) / dtcool
              end do
           end if
        end do

        if(cooling .and. delayed_cooling) then
           cooling_on(1:nleaf)=.true.
           do i=1,nleaf
              if(uold(ind_leaf(i),idelay)/uold(ind_leaf(i),1) .gt. 1d-3) &
                   cooling_on(i)=.false.
           end do
        end if
        if(isothermal)cooling_on(1:nleaf)=.false.
     endif

     if(rt_vc) then ! Do the Lorentz boost. Eqs A4 and A5. in RT15
        do i=1,nleaf
           do ig=1,nGroups
              Npc=Np(ig,i)*rt_c_cgs
              call cmp_Eddington_tensor(Npc,Fp(:,ig,i),tEdd)
              Np_boost(ig,i) = - 2d0/c_cgs/rt_c_cgs * sum(u_gas(:,i)*Fp(:,ig,i))
              do idim=1,ndim
                 Fp_boost(idim,ig,i) =  &
                      -u_gas(idim,i)*Np(ig,i) * rt_c_cgs/c_cgs &
                      -sum(u_gas(:,i)*tEdd(idim,:))*Np(ig,i)*rt_c_cgs/c_cgs
              end do
           end do
           Np(:,i)   = Np(:,i) + Np_boost(:,i)
           Fp(:,:,i) = Fp(:,:,i) + Fp_boost(:,:,i)
        end do
     endif
#endif

     ! grackle tabular cooling
#ifdef grackle
     if(use_grackle==1)then
        gr_rank = 3
        do i = 1, gr_rank
           gr_dimension(i) = 1
           gr_start(i) = 0
           gr_end(i) = 0
        enddo
        gr_dimension(1) = nvector
        gr_end(1) = nleaf - 1

        if(cosmo)then
           my_grackle_units%a_value = aexp
           my_grackle_units%density_units = scale_d
           my_grackle_units%length_units = scale_l
           my_grackle_units%time_units = scale_t
           my_grackle_units%velocity_units = scale_v
        endif

        do i = 1, nleaf
           gr_density(i) = uold(ind_leaf(i),1)
           if(metal)then
              gr_metal_density(i) = uold(ind_leaf(i),imetal)
           else
              gr_metal_density(i) = uold(ind_leaf(i),1)*0.02*z_ave
           endif
           gr_energy(i) = T2(i)/(scale_T2*(gamma-1.0))
           gr_HI_density(i) = X*gr_density(i)
           gr_HeI_density(i) = (1.0-X)*gr_density(i)
           gr_DI_density(i) = 2.0*3.4e-5*gr_density(i)
        enddo
        ! Update grid properties
        my_grackle_fields%grid_rank = gr_rank
        my_grackle_fields%grid_dx = dx_loc

        iresult = solve_chemistry(my_grackle_units, my_grackle_fields, dtnew(ilevel))
        if(iresult.eq.0)then
            write(*,*) 'Grackle: error in solve_chemistry'
#ifndef WITHOUTMPI
            call MPI_ABORT(MPI_COMM_WORLD,1,info)
#else
            stop
#endif
        endif

        do i = 1, nleaf
           T2_new(i) = gr_energy(i)*scale_T2*(gamma-1.0)
        end do
        delta_T2(1:nleaf) = T2_new(1:nleaf) - T2(1:nleaf)
     else
        ! Compute net cooling at constant nH
        if(cooling.and..not.neq_chem)then
#if NCHEM>0
           call solve_cooling(nH,T2,Zsolar,Zchem,fdust,sigma,boost,dtcool,delta_T2,nleaf,ilevel)
#else
           call solve_cooling(nH,T2,Zsolar,fdust,sigma,boost,dtcool,delta_T2,nleaf,ilevel)
#endif
        endif
     endif
#else
     ! Compute net cooling at constant nH
     if(cooling.and..not.neq_chem)then
#if NCHEM>0
           call solve_cooling(nH,T2,Zsolar,Zchem,fdust,sigma,boost,dtcool,delta_T2,nleaf,ilevel,dM_dust_add)
#else
           call solve_cooling(nH,T2,Zsolar,fdust,sigma,boost,dtcool,delta_T2,nleaf,ilevel)
#endif
     endif
#endif
#ifdef RT
     if(neq_chem) then
        T2_new(1:nleaf) = T2(1:nleaf)
        call rt_solve_cooling(T2_new, xion, Np, Fp, p_gas, dNpdt, dFpdt  &
                         , nH, cooling_on, Zsolar, dtcool, aexp_loc,nleaf)
        delta_T2(1:nleaf) = T2_new(1:nleaf) - T2(1:nleaf)
     endif
#endif

#ifdef RT
     if(.not. static) then
        ! Update gas momentum and kinetic energy:
        do i=1,nleaf
           uold(ind_leaf(i),2:1+ndim) = p_gas(:,i) /scale_d /scale_v
        end do
        ! Energy update ==================================================
        ! Calculate NEW pressure from updated momentum
        ekk_new(1:nleaf) = 0d0
        do i=1,nleaf
           do idim=1,ndim
              ekk_new(i) = ekk_new(i) &
                   + 0.5*uold(ind_leaf(i),idim+1)**2 / max(uold(ind_leaf(i),1), smallr)
           end do
        end do
        do i=1,nleaf
           ! Update the pressure variable with the new kinetic energy:
           uold(ind_leaf(i),ndim+2) = uold(ind_leaf(i),ndim+2)           &
                                    - ekk(i) + ekk_new(i)
        end do
        do i=1,nleaf
           ekk(i)=ekk_new(i)
        end do

#if NGROUPS>0
        if(rt_vc) then ! Photon work: subtract from the IR ONLY radiation
           do i=1,nleaf
              Np(iIR,i) = Np(iIR,i) + (ekk(i) - ekk_new(i))              &
                   /scale_d/scale_v**2 / group_egy(iIR) / ev_to_erg
           end do
        endif
#endif
        ! End energy update ==============================================
     endif ! if(.not. static)
#endif

     ! Compute rho
     do i=1,nleaf
        nH(i) = nH(i)/scale_nH
     end do

     ! Deal with cooling
     if(cooling.or.neq_chem)then
        ! Compute net energy sink
        do i=1,nleaf
           delta_T2(i) = delta_T2(i)*nH(i)/scale_T2/(gamma-1.0)
        end do
        ! Compute initial fluid internal energy
        do i=1,nleaf
           T2(i) = T2(i)*nH(i)/scale_T2/(gamma-1.0)
        end do
        ! Turn off cooling in blast wave regions
        if(delayed_cooling)then
           do i=1,nleaf
              cooling_switch = uold(ind_leaf(i),idelay)/max(uold(ind_leaf(i),1),smallr)
              if(cooling_switch > 1d-3)then
                 delta_T2(i) = MAX(delta_T2(i),real(0,kind=dp))
              endif
           end do
        endif
     endif

     ! Compute polytrope internal energy
     do i=1,nleaf
        T2min(i) = T2min(i)*nH(i)/scale_T2/(gamma-1.0)
     end do

     ! Update fluid internal energy
     if(cooling.or.neq_chem)then
        do i=1,nleaf
           T2(i) = T2(i) + delta_T2(i)
        end do
     endif

     ! Update total fluid energy
     if(isothermal)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = T2min(i) + ekk(i) + err(i) + emag(i)
        end do
     else if(cooling .or. neq_chem)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = T2(i) + T2min(i) + ekk(i) + err(i) + emag(i)
        end do
     endif

     ! Update delayed cooling switch
     if(delayed_cooling)then
        t_blast=t_diss*1d6*(365.*24.*3600.)
        damp_factor=exp(-dtcool/t_blast)
        do i=1,nleaf
           uold(ind_leaf(i),idelay)=max(uold(ind_leaf(i),idelay)*damp_factor,0d0)
        end do
     endif

     if(dust)then

        if(.not.dust_dest_within_cool)then !! dust_cooling=false
        do i=1,nleaf !!on all cell, compute temp
           T2(i)=(uold(ind_leaf(i),ndim+2) - ekk(i) - err(i) - emag(i))
           T2(i)=(gamma-1d0)*T2(i)
           T2(i)=T2(i)/nH(i)*scale_T2
        end do
        do ii=1,ndust!!$dust_dev
           t0_des(ii)=t_sputter_ref*year*asize(ii)/0.1d0 !!t0 for sputtering, depend on size
           t0_acc(ii)=t_growth_ref *year*asize(ii)/0.005d0  !Aoyama, 2016
        enddo
        !!$work only for 2 dust bins, but didn't want to hard-code
        t0_sha = t_sha_ref*year*(asize(2)/0.1) !!for v=10km/s and density_grain=3g.cm-3
        t0_coa = t_coa_ref*year*(0.1d0/0.1d0)*(asize(1)/0.005)  !! velocity dispersion=0.1km/s
        do i=1,nleaf

           rhoG0 =uold(ind_leaf(i),1) !!density
           rhoZ0 =uold(ind_leaf(i),imetal) !!metal dens
           d=nH(i)*scale_nH !!gas density in real units
           myT2=MAX(T2(i),T2_star) ! -> replace T2(i) by myT2
           T6=myT2/1d6 !!temp/1e6
           select case (thermal_sputtering) !!Novak or Tsai (did in Granato et al.)
           case('novak')
              do ii=1,ndust
                t_des(ii)=t0_des(ii)/d*(1d0+1d0/T6**3) ! Draine & Salpeter (1979) (see also Novak et al, 2012)
              end do
           case('tsai')
              do ii=1,ndust
                t_des(ii)=1.65d0*t0_des(ii)/d*(1d0+(2d0/T6)**2.5) ! Tsai & Matthews (1998)
              end do
           case default
              do ii=1,ndust
                t_des(ii)=1.65d0*t0_des(ii)/d*(1d0+(2d0/T6)**2.5) ! Tsai & Matthews (1998)
              end do
           end select
           select case (sticking_coef) !!For accretion ->chaabouni
           case('constant')
              do ii=1,ndust
                 t_acc(ii)=t0_acc(ii)* 1d3/d*sqrt(50d0/myT2)*(0.0134*rhoG0/rhoZ0)*(1d0+10.0d0*T6) ! Bekki (2015)
              enddo
           case('step')
              if(myT2.gt.1d3)then
                 do ii=1,ndust
                    t_acc(ii)=1d12*year
                 enddo
              else
                 do ii=1,ndust
                    t_acc(ii)=t0_acc(ii)* 1d3/d*sqrt(50d0/myT2)*(0.0134*rhoG0/rhoZ0)
                 enddo
              endif
           case('chaabouni')
              do ii=1,ndust
                 t_acc(ii)=t0_acc(ii)* 1d3/d*sqrt(50d0/myT2)*(0.0134*rhoG0/rhoZ0) / (0.95d0*(1d0+2.22d0*myT2/56d0)/(1d0+myT2/56d0)**2.22)
              enddo
           case('leitch')
              do ii=1,ndust
                 t_acc(ii)=t0_acc(ii) * 1d3/d*sqrt(50d0/myT2)*(0.0134*rhoG0/rhoZ0) / (1.9d-2*myT2*(1.7d-2+0.4d0)*exp(-7d-3*myT2))
              enddo
           case('bourlot') ! Le Bourlot+ 2012
              do ii=1,ndust
                 t_acc(ii)=t0_acc(ii) * 1d3/d*sqrt(50d0/myT2)*(0.0134*rhoG0/rhoZ0) * (1d0+1d-4*myT2**1.5)
              enddo
           case default ! chaabouni case
              do ii=1,ndust
                 t_acc(ii)=t0_acc(ii) * 1d3/d*sqrt(50d0/myT2)*(0.0134*rhoG0/rhoZ0) / (0.95d0*(1d0+2.22d0*myT2/56d0)/(1d0+myT2/56d0)**2.22)
              enddo
           end select

           !Dust Coagulation
           if(d>1d2) then !!between 1d2 and 1d3
              t_coa=t0_coa*(0.01/(uold(ind_leaf(i),idust)/nH(i)))
           else                                         !!take only dtg ratio of small bins
              t_coa=1d12*year
           endif

           !Dust Shattering
           if(d<1.0d0) then
              t_sha=(t0_sha/d)*(0.01/(uold(ind_leaf(i),idust+1)/nH(i))) !!only dtg ratio of large bins
           elseif (d<1d3) then
              t_sha=(t0_sha/(d**(1.0d0/3.0d0)))*(0.01/(uold(ind_leaf(i),idust+1)/nH(i)))
           else
              t_sha=1d12*year
           endif

           if(rhoZ0.gt.0.0d0)then !!test if rhoZ not null in the cell
              okdust=.false.
           else
              okdust=.true. !!else no rk4
           endif

           rhoDT00=0.0d0
           icount=0
           dtremain=dtcool !!timestep for rk4
           do ii=1,ndust!!$dust_dev
              oneovertdes(ii)=1d0/t_des(ii) !!1/t_sputt
              oneovertacc(ii)=1d0/t_acc(ii) !!1/t_acc
              rhoD00(ii)=uold(ind_leaf(i),idust-1+ii) !!initial dust density
              rhoDT00=rhoDT00+rhoD00(ii)
           enddo
           rhoGZ00=rhoZ0-rhoDT00
           oneovertsha=1d0/t_sha
           oneovertcoa=1d0/t_coa
           t0=1d12*year
           do ii=1,ndust
              if(dust_accretion)   t0(ii,1)=0.1d0*t_acc(ii)
              if(dust_sputtering)  t0(ii,2)=0.1d0*t_des(ii)
              if(dust_shattering)  t0(ii,3)=0.1d0*t_sha
              if(dust_coagulation) t0(ii,4)=0.1d0*t_coa
           enddo

           do while (okdust .eqv. .false.)
              rhoDT0=0.0d0 ! Total dust mass over all dust bins
              do ii=1,ndust
                 rhoDT0=rhoDT0+uold(ind_leaf(i),idust-1+ii)!!total dust dentity over bins
                 rhoD0(ii)=uold(ind_leaf(i),idust-1+ii)    !!current dust density in bins
                 if(icount(ii)==0) dtloc_bin(ii)=MINVAL(t0(ii,:)) !! first timestep -> 10% fastest processus
              enddo
              rhoGZ0=rhoZ0-rhoDT0
              dtloc=MINVAL(dtloc_bin)
              dtloc=MIN(dtloc,dtremain) !!min timestep
              halfdtloc=0.5d0*dtloc !!timestep/2 for rk4

              if(zdmax.gt.0.0d0)then !!to speed-up, test dust-to-metal ratio not too high (zdmax=0.9999d0, amr_parameters)
                 if(rhoDT0/rhoZ0.gt.zdmax)then
                    write(*,*) "Entering fix dust density :"
                    okdust=.true. !!if sup, stop rk4
                    do ii=1,ndust
                       uold(ind_leaf(i),idust-1+ii)=zdmax*rhoZ0*(uold(ind_leaf(i),idust-1+ii)/rhoDT0) !!dust_dens_bin=max_ratio*metal_dens*dust_bin/total_dust
                       write(*,*) uold(ind_leaf(i),idust-1+ii), ii
                    enddo
                    write(*,*)uold(ind_leaf(i),idust)+uold(ind_leaf(i),idust+1)
                 endif
              endif

              !!$ RK4 for all processus
              k1=0
              do ii=1,ndust
                 if(dust_accretion)  k1(ii) = k1(ii) + (rhoGZ0/rhoZ0)*rhoD0(ii)   *oneovertacc(ii)
                 if(dust_sputtering) k1(ii) = k1(ii) -    rhoD0(ii)       *oneovertdes(ii)
                 if(dust_shattering) then
                    if(ii==1)        k1(ii) = k1(ii) +    rhoD0(2)        *oneovertsha
                    if(ii==2)        k1(ii) = k1(ii) -    rhoD0(2)        *oneovertsha
                 endif
                 if(dust_coagulation) then
                    if(ii==1)        k1(ii) = k1(ii) -    rhoD0(1)        *oneovertcoa
                    if(ii==2)        k1(ii) = k1(ii) +    rhoD0(1)        *oneovertcoa
                 endif
                 rhoD0k1(ii)=rhoD0(ii)+halfdtloc*k1(ii)
                 rhoGZ0k1(ii)=rhoGZ0-halfdtloc*k1(ii)
              enddo

              k2=0
              do ii=1,ndust
                 if(dust_accretion)  k2(ii) = k2(ii) + (rhoGZ0k1(ii)/rhoZ0)*rhoD0k1(ii) *oneovertacc(ii)
                 if(dust_sputtering) k2(ii) = k2(ii) -    rhoD0k1(ii)       *oneovertdes(ii)
                 if(dust_shattering) then
                    if(ii==1)        k2(ii) = k2(ii) +    rhoD0k1(2)        *oneovertsha
                    if(ii==2)        k2(ii) = k2(ii) -    rhoD0k1(2)        *oneovertsha
                 endif
                 if(dust_coagulation) then
                    if(ii==1)        k2(ii) = k2(ii) -    rhoD0k1(1)        *oneovertcoa
                    if(ii==2)        k2(ii) = k2(ii) +    rhoD0k1(1)        *oneovertcoa
                 endif
                 rhoD0k2(ii)=rhoD0(ii)+halfdtloc*k2(ii)
                 rhoGZ0k2(ii)=rhoGZ0-halfdtloc*k2(ii)
              enddo

              k3=0
              do ii=1,ndust
                 if(dust_accretion)  k3(ii) = k3(ii) + (rhoGZ0k2(ii)/rhoZ0)*rhoD0k2(ii) *oneovertacc(ii)
                 if(dust_sputtering) k3(ii) = k3(ii) -    rhoD0k2(ii)       *oneovertdes(ii)
                 if(dust_shattering) then
                    if(ii==1)        k3(ii) = k3(ii) +    rhoD0k2(2)        *oneovertsha
                    if(ii==2)        k3(ii) = k3(ii) -    rhoD0k2(2)        *oneovertsha
                 endif
                 if(dust_coagulation) then
                    if(ii==1)        k3(ii) = k3(ii) -    rhoD0k2(1)        *oneovertcoa
                    if(ii==2)        k3(ii) = k3(ii) +    rhoD0k2(1)        *oneovertcoa
                 endif
                 rhoD0k3(ii)=rhoD0(ii)+dtloc*k3(ii)
                 rhoGZ0k3(ii)=rhoGZ0-dtloc*k3(ii)
              enddo

              k4=0
              do ii=1,ndust
                 if(dust_accretion)  k4(ii) = k4(ii) + (rhoGZ0k3(ii)/rhoZ0)*rhoD0k3(ii) *oneovertacc(ii)
                 if(dust_sputtering) k4(ii) = k4(ii) -    rhoD0k3(ii)       *oneovertdes(ii)
                 if(dust_shattering) then
                    if(ii==1)        k4(ii) = k4(ii) +    rhoD0k3(2)        *oneovertsha
                    if(ii==2)        k4(ii) = k4(ii) -    rhoD0k3(2)        *oneovertsha
                 endif
                 if(dust_coagulation) then
                    if(ii==1)        k4(ii) = k4(ii) -    rhoD0k3(1)        *oneovertcoa
                    if(ii==2)        k4(ii) = k4(ii) +    rhoD0k3(1)        *oneovertcoa
                 endif
                 drhoD(ii)=dtloc/6d0*(k1(ii)+2d0*k2(ii)+2d0*k3(ii)+k4(ii))
                 rhoD(ii)=rhoD0(ii)+drhoD(ii) !!new dust density (eqn 24 Trebitsch, 2020)
              enddo

              okdt_bin=.false.
              do ii=1,ndust
                 if(rhoD0(ii)>1d-20)then
                    error_rel1=abs(drhoD(ii))/MIN(rhoD0(ii),rhoD(ii))
                    den0=(1d0-rhoD0(ii)/rhoZ0)*rhoD0(ii)
                    den =(1d0-rhoD(ii) /rhoZ0)*rhoD(ii)
                    if(MIN(den0,den)<1d-10) then
                       error_rel=error_rel1
                    else
                       error_rel2=abs(drhoD(ii))/MIN(den0,den)
                       error_rel=MAX(error_rel1,error_rel2)
                    endif
                 else
                    okdust=.true. !!else stop rk4
                 endif

                 if(.not.okdust)then !!if still in the do while (rk4)
                    if(error_rel.le.errmax.and.error_rel.ge.0.0d0) then
                       uold(ind_leaf(i),idust-1+ii)=rhoD(ii) !!save dust density in bin
                       okdt_bin(ii)=.true.
                       if(error_rel.le.0.5d0*errmax)dtloc_bin(ii)=dtloc*2.0d0 !!new timestep
                    endif
                    if(error_rel.gt.errmax.or.error_rel.lt.0.0d0)dtloc_bin(ii)=dtloc/2.0d0
                    icount(ii)=icount(ii)+1
                 endif !! still in rk4
                 if(icount(ii)>countmax)then !!error, stop do while
                    write(*,*)'stopping in dust processing icount>',countmax
                    write(*,*) 'bin :', ii
                    write(*,*)'rhoG      rhoZ     rhoGZ      rhoD      Temperature'
                    write(*,'(5e10.2)')rhoG0*scale_nH,rhoZ0*scale_nH,rhoGZ0 *scale_nH,rhoD0(ii) *scale_nH,T6*1d6
                    write(*,'(5e10.2)')rhoG0*scale_nH,rhoZ0*scale_nH,rhoGZ00*scale_nH,rhoD00(ii)*scale_nH,T6*1d6
                    write(*,*)'dtloc     dtstep     dtremain'
                    write(*,'(3e10.3,A)')dtloc_bin(ii)/3.15d13,dtcool/3.15d13,dtremain/3.15d13,' Myr'
                    write(*,*)'t_acc     t_des'
                    write(*,'(2e10.3,A)')t_acc(ii)/3.15d13,t_des(ii)/3.15d13,' Myr'
                    write(*,'(e10.3,2I10)')error_rel,i,nleaf
                    stop
                 endif
                 if(uold(ind_leaf(i),idust)<0.)then !!error, stop do while
                    write(*,'(A,e10.2,I10)')'stopping in dust processing rhoD<0',rhoD(ii),icount(ii)
                    write(*,*) 'bin :', ii
                    write(*,*)'rhoG      rhoZ     rhoGZ      rhoD      Temperature'
                    write(*,'(5e10.2)')rhoG0*scale_nH,rhoZ0/rhoG0/0.02,rhoGZ0 /rhoG0/0.02,rhoD0(ii) /rhoG0/0.02,T6*1d6
                    write(*,'(5e10.2)')rhoG0*scale_nH,rhoZ0/rhoG0/0.02,rhoGZ00/rhoG0/0.02,rhoD00(ii)/rhoG0/0.02,T6*1d6
                    write(*,*)'dtloc     dtstep     dtremain'
                    write(*,'(3e10.3,A)')dtloc_bin(ii)/3.15d13,dtcool/3.15d13,dtremain/3.15d13,' Myr'
                    write(*,*)'t_acc     t_des'
                    write(*,'(2e10.3,A)')t_acc(ii)/3.15d13,t_des(ii)/3.15d13,' Myr'
                    write(*,'(e10.3,2I10)')error_rel,i,nleaf
                    stop
                 endif
                 if(.not.okdust .and. i==1) then
                    write(*,*)'BIN',ii, ': icount = ', icount(ii)
                    write(*,'(A,E12.5)') ' fdust =', uold(ind_leaf(i),idust-1+ii)/nH(i)
                    write(*,'(A,E12.5,A,E12.5)') 'rhoGZ0 =', rhoGZ0, ', rhoZ0 =', rhoZ0
                    write(*,'(A,E12.5,A,E12.5,A,E12.5)') ' rhoD0 =', rhoD0(ii), ', drhoD =',drhoD(ii), ', rhoD =', rhoD(ii)
                    write(*,'(A,E12.5)') 'rhoDT0 =', rhoDT0
                    write(*,*)'dtloc     dtstep     dtremain'
                    write(*,'(3e10.3,A)')dtloc_bin(ii)/3.15d13,dtcool/3.15d13,dtremain/3.15d13,' Myr'
                    write(*,*)'t_acc     t_des     t_sha     t_coa'
                    write(*,'(4e10.3,A)')t_acc(ii)/3.15d13,t_des(ii)/3.15d13,t_sha/3.15d13,t_coa/3.15d13,' Myr'
                    write(*,'(A,E12.5,A,E12.5)') 'error rel =', error_rel, ', error max =', errmax
                 endif
              enddo !!on bins
              if(.not.okdust)then !!if still in the do while (rk4)
                 if(ALL(okdt_bin).eqv..true.) dtremain=dtremain-dtloc
                 if(dtremain.le.0.0d0)okdust=.true.
              endif
           enddo !!on do while
           rhoDT=0.0d0
           do ii=1,ndust
              rhoDT=rhoDT+uold(ind_leaf(i),idust-1+ii)
           enddo
        enddo !!on nleaf
        else !!if dust_dest_within_cool=.true. -> compute in solve_cooling
           do i=1,nleaf
              do ii=1,ndust
                if(fdust(i,ii)<0.0)then
                   write(*,*)'END COOLING : fdust<0 on bin :',ii
                   write(*,*) ind_leaf(i),uold(ind_leaf(i),idust-1+ii),nH(i)
                endif
                uold(ind_leaf(i),idust-1+ii)=fdust(i,ii)*uold(ind_leaf(i),1) !!dust_dens=dust2gas_ratio*gas_dens if dust_dens>0
        enddo
           enddo
        endif

        if(zdmax.gt.0d0)then
           do i=1,nleaf
              rhoDT=0.0d0
              do ii=1,ndust
                 rhoDT=rhoDT+uold(ind_leaf(i),idust-1+ii)
           enddo
              do ii=1,ndust !!$dust_dev
                 if(uold(ind_leaf(i),idust-1+ii)/uold(ind_leaf(i),imetal).gt.zdmax)then
                    uold(ind_leaf(i),idust-1+ii)=zdmax*uold(ind_leaf(i),imetal)*uold(ind_leaf(i),idust-1+ii)/rhoDT
        endif
              enddo
           enddo
     endif

!!$        if(dust_chem)then
!!$           do i=1,nleaf
!!$              do ii=1,ndchemtype
!!$                 if(myid==25.and.ind_leaf(i)==504376)write(*,'(A,I2,3es13.5,I4,es13.5)')'(1) C  dust',ii,uold(ind_leaf(i),idust-1+ii)/nH(i)/0.02,nh(i)*scale_nH,T2(i)/nH(i)*(gamma-1.0)*scale_T2,myid,uold(ind_leaf(i),ichem+5-1)/nH(i)/0.02
!!$              enddo
!!$              do ii=ndchemtype+1,ndust
!!$                 if(myid==25.and.ind_leaf(i)==504376)write(*,'(A,I2,3es13.5,I4,2es13.5)')'(1) Si dust',ii,uold(ind_leaf(i),idust-1+ii)/SioverSil/nH(i)/0.02,nh(i)*scale_nH,T2(i)/nH(i)*(gamma-1.0)*scale_T2,myid,uold(ind_leaf(i),ichem+7-1),uold(ind_leaf(i),idust-1+ii)/nH(i)/0.02
!!$              enddo
!!$           enddo
!!$        endif

     endif !!endif(dust)

#ifdef RT
     if(neq_chem) then
        ! Update ionization fraction
        do ii=0,nIons-1
           do i=1,nleaf
              uold(ind_leaf(i),iIons+ii) = xion(1+ii,i)*nH(i)
           end do
        end do
     endif
#if NGROUPS>0
     if(rt) then
        ! Update photon densities and flux magnitudes
        do ig=1,nGroups
           do i=1,nleaf
              rtuold(ind_leaf(i),iGroups(ig)) = (Np(ig,i)-Np_boost(ig,i)) /scale_Np
              rtuold(ind_leaf(i),iGroups(ig)) = &
                   max(rtuold(ind_leaf(i),iGroups(ig)),smallNp)
              rtuold(ind_leaf(i),iGroups(ig)+1:iGroups(ig)+ndim)         &
                               = (Fp(1:ndim,ig,i)-Fp_boost(1:ndim,ig,i)) /scale_Fp
           enddo
        end do
     endif

     ! Split IR photons into trapped and freeflowing
     if(rt_isIRtrap) then
        if(nener .le. 0) then
           print*,'Trying to store E_trapped pressure, but NERAD too small!!'
           STOP
        endif
        iNp=iGroups(iIR)
        unit_tau = 1.5d0 * dx_loc * scale_d * scale_l
        do i=1,nleaf
           il=ind_leaf(i)
           NIRtot =max(rtuold(il,iNp),smallNp)      ! Total photon density
           kIR  = kappaSc(iIR)
           if(is_kIR_T) then                        ! kIR depends on T_rad
              ! For rad. temperature,  weigh the energy in each group by
              ! its opacity over IR opacity (derived from IR temperature)
              E_rad = group_egy(iIR) * ev_to_erg * NIRtot * scale_Np
              TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)     ! IR temp.
              kIR = kappaAbs(iIR) * (TR/10d0)**2
              do ig=1,nGroups
                 if(i .ne. iIR)                                          &
                      E_rad = E_rad + kappaAbs(ig) / kIR                 &
                            * max(rtuold(il,iGroups(ig)),smallNp)        &
                            * ev_to_erg * scale_Np
              end do
              TR = max(0d0,(E_rad*rt_c_fraction/a_r)**0.25)   ! Rad. temp.
              ! Set the IR opacity according to the rad. temperature:
              kIR  = kappaSc(iIR)  * (TR/10d0)**2 * exp(-TR/1d3)
           endif
           f_dust = 1.-xion(ixHII,i)              ! No dust in ionised gas
           tau = nH(i) * Zsolar(i) * f_dust * unit_tau * kIR
           f_trap = 0d0             ! Fraction IR photons that are trapped
           if(tau .gt. 0d0) f_trap = min(max(exp(-1d0/tau), 0d0), 1d0)
           ! Update streaming photons, trapped photons, and tot energy:
           rtuold(il,iNp) = max(smallnp,(1d0-f_trap) * NIRtot) ! Streaming
           ! Limit streaming flux
           rtuold(il,iNp+1:iNp+ndim) = &
                                  rtuold(il,iNp+1:iNp+ndim) * (1d0-f_trap)
           EIR_trapped = max(0d0, NIRtot-rtuold(il,iNp)) * Np2Ep ! Trapped
           ! Update tot energy due to change in trapped radiation energy:
           uold(il,ndim+2)=uold(il,ndim+2)-uold(il,iIRtrapVar)+EIR_trapped
           ! Update the trapped photon energy:
           uold(il,iIRtrapVar) = EIR_trapped

           call reduce_flux(rtuold(il,iNp+1:iNp+ndim),rtuold(il,iNp)*rt_c)
        end do ! i=1,nleaf

     endif  !rt_isIRtrap
#endif
#endif

  end do
  ! End loop over cells

end subroutine coolfine1

#ifdef RT
!************************************************************************
subroutine cmp_Eddington_tensor(Npc,Fp,T_Edd)

! Compute Eddington tensor for given radiation variables
! Npc     => Photon number density times light speed
! Fp     => Photon number flux
! T_Edd  <= Returned Eddington tensor
!------------------------------------------------------------------------
  use amr_commons
  implicit none
  real(dp)::Npc
  real(dp),dimension(1:ndim)::Fp ,u
  real(dp),dimension(1:ndim,1:ndim)::T_Edd
  real(dp)::iterm,oterm,Np_c_sq,Fp_sq,fred_sq,chi
  integer::p,q
!------------------------------------------------------------------------
  if(Npc .le. 0.d0) then
     write(*,*)'negative photon density in cmp_Eddington_tensor. -EXITING-'
     call clean_stop
  endif
  T_Edd(:,:) = 0.d0
  Np_c_sq = Npc**2
  Fp_sq = sum(Fp**2)              !  Sq. photon flux magnitude
  u(:) = 0.d0                           !           Flux unit vector
  if(Fp_sq .gt. 0.d0) u(:) = Fp/sqrt(Fp_sq)
  fred_sq = Fp_sq/Np_c_sq           !      Reduced flux, squared
  chi = max(4.d0-3.d0*fred_sq, 0.d0)   !           Eddington factor
  chi = (3.d0+ 4.d0*fred_sq)/(5.d0 + 2.d0*sqrt(chi))
  iterm = (1.d0-chi)/2.d0               !    Identity term in tensor
  oterm = (3.d0*chi-1.d0)/2.d0          !         Outer product term
  do p = 1, ndim
     do q = 1, ndim
        T_Edd(p,q) = oterm * u(p) * u(q)
     enddo
     T_Edd(p,p) = T_Edd(p,p) + iterm
  enddo

end subroutine cmp_Eddington_tensor
#endif
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_sigma_turb(icell,sigma2,ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel,icell,ncell
  real(dp)::d,d1,d2,d3,d4,d5,d6,ul,ur
  real(dp)::sigma2,sigma2_comp,sigma2_sole
  integer ,dimension(1:nvector)::ind_cell2
  integer ,dimension(1:nvector,0:twondim)::ind_nbor

  ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
  ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor.
  ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation
  ! from neighbouring cell values and differentiate.
  ! Get neighbor cells if they exist, otherwise use straight injection from local cell
  ncell = 1 ! we just want the neighbors of that cell
  ind_cell2(1)=icell
  d=uold(icell,1)
  call getnbor(ind_cell2,ind_nbor,ncell,ilevel)
  d1           = uold(ind_nbor(1,1),1) ; d2 = uold(ind_nbor(1,2),1) ; d3 = uold(ind_nbor(1,3),1)
  d4           = uold(ind_nbor(1,4),1) ; d5 = uold(ind_nbor(1,5),1) ; d6 = uold(ind_nbor(1,6),1)
  sigma2       = 0d0 ; sigma2_comp = 0d0 ; sigma2_sole = 0d0
  !!!!!!!!!!!!!!!!!!
  ! Divergence terms
  !!!!!!!!!!!!!!!!!!
  ul        = (uold(ind_nbor(1,2),2) + uold(icell,2))/(d2+d)
  ur        = (uold(ind_nbor(1,1),2) + uold(icell,2))/(d1+d)
  sigma2_comp = sigma2_comp + (ur-ul)**2
  ul        = (uold(ind_nbor(1,4),3) + uold(icell,3))/(d4+d)
  ur        = (uold(ind_nbor(1,3),3) + uold(icell,3))/(d3+d)
  sigma2_comp = sigma2_comp + (ur-ul)**2
  ul        = (uold(ind_nbor(1,6),4) + uold(icell,4))/(d6+d)
  ur        = (uold(ind_nbor(1,5),4) + uold(icell,4))/(d5+d)
  sigma2_comp = sigma2_comp + (ur-ul)**2
  !!!!!!!!!!!!
  ! Curl terms
  !!!!!!!!!!!!
  ul        = (uold(ind_nbor(1,6),3) + uold(icell,3))/(d6+d)
  ur        = (uold(ind_nbor(1,5),3) + uold(icell,3))/(d5+d)
  sigma2_sole = sigma2_sole + (ur-ul)**2
  ul        = (uold(ind_nbor(1,4),4) + uold(icell,4))/(d4+d)
  ur        = (uold(ind_nbor(1,3),4) + uold(icell,4))/(d3+d)
  sigma2_sole = sigma2_sole + (ur-ul)**2
  ul        = (uold(ind_nbor(1,6),2) + uold(icell,2))/(d6+d)
  ur        = (uold(ind_nbor(1,5),2) + uold(icell,2))/(d5+d)
  sigma2_sole = sigma2_sole + (ur-ul)**2
  ul        = (uold(ind_nbor(1,2),4) + uold(icell,4))/(d2+d)
  ur        = (uold(ind_nbor(1,1),4) + uold(icell,4))/(d1+d)
  sigma2_sole = sigma2_sole + (ur-ul)**2
  ul        = (uold(ind_nbor(1,4),2) + uold(icell,2))/(d4+d)
  ur        = (uold(ind_nbor(1,3),2) + uold(icell,2))/(d3+d)
  sigma2_sole = sigma2_sole + (ur-ul)**2
  ul        = (uold(ind_nbor(1,2),3) + uold(icell,3))/(d2+d)
  ur        = (uold(ind_nbor(1,1),3) + uold(icell,3))/(d1+d)
  sigma2_sole = sigma2_sole + (ur-ul)**2
  sigma2    = sigma2_comp+sigma2_sole

end subroutine cmp_sigma_turb
!###########################################################
!###########################################################
!###########################################################
!###########################################################
!!$subroutine cmp_sigma_turb2(ind_cell,ncell,ilevel)
!!$  implicit none
!!$  use amr_commons
!!$  use pm_commons
!!$  use hydro_commons
!!$  use poisson_commons
!!$  integer::ilevel
!!$  integer ::i,ncell,nx_loc
!!$  real(dp),dimension(1:3)::skip_loc
!!$  real(dp)::d
!!$  real(dp)::ul,ur,fl,fr,trgv
!!$  real(dp)::sigma2,sigma2_comp,sigma2_sole
!!$  real(dp)::divv,divv2,curlv,curlva,curlvb,curlvc,curlv2
!!$  real(dp)::dx,dx_loc,scale,d1,d2,d3,d4,d5,d6
!!$  integer ,dimension(1:ncell)::ind_cell,ind_cell2
!!$  integer ,dimension(1:ncell,0:twondim)::ind_nbor
!!$
!!$  ! Mesh spacing in that level
!!$  dx=0.5D0**ilevel
!!$  nx_loc=(icoarse_max-icoarse_min+1)
!!$  skip_loc=(/0.0d0,0.0d0,0.0d0/)
!!$  if(ndim>0)skip_loc(1)=dble(icoarse_min)
!!$  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
!!$  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
!!$  scale=boxlen/dble(nx_loc)
!!$  dx_loc=dx*scale
!!$
!!$  do i=1,ncell
!!$     ! We need to estimate the norm of the gradient of the velocity field in the cell (tensor of 2nd rank)
!!$     ! i.e. || A ||^2 = trace( A A^T) where A = grad vec(v) is the tensor.
!!$     ! So construct values of velocity field on the 6 faces of the cell using simple linear interpolation
!!$     ! from neighbouring cell values and differentiate.
!!$     ! Get neighbor cells if they exist, otherwise use straight injection from local cell
!!$     ncell2 = 1 ! we just want the neighbors of that cell
!!$     ind_cell2(1) = ind_cell(i)
!!$     call getnbor(ind_cell2,ind_nbor,ncell2,ilevel)
!!$     d1           = uold(ind_nbor(1,1),1) ; d2 = uold(ind_nbor(1,2),1) ; d3 = uold(ind_nbor(1,3),1)
!!$     d4           = uold(ind_nbor(1,4),1) ; d5 = uold(ind_nbor(1,5),1) ; d6 = uold(ind_nbor(1,6),1)
!!$     sigma2       = 0d0 ; sigma2_comp = 0d0 ; sigma2_sole = 0d0
!!$     trgv         = 0d0 ; divv = 0d0 ; curlva = 0d0 ; curlvb = 0d0 ; curlvc = 0d0
!!$     flong        = 0d0
!!$!!!!!!!!!!!!!!!!!!
!!$     ! Divergence terms
!!$!!!!!!!!!!!!!!!!!!
!!$     ul        = (d2*uold(ind_nbor(1,2),2) + d*uold(ind_cell(i),2))/(d2+d)
!!$     ur        = (d1*uold(ind_nbor(1,1),2) + d*uold(ind_cell(i),2))/(d1+d)
!!$     sigma2_comp = sigma2_comp + (ur-ul)**2
!!$     divv      = divv + (ur-ul)
!!$     ul        = (d4*uold(ind_nbor(1,4),3) + d*uold(ind_cell(i),3))/(d4+d)
!!$     ur        = (d3*uold(ind_nbor(1,3),3) + d*uold(ind_cell(i),3))/(d3+d)
!!$     sigma2_comp = sigma2_comp + (ur-ul)**2
!!$     divv      = divv + (ur-ul)
!!$     ul        = (d6*uold(ind_nbor(1,6),4) + d*uold(ind_cell(i),4))/(d6+d)
!!$     ur        = (d5*uold(ind_nbor(1,5),4) + d*uold(ind_cell(i),4))/(d5+d)
!!$     sigma2_comp = sigma2_comp + (ur-ul)**2
!!$     divv      = divv + (ur-ul)
!!$     ftot      = flong
!!$!!!!!!!!!!!!
!!$     ! Curl terms
!!$!!!!!!!!!!!!
!!$     ul        = (d6*uold(ind_nbor(1,6),3) + d*uold(ind_cell(i),3))/(d6+d)
!!$     ur        = (d5*uold(ind_nbor(1,5),3) + d*uold(ind_cell(i),3))/(d5+d)
!!$     sigma2_sole = sigma2_sole + (ur-ul)**2
!!$     curlva    = curlva-(ur-ul)
!!$     ul        = (d4*uold(ind_nbor(1,4),4) + d*uold(ind_cell(i),4))/(d4+d)
!!$     ur        = (d3*uold(ind_nbor(1,3),4) + d*uold(ind_cell(i),4))/(d3+d)
!!$     sigma2_sole = sigma2_sole + (ur-ul)**2
!!$     curlva    = (curlva + (ur-ul))
!!$     ul        = (d6*uold(ind_nbor(1,6),2) + d*uold(ind_cell(i),2))/(d6+d)
!!$     ur        = (d5*uold(ind_nbor(1,5),2) + d*uold(ind_cell(i),2))/(d5+d)
!!$     sigma2_sole = sigma2_sole + (ur-ul)**2
!!$     curlvb    = curlvb+(ur-ul)
!!$     ul        = (d2*uold(ind_nbor(1,2),4) + d*uold(ind_cell(i),4))/(d2+d)
!!$     ur        = (d1*uold(ind_nbor(1,1),4) + d*uold(ind_cell(i),4))/(d1+d)
!!$     sigma2_sole = sigma2_sole + (ur-ul)**2
!!$     curlvb    = (curlvb - (ur-ul))
!!$     ul        = (d4*uold(ind_nbor(1,4),2) + d*uold(ind_cell(i),2))/(d4+d)
!!$     ur        = (d3*uold(ind_nbor(1,3),2) + d*uold(ind_cell(i),2))/(d3+d)
!!$     sigma2_sole = sigma2_sole + (ur-ul)**2
!!$     curlvc    = curlvc-(ur-ul)
!!$     ul        = (d2*uold(ind_nbor(1,2),3) + d*uold(ind_cell(i),3))/(d2+d)
!!$     ur        = (d1*uold(ind_nbor(1,1),3) + d*uold(ind_cell(i),3))/(d1+d)
!!$     sigma2_sole = sigma2_sole + (ur-ul)**2
!!$     curlvc    = (curlvc + (ur-ul))
!!$     sigma2    = sigma2_comp+sigma2_sole
!!$     ! Trace of gradient velocity tensor
!!$     trgv      = sigma2/dx_loc**2
!!$     ! Velocity vector divergence
!!$     divv      = divv/dx_loc
!!$     ! Velocity vector curl
!!$     curlv     = (curlva+curlvb+curlvc)/dx_loc
!!$     divv2     = divv**2
!!$     curlv2    = curlv**2
!!$  enddo
!!$
!!$end subroutine cmp_sigma_turb2
