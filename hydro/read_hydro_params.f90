subroutine read_hydro_params(nml_ok)
  use amr_commons
  use hydro_commons
  use mpi_mod
  implicit none
  logical::nml_ok
  !--------------------------------------------------
  ! Local variables
  !--------------------------------------------------
  integer::i,idim,ifixed,nboundary_true=0,ipar,ich
  integer ,dimension(1:MAXBOUND)::bound_type
  real(dp)::ek_bound,eta_sn_ini
  real(dp)::numtot,fracMgD,fracFeD,fracSiD,fracOD
  logical :: dummy
#ifdef SOLVERmhd
  real(dp)::em_bound
#endif

  !--------------------------------------------------
  ! Namelist definitions
  !--------------------------------------------------

  ! Initial conditions parameters
  namelist/init_params/filetype,initfile,multiple,nregion,region_type &
       & ,x_center,y_center,z_center,aexp_ini &
       & ,length_x,length_y,length_z,exp_region &
       & ,d_region,u_region,v_region,w_region,p_region &
#ifdef SOLVERmhd
       & ,A_region,B_region,C_region,B_ave &
#if NVAR>8+NENER
       & ,var_region &
#endif
#else
#if NVAR>NDIM+2+NENER
       & ,var_region &
#endif
#endif
#if NENER>0
       & ,prad_region &
#endif
       & ,omega_b,galage

  ! Hydro parameters
  namelist/hydro_params/gamma,courant_factor,smallr,smallc &
       & ,niter_riemann,slope_type,difmag &
#if NENER>0
       & ,gamma_rad &
#endif
#ifdef SOLVERmhd
       & ,riemann2d,slope_mag_type,eta_mag &
#endif
       & ,pressure_fix,beta_fix,scheme,riemann,frozen,checkhydro

  ! Refinement parameters
  namelist/refine_params/x_refine,y_refine,z_refine,r_refine &
       & ,a_refine,b_refine,exp_refine,jeans_refine,mass_cut_refine &
       & ,m_refine,mass_sph,err_grad_d,err_grad_p,err_grad_u &
       & ,floor_d,floor_u,floor_p,ivar_refine,var_cut_refine,trans_smooth &
#ifdef SOLVERmhd
       & ,err_grad_A,err_grad_B,err_grad_C,err_grad_B2 &
       & ,floor_A,floor_B,floor_C,floor_B2,interpol_mag_type &
#endif
       & ,interpol_var,interpol_type,sink_refine,del_jeans &
       & ,level_zoom,xzoom,yzoom,zzoom,rzoom,dens_jeans

  ! Boundary parameters
  namelist/boundary_params/nboundary,bound_type &
       & ,ibound_min,ibound_max,jbound_min,jbound_max &
       & ,kbound_min,kbound_max &
#if NENER>0
       & ,prad_bound &
#endif
#ifdef SOLVERmhd
#if NVAR>8+NENER
       & ,var_bound &
#endif
       & ,A_bound,B_bound,C_bound &
#else
#if NVAR>NDIM+2+NENER
       & ,var_bound &
#endif
#endif
       & ,d_bound,u_bound,v_bound,w_bound,p_bound,no_inflow

  ! Feedback parameters
  namelist/feedback_params/eta_sn,eta_ssn,yield,rbubble,f_ek,ndebris &
       & ,f_w,mass_gmc,kappa_IR,delayed_cooling,momentum_feedback &
       & ,A_SN,A_SN_geen,expN_SN,E_SNII,M_SNII,SF_kick_kms,SN_dT2_min &
       & ,log_mfb,log_mfb_mega,mechanical_feedback,mechanical_geen &
       & ,nsn2mass,sn2_real_delay,sn_IC,sn_trelax &
       & ,ir_feedback,ir_eff,t_diss,t_sne,mass_star_max,mass_sne_min &
       & ,stellar_winds,stellar_winds_file,supernovae_II_file,chem_list,snII_freq &
       & ,snIa,E_SNIa,phi_snIa,t_ini_SNIa,t_fin_snIa,SNII_zdep_yield &
       & ,use_initial_mass,no_wind_energy,snyield_model

  ! Cooling / basic chemistry parameters
  namelist/cooling_params/cooling,metal,isothermal,haardt_madau,J21 &
       & ,dust,dust_cooling,sticking_coef,zdmax,metal_gasonly,dust_dest_within_cool,errmax,asize &
       & ,dust_SNdest,dust_cond_eff,dust_cond_eff_Ia,dust_SNdest_eff,thermal_sputtering &
       & ,dust_acc_neglected_large_bin,dust_prod_neglected_small_bin &
       & ,dust_shattering,dust_coagulation &
       & ,dust_sputtering,dust_accretion,dust_turb &
       & ,t_sputter_ref, t_growth_ref, boost_growth,nhboost_growth, Tdust_growth, nh_growth, Sconstant &
       & ,t_sha_ref, t_coa_ref, nh_coa, power_coa &
       & ,a_spec,self_shielding,z_ave,z_reion,T2max,neq_chem &
       & ,muMg,muFe,muSi,muO,nsilMg,nsilFe,nsilSi,nsilO &
       & ,fsmall_ej,dust_chem,dustdebug,sgrain,DTMini,fsmall_ini

  ! Star formation parameters
  namelist/sf_params/m_star,t_star,n_star,T2_star,g_star,del_star &
       & ,eps_star,jeans_ncells,sf_virial,sf_trelax,sf_tdiss,sf_model &
       & ,fstar_min,star_imf,sf_lamjt,write_stellar_densities &
       & ,sf_mach_threshold,sf_log_properties,sf_imf,sf_compressive,tsfr_damp_IC

  ! Sink-SMBH parameters
  namelist/smbh_params/agn,smbh,sink_AGN,bondi,drag,spin_bh,force_exact_mseed &
       & ,bhspinmerge,vrel_merge,random_jet,mad_jet,selfgrav,Mseed,n_sink,ns_sink,ns2_sink &
       & ,eAGN_K,eAGN_T,X_floor,boost_acc,boost_drag,boost_drag_part,T2maxAGN,TAGN,mloadAGN &
       & ,f_bondi,ind_rsink,rAGN,rAGN_dx,r_gal,r_bhr,rmerge,sigmav_max,star_ratio_floor &
       & ,drag_part,DF_ncells,d_boost,no_accretion,fix_smbh_position &
       & ,eddington_limit,jetfrac,maximum_accretion,finestep_AGN &
       & ,point_mass_sink,n_gal,sig_sink,t_que,weighted_drag,adfmax,stellar_velocity_seed &
       & ,ns_sink_scaled

  ! Units parameters
  namelist/units_params/units_density,units_time,units_length

  ! Dummy namelist for physics
  namelist/physics_params/ dummy

#ifdef grackle
   namelist/grackle_params/use_grackle,grackle_with_radiative_cooling,grackle_primordial_chemistry,grackle_metal_cooling &
       & ,grackle_UVbackground,grackle_cmb_temperature_floor,grackle_h2_on_dust,grackle_photoelectric_heating &
       & ,grackle_use_volumetric_heating_rate,grackle_use_specific_heating_rate,grackle_three_body_rate,grackle_cie_cooling &
       & ,grackle_h2_optical_depth_approximation,grackle_ih2co,grackle_ipiht,grackle_NumberOfTemperatureBins,grackle_CaseBRecombination &
       & ,grackle_Compton_xray_heating,grackle_LWbackground_sawtooth_suppression,grackle_NumberOfDustTemperatureBins,grackle_use_radiative_transfer &
       & ,grackle_radiative_transfer_coupled_rate_solver,grackle_radiative_transfer_intermediate_step,grackle_radiative_transfer_hydrogen_only &
       & ,grackle_self_shielding_method,grackle_Gamma,grackle_photoelectric_heating_rate,grackle_HydrogenFractionByMass &
       & ,grackle_DeuteriumToHydrogenRatio,grackle_SolarMetalFractionByMass,grackle_TemperatureStart,grackle_TemperatureEnd &
       & ,grackle_DustTemperatureStart,grackle_DustTemperatureEnd,grackle_LWbackground_intensity,grackle_UVbackground_redshift_on &
       & ,grackle_UVbackground_redshift_off,grackle_UVbackground_redshift_fullon,grackle_UVbackground_redshift_drop &
       & ,grackle_cloudy_electron_fraction_factor,grackle_data_file
#endif

  ! Read namelist file
  rewind(1)
  read(1,NML=init_params,END=101)
  goto 102
101 write(*,*)' You need to set up namelist &INIT_PARAMS in parameter file'
  call clean_stop
102 rewind(1)

  ! Fail if physics params is found
  read(1, NML=physics_params, end=110)
  if (myid == 1) &
       write(*, *) 'ERROR: the namelist contains the old `physics_params`. It has been depreciated in favor of the sections feedback_, cooling_, smbh_, sf_ and units_params.'
  call clean_stop
110 rewind(1)
  if(nlevelmax>levelmin)read(1,NML=refine_params)
  rewind(1)
  if(hydro)read(1,NML=hydro_params)
  rewind(1)
  read(1,NML=boundary_params,END=103)
  simple_boundary=.true.
  goto 104
103 simple_boundary=.false.
104 if(nboundary>MAXBOUND)then
    write(*,*) 'Error: nboundary>MAXBOUND'
    call clean_stop
  end if
  rewind(1)
  read(1,NML=feedback_params,END=105)
105 continue
  rewind(1)
  read(1,NML=cooling_params,END=106)
106 continue
  rewind(1)
  read(1,NML=sf_params,END=107)
107 continue
  rewind(1)
  read(1,NML=smbh_params,END=108)
108 continue
  rewind(1)
  read(1,NML=units_params,END=109)
109 continue
#ifdef grackle
  rewind(1)
  read(1,NML=grackle_params)
#endif
#ifdef ATON
  if(aton)call read_radiation_params(1)
#endif

! MHD sets enums for riemann solvers, hydro just does a string compare
#ifdef SOLVERmhd
  !------------------------------------------------
  ! set ischeme
  !------------------------------------------------
  SELECT CASE (scheme)
  CASE ('muscl')
    ischeme = 0
  CASE ('induction')
    ischeme = 1

  CASE DEFAULT
    write(*,*)'unknown scheme'
    call clean_stop
  END SELECT
  !------------------------------------------------
  ! set iriemann
  !------------------------------------------------
  SELECT CASE (riemann)
  CASE ('llf')
    iriemann = 0
  CASE ('roe')
    iriemann = 1
  CASE ('hll')
    iriemann = 2
  CASE ('hlld')
    iriemann = 3
  CASE ('upwind')
    iriemann = 4
  CASE ('hydro')
    iriemann = 5

  CASE DEFAULT
    write(*,*)'unknown riemann solver'
    call clean_stop
  END SELECT
  !------------------------------------------------
  ! set iriemann
  !------------------------------------------------
  SELECT CASE (riemann2d)
  CASE ('llf')
    iriemann2d = 0
  CASE ('roe')
    iriemann2d = 1
  CASE ('upwind')
    iriemann2d = 2
  CASE ('hll')
    iriemann2d = 3
  CASE ('hlla')
    iriemann2d = 4
  CASE ('hlld')
    iriemann2d = 5
  CASE DEFAULT
    write(*,*)'unknown 2D riemann solver'
    call clean_stop
  END SELECT

  !--------------------------------------------------
  ! Make sure virtual boundaries are expanded to
  ! account for staggered mesh representation
  !--------------------------------------------------
  nexpand_bound=2
#endif

  !--------------------------------------------------
  ! Check for dm only cosmo run
  !--------------------------------------------------
  if(.not.hydro)then
     omega_b = 0.0D0
  endif

  !--------------------------------------------------
  ! Check for star formation
  !--------------------------------------------------
  if(t_star>0)then
     star=.true.
     if (.not. pic) then
        write(*, *) 'Warning: activating PIC because of star formation'
     end if
     pic=.true.
  else if(eps_star>0)then
     t_star=0.1635449*(n_star/0.1)**(-0.5)/eps_star
     if (.not. pic) then
        write(*, *) 'Warning: activating PIC because of star formation'
     end if
     star=.true.
     pic=.true.
  endif

  !--------------------------------------------------
  ! Check for metal
  !--------------------------------------------------
#ifdef SOLVERmhd
  if(metal.and.nvar<(ndim+6))then
#else
  if(metal.and.nvar<(ndim+3))then
#endif
     if(myid==1)write(*,*)'Error: metals need nvar >= ndim+3'
     if(myid==1)write(*,*)'Set METALS = 1 in the Makefile and recompile'
     nml_ok=.false.
  endif

  !--------------------------------------------------
  ! Check for non-thermal energies
  !--------------------------------------------------
#if NENER>0
#ifdef SOLVERmhd
  if(nvar<(8+nener))then
#else
  if(nvar<(ndim+2+nener))then
#endif
     if(myid==1)write(*,*)'Error: non-thermal energy need nvar >= ndim+2+nener'
     if(myid==1)write(*,*)'Modify NENER and recompile'
     nml_ok=.false.
  endif
#endif

  !--------------------------------------------------
  ! Check ind_rsink
  !--------------------------------------------------
  if(ind_rsink<=0.0d0)then
     if(myid==1)write(*,*)'Error in the namelist'
     if(myid==1)write(*,*)'Check ind_rsink'
     nml_ok=.false.
  end if

  !--------------------------------------------------
  ! Check initial mass [TK]
  !--------------------------------------------------
  if(sn2_real_delay.or.snIa) use_initial_mass=.true. ! for continuous SN2 or SNIa
  if(sn2_real_delay)         t_sne=50.             ! SNe explode up to 50 Myr.
  !Note: t_sne<<50 would reduce the SN rate per particle if sn2_real_delay=.true.

  !-------------------------------------------------
  if(dust.and.(.not. metal))then
     if(myid==1)write(*,*)'Error: dust requires metal=.true.'
     nml_ok=.false.
  endif
  if(dust.and.dust_cooling)dust_dest_within_cool=.true.
#ifdef SOLVERmhd
  if(dust.and.nvar<(ndim+7))then
     if(myid==1)write(*,*)'Error: dust needs nvar >= ndim+7'
#else
  if(dust.and.nvar<(ndim+4))then
     if(myid==1)write(*,*)'Error: dust needs nvar >= ndim+4'
#endif
     if(myid==1)write(*,*)'Make sure METALS = 1 and NVAR_EXTRA > 0 in the Makefile and recompile'
     nml_ok=.false.
  endif

  !--------------------------------------------------
  ! IMF-dependent M_SNII and mass loss [TK]
  !--------------------------------------------------
  eta_sn_ini = eta_sn
  if(TRIM(star_imf).eq.'salpeter')then !0.08-100
     M_SNII=18.728254
     eta_sn=0.12761627
  else if(TRIM(star_imf).eq.'kroupa')then
     M_SNII=19.134730
     eta_sn=0.20913717
  else if(TRIM(star_imf).eq.'chabrier')then
     M_SNII=10. !TODO
     eta_sn=0.313706  !TODO
  else if(TRIM(star_imf).ne.'')then
     write(*,*) 'your star_imf seems wrong -->'//TRIM(star_imf)
     stop
  endif
  if (eta_sn_ini .lt. 0) eta_sn = 0.0d0
  if(myid==1) then
     if(SNII_zdep_yield)then
        write(*,*) '>>> SN parameters: SNII_zdep_yield=.true.'
     else
        write(*,*) '>>> SN parameters: M_SNII, eta_sn, t_sne =',sngl(M_SNII), sngl(eta_sn), sngl(t_sne)
     endif
  endif

  !-------------------------------------------------
  ! This section deals with hydro boundary conditions
  !-------------------------------------------------
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  if (simple_boundary)then

     ! Compute new coarse grid boundaries
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           nx=nx+1
           if(ibound_min(i)==-1)then
              icoarse_min=icoarse_min+1
              icoarse_max=icoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ny=ny+1
           if(jbound_min(i)==-1)then
              jcoarse_min=jcoarse_min+1
              jcoarse_max=jcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do
     do i=1,nboundary
        if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           nz=nz+1
           if(kbound_min(i)==-1)then
              kcoarse_min=kcoarse_min+1
              kcoarse_max=kcoarse_max+1
           end if
           nboundary_true=nboundary_true+1
        end if
     end do

     ! Compute boundary geometry
     do i=1,nboundary
        if(ibound_min(i)*ibound_max(i)==1.and.ndim>0.and.bound_type(i)>0)then
           if(ibound_min(i)==-1)then
              ibound_min(i)=icoarse_min+ibound_min(i)
              ibound_max(i)=icoarse_min+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=1
              if(bound_type(i)==2)boundary_type(i)=11
              if(bound_type(i)==3)boundary_type(i)=21
           else
              ibound_min(i)=icoarse_max+ibound_min(i)
              ibound_max(i)=icoarse_max+ibound_max(i)
              if(bound_type(i)==1)boundary_type(i)=2
              if(bound_type(i)==2)boundary_type(i)=12
              if(bound_type(i)==3)boundary_type(i)=22
           end if
           if(ndim>1)jbound_min(i)=jcoarse_min+jbound_min(i)
           if(ndim>1)jbound_max(i)=jcoarse_max+jbound_max(i)
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(jbound_min(i)*jbound_max(i)==1.and.ndim>1.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           if(jbound_min(i)==-1)then
              jbound_min(i)=jcoarse_min+jbound_min(i)
              jbound_max(i)=jcoarse_min+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=3
              if(bound_type(i)==2)boundary_type(i)=13
              if(bound_type(i)==3)boundary_type(i)=23
           else
              jbound_min(i)=jcoarse_max+jbound_min(i)
              jbound_max(i)=jcoarse_max+jbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=4
              if(bound_type(i)==2)boundary_type(i)=14
              if(bound_type(i)==3)boundary_type(i)=24
           end if
           if(ndim>2)kbound_min(i)=kcoarse_min+kbound_min(i)
           if(ndim>2)kbound_max(i)=kcoarse_max+kbound_max(i)
        else if(kbound_min(i)*kbound_max(i)==1.and.ndim>2.and.bound_type(i)>0)then
           ibound_min(i)=icoarse_min+ibound_min(i)
           ibound_max(i)=icoarse_max+ibound_max(i)
           jbound_min(i)=jcoarse_min+jbound_min(i)
           jbound_max(i)=jcoarse_max+jbound_max(i)
           if(kbound_min(i)==-1)then
              kbound_min(i)=kcoarse_min+kbound_min(i)
              kbound_max(i)=kcoarse_min+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=5
              if(bound_type(i)==2)boundary_type(i)=15
              if(bound_type(i)==3)boundary_type(i)=25
           else
              kbound_min(i)=kcoarse_max+kbound_min(i)
              kbound_max(i)=kcoarse_max+kbound_max(i)
              if(bound_type(i)==1)boundary_type(i)=6
              if(bound_type(i)==2)boundary_type(i)=16
              if(bound_type(i)==3)boundary_type(i)=26
           end if
        end if
     end do
     do i=1,nboundary
        ! Check for errors
        if( (ibound_min(i)<0.or.ibound_max(i)>(nx-1)) .and. (ndim>0) .and.bound_type(i)>0 )then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along X direction',i
           nml_ok=.false.
        end if
        if( (jbound_min(i)<0.or.jbound_max(i)>(ny-1)) .and. (ndim>1) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Y direction',i
           nml_ok=.false.
        end if
        if( (kbound_min(i)<0.or.kbound_max(i)>(nz-1)) .and. (ndim>2) .and.bound_type(i)>0)then
           if(myid==1)write(*,*)'Error in the namelist'
           if(myid==1)write(*,*)'Check boundary conditions along Z direction',i
           nml_ok=.false.
        end if
     end do
  end if
  nboundary=nboundary_true
  if(simple_boundary.and.nboundary==0)then
     simple_boundary=.false.
  endif

  !--------------------------------------------------
  ! Compute boundary conservative variables
  !--------------------------------------------------
  do i=1,nboundary
     boundary_var(i,1)=MAX(d_bound(i),smallr)
     boundary_var(i,2)=d_bound(i)*u_bound(i)
#if NDIM>1 || SOLVERmhd
     boundary_var(i,3)=d_bound(i)*v_bound(i)
#endif
#if NDIM>2 || SOLVERmhd
     boundary_var(i,4)=d_bound(i)*w_bound(i)
#endif
     ek_bound=0.0d0
     do idim=1,ndim
        ek_bound=ek_bound+0.5d0*boundary_var(i,idim+1)**2/boundary_var(i,1)
     end do
     boundary_var(i,ndim+2)=ek_bound+P_bound(i)/(gamma-1.0d0)
#ifdef SOLVERmhd
     boundary_var(i,6)=A_bound(i)
     boundary_var(i,7)=B_bound(i)
     boundary_var(i,8)=C_bound(i)
     boundary_var(i,nvar+1)=A_bound(i)
     boundary_var(i,nvar+2)=B_bound(i)
     boundary_var(i,nvar+3)=C_bound(i)
     ek_bound=0.5d0*d_bound(i)*(u_bound(i)**2+v_bound(i)**2+w_bound(i)**2)
     em_bound=0.5d0*(A_bound(i)**2+B_bound(i)**2+C_bound(i)**2)
     boundary_var(i,5)=ek_bound+em_bound+P_bound(i)/(gamma-1.0d0)
#endif
  end do

  !-----------------------------------
  ! Rearrange level dependent arrays
  !-----------------------------------
  do i=nlevelmax,levelmin,-1
     jeans_refine(i)=jeans_refine(i-levelmin+1)
  end do
  do i=1,levelmin-1
     jeans_refine(i)=-1.0
  end do

  !-----------------------------------
  ! Sort out passive variable indices
  !-----------------------------------
#ifdef SOLVERmhd
  ! Hard-coded variables are rho,v*ndim,P,B*ndim
  ! MHD only works in 3D, so ndim=3
  ifixed=8
#else
  ! Hard-coded variables are rho,v*ndim,P
  ifixed=ndim+2
#endif
  ipar=ifixed+1
  if(nener>0)then
     inener=ipar
     ipar=ipar+nener
  endif
  if(metal)then
     imetal=ipar
     ipar=ipar+1
  endif
  if(nchem>0)then
     ichem=ipar
     ipar=ipar+nchem
  end if
  if(dust)then
     idust=ipar
     ipar=ipar+ndust
     if (ndust>1) then
        if (asize(1)>asize(2)) then
           write(*,*) 'ERROR : dust grains sizes must be in increasing order in asize (amr_parameters)'
           write(*,*) 'asize :', asize
           call clean_stop
        endif
     endif
     if (ndust.eq.1 .and. (dust_coagulation .or. dust_shattering)) then
        write(*,*) "ERROR : ndust can't be !=2 with coagulation/shattering on"
        call clean_stop
     endif
     flarge_ej=MIN(MAX(1d0-fsmall_ej,0.0d0),1.0d0)
     sgrain=sgrain/3d0 ! renormalised to 3g.cm^-3 (implicity assumed in all calculations)
  endif
  if(delayed_cooling)then
     idelay=ipar
     ipar=ipar+1
  endif
  if(sf_virial)then
     ivirial1=ipar
     if(sf_compressive)then
        ivirial2=ivirial1+1
        ipar=ipar+1
     endif
     ixion=ivirial2+1
     ipar=ipar+2
  endif
  if(aton)then
     ichem=ipar
     ipar=ipar+1
  endif
  numtot=muMg*nsilMg+muFe*nsilFe+muSi*nsilSi+muO*nsilO
  MgoverSil=muMg*nsilMg/numtot
  FeoverSil=muFe*nsilFe/numtot
  SioverSil=muSi*nsilSi/numtot
  OoverSil =muO *nsilO /numtot
  multMgoverSi=MgoverSil/SioverSil
  multFeoverSi=FeoverSil/SioverSil
  multOoverSi = OoverSil/SioverSil
  if(dust_chem)then
     do ich=1,nchem
        if(TRIM(chem_list(ich))=='H' )ichH =ich
        if(TRIM(chem_list(ich))=='He')ichHe=ich
        if(TRIM(chem_list(ich))=='N' )ichN =ich
        if(TRIM(chem_list(ich))=='S' )ichS =ich
        if(TRIM(chem_list(ich))=='D' )ichD =ich

        if(TRIM(chem_list(ich))=='C' )ichC =ich
        if(TRIM(chem_list(ich))=='Mg')ichMg=ich
        if(TRIM(chem_list(ich))=='Fe')ichFe=ich
        if(TRIM(chem_list(ich))=='Si')ichSi=ich
        if(TRIM(chem_list(ich))=='O' )ichO =ich
     enddo
     ndchemtype=2
     if(ndust==2)dndsize=0
     if(ndust==4)dndsize=1
  elseif(dust)then
     ndchemtype=1
     if(ndust==1)dndsize=0
     if(ndust==2)dndsize=1
  endif
  if(myid==1.and.hydro.and.(nvar>ndim+2)) then
     write(*,'(A50)')"__________________________________________________"
     write(*,*) 'Hydro var indices:'
#if NENER>0
     write(*,*) '   inener   = ',inener
#endif
     if(metal)           write(*,*) '   imetal   = ',imetal
     if(nchem>0)         write(*,*) '   ichems   = ',ichem,'-',ichem+nchem-1
     if(dust) then
        write(*,*) 'cooling  = ',cooling
        write(*,*) '   =============***Dust Parameters***=============   '
        write(*,*) '   idust      = ',idust,',   ndust    = ',ndust
        write(*,*) 'dust_chem     = ',dust_chem
        if(dust_chem)then
           if(ndust.eq.1)then
              write(*,*)'if you decide to use dust_chem, ndust>1. Stopping...'
              stop
           elseif(ndust==2)then
              write(*,'(I3,A,I3,A)')idust,' for carbon, and',idust+1,' for silicates'
           elseif(ndust==4)then
              write(*,'(I3,A,I3,A,I3,A,I3,A)')idust,' and',idust+1,' for carbon, and',idust+2,' and',idust+3,' for silicates'
           endif
        endif
        write(*,*) 'dust_cooling  = ',dust_cooling
        write(*,*) 'dust_acc      = ',dust_accretion  ,',        dust_sput     = ',dust_sputtering
        write(*,*) 'dust_coa      = ',dust_coagulation,',        dust_sha      = ',dust_shattering
        write(*,*) 'dust prod neglected on small grains       = ', dust_prod_neglected_small_bin
        write(*,*) 'dust acc  neglected on large grains       = ', dust_acc_neglected_large_bin
        write(*,'(A,E10.2,A,E10.2)') ' dust_cond_eff =',dust_cond_eff,',dust_SNdest_eff =',dust_SNdest_eff
        if (ndust==4) then
           write(*,'(A,E10.2,A,E10.2)') ' Sizes (µm) carbon   : S =',asize (1),',              L =',asize (2)
           write(*,'(A,E10.2,A,E10.2)') ' Sizes (µm) silicate : S =',asize (3),',              L =',asize (4)
           write(*,'(A,E10.2,A,E10.2)') ' s (g/cm^3) carbon   : S =',sgrain(1)*3d0,',              L =',sgrain(2)*3d0
           write(*,'(A,E10.2,A,E10.2)') ' s (g/cm^3) silicate : S =',sgrain(3)*3d0,',              L =',sgrain(4)*3d0
        else if(ndust==2)then
           write(*,'(A,E10.2,A,E10.2)') ' Sizes (µm): S =',asize (1),',              L =',asize (2)
           write(*,'(A,E10.2,A,E10.2)') ' s (g/cm^3): S =',sgrain(1),',              L =',sgrain(2)
        endif
        write(*,*)'           t_acc_ref t_des_ref t_sha_ref t_coa_ref'
        write(*,'(A,4e10.3,A)')'           ', t_growth_ref/1d6,t_sputter_ref/1d6,t_sha_ref/1d6,t_coa_ref/1d6,' Myr'
        write(*,'(A,4e10.3,A)')'  Model :  ',4d5/1d6,1d5/1d6,5.41d7/1d6,2.71d5/1d6,' Myr'
        write(*,'(A,2f7.4)')'Respective fraction of small and large grains in stellar ejecta:',fsmall_ej,flarge_ej
        if(dust_chem)then
           write(*,*) 'Mass fraction in silicates of:'
           write(*,*) 'Mg    Fe    Si    O'
           write(*,'(4f6.3)') nsilMg,nsilFe,nsilSi,nsilO
           write(*,'(4f6.3)') MgoverSil,FeoverSil,SioverSil,OoverSil
           write(*,*)'Key element to track silicate is Si'
        endif
        write(*,*) '   ===============================================   '
     endif
     if(delayed_cooling) write(*,*) '   idelay   = ',idelay
     if(sf_virial .and. sf_model /= 6)then
           write(*,*) '   ivirial1 = ',ivirial1
        if(sf_compressive)then
           write(*,*) '   ivirial2 = ',ivirial2
        end if
     endif
     if(aton)            write(*,*) '   ixion    = ',ixion
#ifdef RT
     if(rt) write(*,*) '   iIons    = ',ichem
#endif
     if(ivar_refine>0)   write(*,*) '   refmask  = ',ivar_refine
     write(*,'(A50)')"__________________________________________________"
  endif

  !-----------------------------------
  ! Sort out some of the smbh dependencies
  !-----------------------------------
  if(finestep_AGN.and.(.not.sink_AGN))then
     if(myid==1)then
        write(*,*)"sink_AGN=.false. so AGN feedback is turned off"
        write(*,*)"To use AGN feedback on the timestep, set both finestep_AGN=.true. and sink_AGN=.true."
     endif
  endif

  if(sinkprops.and.finestep_AGN)then
     if(myid==1)then
        write(*,*)"sinkprops cannot be used with AGN feedback on the fine timestep"
        write(*,*)"No sinkprop files will be produced"
     endif
     sinkprops=.false.
  endif

  if(maximum_accretion)then
     if(.not.bondi)then
        if(myid==1)then
           write(*,*)"Maximum accretion is built on the Bondi formalism. Setting bondi=.true."
        endif
        bondi=.true.
     endif
     if(sink_AGN)then
        if(myid==1)then
           write(*,*)"WARNING: Maximum accretion cannot be used with AGN feedback."
           write(*,*)"Switch to bondi accretion by setting maximum_accretion=.false. and bondi=.true."
           write(*,*)"Or switch off AGN feedback using sink_AGN=.false."
        endif
        call clean_stop
     endif
  endif


  ! Last variable is ichem

#ifdef SOLVERmhd
  !-----------------------------------
  ! Set magnetic slope limiters
  !-----------------------------------
  if (slope_mag_type == -1) then
    slope_mag_type = slope_type
  endif
  if (interpol_mag_type == -1) then
    interpol_mag_type = interpol_type
  endif
#endif


end subroutine read_hydro_params
