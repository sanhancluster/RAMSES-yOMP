!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godunov_fine(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine is a wrapper to the second order Godunov solver.
  ! Small grids (2x2x2) are gathered from level ilevel and sent to the
  ! hydro solver. On entry, hydro variables are gathered from array uold.
  ! On exit, unew has been updated.
  !--------------------------------------------------------------------------
  integer::i,igrid,ncache,ngrid
  integer,dimension(1:nvector)::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(static)return
  if(verbose)write(*,111)ilevel

  ! Loop over active grids by vector sweeps
  ncache=active(ilevel)%ngrid
!$omp parallel do private(ngrid,ind_grid)
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call godfine1(ind_grid,ngrid,ilevel)
  end do

111 format('   Entering godunov_fine for level ',i2)

end subroutine godunov_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine check_uold_unew(ilevel,check_mynumber)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine checks the values of uold/unew
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip,check_mynumber
  integer::ii,ich,icell,ilow,ihigh
  real(dp)::d,u,v,w,e
  real(dp)::Zsolar,mdustC,mdustSil,mdust,Zd
  real(dp),dimension(1:nchem)::Zchem
#if NENER>0
  integer::irad
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
     do i=1,active(ilevel)%ngrid

        if(dust_chem)then
           icell=active(ilevel)%igrid(i)+iskip

           ! For uold
           Zsolar=uold(icell,imetal)/uold(icell,1)/0.02
           do ich=1,nchem
              Zchem(ich)=uold(icell,ichem+ich-1)/uold(icell,1)/0.02
           enddo
           if(Zsolar<0.0d0)then
              write(*,*)'in check uold Ztot:',Zsolar,check_mynumber
           endif
           do ich=1,nchem
              if(Zchem(ich)<0.0d0)then
                 write(*,'(A,I2,A,2I9,es13.5,2I9)')'in check uold Ztotchem(',ich,'):',icell,myid,Zchem(ich),check_mynumber,ilevel
              endif
           enddo
           mdustC=0.0d0;mdustSil=0.0d0
           ilow=idust;ihigh=ilow+dndsize
           mdustC  =SUM(uold(icell,ilow:ihigh))
           ilow=ihigh+1;ihigh=ilow+dndsize
           mdustSil=SUM(uold(icell,ilow:ihigh))/SioverSil
!!$           ndchemtype=ndust/2
!!$           mdustC=0.0d0;mdustSil=0.0d0
!!$           do ii=1,ndchemtype
!!$              mdustC  =mdustC  +uold(icell,idust-1+ii)
!!$           enddo
!!$           do ii=ndchemtype+1,ndust
!!$              mdustSil=mdustSil+uold(icell,idust-1+ii)/SioverSil
!!$           enddo
           mdust=mdustC+mdustSil
           Zsolar=Zsolar-mdust/uold(icell,1)/0.02 !! gas met
           do ich=1,nchem
              if(TRIM(chem_list(ich))=='C' )Zchem(ich)=Zchem(ich)-mdustC/uold(icell,1)/0.02
              if(TRIM(chem_list(ich))=='Mg')Zchem(ich)=Zchem(ich)-mdustSil*MgoverSil/uold(icell,1)/0.02
              if(TRIM(chem_list(ich))=='Fe')Zchem(ich)=Zchem(ich)-mdustSil*FeoverSil/uold(icell,1)/0.02
              if(TRIM(chem_list(ich))=='Si')Zchem(ich)=Zchem(ich)-mdustSil*SioverSil/uold(icell,1)/0.02
              if(TRIM(chem_list(ich))=='O' )Zchem(ich)=Zchem(ich)-mdustSil* OoverSil/uold(icell,1)/0.02
           enddo
           if(Zsolar<0.0d0)then
              write(*,'(A,2I9,es13.5,2I9)')'in check uold Zgas:',icell,myid,Zsolar,check_mynumber,ilevel
           endif
           do ich=1,nchem
              if(Zchem(ich)<0.0d0)then
                 if(TRIM(chem_list(ich))=='C' )Zd=mdustC/uold(icell,1)/0.02
                 if(TRIM(chem_list(ich))=='Mg')Zd=mdustSil*MgoverSil/uold(icell,1)/0.02
                 if(TRIM(chem_list(ich))=='Fe')Zd=mdustSil*FeoverSil/uold(icell,1)/0.02
                 if(TRIM(chem_list(ich))=='Si')Zd=mdustSil*SioverSil/uold(icell,1)/0.02
                 if(TRIM(chem_list(ich))=='O' )Zd=mdustSil* OoverSil/uold(icell,1)/0.02
                 write(*,'(A,I2,A,2I9,2es13.5,2I9)')'in check uold Zgaschem(',ich,'):',icell,myid,Zchem(ich),Zd,check_mynumber,ilevel
              endif
           enddo

           ! For unew
           Zsolar=unew(icell,imetal)/unew(icell,1)/0.02
           do ich=1,nchem
              Zchem(ich)=unew(icell,ichem+ich-1)/unew(icell,1)/0.02
           enddo
           if(Zsolar<0.0d0)then
              write(*,'(A,2I9,es13.5,2I9)')'in check unew Ztot:',icell,myid,Zsolar,check_mynumber,ilevel
           endif
           do ich=1,nchem
              if(Zchem(ich)<0.0d0)then
                 write(*,'(A,I2,A,2I9,es13.5,2I9)')'in check unew Ztotchem(',ich,'):',icell,myid,Zchem(ich),check_mynumber,ilevel
              endif
           enddo
           mdustC=0.0d0;mdustSil=0.0d0
           ilow=idust;ihigh=ilow+dndsize
           mdustC  =SUM(unew(icell,ilow:ihigh))
           ilow=ihigh+1;ihigh=ilow+dndsize
           mdustSil=SUM(unew(icell,ilow:ihigh))/SioverSil
!!$           ndchemtype=ndust/2
!!$           mdustC=0.0d0;mdustSil=0.0d0
!!$           do ii=1,ndchemtype
!!$              mdustC  =mdustC  +unew(icell,idust-1+ii)
!!$           enddo
!!$           do ii=ndchemtype+1,ndust
!!$              mdustSil=mdustSil+unew(icell,idust-1+ii)/SioverSil
!!$           enddo
           mdust=mdustC+mdustSil
           Zsolar=Zsolar-mdust/unew(icell,1)/0.02 !! gas met
           do ich=1,nchem
              if(TRIM(chem_list(ich))=='C' )Zchem(ich)=Zchem(ich)-mdustC/unew(icell,1)/0.02
              if(TRIM(chem_list(ich))=='Mg')Zchem(ich)=Zchem(ich)-mdustSil*MgoverSil/unew(icell,1)/0.02
              if(TRIM(chem_list(ich))=='Fe')Zchem(ich)=Zchem(ich)-mdustSil*FeoverSil/unew(icell,1)/0.02
              if(TRIM(chem_list(ich))=='Si')Zchem(ich)=Zchem(ich)-mdustSil*SioverSil/unew(icell,1)/0.02
              if(TRIM(chem_list(ich))=='O' )Zchem(ich)=Zchem(ich)-mdustSil* OoverSil/unew(icell,1)/0.02
           enddo
           if(Zsolar<0.0d0)then
              write(*,'(A,2I9,es13.5,2I9)')'in check unew Zgas:',icell,myid,Zsolar,check_mynumber,ilevel
           endif
           do ich=1,nchem
              if(Zchem(ich)<0.0d0)then
                 if(TRIM(chem_list(ich))=='C' )Zd=mdustC/unew(icell,1)/0.02
                 if(TRIM(chem_list(ich))=='Mg')Zd=mdustSil*MgoverSil/unew(icell,1)/0.02
                 if(TRIM(chem_list(ich))=='Fe')Zd=mdustSil*FeoverSil/unew(icell,1)/0.02
                 if(TRIM(chem_list(ich))=='Si')Zd=mdustSil*SioverSil/unew(icell,1)/0.02
                 if(TRIM(chem_list(ich))=='O' )Zd=mdustSil* OoverSil/unew(icell,1)/0.02
                 write(*,'(A,I2,A,2I9,2es13.5,2I9)')'in check unew Zgaschem(',ich,'):',icell,myid,Zchem(ich),Zd,check_mynumber,ilevel
              endif
           enddo
!!$           if(icell==6181.and.myid==26.and.ilevel==14)then
!!$              ich=5
!!$              if(TRIM(chem_list(ich))=='C' )Zd=mdustC/unew(icell,1)/0.02
!!$           endif
!!$           if(icell==1005155.and.myid==37.and.ilevel==14)then
!!$              ich=3
!!$              if(TRIM(chem_list(ich))=='Fe')Zd=mdustSil*FeoverSil/unew(icell,1)/0.02
!!$              write(*,'(A,I2,A,2I9,2es13.5,2I9)')'in check unew Zgaschem(',ich,'):',icell,myid,Zchem(ich),Zd,check_mynumber,ilevel
!!$           endif
        endif


        if(unew(icell,imetal)<0.0d0)then
           write(*,'(A,2I9,e18.9,2I9)')'check unew(metal) in active CPU:',icell,myid,unew(icell,imetal),check_mynumber,ilevel
        endif
        if(uold(icell,imetal)<0.0d0)then
           write(*,'(A,2I9,e18.9,2I9)')'check uold(metal) in active CPU:',icell,myid,uold(icell,imetal),check_mynumber,ilevel
        endif
        if(unew(icell,idust)<0.0d0)then
           write(*,'(A,2I9,e18.9,2I9)')'check unew in active CPU:',icell,myid,unew(icell,idust),check_mynumber,ilevel
        endif
        if(uold(icell,idust)<0.0d0)then
           write(*,'(A,2I9,e18.9,2I9)')'check uold in active CPU:',icell,myid,uold(icell,idust),check_mynumber,ilevel
        endif

        if(unew(icell,imetal)<unew(icell,idust))then
           write(*,'(A,2I9,2e18.9,2I9)')'check unew(metal) in active CPU:',icell,myid,unew(icell,imetal),unew(icell,idust),check_mynumber,ilevel
        endif
        if(uold(icell,imetal)<uold(icell,idust))then
           write(*,'(A,2I9,2e18.9,2I9)')'check uold(metal) in active CPU:',icell,myid,uold(icell,imetal),uold(icell,idust),check_mynumber,ilevel
        endif

!!$        if(icell==1002142)then
!!$           write(*,'(A,2I9,4e18.9)')'***:',myid,icell,unew(icell,imetal),unew(icell,idust),uold(icell,imetal),uold(icell,idust)
!!$        endif

     end do
  end do

!!$  do icpu=1,ncpu
!!$  do ind=1,twotondim
!!$     iskip=ncoarse+(ind-1)*ngridmax
!!$     do i=1,reception(icpu,ilevel)%ngrid
!!$        if(unew(reception(icpu,ilevel)%igrid(i)+iskip,imetal)<0.0d0)then
!!$           write(*,*)'check unew(metal) in virtual cells:',reception(icpu,ilevel)%igrid(i)+iskip,unew(reception(icpu,ilevel)%igrid(i)+iskip,imetal)
!!$        endif
!!$        if(uold(reception(icpu,ilevel)%igrid(i)+iskip,imetal)<0.0d0)then
!!$           write(*,*)'check uold(metal) in virtual cells:',reception(icpu,ilevel)%igrid(i)+iskip,uold(reception(icpu,ilevel)%igrid(i)+iskip,imetal)
!!$        endif
!!$        if(unew(reception(icpu,ilevel)%igrid(i)+iskip,idust)<0.0d0)then
!!$           write(*,*)'check unew in virtual cells:',reception(icpu,ilevel)%igrid(i)+iskip,unew(reception(icpu,ilevel)%igrid(i)+iskip,idust)
!!$        endif
!!$        if(uold(reception(icpu,ilevel)%igrid(i)+iskip,idust)<0.0d0)then
!!$           write(*,*)'check uold in virtual cells:',reception(icpu,ilevel)%igrid(i)+iskip,uold(reception(icpu,ilevel)%igrid(i)+iskip,idust)
!!$        endif
!!$     end do
!!$  end do
!!$  end do

111 format('   Entering set_unew for level ',i2)

end subroutine check_uold_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_unew(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine sets array unew to its initial value uold before calling
  ! the hydro scheme. unew is set to zero in virtual boundaries.
  !--------------------------------------------------------------------------
  integer::i,ivar,ind,icpu,iskip
  real(dp)::d,u,v,w,e
#if NENER>0
  integer::irad
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set unew to uold for myid cells
!$omp parallel private(iskip)
  do ind=1,twotondim
	 iskip=ncoarse+(ind-1)*ngridmax
!$omp do private(d,u,v,w,e) schedule(static)
     do i=1,active(ilevel)%ngrid
        do ivar=1,nvar
           unew(active(ilevel)%igrid(i)+iskip,ivar) = uold(active(ilevel)%igrid(i)+iskip,ivar)
        end do
        if(momentum_feedback)then
              pstarnew(active(ilevel)%igrid(i)+iskip) = 0.0
        endif
        if(pressure_fix)then
           divu(active(ilevel)%igrid(i)+iskip) = 0.0

           d=max(uold(active(ilevel)%igrid(i)+iskip,1),smallr)
           u=0.0; v=0.0; w=0.0
           if(ndim>0)u=uold(active(ilevel)%igrid(i)+iskip,2)/d
           if(ndim>1)v=uold(active(ilevel)%igrid(i)+iskip,3)/d
           if(ndim>2)w=uold(active(ilevel)%igrid(i)+iskip,4)/d
           e=uold(active(ilevel)%igrid(i)+iskip,ndim+2)-0.5*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e=e-uold(active(ilevel)%igrid(i)+iskip,ndim+2+irad)
           end do
#endif
           enew(active(ilevel)%igrid(i)+iskip)=e
        end if
     end do
!$omp end do nowait
  end do
!$omp end parallel
  ! Set unew to 0 for virtual boundary cells
!$omp parallel do private(iskip,d,u,v,w,e) schedule(static)
  do ind=1,twotondim
     do icpu=1,ncpu
        iskip=ncoarse+(ind-1)*ngridmax
        do ivar=1,nvar
           do i=1,reception(icpu,ilevel)%ngrid
              unew(reception(icpu,ilevel)%igrid(i)+iskip,ivar)=0.0
           end do
        end do
        if(momentum_feedback)then
           do i=1,reception(icpu,ilevel)%ngrid
              pstarnew(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
           end do
        endif
        if(pressure_fix)then
           do i=1,reception(icpu,ilevel)%ngrid
              divu(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
              enew(reception(icpu,ilevel)%igrid(i)+iskip) = 0.0
           end do
        end if
     end do
  end do

111 format('   Entering set_unew for level ',i2)

end subroutine set_unew
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine set_uold(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine sets array uold to its new value unew
  ! after the hydro step.
  !---------------------------------------------------------
  integer::i,ivar,ind,iskip,nx_loc,ind_cell
  real(dp)::scale,d,u,v,w
  real(dp)::e_kin,e_cons,e_prim,e_trunc,div,dx
#if NENER>0
  integer::irad
#endif

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel*scale

  ! Add gravity source terms to unew
  if(poisson)then
     call add_gravity_source_terms(ilevel)
  end if

  ! Add non conservative pdV terms to unew
  ! for thermal and/or non-thermal energies
  if(pressure_fix.OR.nener>0)then
     call add_pdv_source_terms(ilevel)
  endif

  ! Set uold to unew for myid cells
!$omp parallel private(iskip)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
!$omp do private(ind_cell,d,u,v,w,e_kin,e_cons,e_prim,div,e_trunc) schedule(static)
     do i=1,active(ilevel)%ngrid
        do ivar=1,nvar
              uold(active(ilevel)%igrid(i)+iskip,ivar) = unew(active(ilevel)%igrid(i)+iskip,ivar)
        end do
        if(momentum_feedback)then
              pstarold(active(ilevel)%igrid(i)+iskip) = pstarnew(active(ilevel)%igrid(i)+iskip)
        endif
        if(pressure_fix)then
           ! Correct total energy if internal energy is too small
           ind_cell=active(ilevel)%igrid(i)+iskip
           d=max(uold(ind_cell,1),smallr)
           u=0.0; v=0.0; w=0.0
           if(ndim>0)u=uold(ind_cell,2)/d
           if(ndim>1)v=uold(ind_cell,3)/d
           if(ndim>2)w=uold(ind_cell,4)/d
           e_kin=0.5*d*(u**2+v**2+w**2)
#if NENER>0
           do irad=1,nener
              e_kin=e_kin+uold(ind_cell,ndim+2+irad)
           end do
#endif
           e_cons=uold(ind_cell,ndim+2)-e_kin
           e_prim=enew(ind_cell)
           ! Note: here divu=-div.u*dt
           div=abs(divu(ind_cell))*dx/dtnew(ilevel)
           e_trunc=beta_fix*d*max(div,3.0*hexp*dx)**2
           if(e_cons<e_trunc)then
              uold(ind_cell,ndim+2)=e_prim+e_kin
           end if
        end if
     end do
!$omp end do nowait
  end do
!$omp end parallel

111 format('   Entering set_uold for level ',i2)

end subroutine set_uold
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_gravity_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel
  !--------------------------------------------------------------------------
  ! This routine adds to unew the gravity source terms
  ! with only half a time step. Only the momentum and the
  ! total energy are modified in array unew.
  !--------------------------------------------------------------------------
  integer::i,ind,iskip,ind_cell
  real(dp)::d,u,v,w,e_kin,e_prim,d_old,fact

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Add gravity source term at time t with half time step
!$omp parallel private(iskip)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
!$omp do private(ind_cell,d,u,v,w,e_kin,d_old,e_prim,fact) schedule(static)
     do i=1,active(ilevel)%ngrid
        ind_cell=active(ilevel)%igrid(i)+iskip
        d=max(unew(ind_cell,1),smallr)
        u=0.0; v=0.0; w=0.0
        if(ndim>0)u=unew(ind_cell,2)/d
        if(ndim>1)v=unew(ind_cell,3)/d
        if(ndim>2)w=unew(ind_cell,4)/d
        e_kin=0.5*d*(u**2+v**2+w**2)
        e_prim=unew(ind_cell,ndim+2)-e_kin
        d_old=max(uold(ind_cell,1),smallr)
        fact=d_old/d*0.5*dtnew(ilevel)
        if(ndim>0)then
           u=u+f(ind_cell,1)*fact
           unew(ind_cell,2)=d*u
        endif
        if(ndim>1)then
           v=v+f(ind_cell,2)*fact
           unew(ind_cell,3)=d*v
        end if
        if(ndim>2)then
           w=w+f(ind_cell,3)*fact
           unew(ind_cell,4)=d*w
        endif
        e_kin=0.5*d*(u**2+v**2+w**2)
        unew(ind_cell,ndim+2)=e_prim+e_kin
     end do
!$omp end do nowait
  end do
!$omp end parallel

111 format('   Entering add_gravity_source_terms for level ',i2)

end subroutine add_gravity_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine add_pdv_source_terms(ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp)::scale,dx,dx_loc
  integer::i,ncache,igrid,ngrid,ilevel,nx_loc
  integer ,dimension(1:nvector)::ind_grid
  common /omp_addpdv/ scale,dx,dx_loc,iii,jjj,nx_loc

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5d0**ilevel
  dx_loc=dx*scale

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ncache=active(ilevel)%ngrid
!$omp parallel do private(ngrid,ind_grid) schedule(static)
  do igrid=1,ncache,nvector
     ngrid = MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call addpdv1(ind_grid,ngrid,ilevel)
  enddo
  ! End loop over grids

111 format('   Entering add_pdv_source_terms for level ',i2)

end subroutine add_pdv_source_terms
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine addpdv1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  implicit none
  integer::ilevel
  !---------------------------------------------------------
  ! This routine adds the pdV source term to the internal
  ! energy equation and to the non-thermal energy equations.
  !---------------------------------------------------------
  integer::i,ind,iskip,nx_loc
  integer::ngrid,idim,id1,ig1,ih1,id2,ig2,ih2
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  real(dp)::scale,dx,dx_loc,d,u,v,w,eold

  integer ,dimension(1:nvector)::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim)::igridn
  integer ,dimension(1:nvector,1:ndim)::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim,1:ndim)::velg,veld
  real(dp),dimension(1:nvector,1:ndim)::dx_g,dx_d
  real(dp),dimension(1:nvector)::divu_loc
#if NENER>0
  integer::irad
#endif
  common /omp_addpdv/ scale,dx,dx_loc,iii,jjj,nx_loc

  ! Gather neighboring grids
  do i=1,ngrid
      igridn(i,0)=ind_grid(i)
  end do
  do idim=1,ndim
      do i=1,ngrid
          ind_left (i,idim)=nbor(ind_grid(i),2*idim-1)
          ind_right(i,idim)=nbor(ind_grid(i),2*idim  )
          igridn(i,2*idim-1)=son(ind_left (i,idim))
          igridn(i,2*idim  )=son(ind_right(i,idim))
      end do
  end do

  ! Loop over cells
  do ind=1,twotondim

      ! Compute central cell index
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,ngrid
          ind_cell(i)=iskip+ind_grid(i)
      end do

      ! Gather all neighboring velocities
      do idim=1,ndim
          id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
          ih1=ncoarse+(id1-1)*ngridmax
          do i=1,ngrid
              if(igridn(i,ig1)>0)then
                  velg(i,idim,1:ndim) = uold(igridn(i,ig1)+ih1,2:ndim+1)/max(uold(igridn(i,ig1)+ih1,1),smallr)
                  dx_g(i,idim) = dx_loc
              else
                  velg(i,idim,1:ndim) = uold(ind_left(i,idim),2:ndim+1)/max(uold(ind_left(i,idim),1),smallr)
                  dx_g(i,idim) = dx_loc*1.5_dp
              end if
          enddo
          id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
          ih2=ncoarse+(id2-1)*ngridmax
          do i=1,ngrid
              if(igridn(i,ig2)>0)then
                  veld(i,idim,1:ndim)= uold(igridn(i,ig2)+ih2,2:ndim+1)/max(uold(igridn(i,ig2)+ih2,1),smallr)
                  dx_d(i,idim)=dx_loc
              else
                  veld(i,idim,1:ndim)= uold(ind_right(i,idim),2:ndim+1)/max(uold(ind_right(i,idim),1),smallr)
                  dx_d(i,idim)=dx_loc*1.5_dp
              end if
          enddo
      end do
      ! End loop over dimensions

      ! Compute divu = Trace G
      divu_loc(1:ngrid)=0.0d0
      do i=1,ngrid
          do idim=1,ndim
              divu_loc(i) = divu_loc(i) + (veld(i,idim,idim)-velg(i,idim,idim)) &
                     &                    / (dx_g(i,idim)     +dx_d(i,idim))
          enddo
      end do

      ! Update thermal internal energy
      if(pressure_fix)then
          do i=1,ngrid
              ! Compute old thermal energy
              d=max(uold(ind_cell(i),1),smallr)
              u=0.0; v=0.0; w=0.0
              if(ndim>0)u=uold(ind_cell(i),2)/d
              if(ndim>1)v=uold(ind_cell(i),3)/d
              if(ndim>2)w=uold(ind_cell(i),4)/d
              eold=uold(ind_cell(i),ndim+2)-0.5*d*(u**2+v**2+w**2)
#if NENER>0
              do irad=1,nener
                  eold=eold-uold(ind_cell(i),ndim+2+irad)
              end do
#endif
              ! Add -pdV term
              enew(ind_cell(i))=enew(ind_cell(i)) &
                     & -(gamma-1.0d0)*eold*divu_loc(i)*dtnew(ilevel)
          end do
      end if

#if NENER>0
      do irad=1,nener
          do i=1,ngrid
              ! Add -pdV term
              unew(ind_cell(i),ndim+2+irad)=unew(ind_cell(i),ndim+2+irad) &
                 & -(gamma_rad(irad)-1.0d0)*uold(ind_cell(i),ndim+2+irad)*divu_loc(i)*dtnew(ilevel)
          end do
      end do
#endif

  enddo
  ! End loop over cells
end subroutine addpdv1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine godfine1(ind_grid,ncache,ilevel)
  use hydro_parameters
  use amr_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel,ncache
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  ! This routine gathers first hydro variables from neighboring grids
  ! to set initial conditions in a 6x6x6 grid. It interpolates from
  ! coarser level missing grid variables. It then calls the
  ! Godunov solver that computes fluxes. These fluxes are zeroed at
  ! coarse-fine boundaries, since contribution from finer levels has
  ! already been taken into account. Conservative variables are updated
  ! and stored in array unew(:), both at the current level and at the
  ! coarser level if necessary.
  !-------------------------------------------------------------------
  integer ,dimension(1:nvector,1:threetondim     )::nbors_father_cells
  integer ,dimension(1:nvector,1:twotondim       )::nbors_father_grids
  integer ,dimension(1:nvector,0:twondim         )::ibuffer_father
  real(dp),dimension(1:nvector,0:twondim  ,1:nvar)::u1
  real(dp),dimension(1:nvector,1:twotondim,1:nvar)::u2

  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:nvar)::uloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2,1:ndim)::gloc
  real(dp),dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ploc
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:nvar,1:ndim)::flux
  real(dp),dimension(1:nvector,if1:if2,jf1:jf2,kf1:kf2,1:2,1:ndim)::tmp
  logical ,dimension(1:nvector,iu1:iu2,ju1:ju2,ku1:ku2)::ok
  real(dp),dimension(1:nvector,1:nvar)::uflow
  real(dp),dimension(1:nvector)::dflow,eflow

  integer,dimension(1:nvector)::igrid_nbor,ind_cell,ind_buffer,ind_exist,ind_nexist

  integer::i,j,ivar,idim,ind_son,ind_father,iskip,nbuffer
  integer::i0,j0,k0,i1,j1,k1,i2,j2,k2,i3,j3,k3,nx_loc,nb_noneigh,nexist
  integer::i1min,i1max,j1min,j1max,k1min,k1max
  integer::i2min,i2max,j2min,j2max,k2min,k2max
  integer::i3min,i3max,j3min,j3max,k3min,k3max
  real(dp)::dx,scale,oneontwotondim

  integer::ich,ilow,ihigh
  real(dp)::mdustC,mdustSil

  oneontwotondim = 1.d0/dble(twotondim)

  ! Mesh spacing in that level
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  dx=0.5D0**ilevel*scale

  ploc=0.0d0; gloc=0.0d0

  ! Integer constants
  i1min=0; i1max=0; i2min=0; i2max=0; i3min=1; i3max=1
  j1min=0; j1max=0; j2min=0; j2max=0; j3min=1; j3max=1
  k1min=0; k1max=0; k2min=0; k2max=0; k3min=1; k3max=1
  if(ndim>0)then
     i1max=2; i2max=1; i3max=2
  end if
  if(ndim>1)then
     j1max=2; j2max=1; j3max=2
  end if
  if(ndim>2)then
     k1max=2; k2max=1; k3max=2
  end if

  !------------------------------------------
  ! Gather 3^ndim neighboring father cells
  !------------------------------------------
  do i=1,ncache
     ind_cell(i)=father(ind_grid(i))
  end do
  call get3cubefather(ind_cell,nbors_father_cells,nbors_father_grids,ncache,ilevel)

  !---------------------------
  ! Gather 6x6x6 cells stencil
  !---------------------------
  ! Loop over 3x3x3 neighboring father cells
  do k1=k1min,k1max
  do j1=j1min,j1max
  do i1=i1min,i1max

     ! Check if neighboring grid exists
     nbuffer=0
     nexist=0
     ind_father=1+i1+3*j1+9*k1
     do i=1,ncache
        igrid_nbor(i)=son(nbors_father_cells(i,ind_father))
        if(igrid_nbor(i)>0) then
           nexist=nexist+1
           ind_exist(nexist)=i
        else
          nbuffer=nbuffer+1
          ind_nexist(nbuffer)=i
          ind_buffer(nbuffer)=nbors_father_cells(i,ind_father)
        end if
     end do

     ! If not, interpolate hydro variables from parent cells
     if(nbuffer>0)then
        call getnborfather(ind_buffer,ibuffer_father,nbuffer,ilevel)
        do j=0,twondim
           do ivar=1,nvar
              do i=1,nbuffer
                 u1(i,j,ivar)=uold(ibuffer_father(i,j),ivar)
              end do
           end do

           if(dust_chem)then
              do i=1,nbuffer
              mdustC=0.0d0;mdustSil=0.0d0
              ilow=idust;ihigh=ilow+dndsize
              mdustC  =SUM(u1(i,j,ilow:ihigh))
              ilow=ihigh+1;ihigh=ilow+dndsize
              mdustSil=SUM(u1(i,j,ilow:ihigh))/SioverSil
              do ich=1,nchem
                 if(TRIM(chem_list(ich))=='C' )u1(i,j,ichem+ich-1)=u1(i,j,ichem+ich-1)-mdustC
                 if(TRIM(chem_list(ich))=='Mg')u1(i,j,ichem+ich-1)=u1(i,j,ichem+ich-1)-mdustSil*MgoverSil
                 if(TRIM(chem_list(ich))=='Fe')u1(i,j,ichem+ich-1)=u1(i,j,ichem+ich-1)-mdustSil*FeoverSil
                 if(TRIM(chem_list(ich))=='Si')u1(i,j,ichem+ich-1)=u1(i,j,ichem+ich-1)-mdustSil*SioverSil
                 if(TRIM(chem_list(ich))=='O' )u1(i,j,ichem+ich-1)=u1(i,j,ichem+ich-1)-mdustSil* OoverSil
              enddo
              enddo
           else
              ilow=idust;ihigh=ilow+dndsize
              do i=1,nbuffer
                 u1(i,j,imetal)=u1(i,j,imetal)-SUM(u1(i,j,ilow:ihigh))
              enddo
           endif

        end do
        call interpol_hydro(u1,u2,nbuffer)
     endif

     ! Loop over 2x2x2 cells
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max

        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,nexist
           ind_cell(i)=iskip+igrid_nbor(ind_exist(i))
        end do

        i3=1; j3=1; k3=1
        if(ndim>0)i3=1+2*(i1-1)+i2
        if(ndim>1)j3=1+2*(j1-1)+j2
        if(ndim>2)k3=1+2*(k1-1)+k2

        ! Gather hydro variables
        do ivar=1,nvar
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,ivar)=uold(ind_cell(i),ivar)
           end do
           do i=1,nbuffer
              uloc(ind_nexist(i),i3,j3,k3,ivar)=u2(i,ind_son,ivar)
           end do
        end do

        if(dust_chem)then
           ! only for the cells w/o interpolations
           !(those with interpolations their values have been manipulated on u1, giving u2)
           do i=1,nexist
              mdustC=0.0d0;mdustSil=0.0d0
              ilow=idust;ihigh=ilow+dndsize
              mdustC  =SUM(uloc(ind_exist(i),i3,j3,k3,ilow:ihigh))
              ilow=ihigh+1;ihigh=ilow+dndsize
              mdustSil=SUM(uloc(ind_exist(i),i3,j3,k3,ilow:ihigh))/SioverSil
              do ich=1,nchem
                 if(TRIM(chem_list(ich))=='C' )uloc(ind_exist(i),i3,j3,k3,ichem+ich-1)=uloc(ind_exist(i),i3,j3,k3,ichem+ich-1)-mdustC
                 if(TRIM(chem_list(ich))=='Mg')uloc(ind_exist(i),i3,j3,k3,ichem+ich-1)=uloc(ind_exist(i),i3,j3,k3,ichem+ich-1)-mdustSil*MgoverSil
                 if(TRIM(chem_list(ich))=='Fe')uloc(ind_exist(i),i3,j3,k3,ichem+ich-1)=uloc(ind_exist(i),i3,j3,k3,ichem+ich-1)-mdustSil*FeoverSil
                 if(TRIM(chem_list(ich))=='Si')uloc(ind_exist(i),i3,j3,k3,ichem+ich-1)=uloc(ind_exist(i),i3,j3,k3,ichem+ich-1)-mdustSil*SioverSil
                 if(TRIM(chem_list(ich))=='O' )uloc(ind_exist(i),i3,j3,k3,ichem+ich-1)=uloc(ind_exist(i),i3,j3,k3,ichem+ich-1)-mdustSil* OoverSil
              enddo
           enddo
        else
           ilow=idust;ihigh=ilow+dndsize
           do i=1,nexist
              uloc(ind_exist(i),i3,j3,k3,imetal)=uloc(ind_exist(i),i3,j3,k3,imetal)-SUM(uloc(ind_exist(i),i3,j3,k3,ilow:ihigh))
           enddo
        endif

        ! Gather gravitational acceleration
        if(poisson)then
           do idim=1,ndim
              do i=1,nexist
                 gloc(ind_exist(i),i3,j3,k3,idim)=f(ind_cell(i),idim)
              end do
              ! Use straight injection for buffer cells
              do i=1,nbuffer
                 gloc(ind_nexist(i),i3,j3,k3,idim)=f(ibuffer_father(i,0),idim)
              end do
           end do
        end if
        ! Gather stellar momentum
        if(momentum_feedback)then
           do i=1,nexist
              ploc(ind_exist(i),i3,j3,k3)=pstarold(ind_cell(i))
           end do
           ! Use straight injection for buffer cells
           do i=1,nbuffer
              ploc(ind_nexist(i),i3,j3,k3)=pstarold(ibuffer_father(i,0))
           end do
        end if

        ! Gather refinement flag
        do i=1,nexist
           ok(ind_exist(i),i3,j3,k3)=son(ind_cell(i))>0
        end do
        do i=1,nbuffer
           ok(ind_nexist(i),i3,j3,k3)=.false.
        end do

     end do
     end do
     end do
     ! End loop over cells

  end do
  end do
  end do

  ! End loop over neighboring grids
  !-----------------------------------------------
  ! Compute flux using second-order Godunov method
  !-----------------------------------------------
  call unsplit(uloc,gloc,ploc,flux,tmp,dx,dx,dx,dtnew(ilevel),ncache)

  !--------------------------------------
  ! Store the fluxes for later use
  !--------------------------------------
  if (MC_tracer) then
     do idim=1,ndim
        i0=0; j0=0; k0=0
        if(idim==1)i0=1
        if(idim==2)j0=1
        if(idim==3)k0=1
        do k2=k2min,k2max
           do j2=j2min,j2max
              do i2=i2min,i2max
                 ind_son=1+i2+2*j2+4*k2
                 iskip=ncoarse+(ind_son-1)*ngridmax
                 do i=1,ncache
                    ind_cell(i)=iskip+ind_grid(i)
                 end do
                 i3=1+i2
                 j3=1+j2
                 k3=1+k2
                 do i=1,ncache
                    ! Copy left flux
                    fluxes(ind_cell(i),(idim-1)*2+1)= flux(i,i3   ,j3   ,k3,   1,idim)&
                         / uold(ind_cell(i), 1)
                    ! Copy right flux
                    fluxes(ind_cell(i),(idim-1)*2+2)=-flux(i,i3+i0,j3+j0,k3+k0,1,idim)&
                         / uold(ind_cell(i), 1)
                 end do
              end do
           end do
        end do
     end do
  end if

  !------------------------------------------------
  ! Reset flux along direction at refined interface
  !------------------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     ! loop equivalent to twotondim
     do k3=k3min,k3max+k0
     do j3=j3min,j3max+j0
     do i3=i3min,i3max+i0
        do ivar=1,nvar
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 flux(i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        if(pressure_fix)then
        do ivar=1,2
           do i=1,ncache
              if(ok(i,i3-i0,j3-j0,k3-k0) .or. ok(i,i3,j3,k3))then
                 tmp (i,i3,j3,k3,ivar,idim)=0.0d0
              end if
           end do
        end do
        end if
     end do
     end do
     end do
  end do

  !! OMP note: minimized the overhead from atomic statements
  !--------------------------------------
  ! Conservative update at level ilevel
  !--------------------------------------
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1
     ! loop equivalent to twotondim
     do k2=k2min,k2max
     do j2=j2min,j2max
     do i2=i2min,i2max
        ind_son=1+i2+2*j2+4*k2
        iskip=ncoarse+(ind_son-1)*ngridmax
        do i=1,ncache
           ind_cell(i)=iskip+ind_grid(i)
        end do
        i3=1+i2
        j3=1+j2
        k3=1+k2
        ! Update conservative variables new state vector
        do ivar=1,nvar
           do i=1,ncache
              unew(ind_cell(i),ivar)=unew(ind_cell(i),ivar)+ &
                   & (flux(i,i3   ,j3   ,k3   ,ivar,idim) &
                   & -flux(i,i3+i0,j3+j0,k3+k0,ivar,idim))
           end do
        end do
        if(dust_chem)then
           do i=1,ncache
              do ich=1,nchem
                 if(TRIM(chem_list(ich))=='C' )then
                    ilow=idust;ihigh=ilow+dndsize
                    unew(ind_cell(i),ichem+ich-1)=unew(ind_cell(i),ichem+ich-1)+ &
                         & (SUM(flux(i,i3   ,j3   ,k3   ,ilow:ihigh,idim)) &
                         & -SUM(flux(i,i3+i0,j3+j0,k3+k0,ilow:ihigh,idim)))
                 endif
                 if(TRIM(chem_list(ich))=='Mg')then
                    ilow=idust+dndsize+1;ihigh=ilow+dndsize
                    unew(ind_cell(i),ichem+ich-1)=unew(ind_cell(i),ichem+ich-1)+ &
                         & (SUM(flux(i,i3   ,j3   ,k3   ,ilow:ihigh,idim)) &
                         & -SUM(flux(i,i3+i0,j3+j0,k3+k0,ilow:ihigh,idim)))*MgoverSil/SioverSil
                 endif
                 if(TRIM(chem_list(ich))=='Fe')then
                    ilow=idust+dndsize+1;ihigh=ilow+dndsize
                    unew(ind_cell(i),ichem+ich-1)=unew(ind_cell(i),ichem+ich-1)+ &
                         & (SUM(flux(i,i3   ,j3   ,k3   ,ilow:ihigh,idim)) &
                         & -SUM(flux(i,i3+i0,j3+j0,k3+k0,ilow:ihigh,idim)))*FeoverSil/SioverSil
                 endif
                 if(TRIM(chem_list(ich))=='Si')then
                    ilow=idust+dndsize+1;ihigh=ilow+dndsize
                    unew(ind_cell(i),ichem+ich-1)=unew(ind_cell(i),ichem+ich-1)+ &
                         & (SUM(flux(i,i3   ,j3   ,k3   ,ilow:ihigh,idim)) &
                         & -SUM(flux(i,i3+i0,j3+j0,k3+k0,ilow:ihigh,idim)))
                 endif
                 if(TRIM(chem_list(ich))=='O' )then
                    ilow=idust+dndsize+1;ihigh=ilow+dndsize
                    unew(ind_cell(i),ichem+ich-1)=unew(ind_cell(i),ichem+ich-1)+ &
                         & (SUM(flux(i,i3   ,j3   ,k3   ,ilow:ihigh,idim)) &
                         & -SUM(flux(i,i3+i0,j3+j0,k3+k0,ilow:ihigh,idim)))* OoverSil/SioverSil
                 endif
              enddo
           enddo
        else
           ilow=idust;ihigh=ilow+dndsize
           do i=1,ncache
              unew(ind_cell(i),imetal)=unew(ind_cell(i),imetal)+ &
                   & (SUM(flux(i,i3   ,j3   ,k3   ,ilow:ihigh,idim)) &
                   & -SUM(flux(i,i3+i0,j3+j0,k3+k0,ilow:ihigh,idim)))
           enddo
        endif

        if(pressure_fix) then
           ! Update velocity divergence
           do i=1,ncache
              divu(ind_cell(i))=divu(ind_cell(i))+ &
                   & (tmp(i,i3   ,j3   ,k3   ,1,idim) &
                   & -tmp(i,i3+i0,j3+j0,k3+k0,1,idim))
           end do
           ! Update internal energy
           do i=1,ncache
              enew(ind_cell(i))=enew(ind_cell(i))+ &
                   & (tmp(i,i3   ,j3   ,k3   ,2,idim) &
                   & -tmp(i,i3+i0,j3+j0,k3+k0,2,idim))
           end do
        end if
     end do
     end do
     end do
  end do

  !--------------------------------------
  ! Conservative update at level ilevel-1
  !--------------------------------------
  ! Loop over dimensions
  do idim=1,ndim
     i0=0; j0=0; k0=0
     if(idim==1)i0=1
     if(idim==2)j0=1
     if(idim==3)k0=1

     !----------------------
     ! Left flux at boundary
     !----------------------
     ! Check if grids sits near left boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim-1))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim-1)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     uflow=0d0; dflow=0d0; eflow=0d0
     ! Loop over boundary cells
     do k3=k3min,k3max-k0
     do j3=j3min,j3max-j0
     do i3=i3min,i3max-i0
        do i=1,nb_noneigh
           do ivar=1,nvar
              uflow(i,ivar)=uflow(i,ivar) &
                    & -flux(ind_cell(i),i3,j3,k3,ivar,idim)
           end do
           if(pressure_fix) then
              ! Update velocity divergence
              dflow(i)=dflow(i)-tmp(ind_cell(i),i3,j3,k3,1,idim)
              ! Update internal energy
              eflow(i)=eflow(i)-tmp(ind_cell(i),i3,j3,k3,2,idim)
           end if
        end do
     end do
     end do
     end do

     if(dust_chem)then
        do i=1,nb_noneigh
           do ich=1,nchem
              if(TRIM(chem_list(ich))=='C' )then
                 ilow=idust;ihigh=ilow+dndsize
                 uflow(i,ichem+ich-1)=uflow(i,ichem+ich-1)+SUM(uflow(i,ilow:ihigh))
              endif
              if(TRIM(chem_list(ich))=='Mg')then
                 ilow=idust+dndsize+1;ihigh=ilow+dndsize
                 uflow(i,ichem+ich-1)=uflow(i,ichem+ich-1)+SUM(uflow(i,ilow:ihigh))*MgoverSil/SioverSil
              endif
              if(TRIM(chem_list(ich))=='Fe')then
                 ilow=idust+dndsize+1;ihigh=ilow+dndsize
                 uflow(i,ichem+ich-1)=uflow(i,ichem+ich-1)+SUM(uflow(i,ilow:ihigh))*FeoverSil/SioverSil
              endif
              if(TRIM(chem_list(ich))=='Si')then
                 ilow=idust+dndsize+1;ihigh=ilow+dndsize
                 uflow(i,ichem+ich-1)=uflow(i,ichem+ich-1)+SUM(uflow(i,ilow:ihigh))
              endif
              if(TRIM(chem_list(ich))=='O' )then
                 ilow=idust+dndsize+1;ihigh=ilow+dndsize
                 uflow(i,ichem+ich-1)=uflow(i,ichem+ich-1)+SUM(uflow(i,ilow:ihigh))*OoverSil/SioverSil
              endif
           enddo
        enddo
     else
        ilow=idust;ihigh=ilow+dndsize
        do i=1,nb_noneigh
           uflow(i,imetal)=uflow(i,imetal)+ &
                & SUM(uflow(i,ilow:ihigh))
        enddo
     endif

     ! Update common arrays
     do i=1,nb_noneigh
        do ivar=1,nvar
!$omp atomic update
        unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar)+uflow(i,ivar)*oneontwotondim
        end do
     end do

     if(pressure_fix) then
        do i=1,nb_noneigh
!$omp atomic update
            divu(ind_buffer(i))=divu(ind_buffer(i))+dflow(i)*oneontwotondim
        end do
        do i=1,nb_noneigh
!$omp atomic update
            enew(ind_buffer(i))=enew(ind_buffer(i))+eflow(i)*oneontwotondim
        end do
     end if

     !-----------------------
     ! Right flux at boundary
     !-----------------------
     ! Check if grids sits near right boundary
     ! and gather neighbor father cells index
     nb_noneigh=0
     do i=1,ncache
        if (son(nbor(ind_grid(i),2*idim))==0) then
           nb_noneigh = nb_noneigh + 1
           ind_buffer(nb_noneigh) = nbor(ind_grid(i),2*idim)
           ind_cell(nb_noneigh) = i
        end if
     end do
     ! Conservative update of new state variables
     uflow=0d0; dflow=0d0; eflow=0d0
     ! Loop over boundary cells
     do k3=k3min+k0,k3max
     do j3=j3min+j0,j3max
     do i3=i3min+i0,i3max
        do i=1,nb_noneigh
           do ivar=1,nvar
              uflow(i,ivar)=uflow(i,ivar) &
                    & +flux(ind_cell(i),i3+i0,j3+j0,k3+k0,ivar,idim)
           end do
           if(pressure_fix) then
              ! Update velocity divergence
              dflow(i)=dflow(i)+tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,1,idim)
              ! Update internal energy
              eflow(i)=eflow(i)+tmp(ind_cell(i),i3+i0,j3+j0,k3+k0,2,idim)
           end if
        end do
     end do
     end do
     end do

     if(dust_chem)then
        do i=1,nb_noneigh
           do ich=1,nchem
              if(TRIM(chem_list(ich))=='C' )then
                 ilow=idust;ihigh=ilow+dndsize
                 uflow(i,ichem+ich-1)=uflow(i,ichem+ich-1)+SUM(uflow(i,ilow:ihigh))
              endif
              if(TRIM(chem_list(ich))=='Mg')then
                 ilow=idust+dndsize+1;ihigh=ilow+dndsize
                 uflow(i,ichem+ich-1)=uflow(i,ichem+ich-1)+SUM(uflow(i,ilow:ihigh))*MgoverSil/SioverSil
              endif
              if(TRIM(chem_list(ich))=='Fe')then
                 ilow=idust+dndsize+1;ihigh=ilow+dndsize
                 uflow(i,ichem+ich-1)=uflow(i,ichem+ich-1)+SUM(uflow(i,ilow:ihigh))*FeoverSil/SioverSil
              endif
              if(TRIM(chem_list(ich))=='Si')then
                 ilow=idust+dndsize+1;ihigh=ilow+dndsize
                 uflow(i,ichem+ich-1)=uflow(i,ichem+ich-1)+SUM(uflow(i,ilow:ihigh))
              endif
              if(TRIM(chem_list(ich))=='O' )then
                 ilow=idust+dndsize+1;ihigh=ilow+dndsize
                 uflow(i,ichem+ich-1)=uflow(i,ichem+ich-1)+SUM(uflow(i,ilow:ihigh))*OoverSil/SioverSil
              endif
           enddo
        enddo
     else
        ilow=idust;ihigh=ilow+dndsize
        do i=1,nb_noneigh
           uflow(i,imetal)=uflow(i,imetal)+ &
                & SUM(uflow(i,ilow:ihigh))
        enddo
     endif

     ! Update common arrays
     do i=1,nb_noneigh
        do ivar=1,nvar
!$omp atomic update
         unew(ind_buffer(i),ivar)=unew(ind_buffer(i),ivar)+uflow(i,ivar)*oneontwotondim
        end do
     end do

     if(pressure_fix) then
        do i=1,nb_noneigh
!$omp atomic update
           divu(ind_buffer(i))=divu(ind_buffer(i))+dflow(i)*oneontwotondim
        end do
        do i=1,nb_noneigh
!$omp atomic update
           enew(ind_buffer(i))=enew(ind_buffer(i))+eflow(i)*oneontwotondim
        end do
     end if
  end do
  ! End loop over dimensions
end subroutine godfine1
