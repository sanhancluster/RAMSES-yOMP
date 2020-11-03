!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine phi_fine_cg(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  use mpi_mod
  implicit none
#ifndef WITHOUTMPI
  integer::info
#endif
  integer::ilevel,icount
  !=========================================================
  ! Iterative Poisson solver with Conjugate Gradient method
  ! to solve A x = b.
  !
  ! The algorithm of the routine is inspired by the two papers
  !
  !   http://www.sciencedirect.com/science/article/pii/S0167819113000719 (layout and pipelining)
  !
  !   http://www.vis.uni-stuttgart.de/~weiskopf/publications/pdp10.pdf (2D preconditioner generalized to 3D in code)
  !
  ! r     : stored in f(i, 1)
  ! u     : stored in f(i, 2)
  ! w     : stored in f(i, 3)
  ! m     : stored in f(i, 4)
  ! n     : stored in f(i, 5)
  ! z     : stored in f(i, 6)
  ! q     : stored in f(i, 7)
  ! s     : stored in f(i, 8)
  ! p     : stored in f(i, 9)
  ! x     : stored in phi(i)
  ! b     : stored in rho(i)
  !
  ! Initial guess for phi: interpolated phi from ilevel-1
  !
  ! Initially :
  !   r = b - A x
  !   u = M^-1 r
  !   w = A u
  !
  ! while (iter<cg_itermax .and. error > error_ini * epsilon)
  !
  !   iter = iter + 1
  !   gamma = (r.u)
  !   delta = (w.u)
  !
  !   error = sqrt(gamma / Ncell)
  !   if (iter=1) error_ini = error
  !   if (iter>1)
  !     alpha = gamma / (delta - beta * gamma / alpha_old)
  !     beta  = gamma / gamma_old
  !   else
  !     alpha = gamma / delta
  !     beta = 0
  !   endif
  !   alpha_old = alpha
  !   gamma_old = gamma
  !
  !   m = M^-1 w
  !   n = A m
  !
  !   z = n + beta z
  !   q = m + beta q
  !   s = w + beta s
  !   p = u + beta p
  !
  !   phi = phi + alpha p
  !   r   = r   - alpha s
  !   u   = u   - alpha q
  !   w   = w   - alpha z
  !
  !  end while
  !
  !=========================================================
  integer::i,ind,iter,iskip,itermax,nx_loc,icpu,off,off2
  integer::idx,addr,nact,ntot,ncell,ncache,idim,ig,ih,j,k
  integer::countsend,countrecv
  integer,dimension(ncpu)::reqsend,reqrecv
  real(dp)::error,error_ini
  real(dp)::dx2,fourpi,scale,oneoversix,fact,fact2,prefac
  real(dp)::alpha_cg,alpha_cg_old,beta_cg,gamma_cg,gamma_cg_old,delta_cg
  real(dp), dimension(1:2) :: local,global
  real(kind=8)::rhs_norm,residu
#ifndef WITHOUTMPI
  real(kind=8) :: rhs_norm_all
#endif

  if(gravity_type>0)return
  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Set constants
  dx2=(0.5D0**ilevel)**2
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  fourpi=4.D0*ACOS(-1.0D0)*scale
  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale
  oneoversix=1.0D0/dble(twondim)
  fact=oneoversix*fourpi*dx2
  fact2 = fact*fact
  prefac = 1d0+0.5d0*oneoversix

  !===============================
  ! Compute initial phi
  !===============================
   if(ilevel>levelmin)then
      call make_initial_phi(ilevel,icount)              ! Interpolate phi down
   else
      call make_multipole_phi(ilevel)            ! Fill up with simple initial guess
   endif
   call make_virtual_fine_dp(phi(1),ilevel)      ! Update boundaries
   call make_boundary_phi(ilevel)                ! Update physical boundaries

  !===============================
  ! Initialize f(:,:) array for active and reception indices
  !===============================
!$omp parallel private(iskip,ncache,idx)
   do ind=1,twotondim
      iskip=ncoarse+(ind-1)*ngridmax
      do icpu=1,ncpu
         if(icpu==myid) then
            ncache=active(ilevel)%ngrid
         else
            ncache=reception(icpu,ilevel)%ngrid
         end if
!$omp do
         do i=1,ncache
            if(icpu==myid) then
               idx=active(ilevel)%igrid(i)+iskip
            else
               idx=reception(icpu,ilevel)%igrid(i)+iskip
            end if
            do j=1,10
               f(idx,j)=0d0
            end do
         end do
!$omp end do nowait
      end do
   end do
!$omp end parallel

  !===============================
  ! Compute right-hand side norm (for error)
  !===============================
  rhs_norm=0d0
!$omp parallel private(iskip,idx) reduction(+:rhs_norm)
  do ind=1,twotondim
     iskip=ncoarse+(ind-1)*ngridmax
!$omp do
	 do i=1,active(ilevel)%ngrid
		idx=active(ilevel)%igrid(i)+iskip
		rhs_norm=rhs_norm+fact2*(rho(idx)-rho_tot)*(rho(idx)-rho_tot)
     end do
!$omp end do nowait
  end do
!$omp end parallel

  ! Compute global norms
#ifndef WITHOUTMPI
  call MPI_ALLREDUCE(rhs_norm,rhs_norm_all,1,MPI_DOUBLE_PRECISION,MPI_SUM,&
       & MPI_COMM_WORLD,info)
  rhs_norm=rhs_norm_all
#endif
  rhs_norm=DSQRT(rhs_norm/dble(twotondim*numbtot(1,ilevel)))

  !==============================================
  ! Compute r = b - Ax and store it into f(i,1) (apply linear addressing)
  ! Interpolate down from paraent grid if no neighbor grid available
  !==============================================
  call cmp_residual_cg(ilevel,icount)

  !==============================================
  ! Update boundaries for r
  !==============================================
  call make_virtual_fine_dp(f(1,1),ilevel)

  !==============================================
  ! Set u = Minv r and store it into f(i, 2)
  !==============================================
  call cmp_Minv_cg(ilevel,1,2)

  !==============================================
  ! Set w = Au and store it into f(i, 3)
  !==============================================
   call cmp_A_cg(ilevel,2,3)

  !==============================================
  ! Post receive of ghostzones for w for first iteration
  !==============================================
  call recv_virtual(ilevel, countrecv, reqrecv)

  !====================================
  ! Main iteration loop
  !====================================
  iter=0; itermax=10000
  gamma_cg=0.0; delta_cg=0.0
  error=1.0D0; error_ini=1.0D0
  !! Main bottleneck

  do while(error>epsilon*error_ini.and.iter<itermax)
#ifndef WITHOUTMPI
     if (iter>0) then
        !==============================================
        ! Post update of ghostzones for w
        !==============================================
        ! Wait for full completion of sends from last iteration before filling emission array again
        if (countsend>0) call MPI_WAITALL(countsend,reqsend,MPI_STATUSES_IGNORE,info)
     end if

     !==============================================
     ! Gather and send emission array for w
     !==============================================
     call send_virtual(f(1,3), ilevel, countsend, reqsend)
#endif

     !====================================
     ! Compute dot products. gamma_cg = r.u, delta_cg = w.u
     !====================================
!$omp parallel private(iskip,idx) reduction(+:gamma_cg,delta_cg)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
!$omp do
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           gamma_cg = gamma_cg + f(idx, 1)*f(idx, 2)
           delta_cg = delta_cg + f(idx, 3)*f(idx, 2)
        end do
!$omp end do nowait
     end do
!$omp end parallel

     iter=iter+1

     !==============================================
     ! Compute global norm
     !==============================================
     local(1) = gamma_cg; local(2) = delta_cg
     global(:) = 0.0
     call MPI_ALLREDUCE(local,global,2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,info)
     gamma_cg = global(1); delta_cg = global(2)

     !==============================================
     ! Wait for full completion of receives for ghostzones of w
     !==============================================
     if (countrecv>0) then
        call MPI_WAITALL(countrecv,reqrecv,MPI_STATUSES_IGNORE,info)
        call recv_post(f(1,3), ilevel)
      end if

     !==============================================
     ! Compute error
     !==============================================
     error=DSQRT(gamma_cg/dble(twotondim*numbtot(1,ilevel)))
     if(iter==1)error_ini=error
     if(myid==1 .and. verbose) write(*,112)iter,error/rhs_norm,error/error_ini,alpha_cg,beta_cg,gamma_cg,delta_cg

     !==============================================
     ! Compute alpha, beta factors
     !==============================================
     if (iter > 1) then
        beta_cg = gamma_cg / gamma_cg_old
        alpha_cg_old = alpha_cg
        alpha_cg = gamma_cg / (delta_cg - beta_cg * gamma_cg / alpha_cg_old)
     else
        beta_cg = 0.
        alpha_cg = gamma_cg / delta_cg
     endif

     gamma_cg_old = gamma_cg
     gamma_cg = 0.0; delta_cg = 0.0

     !==============================================
     ! Find m= Minv w. Do it in ghostzones too, so that p.m does not have to be synced below
     !==============================================
     call cmp_Minv_cg(ilevel,3,4)

      !==============================================
      ! Post receives for ghostzones of w for use in next iteration
      !==============================================
      if (error>epsilon*error_ini.and.iter<itermax) then
         call recv_virtual(ilevel, countrecv, reqrecv)
      end if

     !==============================================
     ! Compute n = A m
     !==============================================
     call cmp_A_cg(ilevel,4,5)

     !====================================
     ! Recurrence relations
     !====================================
!$omp parallel private(iskip,idx)
     do ind=1,twotondim
        iskip=ncoarse+(ind-1)*ngridmax
!$omp do
        do i=1,active(ilevel)%ngrid
           idx=active(ilevel)%igrid(i)+iskip
           f(idx,6) = f(idx,5) + beta_cg * f(idx,6) ! z   = n + beta*z
           f(idx,7) = f(idx,4) + beta_cg * f(idx,7) ! q   = m + beta*q
           f(idx,8) = f(idx,3) + beta_cg * f(idx,8) ! s   = w + beta*s
           f(idx,9) = f(idx,2) + beta_cg * f(idx,9) ! p   = u + beta*p

           phi(idx) = phi(idx) + alpha_cg * f(idx,9) ! phi = phi + alpha*p
           f(idx,1) = f(idx,1) - alpha_cg * f(idx,8) ! r   = r   - alpha*s
           f(idx,2) = f(idx,2) - alpha_cg * f(idx,7) ! u   = u   - alpha*q
           f(idx,3) = f(idx,3) - alpha_cg * f(idx,6) ! w   = w   - alpha*z
        end do
!$omp end do nowait
     end do
!$omp end parallel
  end do
  ! End main iteration loop

  if(myid==1)write(*,115)ilevel,iter,error/rhs_norm,error/error_ini
  if(iter >= itermax)then
     if(myid==1)write(*,*)'Poisson failed to converge...'
  end if
  !deallocate(f,nborl,addrl)

   if (countsend>0) call MPI_WAITALL(countsend,reqsend,MPI_STATUSES_IGNORE,info)

  ! Update boundaries
  call make_virtual_fine_dp(phi(1),ilevel)

111 format('   Entering phi_fine_cg for level ',I2)
112 format('   ==> Step=',i5,' Error=',6(1pe10.3,1x))
115 format('   ==> Level=',i5,' Step=',i5,' Error=',2(1pe10.3,1x))

end subroutine phi_fine_cg
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_residual_cg(ilevel,icount)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !------------------------------------------------------------------
  ! This routine computes the residual for the Conjugate Gradient
  ! Poisson solver. The residual is stored in f(i,1).
  !------------------------------------------------------------------
  integer::i,igrid,ngrid,ncache,nx_loc
  real(dp)::dx2,fourpi,scale,oneoversix,fact
  integer,dimension(1:3,1:2,1:8)::iii,jjj
  integer ,dimension(1:nvector)::ind_grid

  ! Set constants
  dx2=(0.5D0**ilevel)**2
  nx_loc=icoarse_max-icoarse_min+1
  scale=boxlen/dble(nx_loc)
  fourpi=4.D0*ACOS(-1.0D0)*scale
  if(cosmo)fourpi=1.5D0*omega_m*aexp*scale
  oneoversix=1.0D0/dble(twondim)
  fact=oneoversix*fourpi*dx2

  iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
  iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
  iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
  iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid,ind_grid)
  do igrid=1,ncache,nvector
      ! Gather nvector grids
      ngrid=MIN(nvector,ncache-igrid+1)
      ! Loop over myid grids by vector sweeps
      do i=1,ngrid
          ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      call cmprescg1(ilevel,icount,ind_grid,ngrid,iii,jjj,oneoversix,fact)
  enddo
end subroutine cmp_residual_cg
!###########################################################
!###########################################################
subroutine cmprescg1(ilevel,icount,ind_grid,ngrid,iii,jjj,oneoversix,fact)
  use amr_commons
  use pm_commons
  use hydro_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !------------------------------------------------------------------
  ! This routine computes the residual for the Conjugate Gradient
  ! Poisson solver. The residual is stored in f(i,1).
  !------------------------------------------------------------------
  integer::i,idim,ngrid,ind,iskip
  integer::id1,id2,ig1,ig2,ih1,ih2
  real(dp)::oneoversix,fact
  integer,dimension(1:3,1:2,1:8)::iii,jjj

  integer ,dimension(1:nvector)::ind_grid,ind_cell
  integer ,dimension(1:nvector,0:twondim)::igridn
  integer ,dimension(1:nvector,1:ndim)::ind_left,ind_right
  real(dp),dimension(1:nvector,1:ndim)::phig,phid
  real(dp),dimension(1:nvector,1:twotondim,1:ndim)::phi_left,phi_right
  real(dp),dimension(1:nvector)::residu

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

  ! Interpolate potential from upper level
  do idim=1,ndim
	  call interpol_phi(ind_left (1,idim),phi_left (1,1,idim),ngrid,ilevel,icount)
	  call interpol_phi(ind_right(1,idim),phi_right(1,1,idim),ngrid,ilevel,icount)
  end do

  ! Loop over cells
  do ind=1,twotondim
	  ! Gather neighboring potential
	  do idim=1,ndim
		  id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
		  ih1=ncoarse+(id1-1)*ngridmax
		  do i=1,ngrid
			  if(igridn(i,ig1)>0)then
				  phig(i,idim)=phi(igridn(i,ig1)+ih1)
			  else
				  phig(i,idim)=phi_left(i,id1,idim)
			  end if
		  end do
		  id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
		  ih2=ncoarse+(id2-1)*ngridmax
		  do i=1,ngrid
			  if(igridn(i,ig2)>0)then
				  phid(i,idim)=phi(igridn(i,ig2)+ih2)
			  else
				  phid(i,idim)=phi_right(i,id2,idim)
			  end if
		  end do
	  end do

	  ! Compute central cell index
	  iskip=ncoarse+(ind-1)*ngridmax
	  do i=1,ngrid
		  ind_cell(i)=iskip+ind_grid(i)
	  end do

	  ! Compute residual using 6 neighbors potential
	  do i=1,ngrid
		  residu(i)=phi(ind_cell(i))
	  end do
	  do idim=1,ndim
		  do i=1,ngrid
			  residu(i)=residu(i)-oneoversix*(phig(i,idim)+phid(i,idim))
		  end do
	  end do
	  do i=1,ngrid
		  residu(i)=residu(i)+fact*(rho(ind_cell(i))-rho_tot)
	  end do

	  ! Store results in f(i,1)
	  do i=1,ngrid
		  f(ind_cell(i),1)=residu(i)
	  end do
  end do
  ! End loop over cells

end subroutine cmprescg1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_Minv_cg(ilevel,x1,x2)
   use amr_commons
   use pm_commons
   use hydro_commons
   use poisson_commons
   implicit none
   integer::ilevel,x1,x2
   !------------------------------------------------------------------
   ! This routine computes Ap for the Conjugate Gradient
   ! Poisson Solver and store the result into f(i,3).
   !------------------------------------------------------------------
   integer::i,idim,igrid,ngrid,ncache,ind,iskip,icpu
   integer::id1,id2,ig1,ig2,ih1,ih2
   real(dp)::oneoversix
   integer,dimension(1:3,1:2,1:8)::iii,jjj
   integer,dimension(1:nvector)::ind_grid

   ! Set constants
   oneoversix=1.0D0/dble(twondim)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ! Loop over myid grids by vector sweeps
!$omp parallel private(ncache,ngrid,ind_grid)
   do icpu=1,ncpu
      if(icpu==myid) then
         ncache=active(ilevel)%ngrid
      else
         ncache=reception(icpu,ilevel)%ngrid
      end if
!$omp do
      do igrid=1,ncache,nvector
         ! Gather nvector grids
         ngrid=MIN(nvector,ncache-igrid+1)
         do i=1,ngrid
            if(icpu==myid) then
               ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
            else
               ind_grid(i)=reception(icpu,ilevel)%igrid(igrid+i-1)
            end if
         end do
         call cmpminvcg1(ilevel,iii,jjj,ind_grid,ngrid,oneoversix,x1,x2)
      end do
!$omp end do nowait
   end do
!$omp end parallel

end subroutine cmp_Minv_cg
!###########################################################
!###########################################################
subroutine cmpminvcg1(ilevel,iii,jjj,ind_grid,ngrid,oneoversix,x1,x2)
   use amr_commons
   use pm_commons
   use hydro_commons
   use poisson_commons
   implicit none
   integer::ilevel,x1,x2
   !------------------------------------------------------------------
   ! This routine computes Ap for the Conjugate Gradient
   ! Poisson Solver and store the result into f(i,3).
   !------------------------------------------------------------------
   integer::i,idim,igrid,ngrid,ncache,ind,iskip
   integer::id1,id2,ig1,ig2,ih1,ih2
   real(dp)::oneoversix,prefac
   integer,dimension(1:3,1:2,1:8)::iii,jjj

   integer,dimension(1:nvector)::ind_grid,ind_cell
   integer,dimension(1:nvector,0:twondim)::igridn
   real(dp),dimension(1:nvector,1:ndim)::phig,phid
   real(dp),dimension(1:nvector)::residu

   prefac = 1d0+0.5d0*oneoversix

   ! Gather neighboring grids
   do i=1,ngrid
      igridn(i,0)=ind_grid(i)
   end do
   do idim=1,ndim
      do i=1,ngrid
         igridn(i,2*idim-1)=son(nbor(ind_grid(i),2*idim-1))
         igridn(i,2*idim  )=son(nbor(ind_grid(i),2*idim  ))
      end do
   end do

   ! Loop over cells
   do ind=1,twotondim

      ! Gather neighboring potential
      do idim=1,ndim
         id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
         ih1=ncoarse+(id1-1)*ngridmax
         do i=1,ngrid
            if(igridn(i,ig1)>0)then
               phig(i,idim)=f(igridn(i,ig1)+ih1,x1)
            else
               phig(i,idim)=0.
            end if
         end do
         id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
         ih2=ncoarse+(id2-1)*ngridmax
         do i=1,ngrid
            if(igridn(i,ig2)>0)then
               phid(i,idim)=f(igridn(i,ig2)+ih2,x1)
            else
               phid(i,idim)=0.
            end if
         end do
      end do

      ! Compute central cell index
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,ngrid
         ind_cell(i)=iskip+ind_grid(i)
      end do

      ! Compute Ap using neighbors potential
      do i=1,ngrid
         residu(i)=prefac*f(ind_cell(i),x1)
      end do
      do idim=1,ndim
         do i=1,ngrid
            residu(i)=residu(i)+oneoversix*(phig(i,idim)+phid(i,idim))
         end do
      end do
      ! Store results in f(i,3)
      do i=1,ngrid
         f(ind_cell(i),x2)=residu(i)
      end do

   end do
   ! End loop over cells

end subroutine cmpminvcg1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine cmp_A_cg(ilevel,x1,x2)
   use amr_commons
   use pm_commons
   use hydro_commons
   use poisson_commons
   implicit none
   integer::ilevel,x1,x2
   !------------------------------------------------------------------
   ! This routine computes Ap for the Conjugate Gradient
   ! Poisson Solver and store the result into f(i,3).
   !------------------------------------------------------------------
   integer::i,idim,igrid,ngrid,ncache,ind,iskip
   integer::id1,id2,ig1,ig2,ih1,ih2
   real(dp)::oneoversix
   integer,dimension(1:3,1:2,1:8)::iii,jjj
   integer,dimension(1:nvector)::ind_grid

   ! Set constants
   oneoversix=1.0D0/dble(twondim)

   iii(1,1,1:8)=(/1,0,1,0,1,0,1,0/); jjj(1,1,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(1,2,1:8)=(/0,2,0,2,0,2,0,2/); jjj(1,2,1:8)=(/2,1,4,3,6,5,8,7/)
   iii(2,1,1:8)=(/3,3,0,0,3,3,0,0/); jjj(2,1,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(2,2,1:8)=(/0,0,4,4,0,0,4,4/); jjj(2,2,1:8)=(/3,4,1,2,7,8,5,6/)
   iii(3,1,1:8)=(/5,5,5,5,0,0,0,0/); jjj(3,1,1:8)=(/5,6,7,8,1,2,3,4/)
   iii(3,2,1:8)=(/0,0,0,0,6,6,6,6/); jjj(3,2,1:8)=(/5,6,7,8,1,2,3,4/)

   ! Loop over myid grids by vector sweeps
   ncache=active(ilevel)%ngrid
   !$omp parallel do private(igrid,ngrid,ind_grid)
   do igrid=1,ncache,nvector
      ! Gather nvector grids
      ngrid=MIN(nvector,ncache-igrid+1)
      do i=1,ngrid
         ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      call cmpacg1(ilevel,iii,jjj,ind_grid,ngrid,oneoversix,x1,x2)
   enddo

end subroutine cmp_A_cg
!###########################################################
!###########################################################
subroutine cmpacg1(ilevel,iii,jjj,ind_grid,ngrid,oneoversix,x1,x2)
   use amr_commons
   use pm_commons
   use hydro_commons
   use poisson_commons
   implicit none
   integer::ilevel,x1,x2
   !------------------------------------------------------------------
   ! This routine computes Ap for the Conjugate Gradient
   ! Poisson Solver and store the result into f(i,x2).
   !------------------------------------------------------------------
   integer::i,idim,igrid,ngrid,ncache,ind,iskip
   integer::id1,id2,ig1,ig2,ih1,ih2
   real(dp)::oneoversix
   integer,dimension(1:3,1:2,1:8)::iii,jjj

   integer,dimension(1:nvector)::ind_grid,ind_cell
   integer,dimension(1:nvector,0:twondim)::igridn
   real(dp),dimension(1:nvector,1:ndim)::phig,phid
   real(dp),dimension(1:nvector)::residu

   ! Gather neighboring grids
   do i=1,ngrid
      igridn(i,0)=ind_grid(i)
   end do
   do idim=1,ndim
      do i=1,ngrid
         igridn(i,2*idim-1)=son(nbor(ind_grid(i),2*idim-1))
         igridn(i,2*idim  )=son(nbor(ind_grid(i),2*idim  ))
      end do
   end do

   ! Loop over cells
   do ind=1,twotondim

      ! Gather neighboring potential
      do idim=1,ndim
         id1=jjj(idim,1,ind); ig1=iii(idim,1,ind)
         ih1=ncoarse+(id1-1)*ngridmax
         do i=1,ngrid
            if(igridn(i,ig1)>0)then
               phig(i,idim)=f(igridn(i,ig1)+ih1,x1)
            else
               phig(i,idim)=0.
            end if
         end do
         id2=jjj(idim,2,ind); ig2=iii(idim,2,ind)
         ih2=ncoarse+(id2-1)*ngridmax
         do i=1,ngrid
            if(igridn(i,ig2)>0)then
               phid(i,idim)=f(igridn(i,ig2)+ih2,x1)
            else
               phid(i,idim)=0.
            end if
         end do
      end do

      ! Compute central cell index
      iskip=ncoarse+(ind-1)*ngridmax
      do i=1,ngrid
         ind_cell(i)=iskip+ind_grid(i)
      end do

      ! Compute Ap using neighbors potential
      do i=1,ngrid
         residu(i)=-f(ind_cell(i),x1)
      end do
      do idim=1,ndim
         do i=1,ngrid
            residu(i)=residu(i)+oneoversix*(phig(i,idim)+phid(i,idim))
         end do
      end do
      ! Store results in f(i,3)
      do i=1,ngrid
         f(ind_cell(i),x2)=residu(i)
      end do

   end do
   ! End loop over cells

end subroutine cmpacg1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_initial_phi(ilevel,icount)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !
  !
  !
  integer::igrid,ncache,i,ngrid
  integer ,dimension(1:nvector)::ind_grid

  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid,ind_grid)
  do igrid=1,ncache,nvector
      ! Gather nvector grids
      ngrid=MIN(nvector,ncache-igrid+1)
      ! Loop over myid grids by vector sweeps
      do i=1,ngrid
          ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      call makeinitphi1(ilevel,icount,ind_grid,ngrid)
  enddo
end subroutine make_initial_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine makeinitphi1(ilevel,icount,ind_grid,ngrid)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel,icount
  !
  !
  !
  integer::i,ngrid,ind,iskip,idim
  integer ,dimension(1:nvector)::ind_grid,ind_cell,ind_cell_father
  real(dp),dimension(1:nvector,1:twotondim)::phi_int

  if(ilevel==1)then
	  ! Loop over cells
	  do ind=1,twotondim
		  iskip=ncoarse+(ind-1)*ngridmax
		  do i=1,ngrid
			  ind_cell(i)=iskip+ind_grid(i)
		  end do
		  do i=1,ngrid
			  phi(ind_cell(i))=0.0d0
		  end do
		  do idim=1,ndim
			  do i=1,ngrid
				  f(ind_cell(i),idim)=0.0
			  end do
		  end do
	  end do
	  ! End loop over cells
  else
	  ! Compute father cell index
	  do i=1,ngrid
		  ind_cell_father(i)=father(ind_grid(i))
	  end do

	  ! Interpolate
	  call interpol_phi(ind_cell_father,phi_int,ngrid,ilevel,icount)

	  ! Loop over cells
	  do ind=1,twotondim
		  iskip=ncoarse+(ind-1)*ngridmax
		  do i=1,ngrid
			  ind_cell(i)=iskip+ind_grid(i)
		  end do
		  do i=1,ngrid
			  phi(ind_cell(i))=phi_int(i,ind)
		  end do
		  do idim=1,ndim
			  do i=1,ngrid
				  f(ind_cell(i),idim)=0.0
			  end do
		  end do
	  end do
	  ! End loop over cells
  end if
  ! End loop over grids
end subroutine makeinitphi1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine make_multipole_phi(ilevel)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel
  integer::i,ncache,igrid,ngrid,ind
  integer::nx_loc,ix,iy,iz
  integer,dimension(1:nvector)::ind_grid

  real(dp)::dx,dx_loc,scale,fourpi,boxlen2,eps
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc


  ! Mesh size at level ilevel
  dx=0.5D0**ilevel

  ! Rescaling factors
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  fourpi=4.D0*ACOS(-1.0D0)
  boxlen2=boxlen**2
  eps=dx_loc

  ! Set position of cell centers relative to grid center
  do ind=1,twotondim
     iz=(ind-1)/4
     iy=(ind-1-4*iz)/2
     ix=(ind-1-2*iy-4*iz)
     if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
     if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
     if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
  end do

  ! Loop over myid grids by vector sweeps
  ncache=active(ilevel)%ngrid
!$omp parallel do private(igrid,ngrid,ind_grid)
  do igrid=1,ncache,nvector
      ! Gather nvector grids
      ngrid=MIN(nvector,ncache-igrid+1)
      do i=1,ngrid
          ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
      end do
      call makemultiphi1(ilevel,ind_grid,ngrid,scale,eps,xc,skip_loc)
  enddo
end subroutine make_multipole_phi
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine makemultiphi1(ilevel,ind_grid,ngrid,scale,eps,xc,skip_loc)
  use amr_commons
  use pm_commons
  use poisson_commons
  implicit none
  integer::ilevel
  integer::idim
  integer::i,ngrid,ind
  integer::iskip
  integer,dimension(1:nvector)::ind_grid,ind_cell

  real(dp)::scale,eps,r2
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector)::rr,pp
  real(dp),dimension(1:nvector,1:ndim)::xx

  ! Loop over cells
  do ind=1,twotondim
	  iskip=ncoarse+(ind-1)*ngridmax
	  do i=1,ngrid
		  ind_cell(i)=iskip+ind_grid(i)
	  end do

	  if(simple_boundary)then
		  ! Compute cell center in code units
		  do idim=1,ndim
				do i=1,ngrid
					xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
				end do
		  end do

		  ! Rescale position from code units to user units
		  rr(1:ngrid)=0.0d0
		  do idim=1,ndim
				do i=1,ngrid
					xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
					rr(i)=rr(i)+(xx(i,idim)-multipole(idim+1)/multipole(1))**2
				end do
		  end do

		  do i=1,ngrid
				rr(i)=max(eps,sqrt(rr(i)))       ! Cutoff
		  end do

		  if(ngrid>0) call phi_ana(rr,pp,ngrid)

		  ! Scatter variables
		  do i=1,ngrid
				phi(ind_cell(i))=pp(i)/scale
		  end do

	  else
		  do i=1,ngrid
				phi(ind_cell(i))=0d0
		  end do
	  endif

	  ! End loop over cells
  end do
  ! End loop over grids

end subroutine makemultiphi1
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine send_virtual(xx, ilevel, countsend, reqsend)
   use amr_commons
   use poisson_commons
   use mpi_mod
   implicit none
   integer::ilevel
   real(dp),dimension(1:ncoarse+twotondim*ngridmax)::xx
   ! -------------------------------------------------------------------
   ! This routine communicates virtual boundaries among all cpu's.
   ! at level ilevel for any double precision array in the AMR grid.
   ! -------------------------------------------------------------------
#ifndef WITHOUTMPI
   integer::icpu,i,j,ncache,iskip,step,idx
   integer::countsend
   integer::info,tag=101
   integer,dimension(ncpu)::reqsend
#endif

#ifndef WITHOUTMPI
!$omp parallel private(ncache,iskip,step,idx)
   do j=1,twotondim
      do icpu=1,ncpu
         ncache=emission(icpu,ilevel)%ngrid
         if (ncache>0) then
            iskip=ncoarse+(j-1)*ngridmax
            step=(j-1)*ncache
!$omp do
            do i=1,ncache
               idx=emission(icpu,ilevel)%igrid(i)+iskip
               emission(icpu,ilevel)%u(i+step,1)=xx(idx)
            end do
!$omp end do nowait
         end if
      end do
   end do
!$omp end parallel

   countsend=0
   do icpu=1,ncpu
      ncache=emission(icpu,ilevel)%ngrid
      if(ncache>0) then
         countsend=countsend+1
         call MPI_ISEND(emission(icpu,ilevel)%u,ncache*twotondim, &
               & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqsend(countsend),info)
      end if
   end do
#endif
end subroutine send_virtual
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine recv_virtual(ilevel, countrecv, reqrecv)
   use amr_commons
   use mpi_mod
   implicit none
   integer::ilevel
   ! -------------------------------------------------------------------
   ! This routine communicates virtual boundaries among all cpu's.
   ! at level ilevel for any double precision array in the AMR grid.
   ! -------------------------------------------------------------------
#ifndef WITHOUTMPI
   integer::icpu,ncache
   integer::countrecv
   integer::info,tag=101
   integer,dimension(ncpu)::reqrecv
#endif

#ifndef WITHOUTMPI
   countrecv=0
   do icpu=1,ncpu
      ncache=reception(icpu,ilevel)%ngrid
      if(ncache>0) then
         countrecv=countrecv+1
         call MPI_IRECV(reception(icpu,ilevel)%u,ncache*twotondim, &
               & MPI_DOUBLE_PRECISION,icpu-1,tag,MPI_COMM_WORLD,reqrecv(countrecv),info)
      end if
   end do
#endif
end subroutine recv_virtual
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine recv_post(xx, ilevel)
   use amr_commons
   implicit none
   integer::ilevel
   real(dp),dimension(1:ncoarse+twotondim*ngridmax)::xx
   ! -------------------------------------------------------------------
   ! This routine communicates virtual boundaries among all cpu's.
   ! at level ilevel for any double precision array in the AMR grid.
   ! -------------------------------------------------------------------
#ifndef WITHOUTMPI
   integer::icpu,i,j,ncache,iskip,step,idx
   integer::countrecv
#endif

#ifndef WITHOUTMPI
!$omp parallel private(ncache,iskip,step,idx)
   do j=1,twotondim
      iskip = ncoarse+(j-1)*ngridmax
      do icpu=1,ncpu
         ncache = reception(icpu,ilevel)%ngrid
         step=(j-1)*ncache
         if(ncache>0) then
!$omp do
            do i=1,ncache
               idx = reception(icpu,ilevel)%igrid(i)+iskip
               xx(idx) = reception(icpu,ilevel)%u(i+step,1)
            end do
!$omp end do nowait
         end if
      end do
   end do
!$omp end parallel
#endif
end subroutine recv_post