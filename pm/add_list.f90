!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_list(ind_part,ind_grid,ok,np)
  use amr_commons
  use pm_commons
  implicit none
  integer, intent(in )::np
  integer,dimension(1:nvector), intent(in)::ind_part,ind_grid
  logical,dimension(1:nvector), intent(in)::ok

  !
  ! Add particles to their new linked lists
  !
  integer::j
  ! should be called with critical if ind_grid is not single grid.
  do j=1,np
     if(ok(j))then
        if (numbp(ind_grid(j)) > 0) then
           ! Add particle at the tail of its linked list
           nextp(tailp(ind_grid(j))) = ind_part(j)
           prevp(ind_part(j)) = tailp(ind_grid(j))
           nextp(ind_part(j))=0
           tailp(ind_grid(j)) = ind_part(j)
           numbp(ind_grid(j)) = numbp(ind_grid(j)) + 1
        else
           ! Initialise linked list
           headp(ind_grid(j)) = ind_part(j)
           tailp(ind_grid(j)) = ind_part(j)
           prevp(ind_part(j))=0
           nextp(ind_part(j))=0
           numbp(ind_grid(j)) = 1
        end if
     end if
  end do
end subroutine add_list
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_list_single(ind_part,ind_grid,ok,np)
  use amr_commons
  use pm_commons
  implicit none
  integer, intent(in)::np,ind_grid
  integer,dimension(1:nvector),intent(in)::ind_part
  logical,dimension(1:nvector),intent(in)::ok
  integer::j,tailp_old,tailp_new,headp_new,nadd

  nadd = 0
  headp_new = 0
  do j=1,np
     if(ok(j))then
        if(headp_new == 0) then
           headp_new = ind_part(j)
        end if
        tailp_new = ind_part(j)
        nadd = nadd + 1
     end if
  end do

  if(nadd > 0) then
     ! Extend the tail
     tailp_old = tailp(ind_grid)
     tailp(ind_grid) = tailp_new
     if (numbp(ind_grid) <= 0) then
        headp(ind_grid) = headp_new
     end if
     numbp(ind_grid) = numbp(ind_grid) + nadd
     nextp(tailp_new) = 0

     ! Fill the middle
     do j=1,np
        if(ok(j)) then
           ! Add particle at the tail of its linked list
           if(tailp_old /= 0) nextp(tailp_old) = ind_part(j)
           prevp(ind_part(j)) = tailp_old
           tailp_old = ind_part(j)
        end if
     end do
  end if
end subroutine add_list_single
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_list_single_critical(ind_part,ind_grid,ok,np)
  use amr_commons
  use pm_commons
  implicit none
  integer, intent(in)::np,ind_grid
  integer,dimension(1:nvector),intent(in)::ind_part
  logical,dimension(1:nvector),intent(in)::ok
  integer::j,tailp_old,tailp_new,headp_new,nadd

  nadd = 0
  headp_new = 0
  do j=1,np
     if(ok(j))then
        if(headp_new == 0) then
           headp_new = ind_part(j)
        end if
        tailp_new = ind_part(j)
        nadd = nadd + 1
     end if
  end do

  if(nadd > 0) then
     ! Extend the tail
!$omp critical
     tailp_old = tailp(ind_grid)
     tailp(ind_grid) = tailp_new
     if (numbp(ind_grid) <= 0) then
        headp(ind_grid) = headp_new
     end if
     numbp(ind_grid) = numbp(ind_grid) + nadd
     nextp(tailp_new) = 0
!$omp end critical

     ! Fill the middle
     do j=1,np
        if(ok(j)) then
           ! Add particle at the tail of its linked list
           if(tailp_old /= 0) nextp(tailp_old) = ind_part(j)
           prevp(ind_part(j)) = tailp_old
           tailp_old = ind_part(j)
        end if
     end do
  end if
end subroutine add_list_single_critical
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_free(ind_part,np)
  use amr_commons
  use pm_commons
#ifdef DICE
  use dice_commons
#endif
  implicit none
  integer, intent(in)::np
  integer,dimension(1:nvector), intent(in)::ind_part
  !
  ! Add particles to the free memory linked list
  ! and reset all particle variables
  !
  integer::j,idim

  do idim=1,ndim
     do j=1,np
        xp(ind_part(j),idim)=0.0
        vp(ind_part(j),idim)=0.0
     end do
  end do
  do j=1,np
     mp(ind_part(j))=0.0
     if(use_initial_mass)mp0(ind_part(j))=0.0
     idp(ind_part(j))=0
     levelp(ind_part(j))=0
     typep(ind_part(j))%family=FAM_UNDEF
     typep(ind_part(j))%tag=0
  end do
  if(star.or.sink)then
     do j=1,np
        tp(ind_part(j))=0.0
        if(write_stellar_densities) then
           st_n_tp(ind_part(j))=0.0
           st_n_SN(ind_part(j))=0.0
           st_e_SN(ind_part(j))=0.0
        endif
     end do
     if(metal)then
        do j=1,np
           zp(ind_part(j))=0.0
        end do
     end if
  end if
#ifdef DICE
  ! DICE patch
  if(dice_init) then
     do j=1,np
        up(ind_part(j))=0.0
     end do
  endif
#endif

!$omp critical(omp_particle_free)
  do j=1,np
     if(numbp_free>0)then
        ! Add particle at the tail of its linked list
        nextp(tailp_free)=ind_part(j)
        prevp(ind_part(j))=tailp_free
        nextp(ind_part(j))=0
        tailp_free=ind_part(j)
        numbp_free=numbp_free+1
     else
        ! Initialise linked list
        headp_free=ind_part(j)
        tailp_free=ind_part(j)
        prevp(ind_part(j))=0
        nextp(ind_part(j))=0
        numbp_free=1
     end if
  end do
  npart=npartmax-numbp_free
!$omp end critical(omp_particle_free)
end subroutine add_free
!################################################################
!################################################################
!################################################################
!################################################################
subroutine add_free_cond(ind_part,ok,np)
  use amr_commons
  use pm_commons
#ifdef DICE
  use dice_commons
#endif
  implicit none
  integer::np
  integer,dimension(1:nvector)::ind_part
  logical,dimension(1:nvector)::ok
  !
  ! Add particles to the free memory linked list
  ! and reset all particle variables
  !
  integer::j,idim

  do idim=1,ndim
     do j=1,np
        if(ok(j))then
           xp(ind_part(j),idim)=0.0
           vp(ind_part(j),idim)=0.0
        endif
     end do
  end do
  do j=1,np
     if(ok(j))then
        mp(ind_part(j))=0.0
        if(use_initial_mass) mp0(ind_part(j))=0.0
        idp(ind_part(j))=0
        levelp(ind_part(j))=0
        typep(ind_part(j))%family = FAM_UNDEF
        typep(ind_part(j))%tag = 0
     endif
  end do
  if(star.or.sink)then
     do j=1,np
        if(ok(j))then
           tp(ind_part(j))=0.0
           if(write_stellar_densities) then
              st_n_tp(ind_part(j))=0.0
              st_n_SN(ind_part(j))=0.0
              st_e_SN(ind_part(j))=0.0
           endif
        endif
     end do
     if(metal)then
        do j=1,np
           if(ok(j))then
              zp(ind_part(j))=0.0
           endif
        end do
     end if
  end if

#ifdef DICE
  ! DICE patch
  if(dice_init) then
     do j=1,np
        if(ok(j))then
           up(ind_part(j))=0.0
        endif
     end do
  endif
#endif

!$omp critical(omp_particle_free)
  do j=1,np
     if(ok(j))then
        if(numbp_free>0)then
           ! Add particle at the tail of its linked list
           nextp(tailp_free)=ind_part(j)
           prevp(ind_part(j))=tailp_free
           nextp(ind_part(j))=0
           tailp_free=ind_part(j)
           numbp_free=numbp_free+1
        else
           ! Initialise linked list
           headp_free=ind_part(j)
           tailp_free=ind_part(j)
           prevp(ind_part(j))=0
           nextp(ind_part(j))=0
           numbp_free=1
        end if
    endif
  end do
  npart=npartmax-numbp_free
!$omp end critical(omp_particle_free)
end subroutine add_free_cond
