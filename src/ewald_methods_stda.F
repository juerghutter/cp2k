!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!!! JUST TEMPORARY FILE FOR sTDA
!> \brief Calculation of Ewald contributions in DFTB
!> \author JGH
! **************************************************************************************************
MODULE ewald_methods_stda
   USE atprop_types,                    ONLY: atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_fm_types,                     ONLY: cp_fm_get_element,&
                                              cp_fm_p_type,&
                                              cp_fm_set_all,&
                                              cp_fm_set_element
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dgs,                             ONLY: dg_sum_patch,&
                                              dg_sum_patch_force_1d
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE ewald_pw_types,                  ONLY: ewald_pw_get,&
                                              ewald_pw_type
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi,&
                                              pi
   USE particle_types,                  ONLY: particle_type
   USE pme_tools,                       ONLY: get_center,&
                                              set_list
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_grids,                        ONLY: get_pw_grid_info
   USE pw_methods,                      ONLY: pw_transfer
   USE pw_poisson_methods,              ONLY: pw_poisson_rebuild,&
                                              pw_poisson_solve
   USE pw_poisson_types,                ONLY: do_ewald_ewald,&
                                              do_ewald_none,&
                                              do_ewald_pme,&
                                              do_ewald_spme,&
                                              greens_fn_type,&
                                              pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_p_type,&
                                              pw_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind_set,&
                                              qs_kind_type
   USE qs_neighbor_list_types,          ONLY: get_iterator_info,&
                                              neighbor_list_iterate,&
                                              neighbor_list_iterator_create,&
                                              neighbor_list_iterator_p_type,&
                                              neighbor_list_iterator_release,&
                                              neighbor_list_set_p_type
   USE qs_p_env_types,                  ONLY: qs_p_env_type
   USE qs_tddfpt_stda_types,            ONLY: stda_env_type
   USE realspace_grid_types,            ONLY: &
        pw2rs, realspace_grid_desc_type, realspace_grid_p_type, realspace_grid_type, rs2pw, &
        rs_grid_create, rs_grid_release, rs_grid_set_box, rs_grid_zero, rs_pw_transfer
   USE spme,                            ONLY: get_patch
   USE virial_types,                    ONLY: virial_type
   USE basis_set_container_types,       ONLY: get_basis_from_container
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE particle_methods,                ONLY: get_particle_set
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE input_section_types,             ONLY: section_vals_get,& 
                                              section_vals_type,&
                                              section_vals_val_get
   USE input_constants,                 ONLY: tddfpt_singlet,&
                                              tddfpt_triplet
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'ewald_methods_stda'

   PUBLIC :: stda_ewald_overlap, stda_spme_evaluate, ewald_or_direct_sum

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param ewald_env ...
!> \param ewald_pw ...
!> \param particle_set ...
!> \param box ...
! **************************************************************************************************
   SUBROUTINE stda_spme_evaluate(qs_env, p_env, ewald_env, ewald_pw, particle_set, box)
      !  gmcharge, mcharge)!, calculate_forces, virial, use_virial, atprop)
!
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(particle_type), DIMENSION(:), INTENT(IN)      :: particle_set
      TYPE(cell_type), POINTER                           :: box

      CHARACTER(len=*), PARAMETER :: routineN = 'stda_spme_evaluate', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: group, handle, i, ii, ipart, ispin, jj, &
                                                            n, npart, o_spline, output_unit, p1
      INTEGER, ALLOCATABLE, DIMENSION(:, :)              :: center
      INTEGER, DIMENSION(3)                              :: npts
      REAL(KIND=dp)                                      :: alpha, charge, dvols, fint, tmp, vgc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :, :)     :: rhos
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(greens_fn_type), POINTER                      :: green
      TYPE(pw_grid_type), POINTER                        :: grid_spme
      TYPE(pw_p_type), DIMENSION(3)                      :: dphi_g
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: pw_pool
      TYPE(pw_type), POINTER                             :: phi_g, rhob_g, rhob_r
      TYPE(realspace_grid_desc_type), POINTER            :: rs_desc
      TYPE(realspace_grid_type), POINTER                 :: rden, rpot

!      REAL(KIND=dp), DIMENSION(:), INTENT(inout)         :: mcharge
!      LOGICAL, INTENT(in)                                :: calculate_forces
!      TYPE(virial_type), POINTER                         :: virial
!      LOGICAL, INTENT(in)                                :: use_virial
!      TYPE(atprop_type), POINTER                         :: atprop
!     REAL(KIND=dp)                                      :: alpha, dvols, fat(3), ffa, ffb, fint, vgc
!   REAL(KIND=dp), DIMENSION(3, 3)                     :: f_stress, h_stress
!     TYPE(realspace_grid_p_type), DIMENSION(:), POINTER :: drpot

      NULLIFY (logger) !get output_unit
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CALL timeset(routineN, handle)

      CALL get_qs_env(qs_env, dft_control=dft_control)

      !-------------- INITIALISATION ---------------------
      CALL ewald_env_get(ewald_env, alpha=alpha, o_spline=o_spline, group=group, &
                         para_env=para_env)
      NULLIFY (green, poisson_env, pw_pool)
      CALL ewald_pw_get(ewald_pw, pw_big_pool=pw_pool, rs_desc=rs_desc, &
                        poisson_env=poisson_env)
      CALL pw_poisson_rebuild(poisson_env)
      green => poisson_env%green_fft
      grid_spme => pw_pool%pw_grid

      CALL get_pw_grid_info(grid_spme, dvol=dvols, npts=npts)

      npart = SIZE(particle_set)

      n = o_spline
      ALLOCATE (rhos(n, n, n))

      CALL rs_grid_create(rden, rs_desc)
      CALL rs_grid_set_box(grid_spme, rs=rden)
      CALL rs_grid_zero(rden)

      ALLOCATE (center(3, npart))
      CALL get_center(particle_set, box, center, npts, n)

      DO ispin = 1, dft_control%nspins
      DO ii = 1, p_env%n_mo(ispin)
      DO jj = 1, p_env%n_mo(ispin)

         !-------------- DENSITY CALCULATION ----------------
         ipart = 0
         DO
            CALL set_list(particle_set, npart, center, p1, rden, ipart)
            IF (p1 == 0) EXIT

            ! calculate function on small boxes
            CALL get_patch(particle_set, box, green, npts, p1, rhos, is_core=.FALSE., &
                           is_shell=.FALSE., unit_charge=.TRUE.)

            CALL cp_fm_get_element(p_env%stda_env%atom_param_set(p1)%atom_param%charges(ispin)%matrix, &
                                   jj, ii, charge)
            rhos(:, :, :) = rhos(:, :, :)*charge

            ! add boxes to real space grid (big box)
            CALL dg_sum_patch(rden, rhos, center(:, p1))
         END DO

         CALL pw_pool_create_pw(pw_pool, rhob_r, use_data=REALDATA3D, &
                                in_space=REALSPACE)

         CALL rs_pw_transfer(rden, rhob_r, rs2pw)

         ! transform density to G space and add charge function
         CALL pw_pool_create_pw(pw_pool, rhob_g, use_data=COMPLEXDATA1D, &
                                in_space=RECIPROCALSPACE)
         CALL pw_transfer(rhob_r, rhob_g)
         ! update charge function
         rhob_g%cc = rhob_g%cc*green%p3m_charge%cr

         !-------------- ELECTROSTATIC CALCULATION -----------

         ! allocate intermediate arrays
         DO i = 1, 3
            NULLIFY (dphi_g(i)%pw)
            CALL pw_pool_create_pw(pw_pool, dphi_g(i)%pw, &
                                   use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
         END DO
         CALL pw_pool_create_pw(pw_pool, phi_g, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
!      IF (use_virial) THEN
!         CALL pw_poisson_solve(poisson_env, rhob_g, vgc, phi_g, dphi_g, h_stress=h_stress)
!      ELSE
         CALL pw_poisson_solve(poisson_env, rhob_g, vgc, phi_g, dphi_g)
!      END IF

         CALL rs_grid_create(rpot, rs_desc)
         CALL rs_grid_set_box(grid_spme, rs=rpot)

         ! Atomic Stress
!      IF (atprop%stress .AND. use_virial) THEN
!
!         CALL rs_grid_zero(rpot)
!         rhob_g%cc = phi_g%cc*green%p3m_charge%cr
!         CALL pw_transfer(rhob_g, rhob_r)
!         CALL rs_pw_transfer(rpot, rhob_r, pw2rs)
!         ipart = 0
!         DO
!            CALL set_list(particle_set, npart, center, p1, rden, ipart)
!            IF (p1 == 0) EXIT
!            ! calculate function on small boxes
!            CALL get_patch(particle_set, box, green, npts, p1, rhos, is_core=.FALSE., &
!                           is_shell=.FALSE., unit_charge=.TRUE.)
!            CALL dg_sum_patch_force_1d(rpot, rhos, center(:, p1), fint)
!            atprop%atstress(1, 1, p1) = atprop%atstress(1, 1, p1)+0.5_dp*mcharge(p1)*fint*dvols
!            atprop%atstress(2, 2, p1) = atprop%atstress(2, 2, p1)+0.5_dp*mcharge(p1)*fint*dvols
!            atprop%atstress(3, 3, p1) = atprop%atstress(3, 3, p1)+0.5_dp*mcharge(p1)*fint*dvols
!         END DO
!
!         CALL pw_pool_create_pw(pw_pool, phib_g, &
!                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
!         ffa = (0.5_dp/alpha)**2
!         ffb = 1.0_dp/fourpi
!         DO i = 1, 3
!            DO ig = grid_spme%first_gne0, grid_spme%ngpts_cut_local
!               phib_g%cc(ig) = ffb*dphi_g(i)%pw%cc(ig)*(1.0_dp+ffa*grid_spme%gsq(ig))
!               phib_g%cc(ig) = phib_g%cc(ig)*green%influence_fn%cc(ig)
!            END DO
!            IF (grid_spme%have_g0) phib_g%cc(1) = 0.0_dp
!            DO j = 1, i
!               nd = 0
!               nd(j) = 1
!               CALL pw_copy(phib_g, rhob_g)
!               CALL pw_derive(rhob_g, nd)
!               rhob_g%cc = rhob_g%cc*green%p3m_charge%cr
!               CALL pw_transfer(rhob_g, rhob_r)
!               CALL rs_pw_transfer(rpot, rhob_r, pw2rs)
!
!               ipart = 0
!               DO
!                  CALL set_list(particle_set, npart, center, p1, rden, ipart)
!                  IF (p1 == 0) EXIT
!                  ! calculate function on small boxes
!                  CALL get_patch(particle_set, box, green, npts, p1, rhos, &
!                                 is_core=.FALSE., is_shell=.FALSE., unit_charge=.TRUE.)
!                  ! integrate box and potential
!                  CALL dg_sum_patch_force_1d(rpot, rhos, center(:, p1), fint)
!                  atprop%atstress(i, j, p1) = atprop%atstress(i, j, p1)+fint*dvols*mcharge(p1)
!                  IF (i /= j) atprop%atstress(j, i, p1) = atprop%atstress(j, i, p1)+fint*dvols*mcharge(p1)
!               END DO
!
!            END DO
!         END DO
!         CALL pw_pool_give_back_pw(pw_pool, phib_g)
!
!      END IF

         CALL pw_pool_give_back_pw(pw_pool, rhob_g)

         CALL rs_grid_zero(rpot)
         phi_g%cc = phi_g%cc*green%p3m_charge%cr
         CALL pw_transfer(phi_g, rhob_r)
         CALL pw_pool_give_back_pw(pw_pool, phi_g)
         CALL rs_pw_transfer(rpot, rhob_r, pw2rs)

         !---------- END OF ELECTROSTATIC CALCULATION --------

!      !------------- STRESS TENSOR CALCULATION ------------
!
!      IF (use_virial) THEN
!         DO i = 1, 3
!            DO j = i, 3
!               f_stress(i, j) = pw_integral_a2b(dphi_g(i)%pw, dphi_g(j)%pw)
!               f_stress(j, i) = f_stress(i, j)
!            END DO
!         END DO
!         ffa = (1.0_dp/fourpi)*(0.5_dp/alpha)**2
!         virial%pv_virial = virial%pv_virial-(ffa*f_stress-h_stress)/REAL(para_env%num_pe, dp)
!      END IF
!
!      !--------END OF STRESS TENSOR CALCULATION -----------
!
!      IF (calculate_forces) THEN
!         ! move derivative of potential to real space grid and
!         ! multiply by charge function in g-space
!         ALLOCATE (drpot(1:3))
!         DO i = 1, 3
!            CALL rs_grid_create(drpot(i)%rs_grid, rs_desc)
!            CALL rs_grid_set_box(grid_spme, rs=drpot(i)%rs_grid)
!            dphi_g(i)%pw%cc = dphi_g(i)%pw%cc*green%p3m_charge%cr
!            CALL pw_transfer(dphi_g(i)%pw, rhob_r)
!            CALL pw_pool_give_back_pw(pw_pool, dphi_g(i)%pw)
!            CALL rs_pw_transfer(drpot(i)%rs_grid, rhob_r, pw2rs)
!         END DO
!      ELSE
         DO i = 1, 3
            CALL pw_pool_give_back_pw(pw_pool, dphi_g(i)%pw)
         END DO
!      END IF
         CALL pw_pool_give_back_pw(pw_pool, rhob_r)
!
!      !----------------- FORCE CALCULATION ----------------
!
         ipart = 0
         DO
!
            CALL set_list(particle_set, npart, center, p1, rden, ipart)
            IF (p1 == 0) EXIT
!
!         ! calculate function on small boxes
            CALL get_patch(particle_set, box, green, npts, p1, rhos, is_core=.FALSE., &
                           is_shell=.FALSE., unit_charge=.TRUE.)
!
            CALL dg_sum_patch_force_1d(rpot, rhos, center(:, p1), fint)

            CALL cp_fm_get_element(p_env%stda_env%atom_param_set(p1)%atom_param%gmcharge(ispin)%matrix, ii, jj, tmp)
            tmp = tmp+fint*dvols
            CALL cp_fm_set_element(p_env%stda_env%atom_param_set(p1)%atom_param%gmcharge(ispin)%matrix, ii, jj, tmp)
!
!         IF (calculate_forces) THEN
!            CALL dg_sum_patch_force_3d(drpot, rhos, center(:, p1), fat)
!            gmcharge(p1, 2) = gmcharge(p1, 2)-fat(1)*dvols
!            gmcharge(p1, 3) = gmcharge(p1, 3)-fat(2)*dvols
!            gmcharge(p1, 4) = gmcharge(p1, 4)-fat(3)*dvols
!         END IF
!
         END DO
!
!      !--------------END OF FORCE CALCULATION -------------
      ENDDO !jj
      ENDDO !ii
      ENDDO !ispin

      !------------------CLEANING UP ----------------------

      CALL rs_grid_release(rden)
      CALL rs_grid_release(rpot)
!      IF (calculate_forces) THEN
!         IF (ASSOCIATED(drpot)) THEN
!            DO i = 1, 3
!               CALL rs_grid_release(drpot(i)%rs_grid)
!            END DO
!            DEALLOCATE (drpot)
!         END IF
!      END IF
      DEALLOCATE (rhos)
      DEALLOCATE (center)

      CALL timestop(handle)

   END SUBROUTINE stda_spme_evaluate

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param stda_env ...
!> \param alpha ...
!> \param n_list ...
! **************************************************************************************************
   SUBROUTINE stda_ewald_overlap(qs_env, p_env, stda_env, alpha, n_list)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(stda_env_type), POINTER                       :: stda_env
      REAL(KIND=dp), INTENT(in)                          :: alpha
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list

      CHARACTER(LEN=*), PARAMETER :: routineN = 'stda_ewald_overlap', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, ii, ispin, jatom, jj, &
                                                            nspins, output_unit
      REAL(KIND=dp)                                      :: dr, fr, loewdincharge_i, &
                                                            loewdincharge_j, rij(3), tmp
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator

!      REAL(KIND=dp), DIMENSION(:, :,:), INTENT(inout)    :: gmcharge
!, nmat
!      REAL(KIND=dp)                                      :: dfr, dr, fr, pfr, rij(3)

      CALL timeset(routineN, handle)

      NULLIFY (logger) !get output_unit
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

!      nmat = SIZE(gmcharge, 2)
!      WRITE(output_unit,*) 'nmat', nmat

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins
      DO ispin = 1, nspins

         CALL neighbor_list_iterator_create(nl_iterator, n_list)
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0) ! here only i /= j
            CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, r=rij)

            dr = SQRT(SUM(rij(:)**2))
            IF (dr > 1.e-10_dp) THEN
!       IF (dr > 1.e-6_dp) THEN  !!! WITH THIS CANCELLATION OF 1/r TERM
      !         WRITE (output_unit, *) 'iatom, jatom', iatom, jatom !
!         fr = 1/dr !!! WITH THIS CANCELLATION OF 1/r TERM
      !         WRITE (output_unit, *) 'dr', dr, 1/dr
               fr = erfc(alpha*dr)/dr
      !         WRITE (output_unit, *) 'alpha', alpha

               CALL cp_fm_set_all(stda_env%atom_param_set(iatom)%atom_param%gmcharge(ispin)%matrix, 0.0_dp, 0.0_dp)
               CALL cp_fm_set_all(stda_env%atom_param_set(jatom)%atom_param%gmcharge(ispin)%matrix, 0.0_dp, 0.0_dp)

               DO ii = 1, p_env%n_mo(ispin)
                  DO jj = 1, p_env%n_mo(ispin)

                     CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%charges(ispin)%matrix, &
                                            jj, ii, loewdincharge_i)
                     CALL cp_fm_get_element(stda_env%atom_param_set(jatom)%atom_param%charges(ispin)%matrix, &
                                            jj, ii, loewdincharge_j)
                     WRITE (output_unit, *) 'loewdincharges,', loewdincharge_i, loewdincharge_j

                     CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%gmcharge(ispin)%matrix, ii, jj, tmp)
                     CALL cp_fm_set_element(stda_env%atom_param_set(iatom)%atom_param%gmcharge(ispin)%matrix, ii, jj, &
                                            tmp+loewdincharge_j*fr)
                     WRITE (output_unit, *) 'gmcharge', tmp+loewdincharge_j*fr, iatom, ii, jj

                     CALL cp_fm_get_element(stda_env%atom_param_set(jatom)%atom_param%gmcharge(ispin)%matrix, ii, jj, tmp)
                     CALL cp_fm_set_element(stda_env%atom_param_set(jatom)%atom_param%gmcharge(ispin)%matrix, ii, jj, &
                                            tmp+loewdincharge_i*fr)
                     WRITE (output_unit, *) 'gmcharge', tmp+loewdincharge_i*fr, jatom, ii, jj

!             gmcharge(jatom,ii,jj) = gmcharge(jatom,ii,jj)+loewdincharge_i*fr

!             WRITE(output_unit,*) 'iatom,', iatom, 'jatom,', jatom
!             WRITE(output_unit,*) 'dr', dr, 'rij', rij
!             WRITE(output_unit,*) 'alpha', alpha, 'fr', fr

!                IF (nmat > 1) THEN
!                   dfr = -2._dp*alpha*EXP(-alpha*alpha*dr*dr)*oorootpi/dr-fr/dr
!                   dfr = -dfr/dr
!                   DO i = 2, nmat
!                      gmcharge(iatom, i) = gmcharge(iatom, i)-rij(i-1)*loewdincharge_j*dfr
!                      gmcharge(jatom, i) = gmcharge(jatom, i)+rij(i-1)*loewdincharge_i*dfr
!                   END DO
!                END IF
                     !   IF (use_virial) THEN
                     !      IF (iatom == jatom) THEN
                     !         pfr = -0.5_dp*dfr*mcharge(iatom)*mcharge(jatom)
                     !      ELSE
                     !         pfr = -dfr*mcharge(iatom)*mcharge(jatom)
                     !      END IF
                     !      CALL virial_pair_force(virial%pv_virial, -pfr, rij, rij)
                     !      IF (atprop%stress) THEN
                     !         CALL virial_pair_force(atprop%atstress(:, :, iatom), -0.5_dp*pfr, rij, rij)
                     !         CALL virial_pair_force(atprop%atstress(:, :, jatom), -0.5_dp*pfr, rij, rij)
                     !      END IF
                     !   END IF
                  ENDDO !jj
               ENDDO !ii
            END IF

         END DO !neighbor_list
         CALL neighbor_list_iterator_release(nl_iterator)
      ENDDO !ispin

!TEST
!      DO ispin=1,nspins
!      CALL neighbor_list_iterator_create(nl_iterator, n_list)
!      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)  ! here only i /= j
!       CALL get_iterator_info(nl_iterator, iatom=iatom, jatom=jatom, r=rij)
!             dr = SQRT(SUM(rij(:)**2))
!             IF (dr > 1.e-6_dp) THEN
!       CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%gmcharge(ispin)%matrix,1,1,tmp)
!       wRITE(output_unit,*) 'tmp', tmp, iatom, jatom, ispin
!       CALL cp_fm_get_element(stda_env%atom_param_set(jatom)%atom_param%gmcharge(ispin)%matrix,1,1,tmp)
!       wRITE(output_unit,*) 'tmp', tmp, iatom, jatom, ispin
!             ENDIF
!      ENDDO
!          CALL neighbor_list_iterator_release(nl_iterator)
!      ENDDO
      CALL timestop(handle)

   END SUBROUTINE stda_ewald_overlap
! **************************************************************************************************
!> \brief Long-range Coulomb contribution via Ewald (just adjusted vom xtb_coulomb)
!> \param qs_env ...
!> \param p_env ...
!> \param stda_env ...
!> \param res ...
! **************************************************************************************************
   SUBROUTINE ewald_or_direct_sum(qs_env, p_env, stda_env, res)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(stda_env_type), POINTER                       :: stda_env
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: res

      CHARACTER(len=*), PARAMETER :: routineN = 'ewald_or_direct_sum', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: aa, ewald_type, handle, ia, iatom, ii, &
                                                            ikind, ispin, jatom, jj, natom, nsgf, &
                                                            nspins, output_unit, jkind, kkind, bb, &
                                                            res_etype, jspin
      INTEGER, DIMENSION(3)                              :: periodic
      INTEGER, DIMENSION(2)                              :: n_occ
      LOGICAL                                            :: do_ewald, do_mic
      REAL(KIND=dp) :: alpha, charge, deth, dr, gmcharge, loewdincharge_i, loewdincharge_j, tmp, &
         tmp1, tmp_1_r, tmp2, tmp_coul, tmpa, tmpb, xmat, xtransformed_double_i, &
         xtransformed_double_j, tmp_ex, tmpintern, tmp3, tmpintern2
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: r_pbc
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list, sab_tbe
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_iatom, basis_set_jatom
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(section_vals_type), POINTER                   :: input
      INTEGER, DIMENSION(:), POINTER                     :: kind_of

      CALL timeset(routineN, handle)
      CPASSERT(ASSOCIATED(stda_env))

      NULLIFY (logger) 
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      NULLIFY (dft_control)
      NULLIFY (cell)
      NULLIFY (particle_set)
      NULLIFY (qs_kind_set)
      NULLIFY(atomic_kind_set)
      NULLIFY(input)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      cell=cell, &
                      particle_set=particle_set,&
                      qs_kind_set=qs_kind_set,& 
                      atomic_kind_set=atomic_kind_set,&
                      input=input)

      natom = SIZE(particle_set)
      nspins = dft_control%nspins
      n_occ(:) = stda_env%n_occ(:)

      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)
      ALLOCATE (first_sgf(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf)
      ALLOCATE (kind_of(natom))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

      CALL section_vals_val_get(input, "DFT%TDDFPT%RES_ETYPE", i_val=res_etype)
!-------------------------------------------------------------------------------------------------------------
! 1/R contribution
!-------------------------------------------------------------------------------------------------------------

      do_ewald = dft_control%tddfpt_control%stda_control%do_ewald
      do_mic = dft_control%tddfpt_control%stda_control%do_mic
! EWALD SUM
      IF (do_ewald) THEN ! Ewald for Coulomb term
         ! Ewald sum
         NULLIFY (ewald_env, ewald_pw)
         CALL get_qs_env(qs_env=qs_env, ewald_env=ewald_env, ewald_pw=ewald_pw)
         CALL get_cell(cell=cell, periodic=periodic, deth=deth)
         CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
         CALL get_qs_env(qs_env=qs_env, sab_tbe=n_list)
! probably reuse routines from xtb, but for now using special stda routines !!! STILL NEEDS TO BE CHANGED
         CALL stda_ewald_overlap(qs_env, p_env, stda_env, alpha, n_list) ! q^point * phi^screened
         SELECT CASE (ewald_type)
         CASE DEFAULT
            CPABORT("Invalid Ewald type")
         CASE (do_ewald_none)
            CPABORT("Not allowed with sTDA")
         CASE (do_ewald_ewald)
            CPABORT("Standard Ewald not implemented in sTDA")
         CASE (do_ewald_pme)
            CPABORT("PME not implemented in sTDA")
         CASE (do_ewald_spme)
            CALL stda_spme_evaluate(qs_env, p_env, ewald_env, ewald_pw, particle_set, cell) ! q^total * phi^smooth
         END SELECT
!
         DO ispin = 1, nspins
            DO ii = 1, n_occ(ispin)
               DO aa = 1, nsgf

                  DO iatom = 1, natom
                     tmp_1_r = 0.0_dp
                     DO jj = 1, n_occ(ispin)
                        ! add self charge interaction and background charge contribution
                        CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%charges(ispin)%matrix, &
                                               jj, ii, charge)
                        CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%gmcharge(ispin)%matrix, &
                                               jj, ii, gmcharge)
                        gmcharge = gmcharge-2._dp*alpha*oorootpi*charge !self charge interaction ! q^point * phi^charge
                        IF (ANY(periodic(:) == 1)) THEN
                           gmcharge = gmcharge-pi/alpha**2/deth !*charge !background charge contribution! q^back * phi^screened
                        END IF
                        CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%xtransformed_double(ispin)%matrix, &
                                               aa, jj, xmat)
                        !    tmp_1_r = tmp_1_r + 0.5_dp*gmcharge*xmat
                        tmp_1_r = tmp_1_r-gmcharge*xmat  !! Coulomb contribution
!   energy = energy+0.5_dp*SUM(mcharge(:)*gmcharge(:, 1))
                     ENDDO !jj

                     !! ONLY correction to non-diagonal terms
                     CALL cp_fm_get_element(res(ispin)%matrix, aa, ii, tmp)
                     CALL cp_fm_set_element(res(ispin)%matrix, aa, ii, tmp+tmp_1_r)
                  ENDDO !iatom
               ENDDO !aa
            ENDDO ! ii

         ENDDO !ispin
      ENDIF !do_ewald

      IF ((.NOT.do_ewald).OR.(.NOT.do_mic)) THEN ! if not ewald: direct sum for coulomb part
                                                 ! if not mic:   direct sum for exchange part
! DIRECT SUM
         tmp_coul=0.0_dp
         tmp_ex=0.0_dp
            DO iatom = 1, natom
             DO jatom = 1, iatom-1 !natom
!                     CALL get_qs_env(qs_env=qs_env, local_particles=local_particles) !which particles are local to this processor
!                     DO ikind = 1, SIZE(local_particles%n_el)
!                        DO ia = 1, local_particles%n_el(ikind)
!                           iatom = local_particles%list(ikind)%array(ia)
                          WRITE(output_unit,*) 'iatom, jatom,', iatom, jatom

                          rij = particle_set(iatom)%r-particle_set(jatom)%r
                          rij = pbc(rij, cell)
                          dr = SQRT(SUM(rij(:)**2))

                          IF (dr > 1.e-6_dp) THEN
! only non-diagonal terms need to be corrected for
                            NULLIFY (basis_set_iatom)
                            NULLIFY (qs_kind)
                            kkind = kind_of(iatom) 
                            qs_kind => qs_kind_set(kkind)
                            CALL get_basis_from_container(qs_kind%basis_sets, basis_set=basis_set_iatom, basis_type="ORB")
!
                            NULLIFY (basis_set_jatom)
                            NULLIFY (qs_kind)
                            jkind = kind_of(jatom) 
                            qs_kind => qs_kind_set(jkind)
                            CALL get_basis_from_container(qs_kind%basis_sets, basis_set=basis_set_jatom, basis_type="ORB")

         DO ispin = 1, nspins
          DO jspin = 1, nspins
             DO ii = 1, n_occ(ispin)
              DO aa = 1, nsgf
                             tmp_1_r = 0.0_dp
                             tmpintern2=0.0_dp
                             tmpintern=0.0_dp
                             DO jj = 1, n_occ(jspin)
! coulomb-part
                              IF (ispin.EQ.jspin) THEN
                               CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%charges(ispin)%matrix, &
                                                        jj, ii, loewdincharge_i)
                               CALL cp_fm_get_element(stda_env%atom_param_set(jatom)%atom_param%charges(ispin)%matrix, &
                                                        jj, ii, loewdincharge_j)
                               CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%xtransformed_double(ispin)%matrix, &
                                                        aa, jj, xtransformed_double_i)
                               CALL cp_fm_get_element(stda_env%atom_param_set(jatom)%atom_param%xtransformed_double(ispin)%matrix, &
                                                        aa, jj, xtransformed_double_j)

                               tmp1 = xtransformed_double_i*loewdincharge_j
                               tmp2 = xtransformed_double_j*loewdincharge_i

                               tmp_1_r = tmp_1_r+(tmp1+tmp2)
                              ENDIF !ispin.eq.jspin
! exchange part
                               DO bb=first_sgf(jatom),first_sgf(jatom)-1+basis_set_jatom%nsgf
                                 CALL cp_fm_get_element(stda_env%ctransformed(jspin)%matrix,bb,jj,tmp1) 
                                 CALL cp_fm_get_element(stda_env%xtransformed(jspin)%matrix,bb,jj,tmp2)
                                 tmpintern=tmpintern + tmp1*tmp2
                               ENDDO
                               DO bb=first_sgf(iatom),first_sgf(iatom)-1+basis_set_iatom%nsgf
                                 CALL cp_fm_get_element(stda_env%ctransformed(jspin)%matrix,bb,jj,tmp1) 
                                 CALL cp_fm_get_element(stda_env%xtransformed(jspin)%matrix,bb,jj,tmp2)
                                 tmpintern2=tmpintern2 + tmp1*tmp2
                               ENDDO
                             ENDDO !jj
                             CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%ctransformed_double(ispin)%matrix,&
                                                    aa,ii,tmp3)
                             tmpintern=tmpintern*tmp3
                             CALL cp_fm_get_element(stda_env%atom_param_set(jatom)%atom_param%ctransformed_double(ispin)%matrix,&
                                                    aa,ii,tmp3)
                             tmpintern2=tmpintern2*tmp3

                            CALL cp_fm_get_element(res(ispin)%matrix, aa, ii, tmp)
! coulomb-part
                            IF (ispin.EQ.jspin) THEN
                             IF (.NOT.do_ewald) tmp_coul = -stda_env%hfx_fraction*tmp_1_r/dr !direct sum for coulomb part
                            ENDIF
! exchange-part
                            IF ((nspins == 1 .AND. res_etype == tddfpt_singlet ).OR.( nspins == 2 )) THEN
                             IF (.NOT.do_mic) tmp_ex= (tmpintern+ tmpintern2)*(2.0_dp/dr) !direct sum for exchange part
                            ENDIF

                            CALL cp_fm_set_element(res(ispin)%matrix,aa,ii,tmp+tmp_coul+tmp_ex)
               ENDDO !aa
            ENDDO !ii
          ENDDO !jspin
         ENDDO !ispin
                          END IF
                        END DO !jatom
                     END DO !iatom

      END IF !do_ewald , do_mic

     DEALLOCATE(kind_of,first_sgf)

      CALL timestop(handle)

   END SUBROUTINE ewald_or_direct_sum
END MODULE ewald_methods_stda

