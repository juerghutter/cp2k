!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate MAO contributions to Lagrangians and forces
!> \par History
!>      11.2021 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE mao_forces
   USE atomic_kind_types,               ONLY: get_atomic_kind
   USE basis_set_types,                 ONLY: gto_basis_set_p_type,&
                                              gto_basis_set_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_distribution_type,&
                                              dbcsr_norm_maxabsnorm,&
                                              dbcsr_set,&
                                              dbcsr_copy,&
                                              dbcsr_norm,&
                                              dbcsr_add,&
                                              dbcsr_dot,&
                                              dbcsr_multiply,&
                                              dbcsr_type,&
                                              dbcsr_p_type,&
                                              dbcsr_release,&
                                              dbcsr_reserve_diag_blocks,&
                                              dbcsr_type_no_symmetry
   USE kinds,                           ONLY: dp
   USE mao_types,                       ONLY: mao_type
   USE mao_methods,                     ONLY: mao_build_q, mao_project_gradient
   USE mao_optimizer,                   ONLY: mao_optimize
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              qs_kind_type
   USE qs_ks_types,                     ONLY: get_ks_env,&
                                              qs_ks_env_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
   USE qs_neighbor_lists,               ONLY: setup_neighbor_list
   USE qs_overlap,                      ONLY: build_overlap_matrix_simple
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mao_forces'

   PUBLIC ::  mao_lagrangian, mao_rhs_zvector, mao_basis_forces

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mao_env ...
!> \param mao_coef ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE mao_lagrangian(qs_env, mao_env, mao_coef, mao_lagr, ksmat, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mao_type), POINTER                            :: mao_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coef
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_lagr
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: ksmat
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      CHARACTER(len=*), PARAMETER                        :: routineN = 'mao_lagrangian'

      INTEGER                                            :: handle, iw, iolev, nspins, ispin, iter
      REAL(KIND=dp)                                      :: rrk, rrk1, dzk, alpha, beta, rnorm
      TYPE(dbcsr_type)                                   :: cgmat, rhsmat, emat, scrmat, dmat, smat
      TYPE(dbcsr_type), POINTER                          :: matrix_q, matrix_s

      CALL timeset(routineN, handle)

      MARK_USED(qs_env)
      MARK_USED(mao_coef)
      MARK_USED(mao_lagr)
      MARK_USED(ksmat)

      ! output
      IF (PRESENT(unit_nr)) THEN
         iw = unit_nr
      ELSE
         iw = -1
      END IF
      iolev = mao_env%iolevel
      IF (iolev == 0) iw = -1

      nspins = SIZE(mao_coef)
      CALL dbcsr_allocate_matrix_set(mao_lagr, nspins)
      DO ispin = 1, nspins
         ALLOCATE (mao_lagr(ispin)%matrix)
         CALL dbcsr_create(mao_lagr(ispin)%matrix, name="MAO LAGR", &
                           template=mao_coef(ispin)%matrix)
         CALL dbcsr_reserve_diag_blocks(matrix=mao_lagr(ispin)%matrix)
      END DO

      DO ispin = 1, nspins
         CALL dbcsr_set(mao_lagr(ispin)%matrix, 0.0_dp)
         CALL dbcsr_create(cgmat, template=mao_coef(ispin)%matrix)
         CALL dbcsr_create(rhsmat, template=mao_coef(ispin)%matrix)
         CALL dbcsr_reserve_diag_blocks(matrix=rhsmat)
         CALL dbcsr_multiply("N", "N", 1.0_dp, ksmat(ispin, 1)%matrix, mao_coef(ispin)%matrix, &
                             0.0_dp, cgmat)
         CALL dbcsr_multiply("N", "N", -2.0_dp, cgmat, mao_env%matrix_p(ispin, 1)%matrix, &
                             0.0_dp, rhsmat, retain_sparsity=.TRUE.)
         matrix_s => mao_env%matrix_smm(ispin)%matrix
         matrix_q => mao_env%matrix_q(ispin)%matrix
         ! emat = T*Q*T(T)
         IF(mao_env%max_iter < 1) THEN
            CALL dbcsr_create(emat, template=mao_env%matrix_s(1,1)%matrix)
            CALL dbcsr_reserve_diag_blocks(matrix=emat)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_q, mao_coef(ispin)%matrix, 0.0_dp, cgmat)
            CALL dbcsr_multiply("T", "N", 1.0_dp, mao_coef(ispin)%matrix, cgmat, 0.0_dp, emat, retain_sparsity=.TRUE.)
         ELSE
            CALL dbcsr_create(emat, template=mao_env%matrix_s(1,1)%matrix)
            CALL dbcsr_multiply("N", "N", 1.0_dp, matrix_q, mao_coef(ispin)%matrix, 0.0_dp, cgmat)
            CALL dbcsr_multiply("T", "N", 1.0_dp, mao_coef(ispin)%matrix, cgmat, 0.0_dp, emat)
         END IF
         CALL dbcsr_release(cgmat)
         ! CG method
         CALL dbcsr_create(scrmat, template=mao_coef(ispin)%matrix)
         CALL dbcsr_reserve_diag_blocks(matrix=scrmat)
         CALL dbcsr_create(dmat, template=mao_coef(ispin)%matrix)
         CALL dbcsr_reserve_diag_blocks(matrix=dmat)
         CALL dbcsr_create(cgmat, template=mao_coef(ispin)%matrix)
         CALL dbcsr_reserve_diag_blocks(matrix=cgmat)
         IF(mao_env%max_iter < 1) THEN
            CALL dbcsr_create(smat, template=matrix_s)
            CALL dbcsr_reserve_diag_blocks(matrix=smat)
            CALL dbcsr_copy(smat, matrix_s, keep_sparsity=.TRUE.)
            CALL mao_project_gradient(mao_coef(ispin)%matrix, rhsmat, smat)
            CALL dbcsr_copy(dmat, rhsmat, keep_sparsity=.TRUE.)
            DO iter = 1, 250
               ! A*d
               CALL dbcsr_multiply("N", "N", 2.0_dp, matrix_q, dmat, &
                                   0.0_dp, cgmat, retain_sparsity=.TRUE.)
               CALL dbcsr_multiply("N", "N", 1.0_dp, dmat, emat, &
                                   0.0_dp, scrmat, retain_sparsity=.TRUE.)
               CALL dbcsr_multiply("N", "N", -2.0_dp, mao_env%matrix_smm(1)%matrix, scrmat, &
                                   1.0_dp, cgmat, retain_sparsity=.TRUE.)
               CALL mao_project_gradient(mao_coef(ispin)%matrix, cgmat, smat)
               ! alpha
               call dbcsr_dot(rhsmat,rhsmat,rrk)
               call dbcsr_dot(dmat,cgmat,dzk)
               alpha = rrk/dzk
               !
               CALL dbcsr_add(mao_lagr(ispin)%matrix, dmat, 1.0_dp, alpha)
               CALL dbcsr_add(rhsmat, cgmat, 1.0_dp, -alpha)
               CALL dbcsr_norm(rhsmat, dbcsr_norm_maxabsnorm, norm_scalar=rnorm)
               !
               IF(iw > 0) THEN
                   WRITE(iw,*) iter,rnorm
               END IF
               !
               call dbcsr_dot(rhsmat,rhsmat,rrk1)
               beta = rrk1/rrk
               !!call dbcsr_dot(rhsmat,cgmat,rrk1)
               !!beta = rrk1/dzk
               CALL dbcsr_add(dmat, rhsmat, beta, 1.0_dp)
            END DO
            CALL mao_project_gradient(mao_coef(ispin)%matrix, mao_lagr(ispin)%matrix, smat)
            CALL dbcsr_release(smat)
         ELSE
            stop
         END IF
         CALL dbcsr_release(dmat)
         CALL dbcsr_release(scrmat)
         CALL dbcsr_release(emat)
         CALL dbcsr_release(cgmat)
         CALL dbcsr_release(rhsmat)
      END DO

      CALL timestop(handle)

   END SUBROUTINE mao_lagrangian

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mao_env ...
!> \param mao_coef ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE mao_rhs_zvector(qs_env, mao_env, mao_coef, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mao_type), POINTER                            :: mao_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coef
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      CHARACTER(len=*), PARAMETER                        :: routineN = 'mao_rhs_zvector'

      INTEGER                                            :: handle, iw, iolev

      CALL timeset(routineN, handle)

      MARK_USED(qs_env)
      MARK_USED(mao_coef)

      ! output
      IF (PRESENT(unit_nr)) THEN
         iw = unit_nr
      ELSE
         iw = -1
      END IF
      iolev = mao_env%iolevel
      IF (iolev == 0) iw = -1

      CALL timestop(handle)

   END SUBROUTINE mao_rhs_zvector

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param mao_env ...
!> \param mao_coef ...
!> \param unit_nr ...
! **************************************************************************************************
   SUBROUTINE mao_basis_forces(qs_env, mao_env, mao_coef, unit_nr)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(mao_type), POINTER                            :: mao_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: mao_coef
      INTEGER, INTENT(IN), OPTIONAL                      :: unit_nr

      CHARACTER(len=*), PARAMETER                        :: routineN = 'mao_basis_forces'

      INTEGER                                            :: handle, iw, iolev

      CALL timeset(routineN, handle)

      MARK_USED(qs_env)
      MARK_USED(mao_coef)

      ! output
      IF (PRESENT(unit_nr)) THEN
         iw = unit_nr
      ELSE
         iw = -1
      END IF
      iolev = mao_env%iolevel
      IF (iolev == 0) iw = -1

      CALL timestop(handle)

   END SUBROUTINE mao_basis_forces

END MODULE mao_forces
