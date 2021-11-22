!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief The types needed for the calculation of modified atomic orbitals (MAO)
!> \par History
!>      03.2016 created [JGH]
!> \author JGH
! **************************************************************************************************
MODULE mao_types
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_distribution_type,&
                                              dbcsr_p_type,&
                                              dbcsr_reserve_diag_blocks,&
                                              dbcsr_type_no_symmetry
   USE cp_dbcsr_operations,             ONLY: dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE kinds,                           ONLY: dp
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type,&
                                              release_neighbor_list_sets
#include "./base/base_uses.f90"
   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'mao_types'

   PUBLIC :: mao_type, create_mao_type, release_mao_type

! **************************************************************************************************
!> \brief Quantities needed for MAO's
!> \author JGH
! **************************************************************************************************
   TYPE mao_type
      LOGICAL                                            :: molecular
      INTEGER                                            :: max_iter
      REAL(KIND=dp)                                      :: eps_grad
      REAL(KIND=dp)                                      :: eps1_mao
      INTEGER                                            :: iolevel
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_q
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_smm, matrix_smo
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_p, matrix_s
      TYPE(dbcsr_p_type), DIMENSION(:,:), POINTER        :: matrix_ks, matrix_w
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: smm_list, smo_list
   END TYPE mao_type

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param mao_env ...
! **************************************************************************************************
   SUBROUTINE create_mao_type(mao_env, molecular, max_iter, eps_grad, eps1_mao, iolevel)
      TYPE(mao_type), POINTER                            :: mao_env
      LOGICAL, INTENT(IN), OPTIONAL                      :: molecular
      INTEGER, INTENT(IN), OPTIONAL                      :: max_iter
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps_grad
      REAL(KIND=dp), INTENT(IN), OPTIONAL                :: eps1_mao
      INTEGER, INTENT(IN), OPTIONAL                      :: iolevel

      IF (ASSOCIATED(mao_env)) THEN
         CALL release_mao_type(mao_env)
      END IF

      ALLOCATE(mao_env)

      mao_env%molecular = .FALSE.
      IF(PRESENT(molecular)) mao_env%molecular = molecular
      mao_env%max_iter = 0
      IF(PRESENT(max_iter)) mao_env%max_iter = max_iter
      mao_env%eps_grad = 0.00001_dp
      IF(PRESENT(eps_grad)) mao_env%eps_grad = eps_grad
      mao_env%eps1_mao = 1000.0_dp
      IF(PRESENT(eps1_mao)) mao_env%eps1_mao = eps1_mao
      mao_env%iolevel = 1
      IF(PRESENT(iolevel)) mao_env%iolevel = iolevel

      NULLIFY(mao_env%matrix_q)
      NULLIFY(mao_env%matrix_p)
      NULLIFY(mao_env%matrix_s)
      NULLIFY(mao_env%matrix_ks)
      NULLIFY(mao_env%matrix_w)
      NULLIFY(mao_env%matrix_smm)
      NULLIFY(mao_env%matrix_smo)
      NULLIFY(mao_env%smm_list)
      NULLIFY(mao_env%smo_list)

   END SUBROUTINE create_mao_type

! **************************************************************************************************
!> \brief ...
!> \param mao_env ...
! **************************************************************************************************
   SUBROUTINE release_mao_type(mao_env)
      TYPE(mao_type), POINTER                            :: mao_env

      IF (ASSOCIATED(mao_env)) THEN
         ! Deallocate the matrix structures
         IF (ASSOCIATED(mao_env%matrix_smm)) CALL dbcsr_deallocate_matrix_set(mao_env%matrix_smm)
         IF (ASSOCIATED(mao_env%matrix_smo)) CALL dbcsr_deallocate_matrix_set(mao_env%matrix_smo)
         IF (ASSOCIATED(mao_env%matrix_q)) CALL dbcsr_deallocate_matrix_set(mao_env%matrix_q)
         IF (ASSOCIATED(mao_env%matrix_p)) CALL dbcsr_deallocate_matrix_set(mao_env%matrix_p)
         IF (ASSOCIATED(mao_env%matrix_s)) CALL dbcsr_deallocate_matrix_set(mao_env%matrix_s)
         IF (ASSOCIATED(mao_env%matrix_ks)) CALL dbcsr_deallocate_matrix_set(mao_env%matrix_ks)
         IF (ASSOCIATED(mao_env%matrix_w)) CALL dbcsr_deallocate_matrix_set(mao_env%matrix_w)
         ! Deallocate the neighbor list structures
         IF(ASSOCIATED(mao_env%matrix_smm)) THEN
            CALL release_neighbor_list_sets(mao_env%smm_list)
         END IF
         IF(ASSOCIATED(mao_env%matrix_smo)) THEN
            CALL release_neighbor_list_sets(mao_env%smo_list)
         END IF
         DEALLOCATE (mao_env)
      END IF

   END SUBROUTINE release_mao_type

! **************************************************************************************************

END MODULE mao_types
