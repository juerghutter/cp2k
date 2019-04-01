!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Simplified Tamm Dancoff approach (sTDA).
! **************************************************************************************************
MODULE qs_tddfpt_stda_utils

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind_set
   USE atprop_types,                    ONLY: atprop_type
   USE cell_types,                      ONLY: cell_type,&
                                              get_cell,&
                                              pbc
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_plus_fm_fm_t,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_row_scale,&
                                              cp_fm_schur_product
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_to_fm,&
                                              cp_fm_type,&
                                              cp_fm_vectorssum
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add_on_diag, dbcsr_create, dbcsr_distribution_type, dbcsr_filter, dbcsr_finalize, &
        dbcsr_get_block_p, dbcsr_iterator_blocks_left, dbcsr_iterator_next_block, &
        dbcsr_iterator_start, dbcsr_iterator_stop, dbcsr_iterator_type, dbcsr_p_type, &
        dbcsr_release, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry, dbcsr_type_symmetric
   USE distribution_1d_types,           ONLY: distribution_1d_type
   USE ewald_environment_types,         ONLY: ewald_env_get,&
                                              ewald_environment_type
   USE ewald_methods_tb,                ONLY: tb_ewald_overlap,&
                                              tb_spme_evaluate
   USE ewald_pw_types,                  ONLY: ewald_pw_type
   USE input_constants,                 ONLY: tddfpt_singlet
   USE input_section_types,             ONLY: section_vals_type,&
                                              section_vals_val_get
   USE iterate_matrix,                  ONLY: matrix_sqrt_Newton_Schulz
   USE kinds,                           ONLY: dp
   USE mathconstants,                   ONLY: oorootpi
   USE message_passing,                 ONLY: mp_sum
   USE particle_methods,                ONLY: get_particle_set
   USE particle_types,                  ONLY: particle_type
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
   USE scf_control_types,               ONLY: scf_control_type
   USE virial_types,                    ONLY: virial_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt_stda_utils'

   PUBLIC:: stda_init_matrices, stda_calculate_kernel

CONTAINS
! **************************************************************************************************
!> \brief Calculate sTDA matrices (Coulomb- and exchange-type) and Lowdin MO coefficients and
!>        transition density charges; initialization before entering the Davidson iterations
!> \param qs_env ...
!> \param p_env ...
! **************************************************************************************************
   SUBROUTINE stda_init_matrices(qs_env, p_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env

      CHARACTER(len=*), PARAMETER :: routineN = 'stda_init_matrices', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle
      TYPE(stda_env_type), POINTER                       :: stda_env

      CALL timeset(routineN, handle)

      NULLIFY (stda_env)
      stda_env => p_env%stda_env

      ! calculate exchange gamma matrix
      CALL setup_gamma(qs_env, stda_env, stda_env%gamma_exchange)

      ! calculate S_half and Lowdin MO coefficients
      CALL get_lowdin_mo_coefficients(qs_env, p_env, stda_env)

      CALL timestop(handle)

   END SUBROUTINE stda_init_matrices
! **************************************************************************************************
!> \brief Calculate sTDA exchange-type contributions
!> \param qs_env ...
!> \param stda_env ...
!> \param gamma_matrix sTDA exchange-type contributions
!> \note  Note the specific sTDA notation exchange-type integrals (ia|jb) refer to Coulomb interaction
! **************************************************************************************************
   SUBROUTINE setup_gamma(qs_env, stda_env, gamma_matrix)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type), POINTER                       :: stda_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gamma_matrix

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_gamma', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: rsmooth = 1.0_dp

      CHARACTER(len=20)                                  :: cstring
      INTEGER                                            :: handle, iatom, icol, ikind, irow, jatom, &
                                                            jkind, natom, nmat
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes
      LOGICAL                                            :: found
      REAL(KIND=dp)                                      :: dr, eta, fcut, r, ralpha, rcut, rcuta, &
                                                            rcutb, x
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gblock
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      cstring = "GAMMA EXCHANGE MATRIX"
      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, dbcsr_dist=dbcsr_dist, &
                      sab_orb=n_list, cell=cell)
      natom = SIZE(particle_set)

      nmat = 1
      NULLIFY (gamma_matrix)
      CALL dbcsr_allocate_matrix_set(gamma_matrix, nmat)

      ALLOCATE (row_blk_sizes(natom))
      row_blk_sizes(1:natom) = 1
      ALLOCATE (gamma_matrix(1)%matrix)

      CALL dbcsr_create(gamma_matrix(1)%matrix, name=TRIM(ADJUSTL(cstring)), dist=dbcsr_dist, &
                        matrix_type=dbcsr_type_symmetric, row_blk_size=row_blk_sizes, &
                        col_blk_size=row_blk_sizes, nze=0, mutable_work=.TRUE.)

      DEALLOCATE (row_blk_sizes)

      ! setup the matrices using the neighbor list
      CALL cp_dbcsr_alloc_block_from_nbl(gamma_matrix(1)%matrix, n_list)
      CALL dbcsr_set(gamma_matrix(1)%matrix, 0.0_dp)

      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, n_list)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij)

         dr = SQRT(SUM(rij(:)**2)) ! interatomic distance

         eta = (stda_env%kind_param_set(ikind)%kind_param%hardness_param+ &
                stda_env%kind_param_set(jkind)%kind_param%hardness_param)/2.0_dp

         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)

         NULLIFY (gblock)
         CALL dbcsr_get_block_p(matrix=gamma_matrix(1)%matrix, &
                                row=irow, col=icol, BLOCK=gblock, found=found)
         CPASSERT(found)

         ! get rcuta and rcutb
         rcuta = stda_env%kind_param_set(ikind)%kind_param%rcut
         rcutb = stda_env%kind_param_set(jkind)%kind_param%rcut
         rcut = rcuta+rcutb

         !>   Computes the short-range gamma parameter from
         !>   Nataga-Mishimoto-Ohno-Klopman formula equivalently as it is done for xTB
         IF (dr < 1.e-6) THEN
            ! on site terms
            gblock(:, :) = gblock(:, :)+eta
         ELSEIF (dr > rcut) THEN
            ! do nothing
         ELSE
            IF (dr < rcut-rsmooth) THEN
               fcut = 1.0_dp
            ELSE
               r = dr-(rcut-rsmooth)
               x = r/rsmooth
               fcut = -10._dp*x**5+15._dp*x**4-10._dp*x**3+1._dp
            END IF
            ralpha = dr**stda_env%alpha_param
            gblock(:, :) = gblock(:, :)+ &
                           fcut*(1._dp/(dr**(stda_env%alpha_param)+eta**(-stda_env%alpha_param))) &
                           **(1._dp/stda_env%alpha_param)-fcut/dr
         END IF

      END DO

      CALL neighbor_list_iterator_release(nl_iterator)
      CALL dbcsr_finalize(gamma_matrix(1)%matrix)

      CALL timestop(handle)

   END SUBROUTINE setup_gamma
! **************************************************************************************************
!> \brief Calculate Lowdin MO coefficients
!> \param qs_env ...
!> \param p_env ...
!> \param stda_env ...
! **************************************************************************************************
   SUBROUTINE get_lowdin_mo_coefficients(qs_env, p_env, stda_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(stda_env_type), POINTER                       :: stda_env

      CHARACTER(len=*), PARAMETER :: routineN = 'get_lowdin_mo_coefficients', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iounit, ispin, max_iter_lanczos, &
                                                            n_occ, ndep, nsgf, nspins, &
                                                            order_lanczos
      LOGICAL                                            :: converged
      REAL(KIND=dp)                                      :: eps_lanczos, threshold
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: fm_s_half, fm_work1, fm_work2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrixkp_s
      TYPE(dbcsr_type)                                   :: sm_hinv
      TYPE(dbcsr_type), POINTER                          :: sm_h, sm_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (logger) !get output_unit
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CPASSERT(ASSOCIATED(stda_env))

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins

      NULLIFY (qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, &
                      qs_kind_set=qs_kind_set)
      ! Get the total number of contracted spherical Gaussian basis functions
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

      NULLIFY (sm_s, matrixkp_s)
      NULLIFY (fmstruct)
      NULLIFY (fm_work1, fm_s_half)

      NULLIFY (para_env)
      NULLIFY (blacs_env)
      NULLIFY (scf_control)

      CALL get_qs_env(qs_env=qs_env, &
                      matrix_s_kp=matrixkp_s, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      scf_control=scf_control)

      CPASSERT(ASSOCIATED(matrixkp_s))
      CPASSERT(ASSOCIATED(p_env%psi0d))
      IF (SIZE(matrixkp_s, 2) > 1) CPWARN("not implemented for k-points.")

      ! Calculate S^1/2 matrix
      IF (iounit > 0) THEN
         WRITE (iounit, "(T3,A)") "STDA | Create Matrix SQRT(S)"
      END IF
      sm_s => matrixkp_s(1, 1)%matrix
      sm_h => p_env%stda_env%shalf
      CALL dbcsr_create(sm_hinv, template=matrixkp_s(1, 1)%matrix)
      CALL dbcsr_add_on_diag(sm_h, 1.0_dp)
      threshold = 1.0e-8_dp
      order_lanczos = 3
      eps_lanczos = 1.0e-4_dp
      max_iter_lanczos = 40
      CALL matrix_sqrt_Newton_Schulz(sm_h, sm_hinv, sm_s, &
                                     threshold, order_lanczos, eps_lanczos, max_iter_lanczos, converged)
      CALL dbcsr_release(sm_hinv)
      IF (.NOT. converged) THEN
         IF (iounit > 0) THEN
            WRITE (iounit, "(T3,A)") "STDA | Newton-Schulz iteration did not converge"
            WRITE (iounit, "(T3,A)") "STDA | Calculate SQRT(S) from diagonalization"
         END IF
         ! Provide full size work matrices
         CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                  para_env=para_env, &
                                  context=blacs_env, &
                                  nrow_global=nsgf, &
                                  ncol_global=nsgf)
         CALL cp_fm_create(matrix=fm_s_half, matrix_struct=fmstruct, name="S^(1/2) MATRIX")
         CALL cp_fm_create(matrix=fm_work1, matrix_struct=fmstruct, name="TMP MATRIX")
         CALL cp_fm_struct_release(fmstruct=fmstruct)
         CALL copy_dbcsr_to_fm(sm_s, fm_s_half)
         CALL cp_fm_power(fm_s_half, fm_work1, 0.5_dp, scf_control%eps_eigval, ndep)
         IF (ndep /= 0) &
            CALL cp_warn(__LOCATION__, &
                         "Overlap matrix exhibits linear dependencies. At least some "// &
                         "eigenvalues have been quenched.")
         CALL copy_fm_to_dbcsr(fm_s_half, sm_h)
         CALL cp_fm_release(matrix=fm_s_half)
         CALL cp_fm_release(matrix=fm_work1)
         IF (iounit > 0) WRITE (iounit, *)
      END IF

      DO ispin = 1, nspins
         n_occ = stda_env%n_occ(ispin) !p_env%n_mo(ispin)
         ! Build Lowdin coefficients C'= S^1/2 C for each spin
         NULLIFY (fm_work2)
         CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                  para_env=para_env, &
                                  context=blacs_env, &
                                  nrow_global=nsgf, &
                                  ncol_global=n_occ)
         CALL cp_fm_create(matrix=fm_work2, matrix_struct=fmstruct, name="TMP MATRIX 2")
         CALL cp_dbcsr_sm_fm_multiply(p_env%stda_env%shalf, p_env%psi0d(ispin)%matrix, &
                                      fm_work2, n_occ, alpha=1.0_dp, beta=0.0_dp)

         CALL cp_fm_to_fm(fm_work2, stda_env%ctransformed(ispin)%matrix)

         CALL cp_fm_struct_release(fmstruct=fmstruct)
         CALL cp_fm_release(matrix=fm_work2)
         NULLIFY (fm_work2, fmstruct)

      ENDDO !nspins

      CALL timestop(handle)

   END SUBROUTINE get_lowdin_mo_coefficients

! **************************************************************************************************
!> \brief Calculate Lowdin transformed Davidson trial vector X
!> \param qs_env ...
!> \param p_env ...
!> \param X ...
! **************************************************************************************************
   SUBROUTINE get_lowdin_x(qs_env, p_env, X)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: X

      CHARACTER(len=*), PARAMETER :: routineN = 'get_lowdin_x', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nsgf, nspins
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: fm_work2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(stda_env_type), POINTER                       :: stda_env

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(p_env%stda_env))
      NULLIFY (stda_env)
      stda_env => p_env%stda_env

      NULLIFY (dft_control)
      NULLIFY (para_env)
      NULLIFY (blacs_env)
      NULLIFY (qs_kind_set)
      CALL get_qs_env(qs_env, dft_control=dft_control, &
                      para_env=para_env, &
                      blacs_env=blacs_env, &
                      qs_kind_set=qs_kind_set)
      nspins = dft_control%nspins
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

      ! Build Lowdin transformed tilde(X)= S^1/2 X for each spin
      DO ispin = 1, nspins
         NULLIFY (fm_work2, fmstruct)
         CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                  para_env=para_env, &
                                  context=blacs_env, &
                                  nrow_global=nsgf, &
                                  ncol_global=stda_env%n_occ(ispin))
         CALL cp_fm_create(matrix=fm_work2, &
                           matrix_struct=fmstruct, &
                           name="TMP MATRIX 2")
         CALL cp_fm_struct_release(fmstruct=fmstruct)

         CALL cp_dbcsr_sm_fm_multiply(p_env%stda_env%shalf, X(ispin)%matrix, &
                                      fm_work2, stda_env%n_occ(ispin), alpha=1.0_dp, beta=0.0_dp)

         CALL cp_fm_to_fm(fm_work2, stda_env%xtransformed(ispin)%matrix)

         CALL cp_fm_release(matrix=fm_work2)
         NULLIFY (fm_work2, fmstruct)
      ENDDO !nspins

      CALL timestop(handle)

   END SUBROUTINE get_lowdin_x

! **************************************************************************************************
!> \brief ...Calculate the sTDA kernel contribution by contracting the Lowdin MO coefficients --
!>           transition charges with the Coulomb-type or exchange-type integrals
!> \param p_env ...
!> \param qs_env ...
!> \param X ...
!> \param res ... vector AX with A being the sTDA matrix and X the Davidson trial vector of the
!>                eigenvalue problem AX=omegaX
! **************************************************************************************************
   SUBROUTINE stda_calculate_kernel(p_env, qs_env, X, res)

      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(IN)       :: X
      TYPE(cp_fm_p_type), DIMENSION(:), INTENT(INOUT)    :: res

      CHARACTER(len=*), PARAMETER :: routineN = 'stda_calculate_kernel', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: blk, ewald_type, handle, ia, iatom, &
                                                            ikind, iounit, is, ispin, jatom, &
                                                            jkind, natom, nsgf, nspins, res_etype
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, kind_of, last_sgf
      INTEGER, DIMENSION(2)                              :: n_occ
      INTEGER, DIMENSION(3)                              :: periodic
      LOGICAL                                            :: calculate_forces, do_coulomb, do_ewald, &
                                                            use_virial
      REAL(KIND=dp)                                      :: alpha, bp, deth, dr, eta, gabr, hfx, &
                                                            rbeta, spinfac
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: tcharge, tv
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: gtcharge
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gab, pblock
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atprop_type), POINTER                         :: atprop
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: ct, cvec, xt
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_iterator_type)                          :: iter
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dbcsr_type)                                   :: pdens
      TYPE(dbcsr_type), POINTER                          :: tempmat
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(distribution_1d_type), POINTER                :: local_particles
      TYPE(ewald_environment_type), POINTER              :: ewald_env
      TYPE(ewald_pw_type), POINTER                       :: ewald_pw
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: n_list
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(section_vals_type), POINTER                   :: input
      TYPE(stda_env_type), POINTER                       :: stda_env
      TYPE(virial_type), POINTER                         :: virial

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(p_env))
      CPASSERT(ASSOCIATED(p_env%stda_env))

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      stda_env => p_env%stda_env
      n_occ(:) = stda_env%n_occ(:)
      nspins = SIZE(X)
      spinfac = 2.0_dp
      IF (nspins == 2) spinfac = 1.0_dp

      NULLIFY (input)
      CALL get_qs_env(qs_env=qs_env, input=input, para_env=para_env)
      CALL section_vals_val_get(input, "DFT%TDDFPT%RES_ETYPE", i_val=res_etype)
      IF ((nspins == 1 .AND. res_etype == tddfpt_singlet) .OR. (nspins == 2)) THEN
         do_coulomb = .TRUE.
      ELSE
         do_coulomb = .FALSE.
      END IF

      CALL get_qs_env(qs_env, dft_control=dft_control)
      do_ewald = dft_control%tddfpt_control%stda_control%do_ewald

      CALL get_qs_env(qs_env, natom=natom, cell=cell, &
                      particle_set=particle_set, qs_kind_set=qs_kind_set)
      ALLOCATE (first_sgf(natom))
      ALLOCATE (last_sgf(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)

      ! calculate Loewdin transformed Davidson trial vector tilde(X)=S^1/2*X
      ! and tilde(tilde(X))=S^1/2_A*tilde(X)_A
      CALL get_lowdin_x(qs_env, p_env, X)

      ALLOCATE (tcharge(natom), gtcharge(natom, 1))

      DO ispin = 1, nspins
         ct => stda_env%ctransformed(ispin)%matrix
         xt => stda_env%xtransformed(ispin)%matrix
         CALL cp_fm_get_info(ct, matrix_struct=fmstruct, nrow_global=nsgf)
         ALLOCATE (tv(nsgf))
         CALL cp_fm_create(cvec, fmstruct)
         !
         ! *** Coulomb contribution
         !
         IF (do_coulomb) THEN
            ! CV(mu,j) = CT(mu,j)*XT(mu,j)
            CALL cp_fm_schur_product(ct, xt, cvec)
            ! TV(mu) = SUM_j CV(mu,j)
            CALL cp_fm_vectorssum(cvec, tv, "R")
            ! contract charges
            ! TC(a) = SUM_(mu of a) TV(mu)
            DO ia = 1, natom
               tcharge(ia) = 0.0_dp
               DO is = first_sgf(ia), last_sgf(ia)
                  tcharge(ia) = tcharge(ia)+tv(is)
               END DO
            END DO
            ! Apply tcharge*gab -> gtcharge
            ! gT(b) = SUM_a g(a,b)*TC(a)
            ! gab = stda_env%gamma_exchange(1)%matrix
            gtcharge = 0.0_dp
            ! short range contribution
            tempmat => stda_env%gamma_exchange(1)%matrix
            CALL dbcsr_iterator_start(iter, tempmat)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iatom, jatom, gab, blk)
               gtcharge(iatom, 1) = gtcharge(iatom, 1)+gab(1, 1)*tcharge(jatom)
               IF (iatom /= jatom) THEN
                  gtcharge(jatom, 1) = gtcharge(jatom, 1)+gab(1, 1)*tcharge(iatom)
               END IF
            END DO
            CALL dbcsr_iterator_stop(iter)
            ! Ewald long range contribution
            IF (do_ewald) THEN
               NULLIFY (ewald_env, ewald_pw)
               CALL get_qs_env(qs_env=qs_env, ewald_env=ewald_env, ewald_pw=ewald_pw)
               CALL get_cell(cell=cell, periodic=periodic, deth=deth)
               CALL ewald_env_get(ewald_env, alpha=alpha, ewald_type=ewald_type)
               CALL get_qs_env(qs_env=qs_env, sab_tbe=n_list, &
                               virial=virial, atprop=atprop)
               use_virial = .FALSE.
               calculate_forces = .FALSE.
               CALL tb_ewald_overlap(gtcharge, tcharge, alpha, n_list, virial, use_virial, atprop)
               CALL tb_spme_evaluate(ewald_env, ewald_pw, particle_set, cell, &
                                     gtcharge, tcharge, calculate_forces, virial, use_virial, atprop)
               ! add self charge interaction contribution
               IF (para_env%ionode) THEN
                  gtcharge(:, 1) = gtcharge(:, 1)-2._dp*alpha*oorootpi*tcharge(:)
               END IF
            ELSE
               CALL get_qs_env(qs_env=qs_env, local_particles=local_particles)
               DO ikind = 1, SIZE(local_particles%n_el)
                  DO ia = 1, local_particles%n_el(ikind)
                     iatom = local_particles%list(ikind)%array(ia)
                     DO jatom = 1, iatom-1
                        rij = particle_set(iatom)%r-particle_set(jatom)%r
                        rij = pbc(rij, cell)
                        dr = SQRT(SUM(rij(:)**2))
                        IF (dr > 1.e-6_dp) THEN
                           gtcharge(iatom, 1) = gtcharge(iatom, 1)+tcharge(jatom)/dr
                           gtcharge(jatom, 1) = gtcharge(jatom, 1)+tcharge(iatom)/dr
                        END IF
                     END DO
                  END DO
               END DO
            END IF
            CALL mp_sum(gtcharge, para_env%group)
            ! expand charges
            ! TV(mu) = TC(a of mu)
            tv(1:nsgf) = 0.0_dp
            DO ia = 1, natom
               DO is = first_sgf(ia), last_sgf(ia)
                  tv(is) = gtcharge(ia, 1)
               END DO
            END DO
            ! CV(mu,i) = TV(mu)*CV(mu,i)
            CALL cp_fm_to_fm(ct, cvec)
            CALL cp_fm_row_scale(cvec, tv)
            ! rho(nu,i) = rho(nu,i) + Shalf(nu,mu)*CV(mu,i)
            CALL cp_dbcsr_sm_fm_multiply(p_env%stda_env%shalf, cvec, res(ispin)%matrix, &
                                         n_occ(ispin), spinfac, 1.0_dp)
         END IF
         !
         ! *** Exchange contribution
         !
         IF (stda_env%hfx_fraction > 0.0_dp) THEN
            ALLOCATE (kind_of(natom))
            CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
            CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)
            !
            CALL get_qs_env(qs_env=qs_env, matrix_s=matrix_s)
            tempmat => matrix_s(1)%matrix
            CALL dbcsr_create(pdens, template=tempmat, matrix_type=dbcsr_type_no_symmetry)
            ! P(nu,mu) = SUM_j XT(nu,j)*CT(mu,j)
            CALL cp_dbcsr_plus_fm_fm_t(pdens, xt, ct, n_occ(ispin), &
                                       1.0_dp, keep_sparsity=.FALSE.)
            CALL dbcsr_filter(pdens, stda_env%eps_td_filter)
!           pdocc = dbcsr_get_occupation(pdens)
!           IF (iounit > 0)THEN
!              WRITE (iounit, '(T2,A,T65,F16.12)') " Transition density matrix sparsity: ", pdocc
!           END IF
            ! Apply PP*gab -> PP; gab = gamma_coulomb
            ! P(nu,mu) = P(nu,mu)*g(a of nu,b of mu)
            bp = stda_env%beta_param
            hfx = stda_env%hfx_fraction
            CALL dbcsr_iterator_start(iter, pdens)
            DO WHILE (dbcsr_iterator_blocks_left(iter))
               CALL dbcsr_iterator_next_block(iter, iatom, jatom, pblock, blk)
               rij = particle_set(iatom)%r-particle_set(jatom)%r
               rij = pbc(rij, cell)
               dr = SQRT(SUM(rij(:)**2))
               ikind = kind_of(iatom)
               jkind = kind_of(jatom)
               eta = (stda_env%kind_param_set(ikind)%kind_param%hardness_param+ &
                      stda_env%kind_param_set(jkind)%kind_param%hardness_param)/2.0_dp
               rbeta = dr**bp
               gabr = (1._dp/(rbeta+(hfx*eta)**(-bp)))**(1._dp/bp)
               pblock = gabr*pblock
            END DO
            CALL dbcsr_iterator_stop(iter)
            ! CV(mu,i) = P(nu,mu)*CT(mu,i)
            CALL cp_dbcsr_sm_fm_multiply(pdens, ct, cvec, &
                                         n_occ(ispin), 1.0_dp, 0.0_dp)
            ! rho(nu,i) = rho(nu,i) + ShalfP(nu,mu)*CV(mu,i)
            CALL cp_dbcsr_sm_fm_multiply(p_env%stda_env%shalf, cvec, res(ispin)%matrix, &
                                         n_occ(ispin), -1.0_dp, 1.0_dp)
            !
            CALL dbcsr_release(pdens)
            DEALLOCATE (kind_of)
         END IF
         !
         CALL cp_fm_release(cvec)
         DEALLOCATE (tv)
      ENDDO

      DEALLOCATE (tcharge, gtcharge)
      DEALLOCATE (first_sgf, last_sgf)

      CALL timestop(handle)

   END SUBROUTINE stda_calculate_kernel

! **************************************************************************************************

END MODULE qs_tddfpt_stda_utils
