!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      09-JUL-2002, TCH, development started
! **************************************************************************************************
MODULE qs_tddfpt_utils
   USE admm_types,                      ONLY: admm_type
   USE cp_array_utils,                  ONLY: cp_1d_r_p_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt_control_type
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              cp_dbcsr_sm_fm_multiply
   USE cp_fm_basic_linalg,              ONLY: cp_fm_scale,&
                                              cp_fm_scale_and_add,&
                                              cp_fm_trace
   USE cp_fm_cholesky,                  ONLY: cp_fm_cholesky_decompose,&
                                              cp_fm_cholesky_invert
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: &
        cp_fm_create, cp_fm_get_submatrix, cp_fm_init_random, cp_fm_p_type, cp_fm_release, &
        cp_fm_set_all, cp_fm_set_submatrix, cp_fm_to_fm, cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type
   USE kinds,                           ONLY: dp
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_mo_methods,                   ONLY: make_mo_eig
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_p_env_methods,                ONLY: p_env_create,&
                                              p_env_psi0_changed
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_scf_post_gpw,                 ONLY: make_lumo_gpw
   USE qs_scf_types,                    ONLY: qs_scf_env_type
   USE qs_tddfpt_stda_utils,            ONLY: stda_init_matrices
   USE qs_tddfpt_types,                 ONLY: tddfpt_env_allocate,&
                                              tddfpt_env_deallocate,&
                                              tddfpt_env_type
   USE scf_control_types,               ONLY: scf_control_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt_utils'

! **************************************************************************************************
   TYPE simple_solution_sorter
      INTEGER                               :: orbit
      INTEGER                               :: lumo
      REAL(KIND=dp)                         :: value
      TYPE(simple_solution_sorter), POINTER :: next
   END TYPE simple_solution_sorter

   PRIVATE

   ! METHODS
   PUBLIC :: tddfpt_cleanup, &
             tddfpt_init, &
             co_initial_guess, &
             find_contributions, &
             make_orbital_sets, &
             normalize, &
             reorthogonalize

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief Initialize some necessary structures for a tddfpt calculation.
!> \param p_env perturbation environment to be initialized
!> \param t_env tddfpt environment to be initialized
!> \param qs_env Quickstep environment with the results of a
!>                   ground state calcualtion
! **************************************************************************************************
   SUBROUTINE tddfpt_init(p_env, t_env, qs_env)

      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(tddfpt_env_type), INTENT(out)                 :: t_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_init', routineP = moduleN//':'//routineN

      TYPE(dft_control_type), POINTER                    :: dft_control

      !------------------!
      ! create the p_env !
      !------------------!
      CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE.)
      CALL p_env_psi0_changed(p_env, qs_env) ! update the m_epsilon matrix

      !------------------!
      ! create the t_env !
      !------------------!
      CALL tddfpt_env_allocate(t_env, p_env, qs_env)
      CALL tddfpt_env_init(t_env, qs_env)

      !---------------------!
      ! sTDA initialization !
      !---------------------!
      CALL get_qs_env(qs_env, dft_control=dft_control)
      IF (dft_control%tddfpt_control%kernel == 0) THEN
         CALL stda_init_matrices(qs_env, p_env)
      ENDIF

   END SUBROUTINE tddfpt_init

! **************************************************************************************************
!> \brief Initialize t_env with meaningfull values.
!> \param t_env ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE tddfpt_env_init(t_env, qs_env)

      TYPE(tddfpt_env_type), INTENT(inout)               :: t_env
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'tddfpt_env_init', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: n_spins, spin
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control

      NULLIFY (matrix_s, dft_control)

      CALL get_qs_env(qs_env, matrix_s=matrix_s, dft_control=dft_control)
      n_spins = dft_control%nspins
      IF (dft_control%tddfpt_control%invert_S) THEN
         DO spin = 1, n_spins
            CALL copy_dbcsr_to_fm(matrix_s(1)%matrix, t_env%invS(spin)%matrix)
            CALL cp_fm_cholesky_decompose(t_env%invS(spin)%matrix)
            CALL cp_fm_cholesky_invert(t_env%invS(spin)%matrix)
         END DO
      END IF

   END SUBROUTINE tddfpt_env_init

! **************************************************************************************************
!> \brief ...
!> \param t_env ...
!> \param p_env ...
! **************************************************************************************************
   SUBROUTINE tddfpt_cleanup(t_env, p_env)

      TYPE(tddfpt_env_type), INTENT(inout)               :: t_env
      TYPE(qs_p_env_type), POINTER                       :: p_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'tddfpt_cleanup', routineP = moduleN//':'//routineN

      CALL tddfpt_env_deallocate(t_env)
      CALL p_env_release(p_env)

   END SUBROUTINE tddfpt_cleanup

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \par History
!>      03.2019 created [jgh]
!> \author jgh
! **************************************************************************************************
   SUBROUTINE make_orbital_sets(qs_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'make_orbital_sets', &
         routineP = moduleN//':'//routineN

      CHARACTER(LEN=20)                                  :: subspace_selection
      INTEGER                                            :: handle, homo, iounit, ispin, min_lumos, &
                                                            nlumo, nlumos, nmo, nspins
      LOGICAL                                            :: compute_lumos, has_homo, has_lumo
      REAL(KIND=dp)                                      :: gap, homo_lumo(2, 2)
      REAL(KIND=dp), DIMENSION(:), POINTER               :: mo_eigenvalues
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_1d_r_p_type), DIMENSION(:), POINTER        :: unoccupied_evals
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: unoccupied_orbs
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: ks_rmpv, mo_derivs
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(qs_scf_env_type), POINTER                     :: scf_env
      TYPE(scf_control_type), POINTER                    :: scf_control
      TYPE(section_vals_type), POINTER                   :: dft_section, input

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      IF (iounit > 0) THEN
         WRITE (UNIT=iounit, FMT='(T3,A)') "Generating Orbital Sets"
      ENDIF

      ! test input for requested options
      CALL get_qs_env(qs_env, input=input)
      dft_section => section_vals_get_subs_vals(input, "DFT")
      nlumo = section_get_ival(dft_section, "TDDFPT%NLUMO")

      has_homo = .FALSE.
      has_lumo = .FALSE.
      homo_lumo = 0.0_dp

      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins

      ! Select active occupied subspace
      subspace_selection = "ALL"
!     subspace_selection = "ENERGY_WINDOW"
!     subspace_selection = "LOCALIZATION"
      IF (subspace_selection == "ALL" .OR. subspace_selection == "ENERGY_WINDOW") THEN
         ! Makes the MOs eigenstates, computes eigenvalues
         CALL get_qs_env(qs_env, mos=mos, mo_derivs=mo_derivs, matrix_ks=ks_rmpv)
         CALL get_qs_env(qs_env, scf_control=scf_control)
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CALL make_mo_eig(mos, nspins, ks_rmpv, scf_control, mo_derivs, admm_env=admm_env)
         ELSE
            CALL make_mo_eig(mos, nspins, ks_rmpv, scf_control, mo_derivs)
         END IF
         DO ispin = 1, nspins
            NULLIFY (mo_eigenvalues)
            CALL get_mo_set(mo_set=mos(ispin)%mo_set, eigenvalues=mo_eigenvalues, homo=homo)
            homo_lumo(ispin, 1) = mo_eigenvalues(homo)
         END DO
         has_homo = .TRUE.
      END IF

      ! Gets the lumos, and eigenvalues for the lumos
      compute_lumos = nlumo .NE. 0
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, homo=homo, nmo=nmo)
         compute_lumos = compute_lumos .AND. homo == nmo
      END DO
      IF (compute_lumos) THEN
         ALLOCATE (unoccupied_orbs(nspins))
         ALLOCATE (unoccupied_evals(nspins))
         CALL get_qs_env(qs_env, scf_env=scf_env)
         CPASSERT(ASSOCIATED(scf_env))
         CPASSERT(scf_env%ref_count > 0)

         min_lumos = nlumo
         CALL make_lumo_gpw(qs_env, scf_env, unoccupied_orbs, unoccupied_evals, min_lumos, nlumos)
         ALLOCATE (dft_control%tddfpt_control%lumos_eigenvalues(nlumos, nspins))
         DO ispin = 1, nspins
            dft_control%tddfpt_control%lumos_eigenvalues(1:nlumos, ispin) = &
               unoccupied_evals(ispin)%array(1:nlumos)
            homo_lumo(ispin, 2) = unoccupied_evals(ispin)%array(1)
         END DO
         dft_control%tddfpt_control%lumos => unoccupied_orbs
         has_lumo = .TRUE.
      ENDIF
      IF (has_homo .AND. has_lumo) THEN
         DO ispin = 1, nspins
            gap = homo_lumo(ispin, 2)-homo_lumo(ispin, 1)
            IF (iounit > 0) WRITE (iounit, '(T3,A,F12.6)') &
               "HOMO - LUMO gap [eV] :", gap*evolt
         ENDDO
      ENDIF

      IF (compute_lumos) THEN
         DO ispin = 1, nspins
            DEALLOCATE (unoccupied_evals(ispin)%array)
         ENDDO
         DEALLOCATE (unoccupied_evals)
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE make_orbital_sets

! **************************************************************************************************
!> \brief ...
!> \param matrices ...
!> \param energies ...
!> \param n_v ...
!> \param qs_env ...
! **************************************************************************************************
   SUBROUTINE co_initial_guess(matrices, energies, n_v, qs_env)

      TYPE(cp_fm_p_type), DIMENSION(:, :), POINTER       :: matrices
      REAL(kind=DP), DIMENSION(:), INTENT(OUT)           :: energies
      INTEGER, INTENT(IN)                                :: n_v
      TYPE(qs_environment_type), POINTER                 :: qs_env

      CHARACTER(len=*), PARAMETER :: routineN = 'co_initial_guess', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, n_cols, n_lumos, n_orbits, &
                                                            n_rows, n_spins, oo, spin, vo
      REAL(KIND=DP)                                      :: evd
      REAL(KIND=DP), ALLOCATABLE, DIMENSION(:, :)        :: guess, lumos
      REAL(KIND=DP), DIMENSION(:), POINTER               :: orbital_eigenvalues
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(simple_solution_sorter), POINTER              :: sorter_iterator, sorter_pointer, &
                                                            sorter_start
      TYPE(tddfpt_control_type), POINTER                 :: tddfpt_control

      CALL timeset(routineN, handle)

      NULLIFY (dft_control)

      CALL get_qs_env(qs_env, dft_control=dft_control)
      tddfpt_control => dft_control%tddfpt_control
      n_spins = dft_control%nspins
      energies = 0.0_dp

      IF (.NOT. ASSOCIATED(tddfpt_control%lumos)) THEN
         CPABORT("LUMOS missing")
      END IF

      DO spin = 1, n_spins

         n_cols = matrices(1, spin)%matrix%matrix_struct%ncol_global
         n_rows = matrices(1, spin)%matrix%matrix_struct%nrow_global

         DO i = 1, n_v
            CALL cp_fm_set_all(matrices(i, spin)%matrix, 0.0_dp)
         END DO

         CALL get_mo_set(qs_env%mos(spin)%mo_set, eigenvalues=orbital_eigenvalues)

         n_lumos = tddfpt_control%lumos(spin)%matrix%matrix_struct%ncol_global

         n_orbits = SIZE(orbital_eigenvalues)

         !-----------------------------------------!
         ! create a SORTED list of initial guesses !
         !-----------------------------------------!
         ! first element
         evd = tddfpt_control%lumos_eigenvalues(1, spin)-orbital_eigenvalues(n_orbits)
         ALLOCATE (sorter_start)
         sorter_start%orbit = n_orbits
         sorter_start%lumo = 1
         sorter_start%value = evd
         NULLIFY (sorter_start%next)
         ! rest of the elements
         DO oo = n_orbits, 1, -1
            DO vo = 1, n_lumos

               IF (oo == n_orbits .AND. vo == 1) CYCLE ! already in list

               evd = tddfpt_control%lumos_eigenvalues(vo, spin)-orbital_eigenvalues(oo)

               sorter_iterator => sorter_start
               NULLIFY (sorter_pointer)
               DO WHILE (ASSOCIATED(sorter_iterator%next))
                  IF (sorter_iterator%next%value > evd) THEN
                     sorter_pointer => sorter_iterator%next
                     EXIT
                  END IF
                  sorter_iterator => sorter_iterator%next
               END DO

               ALLOCATE (sorter_iterator%next)
               sorter_iterator%next%orbit = oo
               sorter_iterator%next%lumo = vo
               sorter_iterator%next%value = evd
               sorter_iterator%next%next => sorter_pointer

            END DO
         END DO

         ALLOCATE (lumos(n_rows, n_lumos), guess(n_rows, n_orbits))
         CALL cp_fm_get_submatrix(tddfpt_control%lumos(spin)%matrix, lumos, &
                                  start_col=1, n_cols=n_lumos)

         !-------------------!
         ! fill the matrices !
         !-------------------!
         sorter_iterator => sorter_start
         DO i = 1, MIN(n_v, n_orbits*n_lumos)
            CALL dcopy(n_rows*n_orbits, 0.0_dp, 0, guess, 1)
            CALL dcopy(n_rows, lumos(:, sorter_iterator%lumo), 1, &
                       guess(:, sorter_iterator%orbit), 1)
            CALL cp_fm_set_submatrix(matrices(i, spin)%matrix, &
                                     guess(:, 1:n_orbits))
            energies(i) = energies(i)+sorter_iterator%value/REAL(n_spins, dp)
            sorter_iterator => sorter_iterator%next
         END DO
         IF (n_v > n_orbits*n_lumos) THEN
            DO i = n_orbits*n_lumos+1, n_v
               CALL cp_fm_init_random(matrices(i, spin)%matrix, n_orbits)
               energies(i) = 1.0E38_dp
            END DO
         END IF

         !--------------!
         ! some cleanup !
         !--------------!
         DEALLOCATE (lumos, guess)
         sorter_iterator => sorter_start
         DO WHILE (ASSOCIATED(sorter_iterator))
            sorter_pointer => sorter_iterator
            sorter_iterator => sorter_iterator%next
            DEALLOCATE (sorter_pointer)
         END DO

      END DO

      CALL timestop(handle)

   END SUBROUTINE co_initial_guess

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param t_env ...
! **************************************************************************************************
   SUBROUTINE find_contributions(qs_env, t_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(tddfpt_env_type), INTENT(IN)                  :: t_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'find_contributions', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, ii, iounit, jspin, n_ev, &
                                                            n_ev_from_lumos, n_spins, occ, spin, &
                                                            virt
      INTEGER, DIMENSION(2)                              :: nhomos, nlumos, nrows
      REAL(KIND=dp)                                      :: contribution, summed_contributions
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: homo_coeff_col, lumo_coeff_col, ovloc
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: S_lumos
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: ovmat
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(tddfpt_control_type)                          :: t_control

      CALL timeset(routineN, handle)

      NULLIFY (S_lumos, matrix_s, dft_control)
      iounit = cp_logger_get_default_io_unit()
      CALL get_qs_env(qs_env, matrix_s=matrix_s, dft_control=dft_control)
      CALL get_qs_env(qs_env, para_env=para_env, blacs_env=blacs_env)

      IF (iounit > 0) WRITE (iounit, *)
      IF (iounit > 0) WRITE (iounit, *)

      t_control = dft_control%tddfpt_control
      n_ev = t_control%n_ev
      n_spins = dft_control%nspins

      ALLOCATE (S_lumos(n_spins))

      DO spin = 1, n_spins
         nrows(spin) = t_control%lumos(spin)%matrix%matrix_struct%nrow_global
         nhomos(spin) = t_env%evecs(1, spin)%matrix%matrix_struct%ncol_global
         nlumos(spin) = t_control%lumos(spin)%matrix%matrix_struct%ncol_global
         !   nlumos(spin) = nrows(spin)-nhomos(spin)
         CALL cp_fm_create(S_lumos(spin)%matrix, t_control%lumos(spin)%matrix%matrix_struct, &
                           "S times lumos")
         CALL cp_dbcsr_sm_fm_multiply(matrix_s(1)%matrix, t_control%lumos(spin)%matrix, &
                                      S_lumos(spin)%matrix, nlumos(spin), 1.0_dp, 0.0_dp)
      END DO

      ALLOCATE (homo_coeff_col(MAXVAL(nrows(1:n_spins)), 1), &
                lumo_coeff_col(MAXVAL(nrows(1:n_spins)), 1))
      ii = 0
      DO spin = 1, n_spins
         DO i = 1, nhomos(spin)*nlumos(spin) !n_ev
            ii = ii+1
            IF (ii .LE. n_ev) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, '(A,I3,5X,F15.6)') "  excited state : ", ii, t_env%evals(ii)*evolt
               WRITE (iounit, *)
            END IF
            DO jspin = 1, n_spins
               summed_contributions = 0.0_dp
               IF (n_spins == 2) THEN
                  IF (jspin == 1) THEN
                     IF (iounit > 0) WRITE (iounit, *) 'alpha:'
                  ELSE
                     IF (iounit > 0) WRITE (iounit, *) 'beta:'
                  END IF
               END IF
               !
               CALL cp_fm_struct_create(fmstruct, para_env, blacs_env, nhomos(jspin), nlumos(jspin))
               NULLIFY (ovmat)
               CALL cp_fm_create(ovmat, fmstruct)
               CALL cp_fm_struct_release(fmstruct)
               CALL cp_gemm("T", "N", nhomos(jspin), nlumos(jspin), nrows(jspin), 1.0_dp, &
                            t_env%evecs(ii, jspin)%matrix, S_lumos(jspin)%matrix, 0.0_dp, ovmat)
               ALLOCATE (ovloc(nhomos(jspin), nlumos(jspin)))
               CALL cp_fm_get_submatrix(ovmat, ovloc, 1, 1, nhomos(jspin), nlumos(jspin))
               CALL cp_fm_release(ovmat)
               occloop: DO occ = nhomos(jspin), 1, -1
                  DO virt = 1, nlumos(jspin)
                     contribution = ovloc(occ, virt)
                     summed_contributions = summed_contributions+(contribution)**2
                     IF (ABS(contribution) > 5.0e-2_dp) THEN
                        IF (iounit > 0) WRITE (iounit, '(14X,I5,A,I5,10X,F8.3,5X,F8.3)') &
                           occ, " ->", nhomos(spin)+virt, ABS(contribution), summed_contributions
                     END IF
                     IF (ABS(summed_contributions-1.0_dp) < 1.0e-3_dp) EXIT occloop
                  END DO
               END DO occloop
               DEALLOCATE (ovloc)
               !
            END DO
            IF (iounit > 0) WRITE (iounit, *)
            ENDIF
         END DO
      END DO

      n_ev_from_lumos = SUM(nhomos(1:n_spins)*nlumos(1:n_spins))
      IF (n_ev .GT. n_ev_from_lumos) THEN
         IF (iounit > 0) THEN
            WRITE (iounit, '(A)') " No lumos available for the remaining excitation energies."
            WRITE (iounit, '(A,A)') " -> From now on printing only excitation energies", &
               " without further assignment."
            WRITE (iounit, '(A)') " Increase # lumos with keyword NLUMO if desired otherwise."
            WRITE (iounit, *)
         END IF
         DO ii = n_ev_from_lumos+1, n_ev
            IF (iounit > 0) THEN
               WRITE (iounit, '(A,I3,5X,F15.6)') "  excited state : ", ii, t_env%evals(ii)*evolt
               WRITE (iounit, *)
            END IF
         ENDDO
      ENDIF

      !
      ! punch a checksum for the regs
      IF (iounit > 0) THEN
         WRITE (iounit, '(T2,A,E14.6)') ' TDDFPT : CheckSum  =', SQRT(SUM(t_env%evals**2))
      ENDIF

      DO spin = 1, n_spins
         CALL cp_fm_release(S_lumos(spin)%matrix)
      END DO
      DEALLOCATE (homo_coeff_col, lumo_coeff_col)

      DEALLOCATE (S_lumos)

      CALL timestop(handle)

   END SUBROUTINE find_contributions

! **************************************************************************************************
!> \brief ...
!> !---------------------------------------!
!> ! x must not be changed in this routine !
!> ! tmp_vec may be changed                !
!> !---------------------------------------!
!> \param X ...
!> \param tmp_vec ...
!> \param metric ...
! **************************************************************************************************
   SUBROUTINE normalize(X, tmp_vec, metric)

      TYPE(cp_fm_p_type), DIMENSION(:)                   :: x, tmp_vec
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: metric

      CHARACTER(LEN=*), PARAMETER :: routineN = 'normalize', routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, n_spins, spin
      REAL(KIND=dp)                                      :: norm, tmp

      CALL timeset(routineN, handle)

      n_spins = SIZE(x)
      norm = 0.0_dp

      DO spin = 1, n_spins
         tmp = 0.0_dp
         CALL cp_dbcsr_sm_fm_multiply(metric(1)%matrix, X(spin)%matrix, &
                                      tmp_vec(spin)%matrix, &
                                      X(spin)%matrix%matrix_struct%ncol_global, &
                                      1.0_dp, 0.0_dp)
         CALL cp_fm_trace(X(spin)%matrix, tmp_vec(spin)%matrix, tmp)
         norm = norm+tmp
      END DO

      norm = SQRT(norm)
      DO spin = 1, n_spins
         CALL cp_fm_scale((1.0_dp/norm), X(spin)%matrix)
      END DO

      CALL timestop(handle)

   END SUBROUTINE normalize

! **************************************************************************************************
!> \brief ...
!> \param X ...
!> \param V_set ...
!> \param SV_set ...
!> \param work ...
!> \param n ...
! **************************************************************************************************
   SUBROUTINE reorthogonalize(X, V_set, SV_set, work, n)

      TYPE(cp_fm_p_type), DIMENSION(:)                   :: X
      TYPE(cp_fm_p_type), DIMENSION(:, :)                :: V_set, SV_set
      TYPE(cp_fm_p_type), DIMENSION(:)                   :: work
      INTEGER, INTENT(IN)                                :: n

      CHARACTER(LEN=*), PARAMETER :: routineN = 'reorthogonalize', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, i, n_spins, spin
      REAL(KIND=dp)                                      :: dotp, tmp

      CALL timeset(routineN, handle)

      IF (n > 0) THEN

         n_spins = SIZE(X)
         DO spin = 1, n_spins
            CALL cp_fm_to_fm(X(spin)%matrix, work(spin)%matrix)
         END DO

         DO i = 1, n
            dotp = 0.0_dp
            DO spin = 1, n_spins
               CALL cp_fm_trace(SV_set(i, spin)%matrix, work(spin)%matrix, tmp)
               dotp = dotp+tmp
            END DO
            DO spin = 1, n_spins
               CALL cp_fm_scale_and_add(1.0_dp, X(spin)%matrix, &
                                        -1.0_dp*dotp, V_set(i, spin)%matrix)
            END DO
         END DO

      END IF

      CALL timestop(handle)

   END SUBROUTINE reorthogonalize

END MODULE qs_tddfpt_utils
