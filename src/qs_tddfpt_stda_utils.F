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
   USE basis_set_container_types,       ONLY: get_basis_from_container
   USE basis_set_types,                 ONLY: gto_basis_set_type
   USE cell_types,                      ONLY: cell_type
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              tddfpt_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              dbcsr_allocate_matrix_set
   USE cp_fm_basic_linalg,              ONLY: cp_fm_trace
   USE cp_fm_diag,                      ONLY: cp_fm_power
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_submatrix,&
                                              cp_fm_to_fm,&
                                              cp_fm_to_fm_submat,&
                                              cp_fm_type, &
                                              cp_fm_set_element,&
                                              cp_fm_get_element
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_create,&
                                              dbcsr_distribution_type,&
                                              dbcsr_finalize,&
                                              dbcsr_get_block_p,&
                                              dbcsr_p_type,&
                                              dbcsr_set,&
                                              dbcsr_type,&
                                              dbcsr_type_symmetric
   USE ewald_methods_stda,              ONLY: ewald_or_direct_sum
   USE kinds,                           ONLY: dp
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
   USE input_section_types,             ONLY: section_vals_get,& 
                                              section_vals_type,&
                                              section_vals_val_get
   USE input_constants,                 ONLY: tddfpt_singlet,&
                                              tddfpt_triplet
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

! calculate Coulomb and exchange gamma matrices

      CALL setup_gamma(qs_env, stda_env, stda_env%gamma_exchange, .TRUE.)
      CALL setup_gamma(qs_env, stda_env, stda_env%gamma_coulomb, .FALSE.)

! calculate Loewdin transition density charges

      CALL get_lowdin_transition_charges(qs_env, p_env, stda_env)

      CALL timestop(handle)

   END SUBROUTINE stda_init_matrices
! **************************************************************************************************
!> \brief Calculate sTDA Coulomb- and exchange-type contributions
!> \param qs_env ...
!> \param stda_env ...
!> \param gamma_matrix sTDA Coulomb- or exchange-type contributions
!> \param lexchange logical to choose between Coulomb or exchange
!> \note  Note the specific sTDA notation exchange-type integrals (ia|jb) refer to Coulomb interaction
!>                                        Coulomb-type integrals (ij|ab) refer to exchange interaction
! **************************************************************************************************
   SUBROUTINE setup_gamma(qs_env, stda_env, gamma_matrix, lexchange)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type), POINTER                       :: stda_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: gamma_matrix
      LOGICAL, INTENT(IN)                                :: lexchange

      CHARACTER(len=*), PARAMETER :: routineN = 'setup_gamma', routineP = moduleN//':'//routineN
      REAL(KIND=dp), PARAMETER                           :: rsmooth = 1.0_dp

      CHARACTER(len=20)                                  :: cstring
      INTEGER                                            :: handle, iatom, icol, ikind, irow, jatom, &
                                                            jkind, natom, neighbor_list_id, nmat
      INTEGER                                            :: output_unit
      INTEGER, DIMENSION(:), POINTER                     :: row_blk_sizes
      LOGICAL                                            :: do_ewald, do_mic, do_periodic, found
      REAL(KIND=dp)                                      :: dr, eta, fcut, r, ralpha, rbeta, rcut, &
                                                            rcuta, rcutb, x
      REAL(KIND=dp), DIMENSION(3)                        :: rij
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: gblock
      TYPE(cell_type), POINTER                           :: cell
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set

      CALL timeset(routineN, handle)

      NULLIFY (logger) !get output_unit
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      CALL get_qs_env(qs_env=qs_env, particle_set=particle_set, &
                      neighbor_list_id=neighbor_list_id, &
                      dbcsr_dist=dbcsr_dist, &
                      dft_control=dft_control, &
                      cell=cell)
      natom = SIZE(particle_set)

      do_mic = dft_control%tddfpt_control%stda_control%do_mic
      do_ewald = dft_control%tddfpt_control%stda_control%do_ewald
      do_periodic = (.NOT. (ALL(cell%perd == 0)))

      IF (do_mic.and.do_periodic) THEN
         CALL get_qs_env(qs_env=qs_env, sab_orb_mic=sab_orb)
      ELSE
         CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb)
      ENDIF

! setup gamma matrices as in xTB -----------------------------------------------------------------------------
      nmat = 1
      NULLIFY (gamma_matrix)
      CALL dbcsr_allocate_matrix_set(gamma_matrix, nmat)

      ALLOCATE (row_blk_sizes(natom))
      row_blk_sizes(1:natom) = 1
      ALLOCATE (gamma_matrix(1)%matrix)

      IF (lexchange) THEN
         cstring = "GAMMA EXCHANGE MATRIX"
      ELSE
         cstring = "GAMMA COULOMB MATRIX"
      ENDIF

      CALL dbcsr_create(gamma_matrix(1)%matrix, name=TRIM(ADJUSTL(cstring)), dist=dbcsr_dist, &
                        matrix_type=dbcsr_type_symmetric, row_blk_size=row_blk_sizes, col_blk_size=row_blk_sizes, &
                        nze=0, mutable_work=.TRUE.)
      CALL dbcsr_set(gamma_matrix(1)%matrix, 0.0_dp)

      DEALLOCATE (row_blk_sizes)

! setup the matrices using the neighbor list
      CALL cp_dbcsr_alloc_block_from_nbl(gamma_matrix(1)%matrix, sab_orb)

! Calculate gamma matrices------------------------------------------------------------------------------------
      NULLIFY (nl_iterator)
      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
         CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, &
                                iatom=iatom, jatom=jatom, r=rij)

         dr = SQRT(SUM(rij(:)**2)) ! interatomic distance
!
         eta = (stda_env%kind_param_set(ikind)%kind_param%hardness_param+ &
                stda_env%kind_param_set(jkind)%kind_param%hardness_param)/2.0_dp
!
         icol = MAX(iatom, jatom)
         irow = MIN(iatom, jatom)

         NULLIFY (gblock)
         CALL dbcsr_get_block_p(matrix=gamma_matrix(1)%matrix, row=irow, col=icol, BLOCK=gblock, found=found)
         CPASSERT(found)
!
! get rcuta and rcutb
         rcuta = stda_env%kind_param_set(ikind)%kind_param%rcut
         rcutb = stda_env%kind_param_set(jkind)%kind_param%rcut
         rcut = rcuta+rcutb

!>   Computes the derivative of the short-range gamma parameter from
!>   Nataga-Mishimoto-Ohno-Klopman formula equivalently as it is done for xTB
!>   WARNING: The function (gamma - 1/r) has still an l-dependent longrange (1/r^3)
!>            behaviour. We use a cutoff function to smoothly remove this part.
!>            However, this will change energies and effect final results.
         IF (dr < 1.e-6) THEN
            ! on site terms
            IF (lexchange) THEN
               gblock(:, :) = gblock(:, :)+eta
            ELSE
               gblock(:, :) = gblock(:, :)+stda_env%hfx_fraction*eta
            ENDIF
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
            IF (lexchange) THEN
               ralpha = dr**stda_env%alpha_param
               IF (do_periodic.AND.(.NOT.do_mic)) THEN !direct sum
                  gblock(:, :) = gblock(:, :)+ &
                                 fcut*(1._dp/(dr**(stda_env%alpha_param)+eta**(-stda_env%alpha_param))) &
                                 **(1._dp/stda_env%alpha_param)-fcut/dr
               ELSEIF ((do_mic.AND.do_periodic).OR.(.NOT.do_periodic)) THEN !mic or non-periodic 
                  gblock(:, :) = gblock(:, :)+ &
                                 (1._dp/(dr**(stda_env%alpha_param)+eta**(-stda_env%alpha_param))) &
                                 **(1._dp/stda_env%alpha_param)
               ENDIF 
            ELSE
               rbeta = dr**stda_env%beta_param
               IF (stda_env%hfx_fraction > 0._dp) THEN
                  IF (do_periodic) THEN !both for ewald and direct sum 
                     gblock(:, :) = gblock(:, :)+ &
                                    fcut*(1._dp/(rbeta+(stda_env%hfx_fraction*eta)**(-stda_env%beta_param))) &
                                    **(1._dp/stda_env%beta_param) &
                                    -stda_env%hfx_fraction*fcut/dr
                  ELSE !non-periodic
                     gblock(:, :) = gblock(:, :)+ &
                                    (1._dp/(rbeta+(stda_env%hfx_fraction*eta)**(-stda_env%beta_param))) &
                                    **(1._dp/stda_env%beta_param)
                  ENDIF
               ELSE
                  !            gblock(:,:) = gblock(:,:) + &
                  !                          (1._dp/(dr**(stda_env%beta_param)))**(1._dp/stda_env%beta_param)
                  ! to be consistent with sTDA program of Grimme gblock is set to 0.
                  gblock(:, :) = 0._dp
               ENDIF
            ENDIF
         END IF

          WRITE(output_unit,*) 'gblock(:)', gblock(:,:), lexchange
      END DO

      CALL neighbor_list_iterator_release(nl_iterator)
      CALL dbcsr_finalize(gamma_matrix(1)%matrix)

      CALL timestop(handle)

   END SUBROUTINE setup_gamma
! **************************************************************************************************
!> \brief Calculate Lowdin transition density charges for sTDA (from Lowdin MO coefficients)
!> \param qs_env ...
!> \param p_env ...
!> \param stda_env ...
! **************************************************************************************************
   SUBROUTINE get_lowdin_transition_charges(qs_env, p_env, stda_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(stda_env_type), POINTER                       :: stda_env

      CHARACTER(len=*), PARAMETER :: routineN = 'get_lowdin_transition_charges', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ii, ispin, jj,  &
                                                            n_occ, ndep, nsgf, nspins, &
                                                            output_unit
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: values
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: fm_s_half, fm_work1, fm_work2
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrixkp_s
      TYPE(dbcsr_type), POINTER                          :: sm_s
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(scf_control_type), POINTER                    :: scf_control

      CALL timeset(routineN, handle)

      NULLIFY (logger) !get output_unit
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

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

      sm_s => matrixkp_s(1, 1)%matrix ! Overlap matrix in sparse format

      ! Provide full size work matrices
      CALL cp_fm_struct_create(fmstruct=fmstruct, &
                               para_env=para_env, &
                               context=blacs_env, &
                               nrow_global=nsgf, &
                               ncol_global=nsgf)
      CALL cp_fm_create(matrix=fm_s_half, &
                        matrix_struct=fmstruct, &
                        name="S^(1/2) MATRIX")
      CALL cp_fm_create(matrix=fm_work1, &
                        matrix_struct=fmstruct, &
                        name="TMP MATRIX")
      CALL cp_fm_struct_release(fmstruct=fmstruct)

      CALL copy_dbcsr_to_fm(sm_s, fm_s_half)
      CALL cp_fm_power(fm_s_half, fm_work1, 0.5_dp, scf_control%eps_eigval, ndep)

      CALL cp_fm_to_fm(fm_s_half, p_env%stda_env%shalf)

      IF (ndep /= 0) &
         CALL cp_warn(__LOCATION__, &
                      "Overlap matrix exhibits linear dependencies. At least some "// &
                      "eigenvalues have been quenched.")

      DO ispin = 1, nspins

         n_occ = stda_env%n_occ(ispin) !p_env%n_mo(ispin)

! Build Lowdin coefficients C'= S^1/2 C for each spin
         NULLIFY (fm_work2)
         CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                  para_env=para_env, &
                                  context=blacs_env, &
                                  nrow_global=nsgf, &
                                  ncol_global=n_occ)
         CALL cp_fm_create(matrix=fm_work2, &
                           matrix_struct=fmstruct, &
                           name="TMP MATRIX 2")
         CALL cp_gemm(transa="N", &
                      transb="N", &
                      m=nsgf, &
                      n=n_occ, &
                      k=nsgf, &
                      alpha=1.0_dp, &
                      matrix_a=fm_s_half, &
                      matrix_b=p_env%psi0d(ispin)%matrix, & ! row: n_ao, column: n_mo
                      beta=0.0_dp, &
                      matrix_c=fm_work2)
!
         ! Test that sum_mu |C_mup|^2 =1 for all p
         ALLOCATE (values(fmstruct%ncol_global))
         values = 0.0_dp
         DO ii = 1, fmstruct%ncol_global
            DO jj = 1, fmstruct%nrow_global
               values(ii) = values(ii)+(ABS(fm_work2%local_data(jj, ii)))**2
            ENDDO
         ENDDO
!         DO ii = 1, fmstruct%ncol_global
!            WRITE (UNIT=output_unit, FMT="(A,F10.8)") 'check sum lowdin coefficients only occ', values(ii)
!         ENDDO
         DEALLOCATE (values)

         CALL cp_fm_to_fm(fm_work2, stda_env%ctransformed(ispin)%matrix)

         CALL cp_fm_struct_release(fmstruct=fmstruct)
         CALL cp_fm_release(matrix=fm_work2)
         NULLIFY (fm_work2, fmstruct)

      ENDDO !nspins

      CALL cp_fm_release(matrix=fm_s_half)
      CALL cp_fm_release(matrix=fm_work1)
      NULLIFY (fm_work1, fm_s_half)

      CALL lowdin_charges_fm(qs_env, stda_env, stda_env%ctransformed)

      CALL timestop(handle)

   END SUBROUTINE get_lowdin_transition_charges
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

      INTEGER                                            :: handle, iatom, ikind, ispin, n_atoms, &
                                                            n_occ, nsgf, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, last_sgf
      INTEGER, DIMENSION(:), POINTER                     :: kind_of
      TYPE(stda_env_type), POINTER                       :: stda_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: fm_work1, fm_work2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(p_env%stda_env))
      NULLIFY(stda_env)
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

         CALL cp_gemm(transa="N", &
                      transb="N", &
                      m=nsgf, &
                      n=stda_env%n_occ(ispin), &
                      k=nsgf, &
                      alpha=1.0_dp, &
                      matrix_a=stda_env%shalf, &
                      matrix_b=X(ispin)%matrix, &
                      beta=0.0_dp, &
                      matrix_c=fm_work2)
!
         CALL cp_fm_to_fm(fm_work2, stda_env%xtransformed(ispin)%matrix)

         CALL cp_fm_release(matrix=fm_work2)
         NULLIFY (fm_work2, fmstruct)
      ENDDO !nspins

! Build Lowdin transformed tilde(tilde(X))_A=S^1/2_A S^1/2 X for each spin
      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      n_atoms = SIZE(particle_set)
      ALLOCATE (kind_of(n_atoms))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

      ALLOCATE (first_sgf(n_atoms))
      ALLOCATE (last_sgf(n_atoms))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)

      DO iatom = 1, n_atoms

         CPASSERT(ASSOCIATED(stda_env%atom_param_set))
         CPASSERT(ASSOCIATED(stda_env%atom_param_set(iatom)%atom_param))

         DO ispin = 1, nspins

            n_occ = stda_env%n_occ(ispin) 

            NULLIFY (basis_set)
            NULLIFY (qs_kind) ! get orbital basis functions per atom
            ikind = kind_of(iatom)
            qs_kind => qs_kind_set(ikind)
            CALL get_basis_from_container(qs_kind%basis_sets, basis_set=basis_set, basis_type="ORB")

            NULLIFY (fm_work1, fmstruct)
            CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                     para_env=para_env, &
                                     context=blacs_env, &
                                     nrow_global=basis_set%nsgf, &
                                     ncol_global=n_occ)
            CALL cp_fm_create(matrix=fm_work1, &
                              matrix_struct=fmstruct, &
                              name="TMP MATRIX")
            CALL cp_fm_struct_release(fmstruct=fmstruct)

            NULLIFY (fm_work2, fmstruct)
            CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                     para_env=para_env, &
                                     context=blacs_env, &
                                     nrow_global=basis_set%nsgf, &
                                     ncol_global=nsgf)
            CALL cp_fm_create(matrix=fm_work2, &
                              matrix_struct=fmstruct, &
                              name="TMP MATRIX 2")
            CALL cp_fm_struct_release(fmstruct=fmstruct)
!
            CALL cp_fm_to_fm_submat(stda_env%shalf, fm_work2, basis_set%nsgf, nsgf, &
                                    first_sgf(iatom), 1, 1, 1)

            CALL cp_fm_to_fm_submat(stda_env%xtransformed(ispin)%matrix, fm_work1, basis_set%nsgf, n_occ, &
                                    first_sgf(iatom), 1, 1, 1)

            CALL cp_gemm(transa="T", &
                         transb="N", &
                         m=nsgf, &
                         n=n_occ, &
                         k=basis_set%nsgf, &
                         alpha=1.0_dp, &
                         matrix_a=fm_work2, &
                         matrix_b=fm_work1, &
                         beta=0.0_dp, &
                         matrix_c=stda_env%atom_param_set(iatom)%atom_param%xtransformed_double(ispin)%matrix)

            CALL cp_fm_release(matrix=fm_work1)
            NULLIFY (fm_work1)
            CALL cp_fm_release(matrix=fm_work2)
            NULLIFY (fm_work2)
         ENDDO !ispin
      ENDDO !iatom

      DEALLOCATE (first_sgf, last_sgf)
      DEALLOCATE (kind_of)

      CALL timestop(handle)

   END SUBROUTINE get_lowdin_x
! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param stda_env ...
!> \param ctransformed_occ ...
! **************************************************************************************************
   SUBROUTINE lowdin_charges_fm(qs_env, stda_env, ctransformed_occ)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type), POINTER                       :: stda_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: ctransformed_occ

      CHARACTER(len=*), PARAMETER :: routineN = 'lowdin_charges_fm', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, iatom, ikind, ispin, n_atoms, &
                                                            n_occ, nsgf, nspins
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf, last_sgf
      INTEGER, DIMENSION(:), POINTER                     :: kind_of
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fmstruct
      TYPE(cp_fm_type), POINTER                          :: fm_work1, fm_work2
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(stda_env))

      NULLIFY (dft_control, para_env, blacs_env)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      blacs_env=blacs_env)
      nspins = dft_control%nspins

      NULLIFY (qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

      NULLIFY (particle_set, atomic_kind_set, kind_of)
      CALL get_qs_env(qs_env=qs_env, &
                      particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      n_atoms = SIZE(particle_set)
      ALLOCATE (kind_of(n_atoms))
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, &
                               kind_of=kind_of)

      ALLOCATE (first_sgf(n_atoms))
      ALLOCATE (last_sgf(n_atoms))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf, last_sgf=last_sgf)

!!!----------------------------------------------------------------------------------
! Calculate Lowdin transition density charges q^A_ip=sum_mu_A C_mu_Ai C_mu_Ap for every atom A

      DO iatom = 1, n_atoms

         CPASSERT(ASSOCIATED(stda_env%atom_param_set))
         CPASSERT(ASSOCIATED(stda_env%atom_param_set(iatom)%atom_param))

         DO ispin = 1, nspins

            n_occ = stda_env%n_occ(ispin) 

            NULLIFY (basis_set)
            NULLIFY (qs_kind) 
            ikind = kind_of(iatom)
            qs_kind => qs_kind_set(ikind)
            CALL get_basis_from_container(qs_kind%basis_sets, basis_set=basis_set, basis_type="ORB")

            NULLIFY (fm_work1, fmstruct)
            CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                     para_env=para_env, &
                                     context=blacs_env, &
                                     nrow_global=basis_set%nsgf, &
                                     ncol_global=n_occ)
            CALL cp_fm_create(matrix=fm_work1, &
                              matrix_struct=fmstruct, &
                              name="TMP MATRIX")
            CALL cp_fm_struct_release(fmstruct=fmstruct)
!
            CALL cp_fm_to_fm_submat(ctransformed_occ(ispin)%matrix, fm_work1, basis_set%nsgf, n_occ, &
                                    first_sgf(iatom), 1, 1, 1)

            CALL cp_gemm(transa="T", &
                         transb="N", &
                         m=n_occ, &
                         n=n_occ, &
                         k=basis_set%nsgf, &
                         alpha=1.0_dp, &
                         matrix_a=fm_work1, &
                         matrix_b=fm_work1, &
                         beta=0.0_dp, &
                         matrix_c=stda_env%atom_param_set(iatom)%atom_param%charges(ispin)%matrix)

            CALL cp_fm_release(matrix=fm_work1)
            NULLIFY (fm_work1)
!
!!!----------------------------------------------------------------------------------
! Calculating tilde(tilde(C))_A = S^1/2_A (S^1/2 C)
            NULLIFY (fm_work1, fmstruct)
            CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                     para_env=para_env, &
                                     context=blacs_env, &
                                     nrow_global=basis_set%nsgf, &
                                     ncol_global=n_occ)
            CALL cp_fm_create(matrix=fm_work1, &
                              matrix_struct=fmstruct, &
                              name="TMP MATRIX")
            CALL cp_fm_struct_release(fmstruct=fmstruct)
            NULLIFY (fm_work2, fmstruct)
            CALL cp_fm_struct_create(fmstruct=fmstruct, &
                                     para_env=para_env, &
                                     context=blacs_env, &
                                     nrow_global=basis_set%nsgf, &
                                     ncol_global=nsgf) !p_env%n_ao(ispin))
            CALL cp_fm_create(matrix=fm_work2, &
                              matrix_struct=fmstruct, &
                              name="TMP MATRIX 2")
            CALL cp_fm_struct_release(fmstruct=fmstruct)
!
            CALL cp_fm_to_fm_submat(ctransformed_occ(ispin)%matrix, fm_work1, basis_set%nsgf, n_occ, &
                                    first_sgf(iatom), 1, 1, 1)
!
            CALL cp_fm_to_fm_submat(stda_env%shalf, fm_work2, basis_set%nsgf,nsgf, &   ! p_env%n_ao(ispin), &
                                    first_sgf(iatom), 1, 1, 1)

            CALL cp_gemm(transa="T", &
                         transb="N", &
                         m=nsgf, & !p_env%n_ao(ispin), &
                         n=n_occ, &
                         k=basis_set%nsgf, &
                         alpha=1.0_dp, &
                         matrix_a=fm_work2, &
                         matrix_b=fm_work1, &
                         beta=0.0_dp, &
                         matrix_c=stda_env%atom_param_set(iatom)%atom_param%ctransformed_double(ispin)%matrix)

            CALL cp_fm_release(matrix=fm_work1)
            NULLIFY (fm_work1)
            CALL cp_fm_release(matrix=fm_work2)
            NULLIFY (fm_work2)
!!!----------------------------------------------------------------------------------
         ENDDO !ispin
      ENDDO !iatom

      DEALLOCATE (first_sgf, last_sgf)
      DEALLOCATE (kind_of)
      CALL timestop(handle)

   END SUBROUTINE lowdin_charges_fm
! **************************************************************************************************
!> \brief ...Calculate the sTDA kernel contribution by contracting the Lowdin MO coefficients / transition
!>            charges with the Coulomb-type or exchange-type integrals
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

      INTEGER                                            :: handle, iatom, icol, ikind, irow, ispin, &
                                                            jatom, jkind, natom, nsgf, nspins, &
                                                            output_unit, jspin
      INTEGER                                            :: aa,ii,jj, bb, res_etype
      REAL(KIND=dp)                                      :: tmp_coul, tmp_ex, tmp3
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: first_sgf
      INTEGER, DIMENSION(:), POINTER                     :: kind_of
      LOGICAL                                            :: do_ewald, do_mic, do_periodic, found
      REAL(KIND=dp)                                      :: tmp, tmp1, tmp2, spinfac
      REAL(KIND=dp), DIMENSION(:, :), POINTER            :: coulombblock, exchangeblock
      TYPE(stda_env_type), POINTER                       :: stda_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: basis_set_iatom, basis_set_jatom
      TYPE(neighbor_list_iterator_p_type), &
         DIMENSION(:), POINTER                           :: nl_iterator
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      TYPE(section_vals_type), POINTER                   :: input
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      INTEGER, DIMENSION(2)                              :: n_occ

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(p_env))
      CPASSERT(ASSOCIATED(p_env%stda_env))

      NULLIFY (logger) ! Get output_unit
      logger => cp_get_default_logger()
      output_unit = cp_logger_get_default_io_unit(logger)

      NULLIFY(stda_env)
      stda_env => p_env%stda_env
      n_occ(:) = stda_env%n_occ(:)

      NULLIFY (particle_set, dft_control, cell, para_env, blacs_env)
      CALL get_qs_env(qs_env, particle_set=particle_set, &
                      dft_control=dft_control, &
                      cell=cell, &
                      para_env=para_env, &
                      blacs_env=blacs_env)
      natom = SIZE(particle_set)
      nspins = SIZE(X)
      spinfac=2.0_dp
      IF (nspins==2) spinfac=1.0_dp

      do_mic = dft_control%tddfpt_control%stda_control%do_mic
      do_ewald = dft_control%tddfpt_control%stda_control%do_ewald
      do_periodic = (.NOT. (ALL(cell%perd == 0)))

!      IF (do_periodic.AND.((.NOT.do_ewald).OR.(.NOT.do_mic))) THEN
!         CPABORT("Please choose minimum image and ewald for periodic calculations")
!      ENDIF

      IF (do_mic.AND.do_periodic) THEN !minimum image
         CALL get_qs_env(qs_env=qs_env, sab_orb_mic=sab_orb)
      ELSE
         CALL get_qs_env(qs_env=qs_env, sab_orb=sab_orb)
      ENDIF

      NULLIFY (qs_kind_set)
      CALL get_qs_env(qs_env=qs_env, qs_kind_set=qs_kind_set)
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)
      ALLOCATE (first_sgf(natom))
      CALL get_particle_set(particle_set, qs_kind_set, first_sgf=first_sgf)

      NULLIFY (atomic_kind_set)
      NULLIFY (kind_of)
      ALLOCATE (kind_of(natom))
      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set)
      CALL get_atomic_kind_set(atomic_kind_set=atomic_kind_set, kind_of=kind_of)

      NULLIFY(input)
      CALL get_qs_env(qs_env=qs_env, input=input)
      CALL section_vals_val_get(input, "DFT%TDDFPT%RES_ETYPE", i_val=res_etype)
      
! calculate Loewdin transformed Davidson trial vector tilde(X)=S^1/2*X and tilde(tilde(X))=S^1/2_A*tilde(X)_A

      CALL get_lowdin_x(qs_env, p_env, X)

      CALL neighbor_list_iterator_create(nl_iterator, sab_orb)
      DO ispin = 1, nspins
!      DO jspin=1,nspins
         DO WHILE (neighbor_list_iterate(nl_iterator) == 0)
            CALL get_iterator_info(nl_iterator, ikind=ikind, jkind=jkind, iatom=iatom, jatom=jatom)
            icol = MAX(iatom, jatom)
            irow = MIN(iatom, jatom)

! Get exchange-type and Coulomb-type contribution for the kernel
            NULLIFY (exchangeblock)
            CALL dbcsr_get_block_p(matrix=stda_env%gamma_exchange(1)%matrix, row=irow, col=icol, &
                                   BLOCK=exchangeblock, found=found)
            CPASSERT(found)

            NULLIFY (coulombblock)
            CALL dbcsr_get_block_p(matrix=stda_env%gamma_coulomb(1)%matrix, row=irow, col=icol, &
                                   BLOCK=coulombblock, found=found)
            CPASSERT(found)

! Get number of basis functions per atom
            NULLIFY (basis_set_iatom)
            NULLIFY (qs_kind)
            ikind = kind_of(iatom) !iatom
            qs_kind => qs_kind_set(ikind)
            CALL get_basis_from_container(qs_kind%basis_sets, basis_set=basis_set_iatom, basis_type="ORB")

            NULLIFY (basis_set_jatom)
            NULLIFY (qs_kind)
            jkind = kind_of(jatom) !jatom
            qs_kind => qs_kind_set(jkind)
            CALL get_basis_from_container(qs_kind%basis_sets, basis_set=basis_set_jatom, basis_type="ORB")

!1. Contract exchange-type contributions to the kernel with coefficients
!   2* gamma_K(A,B) * Cdoubletilde_B *Ctilde_A *Xtilde_A
!            NULLIFY (Xtilde, fmstruct)
!            CALL cp_fm_struct_create(fmstruct=fmstruct, &
!                                     para_env=para_env, &
!                                     context=blacs_env, &
!                                     nrow_global=basis_set_iatom%nsgf, &
!                                     ncol_global=n_occ(ispin))
!            CALL cp_fm_create(matrix=Xtilde, &
!                              matrix_struct=fmstruct, &
!                              name="MATRIX X tilde")
!            CALL cp_fm_struct_release(fmstruct=fmstruct)
!!
!            NULLIFY (ctilde, fmstruct)
!            CALL cp_fm_struct_create(fmstruct=fmstruct, &
!                                     para_env=para_env, &
!                                     context=blacs_env, &
!                                     nrow_global=basis_set_iatom%nsgf, &
!                                     ncol_global=n_occ(ispin))
!            CALL cp_fm_create(matrix=ctilde, &
!                              matrix_struct=fmstruct, &
!                              name="MATRIX CTRANSFORMED_ATOM")
!            CALL cp_fm_struct_release(fmstruct=fmstruct)
!
!            CALL cp_fm_to_fm_submat(stda_env%xtransformed(ispin)%matrix, Xtilde, &
!                                    basis_set_iatom%nsgf, n_occ(ispin), &
!                                    first_sgf(iatom), 1, 1, 1)
!            CALL cp_fm_to_fm_submat(stda_env%ctransformed(ispin)%matrix, ctilde, &
!                                    basis_set_iatom%nsgf, n_occ(ispin), &
!                                    first_sgf(iatom), 1, 1, 1)
!!     0.5*(Ctransformed*X^T + X*Ctransformed^T)
!            CALL cp_fm_trace(ctilde, Xtilde, tmp1)
!            CALL cp_fm_trace(Xtilde, ctilde, tmp2)
!            tmp = 0.5_dp*(tmp1+tmp2) ! symmetrization of n^1 neccessary ?
!
!            CALL cp_fm_release(matrix=Xtilde)
!            CALL cp_fm_release(matrix=ctilde)
!            NULLIFY (Xtilde, ctilde)
!
!!     Exchange-type contribution to the kernel
!            kernel = tmp*(2.0_dp*exchangeblock(1, 1))
!            alpha = kernel
!            beta = 1.0_dp
!
!!            CALL cp_fm_set_submatrix(res(ispin)%matrix, &
!!                                     stda_env%atom_param_set(jatom)%atom_param%ctransformed_double(ispin)%matrix%local_data, &
!!                                     1, 1, nsgf, n_occ(ispin), alpha=alpha, beta=beta)
!            ! iF (iatom.NE.jatom) ------------------------------------------------------------------------------
!            IF (iatom .NE. jatom) THEN
!               NULLIFY (Xtilde, fmstruct)
!               CALL cp_fm_struct_create(fmstruct=fmstruct, &
!                                        para_env=para_env, &
!                                        context=blacs_env, &
!                                        nrow_global=basis_set_jatom%nsgf, &
!                                        ncol_global=n_occ(ispin))
!               CALL cp_fm_create(matrix=Xtilde, &
!                                 matrix_struct=fmstruct, &
!                                 name="MATRIX X tilde")
!               CALL cp_fm_struct_release(fmstruct=fmstruct)
!!
!               NULLIFY (ctilde, fmstruct)
!               CALL cp_fm_struct_create(fmstruct=fmstruct, &
!                                        para_env=para_env, &
!                                        context=blacs_env, &
!                                        nrow_global=basis_set_jatom%nsgf, &
!                                        ncol_global=n_occ(ispin))
!               CALL cp_fm_create(matrix=ctilde, &
!                                 matrix_struct=fmstruct, &
!                                 name="MATRIX CTRANSFORMED_ATOM")
!               CALL cp_fm_struct_release(fmstruct=fmstruct)
!
!               CALL cp_fm_to_fm_submat(stda_env%xtransformed(ispin)%matrix, Xtilde, &
!                                       basis_set_jatom%nsgf, n_occ(ispin), &
!                                       first_sgf(jatom), 1, 1, 1)
!               CALL cp_fm_to_fm_submat(stda_env%ctransformed(ispin)%matrix, ctilde, &
!                                       basis_set_jatom%nsgf, n_occ(ispin), &
!                                       first_sgf(jatom), 1, 1, 1)
!!       0.5*(Ctransformed*X^T + X*Ctransformed^T)
!               CALL cp_fm_trace(ctilde, Xtilde, tmp1)
!               CALL cp_fm_trace(Xtilde, ctilde, tmp2)
!               tmp = 0.5_dp*(tmp1+tmp2) ! symmetrization of n^1 ?
!
!               CALL cp_fm_release(matrix=Xtilde)
!               CALL cp_fm_release(matrix=ctilde)
!               NULLIFY (Xtilde, ctilde)
!
!!       Exchange-type contribution to the kernel
!               kernel = tmp*(2.0_dp*exchangeblock(1, 1))
!               alpha = kernel
!               beta = 1.0_dp
!!               CALL cp_fm_set_submatrix(res(ispin)%matrix, &
!!                                     stda_env%atom_param_set(iatom)%atom_param%ctransformed_double(ispin)%matrix%local_data, &
!!                                        1, 1, nsgf, n_occ(ispin), alpha=alpha, beta=beta)
!            ENDIF
!            ! iF (iatom.NE.jatom) ------------------------------------------------------------------------------
!
!!2. Contract Coulomb-type contributions of the kernel with coefficients
!!   - gamma_J(A,B) Xdoubletilde_B q_A
!            NULLIFY (Xdoubletilde, fmstruct)
!            CALL cp_fm_struct_create(fmstruct=fmstruct, &
!                                     para_env=para_env, &
!                                     context=blacs_env, &
!                                     nrow_global=nsgf, &
!                                     ncol_global=n_occ(ispin))
!            CALL cp_fm_create(matrix=Xdoubletilde, &
!                              matrix_struct=fmstruct, &
!                              name="MATRIX X doubletilde")
!            CALL cp_fm_struct_release(fmstruct=fmstruct)
!!
!            NULLIFY (tmp_res, fmstruct)
!            CALL cp_fm_struct_create(fmstruct=fmstruct, &
!                                     para_env=para_env, &
!                                     context=blacs_env, &
!                                     nrow_global=nsgf, &
!                                     ncol_global=n_occ(ispin))
!            CALL cp_fm_create(matrix=tmp_res, &
!                              matrix_struct=fmstruct, &
!                              name="MATRIX TMP_RES")
!            CALL cp_fm_struct_release(fmstruct=fmstruct)
!!
!            CALL cp_fm_to_fm_submat(stda_env%atom_param_set(jatom)%atom_param%xtransformed_double(ispin)%matrix, &
!                                    Xdoubletilde, &
!                                    nsgf,n_occ(ispin), &
!                                    1, 1, 1, 1)
!!
!!      kernel= -0.5_dp*(coulombblock(1,1))
!            kernel = -(coulombblock(1, 1))
!            CALL cp_gemm(transa="N", &
!                         transb="T", &
!                         m=nsgf, &
!                         n=n_occ(ispin), &
!                         k=n_occ(ispin), &
!                         alpha=kernel, &
!                         matrix_a=Xdoubletilde, &
!                         matrix_b=stda_env%atom_param_set(iatom)%atom_param%charges(ispin)%matrix, &
!                         beta=0.0_dp, &
!                         matrix_c=tmp_res)
!            !     CALL cp_gemm(transa="N",            &  !symmetrization of n^1 neccessary?
!            !                  transb="T",            &
!            !                  m=p_env%n_mo(ispin),   &
!            !                  n=p_env%n_mo(ispin),   &
!            !                  k=nsgf,                &
!            !                  alpha=kernel,          &
!            !                  matrix_a=stda_env%atom_param_set(iatom)%atom_param%charges(ispin)%matrix, &
!            !                  matrix_b=Xdoubletilde, &
!            !                  beta=1.0_dp,           &
!            !                  matrix_c=tmp_res)
!            alpha = 1.0_dp
!            beta = 1.0_dp
!!            CALL cp_fm_set_submatrix(res(ispin)%matrix, tmp_res%local_data, &
!!                                     1, 1, &
!!                                     nsgf, n_occ(ispin), alpha=alpha, beta=beta)
!
!            CALL cp_fm_release(matrix=Xdoubletilde)
!            CALL cp_fm_release(matrix=tmp_res)
!            NULLIFY (Xdoubletilde, tmp_res)
!
!            !IF (iatom.NE.jatom) -----------------------------------------------------------------
!1            IF (iatom .NE. jatom) THEN
!               NULLIFY (Xdoubletilde, fmstruct)
!               CALL cp_fm_struct_create(fmstruct=fmstruct, &
!                                        para_env=para_env, &
!                                        context=blacs_env, &
!                                        nrow_global=nsgf, &
!                                        ncol_global=n_occ(ispin))
!               CALL cp_fm_create(matrix=Xdoubletilde, &
!                                 matrix_struct=fmstruct, &
!                                 name="MATRIX X doubletilde")
!               CALL cp_fm_struct_release(fmstruct=fmstruct)
!!
!               NULLIFY (tmp_res, fmstruct)
!               CALL cp_fm_struct_create(fmstruct=fmstruct, &
!                                        para_env=para_env, &
!                                        context=blacs_env, &
!                                        nrow_global=nsgf, &
!                                        ncol_global=n_occ(ispin))
!               CALL cp_fm_create(matrix=tmp_res, &
!                                 matrix_struct=fmstruct, &
!                                 name="MATRIX TMP_RES")
!               CALL cp_fm_struct_release(fmstruct=fmstruct)
!!
!               CALL cp_fm_to_fm_submat(stda_env%atom_param_set(iatom)%atom_param%xtransformed_double(ispin)%matrix, &
!                                       Xdoubletilde, &
!                                       nsgf, n_occ(ispin), &
!                                       1, 1, 1, 1)
!!
!!        kernel= -0.5_dp*(coulombblock(1,1))
!               kernel = -(coulombblock(1, 1))
!               CALL cp_gemm(transa="N", &
!                            transb="T", &
!                            m=nsgf, &
!                            n=n_occ(ispin), &
!                            k=n_occ(ispin), &
!                            alpha=kernel, &
!                            matrix_a=Xdoubletilde, &
!                            matrix_b=stda_env%atom_param_set(jatom)%atom_param%charges(ispin)%matrix, &
!                            beta=0.0_dp, &
!                            matrix_c=tmp_res)
!!        CALL cp_gemm(transa="N",            &  !symmetrization of n^1 neccessary?
!!                     transb="T",            &
!!                     m=p_env%n_mo(ispin),   &
!!                     n=p_env%n_mo(ispin),   &
!!                     k=nsgf,                &
!!                     alpha=kernel,          &
!!                     matrix_a=stda_env%atom_param_set(jatom)%atom_param%charges(ispin)%matrix, &
!!                     matrix_b=Xdoubletilde, &
!!                     beta=1.0_dp,           &
!!                     matrix_c=tmp_res)
!               alpha = 1.0_dp
!               beta = 1.0_dp
!!               CALL cp_fm_set_submatrix(res(ispin)%matrix, tmp_res%local_data, &
!!                                        1, 1, &
!!                                        nsgf, n_occ(ispin), alpha=alpha, beta=beta)
!
!               CALL cp_fm_release(matrix=Xdoubletilde)
!               CALL cp_fm_release(matrix=tmp_res)
!               NULLIFY (Xdoubletilde, tmp_res)
!-------------------------------------------------------------------------------------------------------------------!
!!!! TEST explicit loops
      DO ii=1,n_occ(ispin)
       DO aa=1,nsgf

          tmp_ex=0.0_dp
          tmp_coul=0.0_dp
          tmp1=0.0_dp; tmp2=0.0_dp

! Exchange-like ! 
      IF ((nspins == 1 .AND. res_etype == tddfpt_singlet ).OR.( nspins == 2 )) THEN
          tmp=0.0_dp
           DO jspin = 1, nspins
            DO jj=1,n_occ(jspin)
             DO bb=first_sgf(jatom),first_sgf(jatom)-1+basis_set_jatom%nsgf
              CALL cp_fm_get_element(stda_env%ctransformed(jspin)%matrix,bb,jj,tmp1) 
              CALL cp_fm_get_element(stda_env%xtransformed(jspin)%matrix,bb,jj,tmp2) 
              tmp=tmp + tmp1*tmp2  
             ENDDO
            ENDDO 
           ENDDO
         
          CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%ctransformed_double(ispin)%matrix,&
               aa,ii,tmp3) 

          tmp_ex= tmp*(spinfac*exchangeblock(1,1))*tmp3
          CALL cp_fm_get_element(res(ispin)%matrix,aa,ii,tmp)
          CALL cp_fm_set_element(res(ispin)%matrix,aa,ii,tmp+tmp_ex)
         ENDIF

! Coulomb-like
       !  IF (ispin.eq.jspin) THEN
          tmp=0.0_dp
          DO jj=1,n_occ(ispin)
           CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%xtransformed_double(ispin)%matrix,&
                 aa,jj,tmp2) 
           CALL cp_fm_get_element(stda_env%atom_param_set(jatom)%atom_param%charges(ispin)%matrix,&
                 jj,ii,tmp1) 
           tmp=tmp + tmp1*tmp2
          ENDDO

          tmp_coul = -coulombblock(1,1)*tmp

!         WRITE(output_unit,*) '-------------------------------------------------'
!         WRITE(UNIT=output_unit,FMT="(A,F18.6,F18.6,F18.6)") "tmp_coulomb, tmp", tmp_coul ,tmp,coulombblock(1,1)
!         WRITE(output_unit,*) '-------------------------------------------------'

          CALL cp_fm_get_element(res(ispin)%matrix,aa,ii,tmp)
          CALL cp_fm_set_element(res(ispin)%matrix,aa,ii,tmp+tmp_coul)
       !  ENDIF
       ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!IATOM.NE.JATOM
       IF (iatom.NE.jatom) THEN
         DO aa=1,nsgf

          tmp_ex=0.0_dp
          tmp_coul=0.0_dp
          tmp1=0.0_dp; tmp2=0.0_dp

! Exchange-like 
         IF ((nspins == 1 .AND. res_etype == tddfpt_singlet ).OR.( nspins == 2 )) THEN
          tmp=0.0_dp
           DO jspin=1, nspins
            DO jj=1,n_occ(jspin)
             DO bb=first_sgf(iatom),first_sgf(iatom)-1+basis_set_iatom%nsgf
              CALL cp_fm_get_element(stda_env%ctransformed(jspin)%matrix,bb,jj,tmp1) 
              CALL cp_fm_get_element(stda_env%xtransformed(jspin)%matrix,bb,jj,tmp2)
              tmp=tmp + tmp1*tmp2
             ENDDO
            ENDDO
           ENDDO

          CALL cp_fm_get_element(stda_env%atom_param_set(jatom)%atom_param%ctransformed_double(ispin)%matrix,&
               aa,ii,tmp3)

          tmp_ex= tmp*(spinfac*exchangeblock(1,1))*tmp3
          CALL cp_fm_get_element(res(ispin)%matrix,aa,ii,tmp)
          CALL cp_fm_set_element(res(ispin)%matrix,aa,ii,tmp+tmp_ex)
         ENDIF

! Coulomb-like
       !  IF (ispin.eq.jspin) THEN
          tmp=0.0_dp
          DO jj=1,n_occ(ispin)
           CALL cp_fm_get_element(stda_env%atom_param_set(jatom)%atom_param%xtransformed_double(ispin)%matrix,&
                 aa,jj,tmp2)
           CALL cp_fm_get_element(stda_env%atom_param_set(iatom)%atom_param%charges(ispin)%matrix,&
                 jj,ii,tmp1)
           tmp=tmp + tmp1*tmp2
          ENDDO

          tmp_coul = -coulombblock(1,1)*tmp

!         WRITE(output_unit,*) '-------------------------------------------------'
!         WRITE(UNIT=output_unit,FMT="(A,F18.6,F18.6,F18.6)") "tmp_coulomb, tmp", tmp_coul ,tmp,coulombblock(1,1)
!         WRITE(output_unit,*) '-------------------------------------------------'

          CALL cp_fm_get_element(res(ispin)%matrix,aa,ii,tmp)
          CALL cp_fm_set_element(res(ispin)%matrix,aa,ii,tmp+tmp_coul)
     !    ENDIF
      ENDDO !!aa

     ENDIF !iatom.NE.jatom
      ENDDO !!ii
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!IATOM.NE.JATOM
!!!! TEST explicit loops
!-------------------------------------------------------------------------------------------------------------------!
            !IF (iatom.NE.jatom) -----------------------------------------------------------------
         ENDDO !nl_iterator
      ENDDO !ispin
!     ENDDO !jspin

! Ewald sum for Coulomb-like contribution
      IF (do_periodic) THEN
            CALL ewald_or_direct_sum(qs_env, p_env, stda_env, res)
      ENDIF
!
      DEALLOCATE (kind_of)
      DEALLOCATE (first_sgf)

      CALL neighbor_list_iterator_release(nl_iterator)

      CALL timestop(handle)

   END SUBROUTINE stda_calculate_kernel

END MODULE qs_tddfpt_stda_utils
