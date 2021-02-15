!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2021 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Calculate the CPKS equation and the resulting forces
!> \par History
!>       03.2014 created
!>       09.2019 Moved from KG to Kohn-Sham
!>       11.2019 Moved from energy_correction
!> \author JGH
! **************************************************************************************************
MODULE response_solver
   USE admm_types,                      ONLY: admm_type
   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE core_ppl,                        ONLY: build_core_ppl
   USE core_ppnl,                       ONLY: build_core_ppnl
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type
   USE cp_dbcsr_cp2k_link,              ONLY: cp_dbcsr_alloc_block_from_nbl
   USE cp_dbcsr_operations,             ONLY: copy_dbcsr_to_fm,&
                                              copy_fm_to_dbcsr,&
                                              cp_dbcsr_sm_fm_multiply,&
                                              dbcsr_allocate_matrix_set,&
                                              dbcsr_deallocate_matrix_set
   USE cp_fm_struct,                    ONLY: cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_set_all,&
                                              cp_fm_type
   USE cp_gemm_interface,               ONLY: cp_gemm
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_unit_nr,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: &
        dbcsr_add, dbcsr_copy, dbcsr_create, dbcsr_deallocate_matrix, dbcsr_distribution_type, &
        dbcsr_p_type, dbcsr_release, dbcsr_scale, dbcsr_set, dbcsr_type, dbcsr_type_no_symmetry
   USE exstates_types,                  ONLY: excited_energy_type
   USE hfx_derivatives,                 ONLY: derivatives_four_center
   USE hfx_energy_potential,            ONLY: integrate_four_center
   USE hfx_types,                       ONLY: hfx_type
   USE input_constants,                 ONLY: do_admm_aux_exch_func_none,&
                                              kg_tnadd_atomic,&
                                              kg_tnadd_embed,&
                                              ot_precond_full_single_inverse
   USE input_section_types,             ONLY: section_get_ival,&
                                              section_get_rval,&
                                              section_vals_get,&
                                              section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kg_correction,                   ONLY: kg_ekin_subset
   USE kg_environment_types,            ONLY: kg_environment_type
   USE kg_tnadd_mat,                    ONLY: build_tnadd_mat
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_sum
   USE mulliken,                        ONLY: ao_charges
   USE particle_types,                  ONLY: particle_type
   USE pw_env_types,                    ONLY: pw_env_get,&
                                              pw_env_type
   USE pw_grid_types,                   ONLY: pw_grid_type
   USE pw_methods,                      ONLY: pw_axpy,&
                                              pw_scale,&
                                              pw_transfer,&
                                              pw_zero
   USE pw_poisson_methods,              ONLY: pw_poisson_solve
   USE pw_poisson_types,                ONLY: pw_poisson_type
   USE pw_pool_types,                   ONLY: pw_pool_create_pw,&
                                              pw_pool_give_back_pw,&
                                              pw_pool_type
   USE pw_types,                        ONLY: COMPLEXDATA1D,&
                                              REALDATA3D,&
                                              REALSPACE,&
                                              RECIPROCALSPACE,&
                                              pw_create,&
                                              pw_p_type,&
                                              pw_release
   USE qs_collocate_density,            ONLY: calculate_rho_elec
   USE qs_core_energies,                ONLY: calculate_ecore_overlap,&
                                              calculate_ecore_self
   USE qs_density_matrices,             ONLY: calculate_whz_matrix,&
                                              calculate_wz_matrix
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type,&
                                              set_qs_env
   USE qs_force_types,                  ONLY: add_qs_force,&
                                              qs_force_type,&
                                              total_qs_force
   USE qs_integrate_potential,          ONLY: integrate_v_core_rspace,&
                                              integrate_v_rspace
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_kinetic,                      ONLY: build_kinetic_matrix
   USE qs_ks_methods,                   ONLY: calc_rho_tot_gspace
   USE qs_ks_types,                     ONLY: qs_ks_env_type
   USE qs_linres_methods,               ONLY: build_dm_response,&
                                              linres_solver
   USE qs_linres_types,                 ONLY: linres_control_create,&
                                              linres_control_release,&
                                              linres_control_type
   USE qs_mo_types,                     ONLY: get_mo_set,&
                                              mo_set_p_type
   USE qs_neighbor_list_types,          ONLY: neighbor_list_set_p_type
   USE qs_overlap,                      ONLY: build_overlap_force,&
                                              build_overlap_matrix
   USE qs_p_env_methods,                ONLY: p_env_create,&
                                              p_env_psi0_changed
   USE qs_p_env_types,                  ONLY: p_env_release,&
                                              qs_p_env_type
   USE qs_rho_types,                    ONLY: qs_rho_get,&
                                              qs_rho_type
   USE qs_vxc,                          ONLY: qs_vxc_create
   USE task_list_types,                 ONLY: task_list_type
   USE virial_types,                    ONLY: virial_type
   USE xc,                              ONLY: xc_calc_2nd_deriv,&
                                              xc_prep_2nd_deriv
   USE xc_derivative_set_types,         ONLY: xc_derivative_set_type,&
                                              xc_dset_release
   USE xc_derivatives,                  ONLY: xc_functionals_get_needs
   USE xc_rho_cflags_types,             ONLY: xc_rho_cflags_type
   USE xc_rho_set_types,                ONLY: xc_rho_set_create,&
                                              xc_rho_set_release,&
                                              xc_rho_set_type,&
                                              xc_rho_set_update
   USE xtb_ehess,                       ONLY: xtb_coulomb_hessian
   USE xtb_ehess_force,                 ONLY: calc_xtb_ehess_force
   USE xtb_matrices,                    ONLY: xtb_hab_force
   USE xtb_types,                       ONLY: get_xtb_atom_param,&
                                              xtb_atom_type
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

! *** Global parameters ***

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'response_solver'

   LOGICAL, PARAMETER                   :: debug_forces = .TRUE.

   PUBLIC ::  response_equation, response_force, response_force_xtb, &
             ks_ref_potential, admm_projection_derivative

! **************************************************************************************************

CONTAINS

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param cpmos RHS of equation as Ax + b = 0 (sign of b)
!> \param p1_mat ...
!> \param lr_section ...
! **************************************************************************************************
   SUBROUTINE response_equation(qs_env, p_env, cpmos, p1_mat, lr_section)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: cpmos
      TYPE(cp_fm_p_type), DIMENSION(:), OPTIONAL, &
         POINTER                                         :: p1_mat
      TYPE(section_vals_type), OPTIONAL, POINTER         :: lr_section

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'response_equation'

      INTEGER                                            :: handle, iounit, ispin, nao, nao_aux, &
                                                            nspins
      LOGICAL                                            :: should_stop
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER          :: psi0, psi1
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: mo_coeff
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_ks, matrix_s, matrix_s_aux
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(linres_control_type), POINTER                 :: linres_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      IF (logger%para_env%ionode) THEN
         iounit = cp_logger_get_default_unit_nr(logger, local=.TRUE.)
      ELSE
         iounit = -1
      ENDIF

      ! initialized linres_control
      NULLIFY (linres_control)
      CALL linres_control_create(linres_control)
      linres_control%do_kernel = .TRUE.
      linres_control%lr_triplet = .FALSE.
      IF (PRESENT(lr_section)) THEN
         CALL section_vals_val_get(lr_section, "RESTART", l_val=linres_control%linres_restart)
         CALL section_vals_val_get(lr_section, "MAX_ITER", i_val=linres_control%max_iter)
         CALL section_vals_val_get(lr_section, "EPS", r_val=linres_control%eps)
         CALL section_vals_val_get(lr_section, "EPS_FILTER", r_val=linres_control%eps_filter)
         CALL section_vals_val_get(lr_section, "RESTART_EVERY", i_val=linres_control%restart_every)
         CALL section_vals_val_get(lr_section, "PRECONDITIONER", i_val=linres_control%preconditioner_type)
         CALL section_vals_val_get(lr_section, "ENERGY_GAP", r_val=linres_control%energy_gap)
      ELSE
         linres_control%linres_restart = .FALSE.
         linres_control%max_iter = 100
         linres_control%eps = 1.0e-10_dp
         linres_control%eps_filter = 1.0e-15_dp
         linres_control%restart_every = 50
         linres_control%preconditioner_type = ot_precond_full_single_inverse
         linres_control%energy_gap = 0.02_dp
      END IF

      ! initialized p_env
      CALL p_env_release(p_env)
      CALL p_env_create(p_env, qs_env, orthogonal_orbitals=.TRUE., &
                        linres_control=linres_control)
      CALL set_qs_env(qs_env, linres_control=linres_control)
      CALL p_env_psi0_changed(p_env, qs_env)
      p_env%os_valid = .FALSE.
      p_env%new_preconditioner = .TRUE.
      !
      CALL get_qs_env(qs_env, dft_control=dft_control, mos=mos)
      nspins = dft_control%nspins

      ALLOCATE (psi0(nspins), psi1(nspins))
      DO ispin = 1, nspins
         CALL get_mo_set(mo_set=mos(ispin)%mo_set, mo_coeff=mo_coeff)
         psi0(ispin)%matrix => mo_coeff
         CALL cp_fm_get_info(mo_coeff, matrix_struct=fm_struct)
         NULLIFY (psi1(ispin)%matrix)
         CALL cp_fm_create(psi1(ispin)%matrix, fm_struct)
         CALL cp_fm_set_all(psi1(ispin)%matrix, 0.0_dp)
      ENDDO

      CALL get_qs_env(qs_env, matrix_s=matrix_s, sab_orb=sab_orb)
      CALL dbcsr_allocate_matrix_set(p_env%p1, nspins)
      CALL dbcsr_allocate_matrix_set(p_env%w1, nspins)
      DO ispin = 1, nspins
         ALLOCATE (p_env%p1(ispin)%matrix, p_env%w1(ispin)%matrix)
         CALL dbcsr_create(matrix=p_env%p1(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL dbcsr_create(matrix=p_env%w1(ispin)%matrix, template=matrix_s(1)%matrix)
         CALL cp_dbcsr_alloc_block_from_nbl(p_env%p1(ispin)%matrix, sab_orb)
         CALL cp_dbcsr_alloc_block_from_nbl(p_env%w1(ispin)%matrix, sab_orb)
      END DO
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, matrix_s_aux_fit=matrix_s_aux)
         CALL dbcsr_allocate_matrix_set(p_env%p1_admm, nspins)
         DO ispin = 1, nspins
            ALLOCATE (p_env%p1_admm(ispin)%matrix)
            CALL dbcsr_create(p_env%p1_admm(ispin)%matrix, &
                              template=matrix_s_aux(1)%matrix)
            CALL dbcsr_copy(p_env%p1_admm(ispin)%matrix, matrix_s_aux(1)%matrix)
            CALL dbcsr_set(p_env%p1_admm(ispin)%matrix, 0.0_dp)
         END DO
      END IF

      CALL linres_solver(p_env, qs_env, psi1, cpmos, psi0, iounit, should_stop)

      DO ispin = 1, nspins
         CALL dbcsr_copy(p_env%p1(ispin)%matrix, matrix_s(1)%matrix)
      END DO
      CALL build_dm_response(psi0, psi1, p_env%p1)
      DO ispin = 1, nspins
         CALL dbcsr_scale(p_env%p1(ispin)%matrix, 0.5_dp)
      END DO
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         CPASSERT(ASSOCIATED(admm_env%work_orb_orb))
         CPASSERT(ASSOCIATED(admm_env%work_aux_orb))
         CPASSERT(ASSOCIATED(admm_env%work_aux_aux))
         nao = admm_env%nao_orb
         nao_aux = admm_env%nao_aux_fit
         DO ispin = 1, nspins
            CALL copy_dbcsr_to_fm(p_env%p1(ispin)%matrix, admm_env%work_orb_orb)
            CALL cp_gemm('N', 'N', nao_aux, nao, nao, &
                         1.0_dp, admm_env%A, admm_env%work_orb_orb, 0.0_dp, &
                         admm_env%work_aux_orb)
            CALL cp_gemm('N', 'T', nao_aux, nao_aux, nao, &
                         1.0_dp, admm_env%work_aux_orb, admm_env%A, 0.0_dp, &
                         admm_env%work_aux_aux)
            CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, p_env%p1_admm(ispin)%matrix, &
                                  keep_sparsity=.TRUE.)
         END DO
      END IF
      CALL get_qs_env(qs_env, matrix_ks=matrix_ks)
      DO ispin = 1, nspins
         CALL calculate_wz_matrix(mos(ispin)%mo_set, psi1(ispin)%matrix, matrix_ks(ispin)%matrix, &
                                  p_env%w1(ispin)%matrix)
      ENDDO
      IF (PRESENT(p1_mat)) THEN
         CPASSERT(.NOT. ASSOCIATED(p1_mat))
         p1_mat => psi1
      ELSE
         DO ispin = 1, nspins
            CALL cp_fm_release(psi1(ispin)%matrix)
         ENDDO
         DEALLOCATE (psi1)
      END IF
      DEALLOCATE (psi0)
      CALL linres_control_release(linres_control)

      CALL timestop(handle)

   END SUBROUTINE response_equation

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param vh_rspace ...
!> \param vxc_rspace ...
!> \param vtau_rspace ...
!> \param matrix_hz ...
!> \param ex_env ...
! **************************************************************************************************
   SUBROUTINE response_force(qs_env, p_env, vh_rspace, vxc_rspace, vtau_rspace, matrix_hz, &
                             ex_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(pw_p_type), POINTER                           :: vh_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rspace, vtau_rspace
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz
      TYPE(excited_energy_type), OPTIONAL, POINTER       :: ex_env

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'response_force'

      INTEGER                                            :: handle, iounit, ispin, mspin, n_rep_hf, &
                                                            nao, nao_aux, natom, nder, nimages, &
                                                            nspins
      INTEGER, DIMENSION(2, 3)                           :: bo
      INTEGER, DIMENSION(:, :, :), POINTER               :: cell_to_index
      LOGICAL                                            :: distribute_fock_matrix, do_ex, do_hfx, &
                                                            hfx_treat_lsd_in_core, lsd, resp_only, &
                                                            s_mstruct_changed, use_virial
      REAL(KIND=dp)                                      :: eh1, ekin_mol, eps_ppnl, exc, focc, &
                                                            total_rhoz, total_rhoz_aux, zehartree
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: ftot1, ftot2, ftot3
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_distribution_type), POINTER             :: dbcsr_dist
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pd, matrix_pz, matrix_wz, mpa, &
                                                            scrm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_h, matrix_p, mhd, mhx, mhy, mhz, &
                                                            mpd, mpz
      TYPE(dbcsr_type), POINTER                          :: dbwork
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(hfx_type), DIMENSION(:, :), POINTER           :: x_data
      TYPE(kg_environment_type), POINTER                 :: kg_env
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb, sac_ppl, sap_ppnl
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_p_type)                                    :: rhoz_tot_gspace, vhxc_rspace, &
                                                            zv_hartree_gspace, zv_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER :: rho_g, rho_g_aux, rho_r, rho_r_aux, rhoz_g, &
         rhoz_g_aux, rhoz_r, rhoz_r_aux, tau_r, tau_r_aux, v_rspace, v_tau_rspace, v_xc
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho, rho_aux_fit
      TYPE(section_vals_type), POINTER                   :: hfx_section, xc_fun_section, xc_section
      TYPE(task_list_type), POINTER                      :: task_list_aux_fit
      TYPE(virial_type), POINTER                         :: virial
      TYPE(xc_derivative_set_type), POINTER              :: deriv_set
      TYPE(xc_rho_cflags_type)                           :: needs
      TYPE(xc_rho_set_type), POINTER                     :: rho1_set, rho_set

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_unit_nr(logger)

      do_ex = .FALSE.
      IF (PRESENT(ex_env)) do_ex = .TRUE.

      NULLIFY (ks_env, sab_orb, sac_ppl, sap_ppnl)
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, dft_control=dft_control, &
                      sab_orb=sab_orb, sac_ppl=sac_ppl, sap_ppnl=sap_ppnl)
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, force=force)
      nspins = dft_control%nspins

      IF (debug_forces) THEN
         CALL get_qs_env(qs_env, natom=natom, atomic_kind_set=atomic_kind_set)
         ALLOCATE (ftot1(3, natom))
         CALL total_qs_force(ftot1, force, atomic_kind_set)
      END IF

      matrix_pz => p_env%p1
      NULLIFY (mpa)
      IF (do_ex) THEN
         CALL dbcsr_allocate_matrix_set(mpa, nspins)
         DO ispin = 1, nspins
            ALLOCATE (mpa(ispin)%matrix)
            CALL dbcsr_create(mpa(ispin)%matrix, template=matrix_pz(ispin)%matrix)
            CALL dbcsr_copy(mpa(ispin)%matrix, matrix_pz(ispin)%matrix)
            CALL dbcsr_add(mpa(ispin)%matrix, ex_env%matrix_pe(ispin)%matrix, 1.0_dp, 1.0_dp)
            CALL dbcsr_set(matrix_hz(ispin)%matrix, 0.0_dp)
         END DO
      ELSE
         mpa => p_env%p1
      END IF

      !
      ! START OF Tr(P+Z)Hcore
      !
      IF (nspins == 2) THEN
         CALL dbcsr_add(mpa(1)%matrix, mpa(2)%matrix, 1.0_dp, 1.0_dp)
      END IF
      ! Kinetic energy matrix
      ! This routine initializes a matrx scrm, later used for scratch purposes
      NULLIFY (scrm)
      IF (debug_forces) fodeb(1:3) = force(1)%kinetic(1:3, 1)
      CALL build_kinetic_matrix(ks_env, matrix_t=scrm, &
                                matrix_name="KINETIC ENERGY MATRIX", &
                                basis_type="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=mpa(1)%matrix)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%kinetic(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dT      ", fodeb
      END IF
      IF (nspins == 2) THEN
         CALL dbcsr_add(mpa(1)%matrix, mpa(2)%matrix, 1.0_dp, -1.0_dp)
      END IF

      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set, particle_set=particle_set, &
                      atomic_kind_set=atomic_kind_set)
      NULLIFY (cell_to_index, virial)
      use_virial = .FALSE.
      nimages = 1
      ALLOCATE (matrix_p(nspins, 1), matrix_h(nspins, 1))
      DO ispin = 1, nspins
         matrix_p(ispin, 1)%matrix => mpa(ispin)%matrix
         matrix_h(ispin, 1)%matrix => scrm(ispin)%matrix
      END DO
      IF (ASSOCIATED(sac_ppl)) THEN
         nder = 1
         IF (debug_forces) fodeb(1:3) = force(1)%gth_ppl(1:3, 1)
         CALL build_core_ppl(matrix_h, matrix_p, force, &
                             virial, .TRUE., use_virial, nder, &
                             qs_kind_set, atomic_kind_set, particle_set, sab_orb, sac_ppl, &
                             nimages, cell_to_index, "ORB")
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppl(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHppl   ", fodeb
         END IF
      END IF
      eps_ppnl = dft_control%qs_control%eps_ppnl
      IF (ASSOCIATED(sap_ppnl)) THEN
         nder = 1
         IF (debug_forces) fodeb(1:3) = force(1)%gth_ppnl(1:3, 1)
         CALL build_core_ppnl(matrix_h, matrix_p, force, &
                              virial, .TRUE., use_virial, nder, &
                              qs_kind_set, atomic_kind_set, particle_set, sab_orb, sap_ppnl, eps_ppnl, &
                              nimages, cell_to_index, "ORB")
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%gth_ppnl(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHppnl  ", fodeb
         END IF
      END IF
      ! Kim-Gordon subsystem DFT
      ! Atomic potential for nonadditive kinetic energy contribution
      IF (dft_control%qs_control%do_kg) THEN
         IF (qs_env%kg_env%tnadd_method == kg_tnadd_atomic) THEN
            CALL get_qs_env(qs_env=qs_env, kg_env=kg_env, dbcsr_dist=dbcsr_dist)
            IF (debug_forces) fodeb(1:3) = force(1)%kinetic(1:3, 1)
            CALL build_tnadd_mat(kg_env=kg_env, matrix_p=matrix_p, force=force, virial=virial, &
                                 calculate_forces=.TRUE., use_virial=use_virial, &
                                 qs_kind_set=qs_kind_set, atomic_kind_set=atomic_kind_set, &
                                 particle_set=particle_set, sab_orb=sab_orb, dbcsr_dist=dbcsr_dist)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%kinetic(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dTnadd  ", fodeb
            END IF
         END IF
      END IF
      DEALLOCATE (matrix_p, matrix_h)
      IF (debug_forces) THEN
         ALLOCATE (ftot2(3, natom))
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1) - ftot1(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Force Pz*dHcore", fodeb
      END IF
      !
      ! END OF Tr(P+Z)Hcore
      !

      !
      ! Vhxc (KS potentials calculated externally)
      !
      CALL get_qs_env(qs_env, pw_env=pw_env)
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)
      CALL pw_pool_create_pw(auxbas_pw_pool, vhxc_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      DO ispin = 1, nspins
         CALL pw_transfer(vh_rspace%pw, vhxc_rspace%pw)
         CALL pw_axpy(vxc_rspace(ispin)%pw, vhxc_rspace%pw)
         CALL integrate_v_rspace(v_rspace=vhxc_rspace, &
                                 hmat=scrm(1), pmat=mpa(ispin), &
                                 qs_env=qs_env, calculate_forces=.TRUE.)
         IF (ASSOCIATED(vtau_rspace)) THEN
            CALL integrate_v_rspace(v_rspace=vtau_rspace(ispin), &
                                    hmat=scrm(1), pmat=mpa(ispin), &
                                    qs_env=qs_env, calculate_forces=.TRUE., compute_tau=.TRUE.)
         END IF
      END DO
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dVhxc   ", fodeb
      END IF
      CALL pw_pool_give_back_pw(auxbas_pw_pool, vhxc_rspace%pw)

      ! KG Embedding
      ! calculate kinetic energy potential and integrate with response density
      IF (dft_control%qs_control%do_kg) THEN
         IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed) THEN
            CALL get_qs_env(qs_env, kg_env=kg_env)
            ekin_mol = 0.0_dp
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            CALL kg_ekin_subset(qs_env=qs_env, &
                                ks_matrix=scrm, &
                                ekin_mol=ekin_mol, &
                                calc_force=.TRUE., &
                                do_kernel=.FALSE., &
                                pmat_ext=mpa)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dVkg   ", fodeb
            END IF
         END IF
      END IF

      !
      ! Hartree potential of response density
      !
      ALLOCATE (rhoz_r(nspins), rhoz_g(nspins))
      DO ispin = 1, nspins
         NULLIFY (rhoz_r(ispin)%pw, rhoz_g(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_r(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_g(ispin)%pw, &
                                use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      ENDDO
      CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, zv_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, zv_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      CALL pw_zero(rhoz_tot_gspace%pw)
      DO ispin = 1, nspins
         CALL calculate_rho_elec(ks_env=ks_env, matrix_p=mpa(ispin)%matrix, &
                                 rho=rhoz_r(ispin), rho_gspace=rhoz_g(ispin), &
                                 total_rho=total_rhoz)
         CALL pw_axpy(rhoz_g(ispin)%pw, rhoz_tot_gspace%pw)
      END DO
      ! calculate associated hartree potential
      CALL pw_poisson_solve(poisson_env, rhoz_tot_gspace%pw, zehartree, &
                            zv_hartree_gspace%pw)
      CALL pw_transfer(zv_hartree_gspace%pw, zv_hartree_rspace%pw)
      CALL pw_scale(zv_hartree_rspace%pw, zv_hartree_rspace%pw%pw_grid%dvol)
      ! Getting nuclear force contribution from the core charge density
      IF (debug_forces) fodeb(1:3) = force(1)%rho_core(1:3, 1)
      CALL integrate_v_core_rspace(zv_hartree_rspace, qs_env)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_core(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Vh(rhoz)*dncore ", fodeb
      END IF
      !
      CALL get_qs_env(qs_env=qs_env, rho=rho)
      CALL qs_rho_get(rho, rho_r=rho_r, rho_g=rho_g, tau_r=tau_r)
      ALLOCATE (v_xc(nspins))
      DO ispin = 1, nspins
         NULLIFY (v_xc(ispin)%pw)
         CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                use_data=REALDATA3D, in_space=REALSPACE)
         CALL pw_zero(v_xc(ispin)%pw)
      END DO
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      END IF
      lsd = (nspins == 2)
      NULLIFY (deriv_set, rho_set, rho1_set)
      CALL xc_prep_2nd_deriv(deriv_set, rho_set, rho_r, auxbas_pw_pool, xc_section=xc_section)
      bo = rhoz_r(1)%pw%pw_grid%bounds_local
      CALL xc_rho_set_create(rho1_set, bo, &
                             rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                             drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                             tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

      xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
      needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

      ! calculate the arguments needed by the functionals
      CALL xc_rho_set_update(rho1_set, rhoz_r, rhoz_g, tau_r, needs, &
                             section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                             section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                             auxbas_pw_pool)
      CALL xc_calc_2nd_deriv(v_xc, deriv_set, rho_set, &
                             rho1_set, auxbas_pw_pool, xc_section=xc_section)
      CALL xc_dset_release(deriv_set)
      CALL xc_rho_set_release(rho_set)
      CALL xc_rho_set_release(rho1_set)
      !
      !
      CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
      IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
      DO ispin = 1, nspins
         CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
         CALL pw_axpy(zv_hartree_rspace%pw, v_xc(ispin)%pw)
         CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                 hmat=matrix_hz(ispin), &
                                 pmat=matrix_p(ispin, 1), &
                                 calculate_forces=.TRUE.)
      END DO
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dK*rhoz ", fodeb
      END IF

      ! KG Embedding
      ! calculate kinetic energy kernel, folded with response density for partial integration0
      IF (dft_control%qs_control%do_kg) THEN
         IF (qs_env%kg_env%tnadd_method == kg_tnadd_embed) THEN
            CALL get_qs_env(qs_env, kg_env=kg_env)
            ekin_mol = 0.0_dp
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            CALL kg_ekin_subset(qs_env=qs_env, &
                                ks_matrix=matrix_hz, &
                                ekin_mol=ekin_mol, &
                                calc_force=.TRUE., &
                                do_kernel=.TRUE., &
                                pmat_ext=matrix_pz)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*d(Kkg)*rhoz ", fodeb
            END IF
         END IF
      END IF

      CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_tot_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, zv_hartree_rspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, zv_hartree_gspace%pw)
      DO ispin = 1, nspins
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_r(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_g(ispin)%pw)
         CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
      END DO
      DEALLOCATE (rhoz_r, rhoz_g, v_xc)
      IF (debug_forces) THEN
         ALLOCATE (ftot3(3, natom))
         CALL total_qs_force(ftot3, force, atomic_kind_set)
         fodeb(1:3) = ftot3(1:3, 1) - ftot2(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Force Pin*V(rhoz)", fodeb
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)

      IF (dft_control%do_admm) THEN
         IF (dft_control%admm_control%aux_exch_func == do_admm_aux_exch_func_none) THEN
            ! nothing to do
            NULLIFY (mpz, mhx, mhy)
         ELSE
            ! add ADMM xc_section_aux terms: Pz*Vxc + P0*K0[rhoz]
            CALL get_qs_env(qs_env, admm_env=admm_env, rho_aux_fit=rho_aux_fit, &
                            matrix_s_aux_fit=scrm, task_list_aux_fit=task_list_aux_fit)
            !
            NULLIFY (mpz, mhx, mhy)
            CALL dbcsr_allocate_matrix_set(mhx, nspins, 1)
            CALL dbcsr_allocate_matrix_set(mhy, nspins, 1)
            CALL dbcsr_allocate_matrix_set(mpz, nspins, 1)
            DO ispin = 1, nspins
               ALLOCATE (mhx(ispin, 1)%matrix)
               CALL dbcsr_create(mhx(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_copy(mhx(ispin, 1)%matrix, scrm(1)%matrix)
               CALL dbcsr_set(mhx(ispin, 1)%matrix, 0.0_dp)
               ALLOCATE (mhy(ispin, 1)%matrix)
               CALL dbcsr_create(mhy(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_copy(mhy(ispin, 1)%matrix, scrm(1)%matrix)
               CALL dbcsr_set(mhy(ispin, 1)%matrix, 0.0_dp)
               ALLOCATE (mpz(ispin, 1)%matrix)
               CALL dbcsr_create(mpz(ispin, 1)%matrix, template=p_env%p1_admm(ispin)%matrix)
               CALL dbcsr_copy(mpz(ispin, 1)%matrix, p_env%p1_admm(ispin)%matrix)
               IF (do_ex) THEN
                  CALL dbcsr_add(mpz(ispin, 1)%matrix, ex_env%matrix_pe_admm(ispin)%matrix, &
                                 1.0_dp, 1.0_dp)
               END IF
            END DO
            !
            xc_section => admm_env%xc_section_aux
            NULLIFY (v_rspace, v_tau_rspace)
            CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho_aux_fit, xc_section=xc_section, &
                               vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE.)
            CPASSERT(.NOT. ASSOCIATED(v_tau_rspace))
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            DO ispin = 1, nspins
               v_rspace(ispin)%pw%cr3d = v_rspace(ispin)%pw%pw_grid%dvol*v_rspace(ispin)%pw%cr3d
               CALL integrate_v_rspace(v_rspace=v_rspace(ispin), &
                                       hmat=mhx(ispin, 1), pmat=mpz(ispin, 1), &
                                       qs_env=qs_env, calculate_forces=.TRUE., &
                                       basis_type="AUX_FIT", task_list_external=task_list_aux_fit)
            END DO
            DO ispin = 1, nspins
               CALL pw_release(v_rspace(ispin)%pw)
            END DO
            DEALLOCATE (v_rspace)
            NULLIFY (rho_g_aux, rho_r_aux, tau_r_aux, rhoz_g_aux, rhoz_r_aux)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*Vxc(rho_admm)", fodeb
            END IF
            !
            CALL qs_rho_get(rho_aux_fit, rho_r=rho_r_aux, rho_g=rho_g_aux, tau_r=tau_r_aux)
            ALLOCATE (v_xc(nspins))
            DO ispin = 1, nspins
               NULLIFY (v_xc(ispin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, v_xc(ispin)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_zero(v_xc(ispin)%pw)
            END DO
            lsd = (nspins == 2)
            ! rhoz_aux
            ALLOCATE (rhoz_r_aux(nspins), rhoz_g_aux(nspins))
            DO ispin = 1, nspins
               NULLIFY (rhoz_r_aux(ispin)%pw, rhoz_g_aux(ispin)%pw)
               CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_r_aux(ispin)%pw, &
                                      use_data=REALDATA3D, in_space=REALSPACE)
               CALL pw_pool_create_pw(auxbas_pw_pool, rhoz_g_aux(ispin)%pw, &
                                      use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
            ENDDO
            DO ispin = 1, nspins
               CALL calculate_rho_elec(ks_env=ks_env, matrix_p=mpz(ispin, 1)%matrix, &
                                       rho=rhoz_r_aux(ispin), rho_gspace=rhoz_g_aux(ispin), &
                                       total_rho=total_rhoz_aux, basis_type="AUX_FIT", &
                                       task_list_external=task_list_aux_fit)
            END DO
            !
            NULLIFY (deriv_set, rho_set, rho1_set)
            CALL xc_prep_2nd_deriv(deriv_set, rho_set, rho_r_aux, auxbas_pw_pool, xc_section=xc_section)
            bo = rhoz_r_aux(1)%pw%pw_grid%bounds_local
            CALL xc_rho_set_create(rho1_set, bo, &
                                   rho_cutoff=section_get_rval(xc_section, "DENSITY_CUTOFF"), &
                                   drho_cutoff=section_get_rval(xc_section, "GRADIENT_CUTOFF"), &
                                   tau_cutoff=section_get_rval(xc_section, "TAU_CUTOFF"))

            xc_fun_section => section_vals_get_subs_vals(xc_section, "XC_FUNCTIONAL")
            needs = xc_functionals_get_needs(xc_fun_section, lsd, .TRUE.)

            ! calculate the arguments needed by the functionals
            CALL xc_rho_set_update(rho1_set, rhoz_r_aux, rhoz_g_aux, tau_r_aux, needs, &
                                   section_get_ival(xc_section, "XC_GRID%XC_DERIV"), &
                                   section_get_ival(xc_section, "XC_GRID%XC_SMOOTH_RHO"), &
                                   auxbas_pw_pool)
            CALL xc_calc_2nd_deriv(v_xc, deriv_set, rho_set, &
                                   rho1_set, auxbas_pw_pool, xc_section=xc_section)
            CALL xc_dset_release(deriv_set)
            CALL xc_rho_set_release(rho_set)
            CALL xc_rho_set_release(rho1_set)
            !
            CALL qs_rho_get(rho_aux_fit, rho_ao_kp=matrix_p)
            IF (debug_forces) fodeb(1:3) = force(1)%rho_elec(1:3, 1)
            DO ispin = 1, nspins
               CALL dbcsr_set(mhy(ispin, 1)%matrix, 0.0_dp)
               CALL pw_scale(v_xc(ispin)%pw, v_xc(ispin)%pw%pw_grid%dvol)
               CALL integrate_v_rspace(qs_env=qs_env, v_rspace=v_xc(ispin), &
                                       hmat=mhy(ispin, 1), pmat=matrix_p(ispin, 1), &
                                       calculate_forces=.TRUE., basis_type="AUX_FIT", &
                                       task_list_external=task_list_aux_fit)
            END DO
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%rho_elec(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pin*dK*rhoz_admm ", fodeb
            END IF
            DO ispin = 1, nspins
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_xc(ispin)%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_r_aux(ispin)%pw)
               CALL pw_pool_give_back_pw(auxbas_pw_pool, rhoz_g_aux(ispin)%pw)
            END DO
            DEALLOCATE (v_xc, rhoz_r_aux, rhoz_g_aux)
            !
            nao = admm_env%nao_orb
            nao_aux = admm_env%nao_aux_fit
            ALLOCATE (dbwork)
            CALL dbcsr_create(dbwork, template=matrix_hz(1)%matrix)
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(mhy(ispin, 1)%matrix, admm_env%A, &
                                            admm_env%work_aux_orb, nao)
               CALL cp_gemm('T', 'N', nao, nao, nao_aux, &
                            1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                            admm_env%work_orb_orb)
               CALL dbcsr_copy(dbwork, matrix_hz(ispin)%matrix)
               CALL dbcsr_set(dbwork, 0.0_dp)
               CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbwork, keep_sparsity=.TRUE.)
               CALL dbcsr_add(matrix_hz(ispin)%matrix, dbwork, 1.0_dp, 1.0_dp)
            END DO
            CALL dbcsr_release(dbwork)
            DEALLOCATE (dbwork)
            CALL dbcsr_deallocate_matrix_set(mpz)
         END IF
      END IF

      ! HFX
      hfx_section => section_vals_get_subs_vals(xc_section, "HF")
      CALL section_vals_get(hfx_section, explicit=do_hfx)
      IF (do_hfx) THEN
         CALL section_vals_get(hfx_section, n_repetition=n_rep_hf)
         CPASSERT(n_rep_hf == 1)
         CALL section_vals_val_get(hfx_section, "TREAT_LSD_IN_CORE", l_val=hfx_treat_lsd_in_core, &
                                   i_rep_section=1)
         mspin = 1
         IF (hfx_treat_lsd_in_core) mspin = nspins
         !
         CALL get_qs_env(qs_env=qs_env, rho=rho, x_data=x_data, &
                         s_mstruct_changed=s_mstruct_changed)
         distribute_fock_matrix = .TRUE.
         IF (dft_control%do_admm) THEN
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CALL get_qs_env(qs_env=qs_env, matrix_s_aux_fit=scrm)
            CALL get_qs_env(qs_env=qs_env, rho_aux_fit=rho_aux_fit)
            CALL qs_rho_get(rho_aux_fit, rho_ao_kp=matrix_p)
            NULLIFY (mpz, mhz, mpd, mhd)
            ALLOCATE (mpd(nspins, 1))
            CALL dbcsr_allocate_matrix_set(mpz, nspins, 1)
            CALL dbcsr_allocate_matrix_set(mhz, nspins, 1)
            CALL dbcsr_allocate_matrix_set(mhd, nspins, 1)
            DO ispin = 1, nspins
               ALLOCATE (mhz(ispin, 1)%matrix, mhd(ispin, 1)%matrix)
               CALL dbcsr_create(mhz(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_create(mhd(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_copy(mhz(ispin, 1)%matrix, scrm(1)%matrix)
               CALL dbcsr_copy(mhd(ispin, 1)%matrix, scrm(1)%matrix)
               CALL dbcsr_set(mhz(ispin, 1)%matrix, 0.0_dp)
               CALL dbcsr_set(mhd(ispin, 1)%matrix, 0.0_dp)
               ALLOCATE (mpz(ispin, 1)%matrix)
               CALL dbcsr_create(mpz(ispin, 1)%matrix, template=scrm(1)%matrix)
               CALL dbcsr_copy(mpz(ispin, 1)%matrix, p_env%p1_admm(ispin)%matrix)
               IF (do_ex) THEN
                  CALL dbcsr_add(mpz(ispin, 1)%matrix, ex_env%matrix_pe_admm(ispin)%matrix, &
                                 1.0_dp, 1.0_dp)
               END IF
               mpd(ispin, 1)%matrix => matrix_p(ispin, 1)%matrix
            END DO
            !
            DO ispin = 1, mspin
               eh1 = 0.0
               CALL integrate_four_center(qs_env, x_data, mhz, eh1, mpz, hfx_section, &
                                          para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                          ispin=ispin)
            END DO
            DO ispin = 1, mspin
               eh1 = 0.0
               CALL integrate_four_center(qs_env, x_data, mhd, eh1, mpd, hfx_section, &
                                          para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                          ispin=ispin)
            END DO
            !
            CALL get_qs_env(qs_env, admm_env=admm_env)
            CPASSERT(ASSOCIATED(admm_env%work_aux_orb))
            CPASSERT(ASSOCIATED(admm_env%work_orb_orb))
            nao = admm_env%nao_orb
            nao_aux = admm_env%nao_aux_fit
            ALLOCATE (dbwork)
            CALL dbcsr_create(dbwork, template=matrix_hz(1)%matrix)
            DO ispin = 1, nspins
               CALL cp_dbcsr_sm_fm_multiply(mhz(ispin, 1)%matrix, admm_env%A, &
                                            admm_env%work_aux_orb, nao)
               CALL cp_gemm('T', 'N', nao, nao, nao_aux, &
                            1.0_dp, admm_env%A, admm_env%work_aux_orb, 0.0_dp, &
                            admm_env%work_orb_orb)
               CALL dbcsr_copy(dbwork, matrix_hz(ispin)%matrix)
               CALL dbcsr_set(dbwork, 0.0_dp)
               CALL copy_fm_to_dbcsr(admm_env%work_orb_orb, dbwork, keep_sparsity=.TRUE.)
               CALL dbcsr_add(matrix_hz(ispin)%matrix, dbwork, 1.0_dp, 1.0_dp)
            END DO
            CALL dbcsr_release(dbwork)
            DEALLOCATE (dbwork)
            ! derivatives Tr (Pz [A(T)H dA/dR])
            IF (debug_forces) fodeb(1:3) = force(1)%overlap_admm(1:3, 1)
            IF (ASSOCIATED(mhx) .AND. ASSOCIATED(mhy)) THEN
               DO ispin = 1, nspins
                  CALL dbcsr_add(mhd(ispin, 1)%matrix, mhx(ispin, 1)%matrix, 1.0_dp, 1.0_dp)
                  CALL dbcsr_add(mhz(ispin, 1)%matrix, mhy(ispin, 1)%matrix, 1.0_dp, 1.0_dp)
               END DO
            END IF
            CALL qs_rho_get(rho, rho_ao=matrix_pd)
            CALL admm_projection_derivative(qs_env, admm_env, mhd, mhz, matrix_pd, mpa)
            IF (debug_forces) THEN
               fodeb(1:3) = force(1)%overlap_admm(1:3, 1) - fodeb(1:3)
               CALL mp_sum(fodeb, para_env%group)
               IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*hfx*S' ", fodeb
            END IF
            CALL dbcsr_deallocate_matrix_set(mpz)
            CALL dbcsr_deallocate_matrix_set(mhz)
            CALL dbcsr_deallocate_matrix_set(mhd)
            CALL dbcsr_deallocate_matrix_set(mhx)
            CALL dbcsr_deallocate_matrix_set(mhy)
            DEALLOCATE (mpd)
         ELSE
            ALLOCATE (mpz(nspins, 1), mhz(nspins, 1))
            DO ispin = 1, nspins
               mhz(ispin, 1)%matrix => matrix_hz(ispin)%matrix
               mpz(ispin, 1)%matrix => mpa(ispin)%matrix
            END DO
            DO ispin = 1, mspin
               eh1 = 0.0
               CALL integrate_four_center(qs_env, x_data, mhz, eh1, mpz, hfx_section, &
                                          para_env, s_mstruct_changed, 1, distribute_fock_matrix, &
                                          ispin=ispin)
            END DO
            DEALLOCATE (mhz, mpz)
         END IF

         resp_only = .TRUE.
         IF (debug_forces) fodeb(1:3) = force(1)%fock_4c(1:3, 1)
         IF (dft_control%do_admm) THEN
            CALL qs_rho_get(rho_aux_fit, rho_ao_kp=matrix_p)
            NULLIFY (matrix_pz)
            CALL dbcsr_allocate_matrix_set(matrix_pz, nspins)
            DO ispin = 1, nspins
               ALLOCATE (matrix_pz(ispin)%matrix)
               CALL dbcsr_create(matrix_pz(ispin)%matrix, template=p_env%p1_admm(ispin)%matrix)
               CALL dbcsr_copy(matrix_pz(ispin)%matrix, p_env%p1_admm(ispin)%matrix)
               IF (do_ex) THEN
                  CALL dbcsr_add(matrix_pz(ispin)%matrix, ex_env%matrix_pe_admm(ispin)%matrix, &
                                 1.0_dp, 1.0_dp)
               END IF
            END DO
            CALL derivatives_four_center(qs_env, matrix_p, matrix_pz, hfx_section, para_env, &
                                         1, use_virial, resp_only=resp_only)
            CALL dbcsr_deallocate_matrix_set(matrix_pz)
         ELSE
            CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
            CALL derivatives_four_center(qs_env, matrix_p, mpa, hfx_section, para_env, &
                                         1, use_virial, resp_only=resp_only)
         END IF
         IF (debug_forces) THEN
            fodeb(1:3) = force(1)%fock_4c(1:3, 1) - fodeb(1:3)
            CALL mp_sum(fodeb, para_env%group)
            IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*hfx ", fodeb
         END IF
      END IF

      ! Overlap matrix
      ! H(drho+dz) + Wz
      matrix_wz => p_env%w1
      focc = 1.0_dp
      IF (nspins == 1) focc = 2.0_dp
      CALL get_qs_env(qs_env, mos=mos)
      DO ispin = 1, nspins
         CALL calculate_whz_matrix(mos(ispin)%mo_set%mo_coeff, matrix_hz(ispin)%matrix, &
                                   matrix_wz(ispin)%matrix, focc)
      END DO
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_wz(1)%matrix, matrix_wz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF
      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      NULLIFY (scrm)
      CALL build_overlap_matrix(ks_env, matrix_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_wz(1)%matrix)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%overlap(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wz*dS ", fodeb
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)

      IF (debug_forces) THEN
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1) - ftot1(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Response Force", fodeb
         DEALLOCATE (ftot1, ftot2, ftot3)
      END IF

      IF (do_ex) THEN
         CALL dbcsr_deallocate_matrix_set(mpa)
      END IF

      CALL timestop(handle)

   END SUBROUTINE response_force

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param p_env ...
!> \param matrix_hz ...
!> \param ex_env ...
! **************************************************************************************************
   SUBROUTINE response_force_xtb(qs_env, p_env, matrix_hz, ex_env)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(qs_p_env_type), POINTER                       :: p_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_hz
      TYPE(excited_energy_type), OPTIONAL, POINTER       :: ex_env

      CHARACTER(LEN=*), PARAMETER :: routineN = 'response_force_xtb'

      INTEGER                                            :: atom_a, handle, iatom, ikind, iounit, &
                                                            is, ispin, na, natom, natorb, nimages, &
                                                            nkind, ns, nsgf, nspins
      INTEGER, DIMENSION(25)                             :: lao
      INTEGER, DIMENSION(5)                              :: occ
      LOGICAL                                            :: do_ex, use_virial
      REAL(KIND=dp)                                      :: focc
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:)           :: mcharge, mcharge1
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: aocg, aocg1, charges, charges1, ftot1, &
                                                            ftot2
      REAL(KIND=dp), DIMENSION(3)                        :: fodeb
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pz, matrix_wz, mpa, p_matrix, scrm
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_p, matrix_s
      TYPE(dbcsr_type), POINTER                          :: s_matrix
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(mo_set_p_type), DIMENSION(:), POINTER         :: mos
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_orb
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(xtb_atom_type), POINTER                       :: xtb_kind

      CALL timeset(routineN, handle)

      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_unit_nr(logger)

      do_ex = .FALSE.
      IF (PRESENT(ex_env)) do_ex = .TRUE.

      NULLIFY (ks_env, sab_orb)
      CALL get_qs_env(qs_env=qs_env, ks_env=ks_env, dft_control=dft_control, &
                      sab_orb=sab_orb)
      CALL get_qs_env(qs_env=qs_env, para_env=para_env, force=force)
      nspins = dft_control%nspins

      IF (debug_forces) THEN
         CALL get_qs_env(qs_env, natom=natom, atomic_kind_set=atomic_kind_set)
         ALLOCATE (ftot1(3, natom))
         ALLOCATE (ftot2(3, natom))
         CALL total_qs_force(ftot1, force, atomic_kind_set)
      END IF

      matrix_pz => p_env%p1
      NULLIFY (mpa)
      IF (do_ex) THEN
         CALL dbcsr_allocate_matrix_set(mpa, nspins)
         DO ispin = 1, nspins
            ALLOCATE (mpa(ispin)%matrix)
            CALL dbcsr_create(mpa(ispin)%matrix, template=matrix_pz(ispin)%matrix)
            CALL dbcsr_copy(mpa(ispin)%matrix, matrix_pz(ispin)%matrix)
            CALL dbcsr_add(mpa(ispin)%matrix, ex_env%matrix_pe(ispin)%matrix, 1.0_dp, 1.0_dp)
            CALL dbcsr_set(matrix_hz(ispin)%matrix, 0.0_dp)
         END DO
      ELSE
         mpa => p_env%p1
      END IF
      !
      ! START OF Tr(P+Z)Hcore
      !
      IF (nspins == 2) THEN
         CALL dbcsr_add(mpa(1)%matrix, mpa(2)%matrix, 1.0_dp, 1.0_dp)
      END IF
      ! Hcore  matrix
      IF (debug_forces) fodeb(1:3) = force(1)%all_potential(1:3, 1)
      CALL xtb_hab_force(qs_env, mpa(1)%matrix)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%all_potential(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Pz*dHcore  ", fodeb
      END IF
      IF (nspins == 2) THEN
         CALL dbcsr_add(mpa(1)%matrix, mpa(2)%matrix, 1.0_dp, -1.0_dp)
      END IF
      !
      ! END OF Tr(P+Z)Hcore
      !
      use_virial = .FALSE.
      nimages = 1
      !
      ! Hartree potential of response density
      !
      IF (dft_control%qs_control%xtb_control%coulomb_interaction) THEN
         ! Mulliken charges
         CALL get_qs_env(qs_env, rho=rho, particle_set=particle_set, matrix_s_kp=matrix_s)
         natom = SIZE(particle_set)
         CALL qs_rho_get(rho, rho_ao_kp=matrix_p)
         ALLOCATE (mcharge(natom), charges(natom, 5))
         ALLOCATE (mcharge1(natom), charges1(natom, 5))
         charges = 0.0_dp
         charges1 = 0.0_dp
         CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
         nkind = SIZE(atomic_kind_set)
         CALL get_qs_kind_set(qs_kind_set, maxsgf=nsgf)
         ALLOCATE (aocg(nsgf, natom))
         aocg = 0.0_dp
         ALLOCATE (aocg1(nsgf, natom))
         aocg1 = 0.0_dp
         p_matrix => matrix_p(:, 1)
         s_matrix => matrix_s(1, 1)%matrix
         CALL ao_charges(p_matrix, s_matrix, aocg, para_env)
         CALL ao_charges(mpa, s_matrix, aocg1, para_env)
         DO ikind = 1, nkind
            CALL get_atomic_kind(atomic_kind_set(ikind), natom=na)
            CALL get_qs_kind(qs_kind_set(ikind), xtb_parameter=xtb_kind)
            CALL get_xtb_atom_param(xtb_kind, natorb=natorb, lao=lao, occupation=occ)
            DO iatom = 1, na
               atom_a = atomic_kind_set(ikind)%atom_list(iatom)
               charges(atom_a, :) = REAL(occ(:), KIND=dp)
               DO is = 1, natorb
                  ns = lao(is) + 1
                  charges(atom_a, ns) = charges(atom_a, ns) - aocg(is, atom_a)
                  charges1(atom_a, ns) = charges1(atom_a, ns) - aocg1(is, atom_a)
               END DO
            END DO
         END DO
         DEALLOCATE (aocg, aocg1)
         DO iatom = 1, natom
            mcharge(iatom) = SUM(charges(iatom, :))
            mcharge1(iatom) = SUM(charges1(iatom, :))
         END DO
         ! Coulomb Kernel
         CALL xtb_coulomb_hessian(qs_env, matrix_hz, charges1, mcharge1, mcharge)
         CALL calc_xtb_ehess_force(qs_env, p_matrix, mpa, charges, mcharge, charges1, &
                                   mcharge1, debug_forces)
         !
         DEALLOCATE (charges, mcharge, charges1, mcharge1)
      END IF
      ! Overlap matrix
      ! H(drho+dz) + Wz
      matrix_wz => p_env%w1
      focc = 0.5_dp
      IF (nspins == 1) focc = 1.0_dp
      CALL get_qs_env(qs_env, mos=mos)
      DO ispin = 1, nspins
         CALL calculate_whz_matrix(mos(ispin)%mo_set%mo_coeff, matrix_hz(ispin)%matrix, &
                                   matrix_wz(ispin)%matrix, focc)
      END DO
      IF (nspins == 2) THEN
         CALL dbcsr_add(matrix_wz(1)%matrix, matrix_wz(2)%matrix, &
                        alpha_scalar=1.0_dp, beta_scalar=1.0_dp)
      END IF
      IF (debug_forces) fodeb(1:3) = force(1)%overlap(1:3, 1)
      NULLIFY (scrm)
      CALL build_overlap_matrix(ks_env, matrix_s=scrm, &
                                matrix_name="OVERLAP MATRIX", &
                                basis_type_a="ORB", basis_type_b="ORB", &
                                sab_nl=sab_orb, calculate_forces=.TRUE., &
                                matrix_p=matrix_wz(1)%matrix)
      IF (debug_forces) THEN
         fodeb(1:3) = force(1)%overlap(1:3, 1) - fodeb(1:3)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T33,3F16.8)") "DEBUG:: Wz*dS ", fodeb
      END IF
      CALL dbcsr_deallocate_matrix_set(scrm)

      IF (debug_forces) THEN
         CALL total_qs_force(ftot2, force, atomic_kind_set)
         fodeb(1:3) = ftot2(1:3, 1) - ftot1(1:3, 1)
         CALL mp_sum(fodeb, para_env%group)
         IF (iounit > 0) WRITE (iounit, "(T3,A,T30,3F16.8)") "DEBUG:: Response Force", fodeb
         DEALLOCATE (ftot1, ftot2)
      END IF

      IF (do_ex) THEN
         CALL dbcsr_deallocate_matrix_set(mpa)
      END IF

      CALL timestop(handle)

   END SUBROUTINE response_force_xtb

! **************************************************************************************************
!> \brief ...
!> \param qs_env ...
!> \param admm_env ...
!> \param matrix_hd ...
!> \param matrix_hz ...
!> \param matrix_pd ...
!> \param matrix_pz ...
! **************************************************************************************************
   SUBROUTINE admm_projection_derivative(qs_env, admm_env, matrix_hd, matrix_hz, matrix_pd, matrix_pz)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(dbcsr_p_type), DIMENSION(:, :), POINTER       :: matrix_hd, matrix_hz
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_pd, matrix_pz

      CHARACTER(LEN=*), PARAMETER :: routineN = 'admm_projection_derivative', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ispin, nao, natom, naux, nspins
      REAL(KIND=dp), ALLOCATABLE, DIMENSION(:, :)        :: admm_force
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER          :: matrix_s_aux_fit, matrix_s_aux_fit_vs_orb
      TYPE(dbcsr_type), POINTER                          :: matrix_w_q, matrix_w_s
      TYPE(neighbor_list_set_p_type), DIMENSION(:), &
         POINTER                                         :: sab_aux_fit_asymm, sab_aux_fit_vs_orb
      TYPE(qs_force_type), DIMENSION(:), POINTER         :: force
      TYPE(qs_ks_env_type), POINTER                      :: ks_env

      CALL timeset(routineN, handle)

      CPASSERT(ASSOCIATED(qs_env))
      CPASSERT(ASSOCIATED(admm_env))
      CPASSERT(ASSOCIATED(matrix_hd))
      CPASSERT(ASSOCIATED(matrix_pd))
      CPASSERT(ASSOCIATED(matrix_hz))
      CPASSERT(ASSOCIATED(matrix_pz))

      CALL get_qs_env(qs_env, &
                      ks_env=ks_env, &
                      matrix_s_aux_fit=matrix_s_aux_fit, &
                      matrix_s_aux_fit_vs_orb=matrix_s_aux_fit_vs_orb, &
                      sab_aux_fit_vs_orb=sab_aux_fit_vs_orb, &
                      sab_aux_fit_asymm=sab_aux_fit_asymm)

      ALLOCATE (matrix_w_q)
      CALL dbcsr_copy(matrix_w_q, matrix_s_aux_fit_vs_orb(1)%matrix, &
                      "W MATRIX AUX Q")
      CALL cp_dbcsr_alloc_block_from_nbl(matrix_w_q, sab_aux_fit_vs_orb)
      ALLOCATE (matrix_w_s)
      CALL dbcsr_create(matrix_w_s, template=matrix_s_aux_fit(1)%matrix, &
                        name='W MATRIX AUX S', &
                        matrix_type=dbcsr_type_no_symmetry)
      CALL cp_dbcsr_alloc_block_from_nbl(matrix_w_s, sab_aux_fit_asymm)

      CALL get_qs_env(qs_env=qs_env, atomic_kind_set=atomic_kind_set, &
                      natom=natom, force=force)
      ALLOCATE (admm_force(3, natom))

      nspins = SIZE(matrix_pz)
      nao = admm_env%nao_orb
      naux = admm_env%nao_aux_fit
      DO ispin = 1, nspins
         CALL copy_dbcsr_to_fm(matrix_hd(ispin, 1)%matrix, admm_env%work_aux_aux)
         CALL cp_gemm("N", "N", naux, naux, naux, 1.0_dp, admm_env%s_inv, &
                      admm_env%work_aux_aux, 0.0_dp, admm_env%work_aux_aux2)
         CALL cp_gemm("N", "N", naux, nao, naux, 1.0_dp, admm_env%work_aux_aux2, &
                      admm_env%A, 0.0_dp, admm_env%work_aux_orb)
         CALL copy_dbcsr_to_fm(matrix_pz(ispin)%matrix, admm_env%work_orb_orb)
         ! admm_env%work_aux_orb2 = S-1*H(D)*A*PZ
         CALL cp_gemm("N", "N", naux, nao, nao, 1.0_dp, admm_env%work_aux_orb, &
                      admm_env%work_orb_orb, 0.0_dp, admm_env%work_aux_orb2)

         CALL copy_dbcsr_to_fm(matrix_hz(ispin, 1)%matrix, admm_env%work_aux_aux)
         CALL cp_gemm("N", "N", naux, naux, naux, 1.0_dp, admm_env%s_inv, &
                      admm_env%work_aux_aux, 0.0_dp, admm_env%work_aux_aux2)
         CALL cp_gemm("N", "N", naux, nao, naux, 1.0_dp, admm_env%work_aux_aux2, &
                      admm_env%A, 0.0_dp, admm_env%work_aux_orb)
         CALL copy_dbcsr_to_fm(matrix_pd(ispin)%matrix, admm_env%work_orb_orb)
         ! admm_env%work_aux_orb2 = S-1*H(Z)*A*PD
         CALL cp_gemm("N", "N", naux, nao, nao, 1.0_dp, admm_env%work_aux_orb, &
                      admm_env%work_orb_orb, 1.0_dp, admm_env%work_aux_orb2)

         ! admm_env%work_aux_aux = [S-1*H(D)*A*PZ+S-1*H(Z)*A*PD]*A(T)
         CALL cp_gemm("N", "T", naux, naux, nao, 1.0_dp, admm_env%work_aux_orb2, &
                      admm_env%A, 0.0_dp, admm_env%work_aux_aux)
         !
         CALL copy_fm_to_dbcsr(admm_env%work_aux_orb2, matrix_w_q, keep_sparsity=.TRUE.)
         CALL copy_fm_to_dbcsr(admm_env%work_aux_aux, matrix_w_s, keep_sparsity=.TRUE.)
! focc = 2?
         CALL dbcsr_scale(matrix_w_q, -2.0_dp)
         CALL dbcsr_scale(matrix_w_s, 2.0_dp)
         !
         admm_force = 0.0_dp
         CALL build_overlap_force(ks_env, admm_force, &
                                  basis_type_a="AUX_FIT", basis_type_b="AUX_FIT", &
                                  sab_nl=sab_aux_fit_asymm, matrix_p=matrix_w_s)
         CALL build_overlap_force(ks_env, admm_force, &
                                  basis_type_a="AUX_FIT", basis_type_b="ORB", &
                                  sab_nl=sab_aux_fit_vs_orb, matrix_p=matrix_w_q)
         ! add forces
         CALL add_qs_force(admm_force, force, "overlap_admm", atomic_kind_set)
      END DO

      DEALLOCATE (admm_force)
      CALL dbcsr_deallocate_matrix(matrix_w_s)
      CALL dbcsr_deallocate_matrix(matrix_w_q)

      CALL timestop(handle)

   END SUBROUTINE admm_projection_derivative

! **************************************************************************************************
!> \brief calculate the Kohn-Sham reference potential
!> \param qs_env ...
!> \param vh_rspace ...
!> \param vxc_rspace ...
!> \param vtau_rspace ...
!> \param ehartree ...
!> \param exc ...
!> \par History
!>      10.2019 created [JGH]
!> \author JGH
! **************************************************************************************************
   SUBROUTINE ks_ref_potential(qs_env, vh_rspace, vxc_rspace, vtau_rspace, ehartree, exc)
      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(pw_p_type), POINTER                           :: vh_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: vxc_rspace, vtau_rspace
      REAL(KIND=dp), INTENT(OUT)                         :: ehartree, exc

      CHARACTER(LEN=*), PARAMETER                        :: routineN = 'ks_ref_potential'

      INTEGER                                            :: handle, iab, ispin, nspins
      REAL(dp)                                           :: eovrl, eself
      TYPE(admm_type), POINTER                           :: admm_env
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(pw_env_type), POINTER                         :: pw_env
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      TYPE(pw_p_type)                                    :: rho_tot_gspace, v_hartree_gspace, &
                                                            v_hartree_rspace
      TYPE(pw_p_type), DIMENSION(:), POINTER             :: v_rspace, v_tau_rspace
      TYPE(pw_poisson_type), POINTER                     :: poisson_env
      TYPE(pw_pool_type), POINTER                        :: auxbas_pw_pool
      TYPE(qs_ks_env_type), POINTER                      :: ks_env
      TYPE(qs_rho_type), POINTER                         :: rho
      TYPE(section_vals_type), POINTER                   :: xc_section

      CALL timeset(routineN, handle)

      ! get all information on the electronic density
      NULLIFY (rho, ks_env)
      CALL get_qs_env(qs_env=qs_env, rho=rho, dft_control=dft_control, &
                      para_env=para_env, blacs_env=blacs_env, ks_env=ks_env)

      nspins = dft_control%nspins

      NULLIFY (pw_env)
      CALL get_qs_env(qs_env=qs_env, pw_env=pw_env)
      CPASSERT(ASSOCIATED(pw_env))

      NULLIFY (auxbas_pw_pool, poisson_env)
      ! gets the tmp grids
      CALL pw_env_get(pw_env, auxbas_pw_pool=auxbas_pw_pool, &
                      poisson_env=poisson_env)

      ! Calculate the Hartree potential
      NULLIFY (v_hartree_gspace%pw, rho_tot_gspace%pw, v_hartree_rspace%pw)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, v_hartree_rspace%pw, &
                             use_data=REALDATA3D, in_space=REALSPACE)
      CALL pw_pool_create_pw(auxbas_pw_pool, rho_tot_gspace%pw, &
                             use_data=COMPLEXDATA1D, in_space=RECIPROCALSPACE)

      ! Get the total density in g-space [ions + electrons]
      CALL calc_rho_tot_gspace(rho_tot_gspace, qs_env, rho)

      CALL pw_poisson_solve(poisson_env, rho_tot_gspace%pw, ehartree, &
                            v_hartree_gspace%pw)
      CALL pw_transfer(v_hartree_gspace%pw, v_hartree_rspace%pw)
      CALL pw_scale(v_hartree_rspace%pw, v_hartree_rspace%pw%pw_grid%dvol)

      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_gspace%pw)
      CALL pw_pool_give_back_pw(auxbas_pw_pool, rho_tot_gspace%pw)
      !
      CALL calculate_ecore_self(qs_env, E_self_core=eself)
      CALL calculate_ecore_overlap(qs_env, para_env, .FALSE., E_overlap_core=eovrl)
      ehartree = ehartree + eovrl + eself

      ! v_rspace and v_tau_rspace are generated from the auxbas pool
      IF (dft_control%do_admm) THEN
         CALL get_qs_env(qs_env, admm_env=admm_env)
         xc_section => admm_env%xc_section_primary
      ELSE
         xc_section => section_vals_get_subs_vals(qs_env%input, "DFT%XC")
      END IF
      NULLIFY (v_rspace, v_tau_rspace)
      CALL qs_vxc_create(ks_env=ks_env, rho_struct=rho, xc_section=xc_section, &
                         vxc_rho=v_rspace, vxc_tau=v_tau_rspace, exc=exc, just_energy=.FALSE.)

      ! allocate potentials
      IF (ASSOCIATED(vh_rspace)) THEN
         CALL pw_release(vh_rspace%pw)
      ELSE
         ALLOCATE (vh_rspace)
         NULLIFY (vh_rspace%pw)
      END IF
      IF (ASSOCIATED(vxc_rspace)) THEN
         DO iab = 1, SIZE(vxc_rspace)
            CALL pw_release(vxc_rspace(iab)%pw)
         END DO
      ELSE
         ALLOCATE (vxc_rspace(nspins))
         DO iab = 1, nspins
            NULLIFY (vxc_rspace(iab)%pw)
         END DO
      END IF
      IF (ASSOCIATED(v_tau_rspace)) THEN
         IF (ASSOCIATED(vtau_rspace)) THEN
            DO iab = 1, SIZE(vtau_rspace)
               CALL pw_release(vtau_rspace(iab)%pw)
            END DO
         ELSE
            ALLOCATE (vtau_rspace(nspins))
            DO iab = 1, nspins
               NULLIFY (vtau_rspace(iab)%pw)
            END DO
         END IF
      ELSE
         NULLIFY (vtau_rspace)
      END IF

      pw_grid => v_hartree_rspace%pw%pw_grid
      CALL pw_create(vh_rspace%pw, pw_grid, use_data=REALDATA3D, in_space=REALSPACE)
      DO ispin = 1, nspins
         NULLIFY (vxc_rspace(ispin)%pw)
         CALL pw_create(vxc_rspace(ispin)%pw, pw_grid, &
                        use_data=REALDATA3D, in_space=REALSPACE)
         IF (ASSOCIATED(vtau_rspace)) THEN
            NULLIFY (vtau_rspace(ispin)%pw)
            CALL pw_create(vtau_rspace(ispin)%pw, pw_grid, &
                           use_data=REALDATA3D, in_space=REALSPACE)
         END IF
      END DO
      !
      CALL pw_transfer(v_hartree_rspace%pw, vh_rspace%pw)
      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_transfer(v_rspace(ispin)%pw, vxc_rspace(ispin)%pw)
            CALL pw_scale(vxc_rspace(ispin)%pw, v_rspace(ispin)%pw%pw_grid%dvol)
            IF (ASSOCIATED(v_tau_rspace)) THEN
               CALL pw_transfer(v_tau_rspace(ispin)%pw, vtau_rspace(ispin)%pw)
               CALL pw_scale(vtau_rspace(ispin)%pw, v_tau_rspace(ispin)%pw%pw_grid%dvol)
            END IF
         END DO
      ELSE
         DO ispin = 1, nspins
            CALL pw_zero(vxc_rspace(ispin)%pw)
         END DO
      END IF

      ! return pw grids
      CALL pw_pool_give_back_pw(auxbas_pw_pool, v_hartree_rspace%pw)
      IF (ASSOCIATED(v_rspace)) THEN
         DO ispin = 1, nspins
            CALL pw_pool_give_back_pw(auxbas_pw_pool, v_rspace(ispin)%pw)
            IF (ASSOCIATED(v_tau_rspace)) THEN
               CALL pw_pool_give_back_pw(auxbas_pw_pool, v_tau_rspace(ispin)%pw)
            END IF
         ENDDO
         DEALLOCATE (v_rspace)
      END IF
      IF (ASSOCIATED(v_tau_rspace)) DEALLOCATE (v_tau_rspace)

      CALL timestop(handle)

   END SUBROUTINE ks_ref_potential

! **************************************************************************************************

END MODULE response_solver
