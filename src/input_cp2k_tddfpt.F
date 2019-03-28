!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2019  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief function that build the dft section of the input
!> \par History
!>      10.2005 moved out of input_cp2k [fawzi]
!> \author fawzi
! **************************************************************************************************
MODULE input_cp2k_tddfpt
   USE bibliography,                    ONLY: Iannuzzi2005
   USE input_constants,                 ONLY: &
        oe_gllb, oe_lb, oe_none, oe_saop, oe_sic, tddfpt_davidson, tddfpt_kernel_full, &
        tddfpt_kernel_stda, tddfpt_lanczos, tddfpt_singlet, tddfpt_triplet
   USE input_cp2k_xc,                   ONLY: create_sic_section,&
                                              create_xc_section
   USE input_keyword_types,             ONLY: keyword_create,&
                                              keyword_release,&
                                              keyword_type
   USE input_section_types,             ONLY: section_add_keyword,&
                                              section_add_subsection,&
                                              section_create,&
                                              section_release,&
                                              section_type
   USE input_val_types,                 ONLY: integer_t,&
                                              real_t
   USE kinds,                           ONLY: dp
   USE string_utilities,                ONLY: s2a
#include "./base/base_uses.f90"

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'input_cp2k_tddfpt'

   PUBLIC :: create_tddfpt_section

CONTAINS

! **************************************************************************************************
!> \brief creates the input section for the tddfpt part
!> \param section the section to create
!> \author teo
! **************************************************************************************************
   SUBROUTINE create_tddfpt_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_tddfpt_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword
      TYPE(section_type), POINTER                        :: subsection

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "tddfpt", &
                          description="parameters needed to set up the Time Dependent Density Functional PT", &
                          n_keywords=14, n_subsections=1, repeats=.FALSE., &
                          citations=(/Iannuzzi2005/))

      NULLIFY (subsection, keyword)

      ! Integer
      CALL keyword_create(keyword, name="MAX_KV", &
                          variants=s2a("MAX_VECTORS"), &
                          description=" maximal number of Krylov space vectors", &
                          usage="MAX_KV someInteger>0", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=60)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="RESTARTS", &
                          variants=s2a("N_RESTARTS"), &
                          description=" maximal number subspace search restarts", &
                          usage="RESTARTS someInteger>0", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=5)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="NEV", &
                          variants=s2a("N_EV", "EV"), &
                          description=" number of excitations to calculate", &
                          usage="NEV someInteger>0", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=1)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="NLUMO", &
                          description=" number of additional unoccupied orbitals ", &
                          usage="NLUMO 10", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=5)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="NREORTHO", &
                          variants=s2a("N_REORTHO", "REORTHO", "REORTHOGONALITAZIONS"), &
                          description=" number of reorthogonalization steps", &
                          usage="NREORTHO someInteger>0", &
                          n_var=1, type_of_var=integer_t, &
                          default_i_val=2)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Logical
      CALL keyword_create(keyword, name="CALCULATE_KERNEL", &
                          variants=s2a("DO_KERNEL"), &
                          description="compute the kernel (debug purpose only)", &
                          usage="CALCULATE_KERNEL logical_value", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="LSD_SINGLETS", &
                          description="compute singlets using lsd vxc kernel", &
                          usage="LSD_SINGLETS logical_value", &
                          default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="INVERT_S", &
                          variants=s2a("INVERT_OVERLAP"), &
                          description="use the inverse of the overlap matrix", &
                          usage="INVERT_S logical_value", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="PRECONDITIONER", &
                          variants=s2a("PRECOND"), &
                          description="use the preconditioner (only for Davidson)", &
                          usage="PRECONDITIONER logical_value", &
                          default_l_val=.TRUE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Character
      CALL keyword_create(keyword, name="KERNEL", &
                          description="options to compute the kernel", &
                          usage="KERNEL FULL", &
                          enum_c_vals=s2a("FULL", "sTDA"), &
                          enum_i_vals=(/tddfpt_kernel_full, tddfpt_kernel_stda/), &
                          default_i_val=tddfpt_kernel_full)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="RES_ETYPE", &
                          variants=s2a("RESTRICTED_EXCITATIONS_TYPE", "RES_E_TYPE"), &
                          description="(singlets/triplets) for restricted calculation", &
                          usage="RES_ETYPE T", &
                          enum_c_vals=s2a("S", "SINGLET", "SINGLETS", "T", "TRIPLET", "TRIPLETS"), &
                          enum_i_vals=(/tddfpt_singlet, tddfpt_singlet, tddfpt_singlet, &
                                        tddfpt_triplet, tddfpt_triplet, tddfpt_triplet/), &
                          default_i_val=tddfpt_singlet)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="DIAG_METHOD", &
                          variants=s2a("DIAGONALIZATION_METHOD", "METHOD"), &
                          description="Diagonalization method used in tddfpt", &
                          usage="DIAG_METHOD DAVIDSON", &
                          enum_c_vals=s2a("DAVIDSON", "LANCZOS"), &
                          enum_i_vals=(/tddfpt_davidson, tddfpt_lanczos/), &
                          default_i_val=tddfpt_davidson)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="OE_CORR", &
                          variants=s2a("ORBITAL_EIGENVALUES_CORRECTION"), &
                          description="Which type of orbital eigenvalue correction to use\n"// &
                          "(to yield better HOMO-LUMO energies)", &
                          usage="OE_CORR SAOP", &
                          enum_c_vals=s2a("NONE", "LB", "LB_ALPHA", "LB94", "GLLB", "GLB", "SAOP", "SIC"), &
                          enum_i_vals=(/oe_none, oe_lb, oe_lb, oe_lb, oe_gllb, oe_gllb, oe_saop, oe_sic/), &
                          default_i_val=oe_none)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      ! Real
      CALL keyword_create(keyword, name="CONVERGENCE", &
                          variants=s2a("CONV"), &
                          description="The convergence of the eigenvalues", &
                          usage="CONVERGENCE 1.0E-6 ", &
                          n_var=1, type_of_var=real_t, &
                          default_r_val=1.0e-5_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL create_xc_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_sic_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

      CALL create_stda_section(subsection)
      CALL section_add_subsection(section, subsection)
      CALL section_release(subsection)

   END SUBROUTINE create_tddfpt_section
! **************************************************************************************************
!> \brief creates the stda input section (simplified Tamm Dancoff Approximation)
!> \param section the section to create
! **************************************************************************************************
   SUBROUTINE create_stda_section(section)
      TYPE(section_type), POINTER                        :: section

      CHARACTER(len=*), PARAMETER :: routineN = 'create_stda_section', &
         routineP = moduleN//':'//routineN

      TYPE(keyword_type), POINTER                        :: keyword

      CPASSERT(.NOT. ASSOCIATED(section))
      CALL section_create(section, "sTDA", &
                          description="parameters needed and setup for sTDA calculations", &
                          n_keywords=4, n_subsections=0, repeats=.FALSE.)
      NULLIFY (keyword)

      CALL keyword_create(keyword, name="FRACTION", &
                          variants=(/"HFX_FRACTION"/), &
                          description="The fraction of Hartree-Fock to add to the total energy. "// &
                          "1.0 implies standard Hartree-Fock if used with XC_FUNCTIONAL NONE. "// &
                          "NOTE: In a mixed potential calculation this should be set to 1.0, otherwise "// &
                          "all parts are multiplied with this factor. ", &
                          usage="FRACTION 0.0", default_r_val=0.0_dp)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

      CALL keyword_create(keyword, name="DO_EWALD", &
                          description="Use Ewald type method for Coulomb interaction", &
                          usage="DO_EWALD", default_l_val=.FALSE., lone_keyword_l_val=.TRUE.)
      CALL section_add_keyword(section, keyword)
      CALL keyword_release(keyword)

   END SUBROUTINE create_stda_section

END MODULE input_cp2k_tddfpt
