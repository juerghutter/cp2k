!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright (C) 2000 - 2018  CP2K developers group                                               !
!--------------------------------------------------------------------------------------------------!
! **************************************************************************************************
!> \brief Simplified Tamm Dancoff approach (sTDA).
! **************************************************************************************************
MODULE qs_tddfpt_stda_types

   USE atomic_kind_types,               ONLY: atomic_kind_type,&
                                              get_atomic_kind
   USE basis_set_types,                 ONLY: get_gto_basis_set,&
                                              gto_basis_set_type
   USE cell_types,                      ONLY: cell_type,&
                                              plane_distance
   USE cp_blacs_env,                    ONLY: cp_blacs_env_type
   USE cp_control_types,                ONLY: dft_control_type,&
                                              stda_control_type,&
                                              tddfpt_control_type
   USE cp_dbcsr_operations,             ONLY: dbcsr_deallocate_matrix_set
   USE cp_fm_pool_types,                ONLY: cp_fm_pool_p_type,&
                                              fm_pools_create_fm_vect
   USE cp_fm_struct,                    ONLY: cp_fm_struct_create,&
                                              cp_fm_struct_release,&
                                              cp_fm_struct_type
   USE cp_fm_types,                     ONLY: cp_fm_create,&
                                              cp_fm_get_info,&
                                              cp_fm_p_type,&
                                              cp_fm_release,&
                                              cp_fm_type
   USE cp_fm_vect,                      ONLY: cp_fm_vect_dealloc
   USE cp_log_handling,                 ONLY: cp_get_default_logger,&
                                              cp_logger_get_default_io_unit,&
                                              cp_logger_type
   USE cp_para_types,                   ONLY: cp_para_env_type
   USE dbcsr_api,                       ONLY: dbcsr_p_type
   USE input_section_types,             ONLY: section_vals_get_subs_vals,&
                                              section_vals_type,&
                                              section_vals_val_get
   USE kinds,                           ONLY: dp
   USE particle_types,                  ONLY: particle_type
   USE physcon,                         ONLY: evolt
   USE qs_environment_types,            ONLY: get_qs_env,&
                                              qs_environment_type
   USE qs_kind_types,                   ONLY: get_qs_kind,&
                                              get_qs_kind_set,&
                                              qs_kind_type
   USE qs_matrix_pools,                 ONLY: mpools_get
   USE qs_mo_types,                     ONLY: get_mo_set
#include "./base/base_uses.f90"

   IMPLICIT NONE

   PRIVATE

   TYPE stda_kind_type
      CHARACTER(LEN=2)                            :: symbol !element symbol
      INTEGER                                     :: z !atomic_number
      INTEGER                                     :: kind_number !kind number
      REAL(KIND=dp)                               :: hardness_param !hardness_parameter eta
      REAL(KIND=dp)                               :: rcut !cutoff radius for short-range Coulomb
   END TYPE

   TYPE stda_atom_type
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER    :: charges => Null() ! Lowdin transition density charges   !NOCC*NOCC
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER    :: xtransformed_double => Null() !NAO*NOCC
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER    :: ctransformed_double => Null() !NAO*NOCC
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER    :: gmcharge => Null()
   END TYPE

   TYPE stda_atom_p_type
      TYPE(stda_atom_type), POINTER               :: atom_param
   END TYPE

   TYPE stda_kind_p_type
      TYPE(stda_kind_type), POINTER               :: kind_param
   END TYPE

   TYPE stda_env_type
      REAL(KIND=dp)                                :: hfx_fraction !fraction of non-local HF exchange !! REDUNDANT also in stda_control
      REAL(KIND=dp)                                :: alpha_param ! empirical parameters
      REAL(KIND=dp)                                :: beta_param ! empirical parameters
      TYPE(stda_atom_p_type), DIMENSION(:), POINTER:: atom_param_set ! for each atom
      TYPE(stda_kind_p_type), DIMENSION(:), POINTER:: kind_param_set ! for each kind of atoms
      TYPE(dbcsr_p_type), DIMENSION(:), POINTER    :: gamma_coulomb, & !approximated Coulomb-type integrals (Mataga-Nishimoto-Ohno-Klopman)
                                                      gamma_exchange !approximated exchange-type integrals (Mataga-Nishimoto-Ohno-Klopman)
      ! BOTH NATOM*NATOM
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER    :: ctransformed => Null() !Lowdin MO coefficients !NAO*NOCC
      TYPE(cp_fm_p_type), DIMENSION(:), POINTER    :: xtransformed => Null() !Lowdin transformed Davidson trial vector X  !NAO*NOCC
      TYPE(cp_fm_type), POINTER                    :: shalf => Null() !Lowdin MO coefficients !for debugging !NAO*NAO
      INTEGER, DIMENSION(2)                         :: n_ao !number of atomic orbitals
      INTEGER, DIMENSION(2)                         :: n_occ !number of occupied orbitals
      !  INTEGER,DIMENSION(2)                         :: n_unocc !number of unoccupied orbitals
   END TYPE

!PARAMETERS
   INTEGER, PARAMETER, PRIVATE :: nelem = 103
   !   H                                                                      He
   !   Li Be                                                 B  C  N  O  F    Ne
   !   Na Mg                                                 Al Si P  S  Cl   Ar
   !   K  Ca Sc                Ti V  Cr Mn Fe Co Ni Cu Zn    Ga Ge As Se Br   Kr
   !   Rb Sr Y                 Zr Nb Mo Tc Ru Rh Pd Ag Cd    In Sn Sb Te I    Xe
   !   Cs Ba La Ce-Lu          Hf Ta W  Re Os Ir Pt Au Hg    Tl Pb Bi Po At   Rn
   !   Fr Ra Ac Th Pa U        Np Pu Am Cm Bk Cf Es Fm Md    No Lr  103

   REAL(KIND=dp), DIMENSION(1:nelem), PARAMETER, PRIVATE:: hardness = &
                                                           (/6.4299544220_dp, 12.544911890_dp, & ! 2 H-He
                                                             2.3745866560_dp, 3.4967633530_dp, 4.6190089720_dp, 5.7409789220_dp, &
                                                             6.8624665290_dp, 7.9854357010_dp, 9.1064753720_dp, 10.23034050_dp, & ! 8 Li-Ne
                                                             2.444141360_dp, 3.0146513830_dp, 3.5849070740_dp, 4.15513090_dp, &
                                                             4.7258039740_dp, 5.2959792410_dp, 5.8661864840_dp, 6.4366187140_dp, & ! 8 Na-Ar
                                                             2.3273178360_dp, 2.7587238140_dp, 2.8581921140_dp, 2.9578300430_dp, &
                                                             3.0573410060_dp, 3.1567254290_dp, 3.2563827230_dp, 3.3559314050_dp, &
                                                             3.4556091170_dp, 3.5550133130_dp, 3.6544183480_dp, 3.7541601450_dp, &
                                                             4.1855197930_dp, 4.6166272460_dp, 5.0662145070_dp, 5.4794960970_dp, &
                                                             5.9110996450_dp, 6.3418467680_dp, & ! 18 K-Kr
                                                             2.1204582570_dp, 2.5373700480_dp, 2.6335468980_dp, 2.7297528930_dp, &
                                                             2.8259738860_dp, 2.9221296040_dp, 3.0183708780_dp, 3.1145981770_dp, &
                                                             3.210756280_dp, 3.3069474480_dp, 3.4031948570_dp, 3.4993761390_dp, &
                                                             3.9163692460_dp, 4.3332332190_dp, 4.7500787860_dp, 5.1669793270_dp, &
                                                             5.5838871020_dp, 6.000897330_dp, & ! 18 Rb-Xe
                                                             0.6829150240_dp, 0.9200946840_dp, 1.1570887860_dp, 1.39427570_dp, &
                                                             1.6314731730_dp, 1.8684389980_dp, 2.1056577930_dp, 2.3426646420_dp, &
                                                             2.5798149820_dp, 2.8170264230_dp, 3.0540365330_dp, 3.2911692310_dp, &
                                                             3.5282971610_dp, 3.7655249290_dp, 4.0025547030_dp, 4.2394783410_dp, &
                                                             4.4765830210_dp, 4.7065224490_dp, 4.9508466940_dp, 5.1879311720_dp, &
                                                             5.4256076210_dp, 5.6619144310_dp, 5.900042920_dp, 6.1367145320_dp, &
                                                             6.3741299770_dp, 6.6102656130_dp, 1.7043485810_dp, 1.9413526120_dp, &
                                                             2.178491510_dp, 2.4158121060_dp, 2.6527780840_dp, 2.8899554570_dp, & ! 32 Cs-Rn
                                                             0.9882529880_dp, 1.2819499970_dp, 1.3497250380_dp, 1.4175257380_dp, &
                                                             1.9368567520_dp, 2.2305576050_dp, 2.5241204960_dp, 3.0436128480_dp, &
                                                             3.4168675260_dp, 3.4049844440_dp, 3.9244199680_dp, 4.2180813280_dp, &
                                                             4.5115926320_dp, 4.8050928950_dp, 5.0989816210_dp, 5.3926054620_dp, &
                                                             5.4606987930_dp/) ! 17 Fr-Lr

   ! Shell occupation
   INTEGER, DIMENSION(1:5, 0:nelem), PARAMETER, PRIVATE :: occupation = RESHAPE((/0, 0, 0, 0, 0, & ! 0
                                                                                  1, 0, 0, 0, 0, 2, 0, 0, 0, 0, & ! 2
      1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 0, 0, 0, 2, 2, 0, 0, 0, 2, 3, 0, 0, 0, 2, 4, 0, 0, 0, 2, 5, 0, 0, 0, 2, 6, 0, 0, 0, & ! 10
      1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 1, 0, 0, 0, 2, 2, 0, 0, 0, 2, 3, 0, 0, 0, 2, 4, 0, 0, 0, 2, 5, 0, 0, 0, 2, 6, 0, 0, 0, & ! 18
           1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 1, 0, 0, 2, 0, 2, 0, 0, 2, 0, 3, 0, 0, 2, 0, 4, 0, 0, 2, 0, 5, 0, 0, 2, 0, 6, 0, 0, &
      2,0,7,0,0,  2,0,8,0,0,  2,0,9,0,0,  2,0,10,0,0, 2,1,0,0,0,  2,2,0,0,0,  2,3,0,0,0,  2,4,0,0,0,  2,5,0,0,0, 2,6,0,0,0, & ! 36
         1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 1, 0, 0, 2, 0, 2, 0, 0, 2, 0, 3, 0, 0, 2, 0, 4, 0, 0, 2, 0, 5, 0, 0, 2, 0, 6, 0, 0, & !
      2,0,7,0,0,  2,0,8,0,0,  2,0,9,0,0,  2,0,10,0,0, 2,1,0,0,0,  2,2,0,0,0,  2,3,0,0,0,  2,4,0,0,0,  2,5,0,0,0, 2,6,0,0,0, & ! 54
                                                                                  1, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 1, 0, 0, &
           2, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 1, 0, 0, &
                           2, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 1, 0, 0, 2, 0, 1, 0, 0, & ! Lanthanides
           2, 0, 2, 0, 0, 2, 0, 3, 0, 0, 2, 0, 4, 0, 0, 2, 0, 5, 0, 0, 2, 0, 6, 0, 0, 2, 0, 7, 0, 0, 2, 0, 8, 0, 0, 2, 0, 9, 0, 0, &
2, 0, 10, 0, 0, 2, 1, 0, 0, 0, 2, 2, 0, 0, 0, 2, 3, 0, 0, 0, 2, 4, 0, 0, 0, 2, 5, 0, 0, 0, 2, 6, 0, 0, 0, & ! 86 (last element defined)
                                                                                  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & !
           0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
                             0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, & ! Actinides
                                                                      0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/), (/5, nelem+1/))

   REAL(KIND=dp), DIMENSION(2), PARAMETER, PRIVATE:: alpha = (/1.420_dp, 0.480_dp/)
   REAL(KIND=dp), DIMENSION(2), PARAMETER, PRIVATE:: beta = (/0.200_dp, 1.830_dp/)

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'qs_tddfpt_stda_types'

   PUBLIC:: stda_env_type, allocate_stda_atom_param, deallocate_stda_atom_param, &
            allocate_stda_env, deallocate_stda_env, stda_init_param

CONTAINS
! **************************************************************************************************
!> \brief Get the parameters needed for an sTDA calculation
!> \param qs_env ...
!> \param stda_env ...
! **************************************************************************************************
   SUBROUTINE stda_init_param(qs_env, stda_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type), POINTER                       :: stda_env

      CHARACTER(len=*), PARAMETER :: routineN = 'stda_init_param', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ikind, iounit, nkind
      LOGICAL                                            :: explicit
      REAL(KIND=dp)                                      :: rc
      TYPE(atomic_kind_type), DIMENSION(:), POINTER      :: atomic_kind_set
      TYPE(atomic_kind_type), POINTER                    :: atomic_kind
      TYPE(cell_type), POINTER                           :: cell
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(gto_basis_set_type), POINTER                  :: tmp_basis_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(qs_kind_type), POINTER                        :: qs_kind
      TYPE(section_vals_type), POINTER                   :: input, stda_section
      TYPE(stda_control_type), POINTER                   :: stda_control
      TYPE(stda_kind_type), POINTER                      :: kind_param

      CALL timeset(routineN, handle)

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      CPASSERT(ASSOCIATED(stda_env))
      CPASSERT(ASSOCIATED(stda_env%kind_param_set))

      NULLIFY (atomic_kind_set, qs_kind_set) !get element symbol and atomic number
      CALL get_qs_env(qs_env, atomic_kind_set=atomic_kind_set, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)

      DO ikind = 1, nkind
         atomic_kind => atomic_kind_set(ikind)
         kind_param => stda_env%kind_param_set(ikind)%kind_param
         CALL get_atomic_kind(atomic_kind, &
                              element_symbol=stda_env%kind_param_set(ikind)%kind_param%symbol, & !element symbol
                              kind_number=stda_env%kind_param_set(ikind)%kind_param%kind_number, & !kind_number
                              z=stda_env%kind_param_set(ikind)%kind_param%z) !atomic number
         NULLIFY (kind_param)
      ENDDO !ikind

! HFx_fraction
      NULLIFY (input)
      CALL get_qs_env(qs_env, input=input)
      stda_section => section_vals_get_subs_vals(input, "DFT%TDDFPT%sTDA")
      CALL section_vals_val_get(stda_section, "FRACTION", r_val=stda_env%hfx_fraction, explicit=explicit)
!      IF (explicit) THEN
!         xc_section => section_vals_get_subs_vals(input, "DFT%XC")
!         CALL section_vals_val_get(input, "DFT%EXCITATIONS", &
!                                   i_val=excitations)
!         IF (excitations == tddfpt_excitations) THEN
!            xc_section => section_vals_get_subs_vals(input, "DFT%TDDFPT%XC")
!            CALL section_vals_get(xc_section, explicit=explicit)
!            IF (.NOT. explicit) THEN
!               xc_section => section_vals_get_subs_vals(input, "DFT%XC")
!            END IF
!         END IF
!         hfx_section => section_vals_get_subs_vals(xc_section, "HF")
!         CALL section_vals_val_get(hfx_section, "FRACTION", r_val=stda_env%hfx_fraction)
      IF (iounit > 0) THEN
         WRITE (iounit, *)
         WRITE (iounit, '(A,T71,F10.3)') " STDA | Fraction a_x [a.u.]", &
            stda_env%hfx_fraction
      END IF
!      ENDIF

! alpha and beta parameter
      stda_env%alpha_param = alpha(1)+stda_env%hfx_fraction*alpha(2)
      stda_env%beta_param = beta(1)+stda_env%hfx_fraction*beta(2)

      IF (iounit > 0) THEN
         WRITE (iounit, '(A,T71,F10.3)') " STDA | Alpha parameter [a.u.]", &
            stda_env%alpha_param
         WRITE (iounit, '(A,T71,F10.3)') " STDA | Beta parameter [a.u.]", &
            stda_env%beta_param
      END IF

      DO ikind = 1, nkind
! hardness parameter
         stda_env%kind_param_set(ikind)%kind_param%hardness_param = &
            hardness(stda_env%kind_param_set(ikind)%kind_param%z)*2.0_dp/evolt
! rcut parameter
         qs_kind => qs_kind_set(ikind)
         CALL get_qs_kind(qs_kind, basis_set=tmp_basis_set)
         CALL get_gto_basis_set(tmp_basis_set, kind_radius=stda_env%kind_param_set(ikind)%kind_param%rcut)
      ENDDO ! nkind

! set coulomb cutoff
      CALL section_vals_val_get(stda_section, "COULOMB_CUTOFF", explicit=explicit)
      CALL get_qs_env(qs_env, cell=cell, dft_control=dft_control)
      stda_control => dft_control%tddfpt_control%stda_control
      IF ((.NOT. explicit) .AND. stda_control%do_mic) THEN
         rc = MAX(0.5*plane_distance(1, 0, 0, cell), &
                  0.5*plane_distance(0, 1, 0, cell), &
                  0.5*plane_distance(0, 0, 1, cell))
         IF (rc /= stda_control%cutoff_cou) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, *)
               WRITE (iounit, '(A,T43,A)') " STDA |", &
                  " Coulomb Integral cutoff was redefined"
               WRITE (iounit, '(A,T71,F10.3)') " STDA | Old value [a.u.]", &
                  stda_control%cutoff_cou
               WRITE (iounit, '(A,T71,F10.3)') " STDA | New value [a.u.]", rc
               WRITE (iounit, *)
            END IF
         END IF
         stda_control%cutoff_cou = rc
      ELSE IF ((.NOT. explicit) .AND. (ALL(cell%perd == 0))) THEN
         rc = MAX(plane_distance(1, 0, 0, cell), &
                  plane_distance(0, 1, 0, cell), &
                  plane_distance(0, 0, 1, cell))
         IF (rc /= stda_control%cutoff_cou) THEN
            IF (iounit > 0) THEN
               WRITE (iounit, *)
               WRITE (iounit, '(A,T37,A)') " STDA |", &
                  " Coulomb Integral cutoff was redefined"
               WRITE (iounit, '(A,T71,F10.3)') " STDA | Old value [a.u.]", &
                  stda_control%cutoff_cou
               WRITE (iounit, '(A,T71,F10.3)') " STDA | New value [a.u.]", rc
               WRITE (iounit, *)
            END IF
         END IF
         stda_control%cutoff_cou = rc
      ELSE
         IF (iounit > 0) THEN
            WRITE (iounit, *)
            WRITE (iounit, '(A,T50,A)') " STDA |", &
               " Coulomb Integral cutoff values"
            WRITE (iounit, '(A,T71,F10.3)') " STDA | Cutoff [a.u.]", &
               stda_control%cutoff_cou
            WRITE (iounit, *)
         END IF
      ENDIF

      CALL timestop(handle)

   END SUBROUTINE stda_init_param
! **************************************************************************************************
!> \brief Allocate the sTDA environment
!> \param qs_env ...
!> \param stda_env ...
! **************************************************************************************************
   SUBROUTINE allocate_stda_env(qs_env, stda_env)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type), POINTER                       :: stda_env

      CHARACTER(len=*), PARAMETER :: routineN = 'allocate_stda_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ii, ispin, n_ao, n_occ, natoms, &
                                                            nkind, nsgf, nspins
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_pool_p_type), DIMENSION(:), POINTER     :: ao_mo_fm_pools
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_fm_type), POINTER                          :: qs_env_c
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(particle_type), DIMENSION(:), POINTER         :: particle_set
      TYPE(qs_kind_type), DIMENSION(:), POINTER          :: qs_kind_set
      TYPE(stda_atom_type), POINTER                      :: atom_param
      TYPE(stda_kind_type), POINTER                      :: kind_param

      CALL timeset(routineN, handle)

      ALLOCATE (stda_env)

      NULLIFY (dft_control)
      CALL get_qs_env(qs_env, dft_control=dft_control)
      nspins = dft_control%nspins

      stda_env%hfx_fraction = 0.0d0
      stda_env%alpha_param = 0.0d0
      stda_env%beta_param = 0.0d0
      stda_env%n_occ = 0

! get number of occ orbitals / atomic orbitals
      DO ispin = 1, nspins
         CALL get_mo_set(qs_env%mos(ispin)%mo_set, mo_coeff=qs_env_c)
         CALL cp_fm_get_info(qs_env_c, ncol_global=n_occ, nrow_global=n_ao)
         stda_env%n_occ(ispin) = n_occ
         stda_env%n_ao(ispin) = n_ao
      ENDDO

! initialize stda_kind_parameters
      NULLIFY (qs_kind_set) !get different kinds
      CALL get_qs_env(qs_env, qs_kind_set=qs_kind_set)
      nkind = SIZE(qs_kind_set)

      ALLOCATE (stda_env%kind_param_set(nkind))
      NULLIFY (kind_param)

      DO ii = 1, nkind
         ALLOCATE (stda_env%kind_param_set(ii)%kind_param)
         kind_param => stda_env%kind_param_set(ii)%kind_param
         CALL allocate_stda_kind_param(kind_param)
      ENDDO

! initialize stda_atom_parameters
      NULLIFY (particle_set) ! get total number of atoms
      CALL get_qs_env(qs_env, particle_set=particle_set)
      natoms = SIZE(particle_set)

      ALLOCATE (stda_env%atom_param_set(natoms))
      NULLIFY (atom_param)

      DO ii = 1, natoms
         ALLOCATE (stda_env%atom_param_set(ii)%atom_param)
         atom_param => stda_env%atom_param_set(ii)%atom_param
         CALL allocate_stda_atom_param(qs_env, stda_env, atom_param)
      ENDDO

      NULLIFY (stda_env%gamma_coulomb, stda_env%gamma_exchange)
      NULLIFY (stda_env%ctransformed)
      NULLIFY (stda_env%xtransformed)
      NULLIFY (ao_mo_fm_pools)
      CALL mpools_get(qs_env%mpools, ao_mo_fm_pools=ao_mo_fm_pools)
      CALL fm_pools_create_fm_vect(ao_mo_fm_pools, stda_env%ctransformed, name=routineP//"ctransformed")
      CALL fm_pools_create_fm_vect(ao_mo_fm_pools, stda_env%xtransformed, name=routineP//"xtransformed")

      NULLIFY (para_env, blacs_env)
      CALL get_qs_env(qs_env, &
                      para_env=para_env, &
                      blacs_env=blacs_env)

      NULLIFY (stda_env%shalf)
      CALL get_qs_kind_set(qs_kind_set, nsgf=nsgf)

      NULLIFY (fm_struct)
      CALL cp_fm_struct_create(fm_struct, para_env, blacs_env, nsgf, nsgf)
      CALL cp_fm_create(stda_env%shalf, fm_struct, routineP//"shalf")
      CALL cp_fm_struct_release(fm_struct)

      CALL timestop(handle)

   END SUBROUTINE allocate_stda_env
! **************************************************************************************************
!> \brief Deallocate the sTDA environment
!> \param stda_env ...
! **************************************************************************************************
   SUBROUTINE deallocate_stda_env(stda_env)

      TYPE(stda_env_type), POINTER                       :: stda_env

      CHARACTER(len=*), PARAMETER :: routineN = 'deallocate_stda_env', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: handle, ii
      TYPE(stda_atom_type), POINTER                      :: atom_param
      TYPE(stda_kind_type), POINTER                      :: kind_param

      CALL timeset(routineN, handle)

      IF (ASSOCIATED(stda_env)) THEN
         IF (ASSOCIATED(stda_env%gamma_coulomb)) THEN
            CALL dbcsr_deallocate_matrix_set(stda_env%gamma_coulomb)
         ENDIF
         IF (ASSOCIATED(stda_env%gamma_exchange)) THEN
            CALL dbcsr_deallocate_matrix_set(stda_env%gamma_exchange)
         ENDIF
         IF (ASSOCIATED(stda_env%ctransformed)) THEN
            CALL cp_fm_vect_dealloc(stda_env%ctransformed)
            NULLIFY (stda_env%ctransformed)
         ENDIF
         IF (ASSOCIATED(stda_env%xtransformed)) THEN
            CALL cp_fm_vect_dealloc(stda_env%xtransformed)
            NULLIFY (stda_env%xtransformed)
         ENDIF
         IF (ASSOCIATED(stda_env%shalf)) THEN
            CALL cp_fm_release(stda_env%shalf)
            NULLIFY (stda_env%shalf)
         ENDIF
! deallocate stda_kind_parameters
         IF (ASSOCIATED(stda_env%kind_param_set)) THEN
            DO ii = 1, SIZE(stda_env%kind_param_set)
               NULLIFY (kind_param)
               kind_param => stda_env%kind_param_set(ii)%kind_param
               CALL deallocate_stda_kind_param(kind_param)
            ENDDO
            DEALLOCATE (stda_env%kind_param_set)
            NULLIFY (stda_env%kind_param_set)
         ENDIF
! deallocate stda_atom_parameters
         IF (ASSOCIATED(stda_env%atom_param_set)) THEN
            DO ii = 1, SIZE(stda_env%atom_param_set)
               NULLIFY (atom_param)
               atom_param => stda_env%atom_param_set(ii)%atom_param
               CALL deallocate_stda_atom_param(atom_param)
            ENDDO
            DEALLOCATE (stda_env%atom_param_set)
            NULLIFY (stda_env%atom_param_set)
         ENDIF
         DEALLOCATE (stda_env)
      ENDIF
      NULLIFY (stda_env)

      CALL timestop(handle)

   END SUBROUTINE deallocate_stda_env
! **************************************************************************************************
!> \brief Allocate sTDA kind parameter
!> \param kind_param ...
! **************************************************************************************************
   SUBROUTINE allocate_stda_kind_param(kind_param)

      TYPE(stda_kind_type), POINTER                      :: kind_param

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_stda_kind_param', &
         routineP = moduleN//':'//routineN

      IF (ASSOCIATED(kind_param)) &
         CALL deallocate_stda_kind_param(kind_param)

      ALLOCATE (kind_param)

      kind_param%symbol = ""
      kind_param%z = 0
      kind_param%kind_number = 0
      kind_param%hardness_param = 0.0_dp
      kind_param%rcut = 0.0_dp

   END SUBROUTINE allocate_stda_kind_param
! **************************************************************************************************
!> \brief Deallocate sTDA kind parameter
!> \param kind_param ...
! **************************************************************************************************
   SUBROUTINE deallocate_stda_kind_param(kind_param)

      TYPE(stda_kind_type), POINTER                      :: kind_param

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_stda_kind_param', &
         routineP = moduleN//':'//routineN

      CPASSERT(ASSOCIATED(kind_param))
      DEALLOCATE (kind_param)
      NULLIFY (kind_param)

   END SUBROUTINE deallocate_stda_kind_param
! **************************************************************************************************
!> \brief Allocate sTDA atom parameter
!> \param qs_env ...
!> \param stda_env ...
!> \param atom_param ...
! **************************************************************************************************
   SUBROUTINE allocate_stda_atom_param(qs_env, stda_env, atom_param)

      TYPE(qs_environment_type), POINTER                 :: qs_env
      TYPE(stda_env_type), POINTER                       :: stda_env
      TYPE(stda_atom_type), POINTER                      :: atom_param

      CHARACTER(LEN=*), PARAMETER :: routineN = 'allocate_stda_atom_param', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: iounit, ispin, n_ao, n_occ, nspins
      TYPE(cp_blacs_env_type), POINTER                   :: blacs_env
      TYPE(cp_fm_struct_type), POINTER                   :: fm_struct
      TYPE(cp_logger_type), POINTER                      :: logger
      TYPE(cp_para_env_type), POINTER                    :: para_env
      TYPE(dft_control_type), POINTER                    :: dft_control
      TYPE(tddfpt_control_type), POINTER                 :: tddfpt_control

! n_unocc

      NULLIFY (logger)
      logger => cp_get_default_logger()
      iounit = cp_logger_get_default_io_unit(logger)

      IF (ASSOCIATED(atom_param)) &
         CALL deallocate_stda_atom_param(atom_param)

      NULLIFY (dft_control, para_env, blacs_env)
      CALL get_qs_env(qs_env, &
                      dft_control=dft_control, &
                      para_env=para_env, &
                      blacs_env=blacs_env)

      NULLIFY (tddfpt_control)
      tddfpt_control => dft_control%tddfpt_control

      nspins = dft_control%nspins
      ALLOCATE (atom_param)
      ALLOCATE (atom_param%charges(nspins))
      ALLOCATE (atom_param%xtransformed_double(nspins))
      ALLOCATE (atom_param%ctransformed_double(nspins))
      ALLOCATE (atom_param%gmcharge(nspins))

      DO ispin = 1, nspins
         n_occ = stda_env%n_occ(ispin)
         n_ao = stda_env%n_ao(ispin)
!charges
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, para_env, blacs_env, n_occ, n_occ)
         CALL cp_fm_create(atom_param%charges(ispin)%matrix, fm_struct, routineP//"charges")
         CALL cp_fm_struct_release(fm_struct)
!xtransformed_double
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, para_env, blacs_env, n_ao, n_occ)
         CALL cp_fm_create(atom_param%xtransformed_double(ispin)%matrix, fm_struct, routineP//"xtransformed_double")
         CALL cp_fm_struct_release(fm_struct)
!ctransformed_double
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, para_env, blacs_env, n_ao, n_occ)
         CALL cp_fm_create(atom_param%ctransformed_double(ispin)%matrix, fm_struct, routineP//"ctransformed_double")
         CALL cp_fm_struct_release(fm_struct)
!gmcharge
         NULLIFY (fm_struct)
         CALL cp_fm_struct_create(fm_struct, para_env, blacs_env, n_occ, n_occ)
         CALL cp_fm_create(atom_param%gmcharge(ispin)%matrix, fm_struct, routineP//"gmcharge")
         CALL cp_fm_struct_release(fm_struct)
      ENDDO

   END SUBROUTINE allocate_stda_atom_param
! **************************************************************************************************
!> \brief Deallocate sTDA atom parameter
!> \param atom_param ...
! **************************************************************************************************
   SUBROUTINE deallocate_stda_atom_param(atom_param)

      TYPE(stda_atom_type), POINTER                      :: atom_param

      CHARACTER(LEN=*), PARAMETER :: routineN = 'deallocate_stda_atom_param', &
         routineP = moduleN//':'//routineN

      INTEGER                                            :: ispin

      CPASSERT(ASSOCIATED(atom_param))
      IF (ASSOCIATED(atom_param)) THEN
!charges
         IF (ASSOCIATED(atom_param%charges)) THEN
            DO ispin = 1, SIZE(atom_param%charges)
               IF (ASSOCIATED(atom_param%charges(ispin)%matrix)) &
                  CALL cp_fm_release(atom_param%charges(ispin)%matrix)
            ENDDO
            DEALLOCATE (atom_param%charges)
            NULLIFY (atom_param%charges)
         ENDIF
!ctransformed_double
         IF (ASSOCIATED(atom_param%ctransformed_double)) THEN
            DO ispin = 1, SIZE(atom_param%ctransformed_double)
               IF (ASSOCIATED(atom_param%ctransformed_double(ispin)%matrix)) &
                  CALL cp_fm_release(atom_param%ctransformed_double(ispin)%matrix)
            ENDDO
            DEALLOCATE (atom_param%ctransformed_double)
            NULLIFY (atom_param%ctransformed_double)
         ENDIF
!gmcharge
         IF (ASSOCIATED(atom_param%gmcharge)) THEN
            DO ispin = 1, SIZE(atom_param%gmcharge)
               IF (ASSOCIATED(atom_param%gmcharge(ispin)%matrix)) &
                  CALL cp_fm_release(atom_param%gmcharge(ispin)%matrix)
            ENDDO
            DEALLOCATE (atom_param%gmcharge)
            NULLIFY (atom_param%gmcharge)
         ENDIF
!xtransformed_double
         IF (ASSOCIATED(atom_param%xtransformed_double)) THEN
            DO ispin = 1, SIZE(atom_param%xtransformed_double)
               IF (ASSOCIATED(atom_param%xtransformed_double(ispin)%matrix)) &
                  CALL cp_fm_release(atom_param%xtransformed_double(ispin)%matrix)
            ENDDO
            DEALLOCATE (atom_param%xtransformed_double)
            NULLIFY (atom_param%xtransformed_double)
         ENDIF
      ENDIF
      DEALLOCATE (atom_param)
      NULLIFY (atom_param)
   END SUBROUTINE deallocate_stda_atom_param
END MODULE qs_tddfpt_stda_types
