! Author: Rolf David
! Date: 03/2021
! License: GNU AGPLv3
! UTF-8, LF, Fortran2003

PROGRAM main
USE read_input

IMPLICIT NONE
!   ------------------------------------------------- Set Double precision
INTEGER, PARAMETER              :: dp=KIND(0.0d0)
!   -------------------------------------------------
REAL(dp), PARAMETER             :: Ha_to_eV=27.211386245988_dp
REAL(dp), PARAMETER             :: Bohr_to_A=0.529177210903_dp
REAL(dp), PARAMETER             :: au_to_eV_per_A=Ha_to_eV/Bohr_to_A
REAL(dp), PARAMETER             :: eV_per_A3_to_GPa=160.21766208_dp

!   -------------------------------------------------
CHARACTER(LEN=100)              :: in_file
!   -------------------------------------------------
INTEGER                         :: nb_argument
!   -------------------------------------------------
REAL(dp), ALLOCATABLE           :: cell_mat(:,:), coord_mat(:,:,:), force_mat(:,:,:)
REAL(dp), ALLOCATABLE           :: energy_mat(:), stress_mat(:,:), cell_volume(:)
REAL(dp), ALLOCATABLE           :: nb_atm_from_coord(:), pot_energy(:)
CHARACTER(LEN=3), ALLOCATABLE   :: atm_name_from_coord(:,:), atm_name_from_force(:,:)
INTEGER, ALLOCATABLE            :: atm_type(:,:)
TYPE t_unique_atm
    CHARACTER(LEN=3),ALLOCATABLE:: u_name(:,:)
    INTEGER,ALLOCATABLE         :: u_type_from_zero(:,:)
END TYPE t_unique_atm
type(t_unique_atm)              :: u_atm
!   -------------------------------------------------
INTEGER                         :: s, i, j, o
CHARACTER(LEN=64)               :: DUMMY

!   ----------------------------------------------- Get arguments (filenames, choices)
nb_argument = COMMAND_ARGUMENT_COUNT()

IF ( nb_argument .EQ. 0 ) THEN
    PRINT*, "No input file, exiting..."
    STOP
END IF

CALL GET_COMMAND_ARGUMENT(1, in_file)
in_file=TRIM(in_file)
call sb_read_input(in_file)

IF ( ( nb_atm .EQ. 0 ) .OR. (nb_step .EQ. 0 ) ) THEN
    PRINT*, "Either the number of atoms ", nb_atm&
    , " or the number of steps ", nb_step, " is not defined (or both), exiting..."
    STOP
END IF

IF ( in_cell_file .NE. '0' ) in_cell_file = TRIM( in_cell_file )
IF ( in_coord_file .NE. '0' ) in_coord_file = TRIM( in_coord_file )
IF ( in_force_file .NE. '0' ) in_force_file = TRIM( in_force_file )
IF ( in_energy_file .NE. '0' ) in_energy_file = TRIM( in_energy_file )
IF ( in_stresstensor_file .NE. '0' ) in_stresstensor_file = TRIM( in_stresstensor_file )
IF ( in_forceeval_file .NE. '0' ) in_forceeval_file = TRIM( in_forceeval_file )

IF ( in_box(4) .EQ. 1 ) THEN
    PRINT*, "Writing box.raw with provided box size (in Å)..."
    OPEN(UNIT=30, FILE='box.raw')
        DO s = 1, nb_step, nb_stride
            WRITE(30,'(F22.10,I2,I2,I2,F22.10,I2,I2,I2,F22.10)') in_box(1), 0, 0, 0, in_box(2), 0, 0, 0, in_box(3)
        END DO
    CLOSE(UNIT=30)
    PRINT*, "Done writing box.raw with provided box size (in Å)."
ELSE IF ( in_cell_file .NE. '0' ) THEN
    PRINT*, "Reading cell file..."
    ALLOCATE(cell_mat(9,nb_step))
    cell_mat(:,:) = 0.0_dp
    OPEN(UNIT=20, FILE=in_cell_file, STATUS='old', FORM='formatted', ACTION='READ')
        s = 1
        DO WHILE( s .LE. nb_step )
            READ(20,*) DUMMY
            IF ( DUMMY .NE. "#" ) THEN
                BACKSPACE(20)
                READ(20,*) DUMMY, DUMMY, cell_mat(1,s), cell_mat(2,s), cell_mat(3,s) &
                , cell_mat(4,s), cell_mat(5,s), cell_mat(6,s) &
                , cell_mat(7,s), cell_mat(8,s), cell_mat(9,s), DUMMY
                s = s + 1
            END IF
        END DO
    CLOSE(UNIT=20)
    PRINT*, "Done reading cell file."
    PRINT*, "Writing box.raw (Å to Å)..."
    OPEN(UNIT=30, FILE='box.raw')
        DO s = 1, nb_step, nb_stride
            WRITE(30,'(*(F22.10))', ADVANCE='no') cell_mat(:,s)
            WRITE(30,'()')
        END DO
    CLOSE(UNIT=30)
    IF (in_stresstensor_file .EQ. '0') DEALLOCATE(cell_mat)
    PRINT*, "Done writing box.raw (Å to Å)."
END IF
IF ( in_stresstensor_file .NE. '0' ) THEN
    IF  ( (in_box(4) .NE. 1 ) .AND. (in_cell_file .EQ. '0') ) THEN
        PRINT*, 'The virial file cannot be written if cell information is not present (volume), exiting...'
        STOP
    ELSE
        IF ( (in_cell_file .NE. '0') .AND. (in_box(4) .NE. 1 ) ) THEN
            IF ( ( ANY( cell_mat(2,:) .NE. 0 ) ) .OR. ( ANY( cell_mat(3,:) .NE. 0 ) ) .OR.&
            ( ANY( cell_mat(4,:) .NE. 0 ) ) .OR. ( ANY( cell_mat(6,:) .NE. 0 ) ) .OR. &
            ( ANY( cell_mat(7,:) .NE. 0 ) ) .OR. ( ANY( cell_mat(8,:) .NE. 0 ) ) ) THEN
                PRINT*, "Virial from stress tensor only for orthorombic cell (for now...), exiting..."
                STOP
            END IF
        END IF
        ALLOCATE(cell_volume(nb_step))
        cell_volume(:) = 0.0_dp
        DO s = 1, nb_step, nb_stride
            IF ( in_box(4) .EQ. 1 ) THEN
                cell_volume(s) = in_box(1) * in_box(2) * in_box(3)
            ELSE IF ( in_cell_file .NE. '0' ) THEN
                cell_volume(s) = cell_mat(1,s) * cell_mat(5,s) * cell_mat(9,s)
            END IF
        END DO
        IF ( in_cell_file .NE. '0' ) DEALLOCATE(cell_mat)
        PRINT*, "Reading stress tensor file..."
        ALLOCATE(stress_mat(9,nb_step))
        stress_mat(:,:) = 0.0_dp
        OPEN(UNIT=26, FILE=in_stresstensor_file, STATUS='OLD', FORM='formatted', ACTION='READ')
        s = 1
        DO WHILE( s .LE. nb_step )
            READ(26,'(A)') DUMMY
            ! CP2K (< 8.1)
            IF ( ( INDEX(TRIM(ADJUSTL(DUMMY)),"STRESS TENSOR") .EQ. 1 ) .OR. &
            ( INDEX(TRIM(ADJUSTL(DUMMY)),"NUMERICAL STRESS TENSOR") .EQ. 1 ) ) THEN
                READ(26,*) DUMMY
                READ(26,*) DUMMY, stress_mat(1,s), stress_mat(2,s), stress_mat(3,s)
                READ(26,*) DUMMY, stress_mat(4,s), stress_mat(5,s), stress_mat(6,s)
                READ(26,*) DUMMY, stress_mat(7,s), stress_mat(8,s), stress_mat(9,s)
                s = s + 1
            ! CP2K (>= 8.1) (https://github.com/cp2k/cp2k/issues/583)
            ELSE IF ( ( INDEX(TRIM(ADJUSTL(DUMMY)),"STRESS| Analytical") .EQ. 1 ) .OR. &
                ( INDEX(TRIM(ADJUSTL(DUMMY))," STRESS| Numerical") .EQ. 1 ) ) THEN
                READ(26,*) DUMMY
                READ(26,*) DUMMY, DUMMY, stress_mat(1,s), stress_mat(2,s), stress_mat(3,s)
                READ(26,*) DUMMY, DUMMY, stress_mat(4,s), stress_mat(5,s), stress_mat(6,s)
                READ(26,*) DUMMY, DUMMY, stress_mat(7,s), stress_mat(8,s), stress_mat(9,s)
                s = s + 1
            END IF
        END DO
        CLOSE(UNIT=26)
        PRINT*, "Done reading stress tensor file..."
        PRINT*, "Writing virial.raw (GPa * A3 to eV)..."
        OPEN(UNIT=36, FILE='virial.raw')
        DO s = 1, nb_step, nb_stride
            WRITE(36,'(*(F22.10))', ADVANCE='no') ( stress_mat(:,s) *  cell_volume(s) ) / eV_per_A3_to_GPa
            WRITE(36,'()')
        END DO
        CLOSE(UNIT=36)
        PRINT*, "Done virial.raw (GPa * A3 to eV)..."
        DEALLOCATE(stress_mat,cell_volume)
    END IF
END IF
IF ( in_coord_file .NE. '0' ) THEN
    PRINT*, "Reading coord file..."
    ALLOCATE(coord_mat(3,nb_atm,nb_step))
    coord_mat(:,:,:) = 0.0_dp
    ALLOCATE(nb_atm_from_coord(nb_step))
    nb_atm_from_coord(:) = 0
    ALLOCATE(atm_name_from_coord(nb_atm,nb_step))
    atm_name_from_coord(:,:) = '000'
    ALLOCATE(atm_type(nb_atm,nb_step))
    atm_type(:,:) = 0
    ALLOCATE(pot_energy(nb_step))
    pot_energy(:) = 0.0_dp
    OPEN(UNIT=21, FILE=in_coord_file, STATUS='old', FORM='formatted', ACTION='READ')
        DO s = 1, nb_step
            READ(21,*) nb_atm_from_coord(s)
            IF ( nb_atm_from_coord(s) .NE. nb_atm ) THEN
                PRINT*, 'Number of atom mismatch between input and coord file, exiting...', s, nb_atm, nb_atm_from_coord(s)
                STOP
            END IF
            IF ( get_energy_from .EQ. 'CO' ) THEN
                READ(21,*) DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, pot_energy(s)
            ELSE
                READ(21,*) DUMMY
            END IF
            DO i = 1, nb_atm
                READ(21,*) atm_name_from_coord(i,s), coord_mat(1,i,s), coord_mat(2,i,s), coord_mat(3,i,s)
            END DO
        END DO
    CLOSE(UNIT=21)
    PRINT*, "Done reading coord file..."
    PRINT*, "Writing coord.raw (Å to Å)..."
    DEALLOCATE(nb_atm_from_coord)
    OPEN(UNIT=31, FILE='coord.raw')
        DO s = 1, nb_step, nb_stride
            WRITE(31,'(*(F22.10))', ADVANCE='no') coord_mat(:,:,s)
            WRITE(31,'()')
        END DO
    CLOSE(UNIT=31)
    PRINT*, "Done writing coord.raw (Å to Å)."
    IF ( reprint_xyz .EQ. 'Y') THEN
        OPEN(UNIT=311, FILE='coord.xyz')
        PRINT*, "Writing coord.xyz (Å to Å)..."
        DO s = 1, nb_step, nb_stride
            WRITE(311,'(I10)') nb_atm
            WRITE(311,'(A10,I10)') "Step nb:", s
            DO i = 1, nb_atm
                WRITE(311,'(A3,1X,F22.10,1X,F22.10,1X,F22.10)') ADJUSTL( atm_name_from_coord(i,s) ) &
                , coord_mat(1,i,s), coord_mat(2,i,s), coord_mat(3,i,s)
            END DO
        END DO
        CLOSE(UNIT=311)
        PRINT*, "Done writing coord.xyz (Å to Å)."
    END IF
    DEALLOCATE(coord_mat)
    PRINT*, "Writing type.raw..."
    ALLOCATE(u_atm%u_type_from_zero(nb_atm,nb_step))
    ALLOCATE(u_atm%u_name(nb_atm,nb_step))
    u_atm%u_name(:,:) = '0'
    u_atm%u_type_from_zero(:,:) = -1
    OPEN(UNIT=32, FILE='type.raw')
        DO s = 1, nb_step, nb_stride
            j = 0
            DO i = 1, nb_atm
                IF ( .NOT. ANY( u_atm%u_name(:,s) .EQ. atm_name_from_coord(i,s) ) ) THEN
                    j = j + 1
                    u_atm%u_name(j,s) = atm_name_from_coord(i,s)
                    u_atm%u_type_from_zero(j,s) = j - 1
                    atm_type(i,s) = u_atm%u_type_from_zero(j,s)
                ELSE
                    DO o = 1, nb_atm
                        IF ( u_atm%u_name(o,s) .EQ. atm_name_from_coord(i,s) ) THEN
                            atm_type(i,s) = u_atm%u_type_from_zero(o,s)
                            EXIT
                        END IF
                    END DO
                END IF
            END DO
        END DO
        IF ( print_type_raw_every_step .EQ. 'Y' ) THEN
            DO s = 1, nb_step
                WRITE(32,'(*(I3))', ADVANCE='no') atm_type(:,s)
                WRITE(32,'()')
            END DO
        ELSE
            WRITE(32,'(*(I3))', ADVANCE='no') atm_type(:,1)
            WRITE(32,'()')
        END IF
        DO s = 1, nb_step, nb_stride
            IF( MAXVAL(u_atm%u_type_from_zero(:,s)) .NE. MAXVAL(u_atm%u_type_from_zero) ) THEN
                PRINT*, "Mismatch between of unique atom type", s &
                , MAXVAL(u_atm%u_type_from_zero(:,s)), MAXVAL(u_atm%u_type_from_zero)
                STOP
            END IF
        END DO
    CLOSE(UNIT=32)
    OPEN(UNIT=321, FILE='type_eq.raw')
        DO i = 1, nb_atm
            IF ( u_atm%u_type_from_zero(i,1) .NE. -1 ) THEN
                WRITE(321,'(A3,I3)') u_atm%u_name(i,1), u_atm%u_type_from_zero(i,1)
            ELSE
                EXIT
            END IF
        END DO
    CLOSE(UNIT=321)
    DEALLOCATE(atm_type,u_atm%u_name,u_atm%u_type_from_zero)
    IF ( in_force_file .EQ. '0' ) DEALLOCATE(atm_name_from_coord)
    PRINT*, "Done writing type.raw."
    IF ( get_energy_from .EQ. 'CO' ) THEN
        PRINT*, "Writing energy.raw (from Ha to eV) [from coord file]..."
        OPEN(UNIT=33, FILE='energy.raw')
            DO s = 1, nb_step, nb_stride
                WRITE(33,'(F22.10)') Ha_to_eV*pot_energy(s)
            END DO
        CLOSE(UNIT=33)
        DEALLOCATE(pot_energy)
        PRINT*, "Done writing energy.raw (from Ha to eV) [from coord file]."
    END IF
END IF
IF ( in_force_file .NE. '0' ) THEN
    PRINT*, "Reading force file..."
    ALLOCATE(force_mat(3,nb_atm,nb_step))
    ALLOCATE(atm_name_from_force(nb_atm,nb_step))
    force_mat(:,:,:) = 0.0_dp
    atm_name_from_force(:,:) = '000'
    OPEN(UNIT=24, FILE=in_force_file, STATUS='old', FORM='formatted', ACTION='READ')
    s = 1
    DO WHILE( s .LE. nb_step )
        READ(24,'(A)') DUMMY
        IF  ( INDEX(TRIM(ADJUSTL(DUMMY)),"ATOMIC FORCES in") .EQ. 1 ) THEN
            DO j = 1, 2
                READ(24,*)
            END DO
            DO i = 1, nb_atm
                READ(24,*) DUMMY, DUMMY, atm_name_from_force(i,s), force_mat(1,i,s), force_mat(2,i,s), force_mat(3,i,s)
                IF ( in_coord_file .NE. '0' ) THEN
                    IF ( atm_name_from_force(i,s) .NE. atm_name_from_coord(i,s) ) THEN
                        PRINT*, "Atom name mismatch between coord and force files, exiting...", s &
                        , atm_name_from_coord(i,s), atm_name_from_force(i,s)
                    END IF
                END IF
            END DO
            s = s + 1
        END IF
    END DO
    CLOSE(UNIT=24)
    PRINT*, "Done reading force file."
    DEALLOCATE(atm_name_from_force)
    PRINT*, "Writing force.raw (from a.u. to eV/Å)..."
    OPEN(UNIT=34, FILE='force.raw')
        DO s = 1, nb_step, nb_stride
            WRITE(34,'(*(F15.8))', ADVANCE='no') au_to_eV_per_A*force_mat(:,:,s)
            WRITE(34,'()')
        END DO
    CLOSE(UNIT=34)
    DEALLOCATE(force_mat)
    PRINT*, "Done force.raw (from au to eV/Å)."
END IF
IF ( ( in_forceeval_file .NE. '0') .AND. ( get_energy_from .EQ. 'FE' ) ) THEN
    PRINT*, "Reading force eval file..."
    ALLOCATE(energy_mat(nb_step))
    energy_mat(:) = 0.0_dp
    OPEN(UNIT=25, FILE=in_forceeval_file, STATUS='old', FORM='formatted', ACTION='READ')
    s = 1
    DO WHILE( s .LE. nb_step )
        READ(25,'(A)') DUMMY
        IF  ( INDEX(TRIM(ADJUSTL(DUMMY)),"ENERGY|") .EQ. 1 ) THEN
            BACKSPACE(25)
            READ(25,*) DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, energy_mat(s)
            s = s + 1
        END IF
    END DO
    CLOSE(UNIT=25)
    PRINT*, "Writing energy.raw (from Ha to eV) [from force eval file]..."
    OPEN(UNIT=35, FILE='energy.raw')
        DO s = 1, nb_step, nb_stride
            WRITE(35,'(F30.15)') Ha_to_eV*energy_mat(s)
        END DO
    CLOSE(UNIT=35)
    DEALLOCATE(energy_mat)
    PRINT*, "Done writing energy.raw (from Ha to eV) [from force eval file]."
END IF

END PROGRAM main