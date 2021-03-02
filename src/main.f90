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
REAL(dp), PARAMETER             :: Ha_to_eV=27.21138386_dp
REAL(dp), PARAMETER             :: au_to_eV_per_A=51.42208619083232_dp
!   -------------------------------------------------
CHARACTER(LEN=100)              :: in_file
!   -------------------------------------------------
INTEGER                         :: nb_argument
!   -------------------------------------------------
REAL(dp), ALLOCATABLE           :: cell_mat(:,:), coord_mat(:,:,:), force_mat(:,:,:)
!REAL(dp), ALLOCATABLE           :: energy_mat(:,:), virial_mat(:,:)
REAL(dp), ALLOCATABLE           :: nb_atm_from_coord(:), pot_energy(:)
CHARACTER(LEN=3), ALLOCATABLE   :: atm_name_from_coord(:,:), atm_name_from_force(:,:)
INTEGER, ALLOCATABLE            :: atm_type(:,:)
!   -------------------------------------------------
INTEGER                         :: s, i, j
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
END IF

IF ( in_cell_file .NE. '0' ) in_cell_file = TRIM( in_cell_file )
IF ( in_coord_file .NE. '0' ) in_coord_file = TRIM( in_coord_file )
IF ( in_force_file .NE. '0' ) in_force_file = TRIM( in_force_file )
IF ( in_energy_file .NE. '0' ) in_energy_file = TRIM( in_energy_file )
IF ( in_virial_file .NE. '0' ) in_virial_file = TRIM( in_virial_file )


IF ( in_cell_file .NE. '0' ) THEN
    PRINT*, "Reading cell file..."
    ALLOCATE(cell_mat(3,nb_step))
    OPEN(UNIT=20, FILE=in_cell_file, STATUS='old', FORM='formatted', ACTION='READ')
        READ(20,*)
        DO s = 1, nb_step
            READ(20,*) DUMMY, DUMMY, cell_mat(1,s), DUMMY, DUMMY, DUMMY, cell_mat(2,s), DUMMY, DUMMY, DUMMY, cell_mat(3,s), DUMMY
        END DO
    CLOSE(UNIT=20)
    PRINT*, "Done reading cell file."
    PRINT*, "Writing cell.raw (Å to Å)..."
    OPEN(UNIT=30, FILE='cell.raw')
        DO s = 1, nb_step
            WRITE(30,'(F22.10,F22.10,F22.10)') cell_mat(1,s), cell_mat(2,s), cell_mat(3,s)
        END DO
    CLOSE(UNIT=30)
    DEALLOCATE(cell_mat)
    PRINT*, "Done writing cell.raw (Å to Å)."
END IF

IF ( in_coord_file .NE. '0' ) THEN
    PRINT*, "Reading coord file..."
    ALLOCATE(coord_mat(3,nb_atm,nb_step))
    ALLOCATE(nb_atm_from_coord(nb_step))
    ALLOCATE(atm_name_from_coord(nb_atm,nb_step))
    ALLOCATE(atm_type(nb_atm,nb_step))
    ALLOCATE(pot_energy(nb_step))
    OPEN(UNIT=21, FILE=in_coord_file, STATUS='old', FORM='formatted', ACTION='READ')
        DO s = 1, nb_step
            READ(21,*) nb_atm_from_coord(s)
            IF ( nb_atm_from_coord(s) .NE. nb_atm ) THEN
                PRINT*, 'Number of atom mismatch between input and coord file, exiting...', s, nb_atm, nb_atm_from_coord(s)
                STOP
            END IF
            READ(21,*) DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, DUMMY, pot_energy(s)
            DO i = 1, nb_atm
                READ(21,*) atm_name_from_coord(i,s), coord_mat(1,i,s), coord_mat(2,i,s), coord_mat(3,i,s)
            END DO
        END DO
    CLOSE(UNIT=21)
    PRINT*, "Done coord file..."
    PRINT*, "Writing coord.raw (Å to Å)..."
    DEALLOCATE(nb_atm_from_coord)
    OPEN(UNIT=31, FILE='coord.raw')
        DO s = 1, nb_step
            WRITE(31,'(*(F22.10))', ADVANCE='no') coord_mat(:,:,s) ! NEED A GOOD FORMATING
            WRITE(31,'()')
        END DO
    CLOSE(UNIT=31)
    DEALLOCATE(coord_mat)
    PRINT*, "Done writing coord.raw (Å to Å)."
    PRINT*, "Writing type.raw..."
    OPEN(UNIT=32, FILE='type.raw')
        DO s = 1, nb_step
            DO i = 1, nb_atm
                IF ( atm_name_from_coord(i,s) .EQ. 'C' ) THEN
                    atm_type(i,s) = 6
                ELSE IF ( atm_name_from_coord(i,s) .EQ. 'O' ) THEN
                    atm_type(i,s) = 8
                ELSE IF ( atm_name_from_coord(i,s) .EQ. 'N' ) THEN
                    atm_type(i,s) = 7
                ELSE IF ( atm_name_from_coord(i,s) .EQ. 'H' ) THEN
                    atm_type(i,s) = 1
                ELSE
                    PRINT*, 'Atom number unknown (Only C/O/N/H), exiting...', s, i, atm_name_from_coord(i,s)
                    STOP
                END IF
            END DO
            WRITE(32,'(*(I3))', ADVANCE='no') atm_type(:,s)
            WRITE(32,'()')
        END DO
    CLOSE(UNIT=32)
    DEALLOCATE(atm_type)
    IF ( 'in_force_file' .EQ. '0' ) DEALLOCATE(atm_name_from_coord)
    PRINT*, "Done writing type.raw."
    PRINT*, "Writing energy.raw (from Ha to eV)..."
    OPEN(UNIT=33, FILE='energy.raw')
        DO s = 1, nb_step
            WRITE(33,'(F22.10)') Ha_to_eV*pot_energy(s)
        END DO
    CLOSE(UNIT=33)
    DEALLOCATE(pot_energy)
    PRINT*, "Done writing energy.raw (from Ha to eV)."
END IF

IF ( in_force_file .NE. '0' ) THEN
    PRINT*, "Reading force file..."
    ALLOCATE(force_mat(3,nb_atm,nb_step))
    ALLOCATE(atm_name_from_force(i,s))
    OPEN(UNIT=24, FILE=in_force_file, STATUS='old', FORM='formatted', ACTION='READ')
        DO s = 1, nb_step
            DO j = 1, 4
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
            READ(24,*)
        END DO
    CLOSE(UNIT=24)
    PRINT*, "Done reading force file."
    DEALLOCATE(atm_name_from_force)
    PRINT*, "Writing force.raw (from a.u. to eV/Å)..."
    OPEN(UNIT=34, FILE='force.raw')
        DO s = 1, nb_step
            WRITE(34,'(*(F15.8))', ADVANCE='no') au_to_eV_per_A*force_mat(:,:,s) ! NEED A GOOD FORMATING
            WRITE(34,'()')
        END DO
    CLOSE(UNIT=34)
    DEALLOCATE(force_mat)
    PRINT*, "Done force.raw (from au to eV/Å)."
END IF

END PROGRAM main