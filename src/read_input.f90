! Author: Rolf David
! Date: 03/2021
! License: GNU AGPLv3
! UTF-8, LF, Fortran2003

MODULE read_input

IMPLICIT NONE
!   ------------------------------------------------- Set Double precision
INTEGER, PARAMETER, PRIVATE :: dp=KIND(0.0d0)
!   ------------------------------------------------- Reading related variables
CHARACTER(LEN=100)          :: label, value
INTEGER                     :: iostatus=0
INTEGER                     :: line_c=0
INTEGER                     :: delim
!   ------------------------------------------------- Variables
CHARACTER(LEN=64)           :: in_cell_file='0', in_coord_file='0'
CHARACTER(LEN=64)           :: in_force_file='0', in_energy_file='0'
CHARACTER(LEN=64)           :: in_virial_file='0', in_forceeval_file='0'
INTEGER                     :: nb_atm=0, nb_step=0, nb_stride=1
REAL(dp)                    :: in_box(4)=0
CHARACTER(LEN=1)            :: print_type_raw_every_step='N'
CHARACTER(LEN=2)            :: get_energy_from="0"

!   -------------------------------------------------

CONTAINS

SUBROUTINE sb_read_input(in_file)

    CHARACTER(LEN=100)      :: in_file

    OPEN(99, file=in_file)
    PRINT'(A100)', 'Input parameters'
    PRINT'(A100)', '--------------------------------------------------' &
    , '--------------------------------------------------'
    RL:DO WHILE ( iostatus == 0 )
        READ(99, '(A)', IOSTAT=iostatus ) value
        IF ( iostatus == 0 )then
            line_c=line_c + 1
            delim=SCAN( value, '    ' )
            label=value(1:delim)
            value=value(delim + 1:)
            IF ( label(1:1) == '!' )THEN
                CYCLE RL
            END IF
            SELECT CASE (label)
                CASE ('in_cell_file')
                    READ(value, * , IOSTAT=iostatus) in_cell_file
                    PRINT'(A50,A64)', 'in_cell_file:', ADJUSTR( in_cell_file )
                CASE ('in_coord_file')
                    READ(value, * , IOSTAT=iostatus) in_coord_file
                    PRINT'(A50,A64)', 'in_coord_file:', ADJUSTR( in_coord_file )
                CASE ('in_force_file')
                    READ(value, * , IOSTAT=iostatus) in_force_file
                    PRINT'(A50,A64)', 'in_force_file:', ADJUSTR( in_force_file )
                CASE ('in_energy_file')
                    READ(value, * , IOSTAT=iostatus) in_energy_file
                    PRINT'(A50,A65)', 'in_energy_file:', ADJUSTR( in_energy_file )
                CASE ('in_forceeval_file')
                    READ(value, * , IOSTAT=iostatus) in_forceeval_file
                    PRINT'(A50,A64)', 'in_forceeval_file:', ADJUSTR( in_forceeval_file )
                CASE ('in_virial_file')
                    READ(value, * , IOSTAT=iostatus) in_virial_file
                    PRINT'(A50,A64)', 'in_virial_file:', ADJUSTR( in_virial_file )
                CASE ('get_energy_from')
                    READ(value, * , IOSTAT=iostatus) get_energy_from
                    PRINT'(A50,A64)', 'get_energy_from:', ADJUSTR( get_energy_from )
                CASE ('nb_atm')
                    READ(value, * , IOSTAT=iostatus) nb_atm
                    PRINT'(A50,I64)', 'nb_atm:', nb_atm
                CASE ('nb_step')
                    READ(value, * , IOSTAT=iostatus) nb_step
                    PRINT'(A50,I64)', 'nb_step:', nb_step
                CASE ('nb_stride')
                    READ(value, * , IOSTAT=iostatus) nb_stride
                    PRINT'(A50,I64)', 'nb_stride:', nb_stride
                CASE ('in_box')
                    READ(value, *, IOSTAT=iostatus) in_box(1), in_box(2), in_box(3)
                    in_box(4) = 1
                    PRINT'(A50,E64.10)', 'in_box_x:', in_box(1)
                    PRINT'(A50,E64.10)', 'in_box_y:', in_box(2)
                    PRINT'(A50,E64.10)', 'in_box_z:', in_box(3)
                CASE DEFAULT
                    PRINT'(A50,I64)', 'Invalid label, line:', line_c
            END SELECT
        END IF
    END DO RL
    PRINT'(A100)', '--------------------------------------------------' &
    , '--------------------------------------------------'

END SUBROUTINE sb_read_input

END MODULE read_input