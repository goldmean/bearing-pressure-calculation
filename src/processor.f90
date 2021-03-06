MODULE PROCESSOR

CONTAINS

    SUBROUTINE GET_TIME(TIME)
        IMPLICIT NONE

        INTEGER :: CLOCK_READING, CLOCK_RATE
        DOUBLE PRECISION, INTENT(OUT) :: TIME

        CALL SYSTEM_CLOCK(CLOCK_READING, CLOCK_RATE)

        TIME = FLOAT(CLOCK_READING) / FLOAT(CLOCK_RATE)

    END SUBROUTINE GET_TIME

END MODULE PROCESSOR
