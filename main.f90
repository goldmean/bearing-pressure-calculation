PROGRAM SOLUTION

    USE SUPPLEMENT
    USE PROCESSOR

    IMPLICIT NONE

    INTEGER :: NI, NJ, UNIT, LIMIT, S
    DOUBLE PRECISION :: REGIME, GROOVE_DEPTH, GROOVE_RADIUS, EXTERNAL_RADIUS, &
        EXTERNAL_ANGLE, ATMOSPHERIC_PRESSURE, RADIUS_STEP, ANGLE_STEP, EPS, &
        REGIME_MIN, REGIME_MAX, REGIME_STEP, CLEARANCE, CLEARANCE_STEP, &
        GROOVE_PRESSURE_DROP, LAYER_PRESSURE_DROP, PRESSURE_RATIO_DROP, &
        START_TIME, STOP_TIME

    LOGICAL :: PRESSURE_LOGGING, RIGIDITY_LOGGING, ANALYTIC_LOGGING
    CHARACTER(16) :: ACCURACY, RESIDUAL_, PRESSURE_, RIGIDITY_, ANALYTIC_
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:, :) :: RADIUS, ANGLE, PRESSURE, PREVIOUS

    NAMELIST /PARAMETERS/ NI, NJ, REGIME, GROOVE_DEPTH, GROOVE_RADIUS, &
        EXTERNAL_RADIUS, EXTERNAL_ANGLE, ATMOSPHERIC_PRESSURE, LIMIT, EPS, &
        REGIME_MIN, REGIME_MAX, REGIME_STEP, CLEARANCE, CLEARANCE_STEP, &
        PRESSURE_LOGGING, RIGIDITY_LOGGING, ANALYTIC_LOGGING, ACCURACY, &
        RESIDUAL_, PRESSURE_, RIGIDITY_, ANALYTIC_

    OPEN(UNIT=UNIT, FILE="PARAMETERS.NML", ACTION="READ")

    READ(UNIT=UNIT, NML=PARAMETERS)

    CLOSE(UNIT=UNIT)

    EXTERNAL_ANGLE = 4.0 * ATAN(1.0) * EXTERNAL_ANGLE / 180.

    S = NINT(GROOVE_RADIUS / EXTERNAL_RADIUS * NI)

    ALLOCATE(RADIUS(NI, NJ), ANGLE(NI, NJ), PRESSURE(NI, NJ), PREVIOUS(NI, NJ))

    CALL CREATE_MESH(NI, NJ, EXTERNAL_RADIUS, EXTERNAL_ANGLE, RADIUS_STEP, &
        ANGLE_STEP, RADIUS, ANGLE)

    CALL GET_TIME(START_TIME)

    IF (PRESSURE_LOGGING) CALL CALCULATE_PRESSURE(NI, NJ, REGIME, GROOVE_DEPTH, &
        GROOVE_RADIUS, EXTERNAL_RADIUS, ATMOSPHERIC_PRESSURE, RADIUS_STEP, ANGLE_STEP, &
        RADIUS, PRESSURE, PREVIOUS, .TRUE., .FALSE., ACCURACY, RESIDUAL_, LIMIT, EPS)

    CALL CREATE_FILES(NI, NJ, RADIUS, ANGLE, PRESSURE, PRESSURE_)

    IF (ANALYTIC_LOGGING) CALL CALCULATE_ANALYTICAL_SOLUTION(NI, NJ, REGIME, &
        GROOVE_RADIUS, EXTERNAL_RADIUS, ATMOSPHERIC_PRESSURE, RADIUS_STEP, ANGLE_STEP, &
        LIMIT, EPS, RADIUS, PRESSURE, PREVIOUS, ACCURACY, RESIDUAL_, ANALYTIC_)

    IF (RIGIDITY_LOGGING) CALL CALCULATE_RIGIDITY(NI, NJ, REGIME, GROOVE_DEPTH, &
        GROOVE_RADIUS, EXTERNAL_RADIUS, ATMOSPHERIC_PRESSURE, RADIUS_STEP, ANGLE_STEP, &
        RADIUS, PRESSURE, PREVIOUS, ACCURACY, RESIDUAL_, RIGIDITY_, LIMIT, EPS, &
        REGIME_MIN, REGIME_MAX, REGIME_STEP, CLEARANCE, CLEARANCE_STEP)

    CALL GET_TIME(STOP_TIME)

    IF (PRESSURE_LOGGING) THEN
        GROOVE_PRESSURE_DROP = MAXVAL(PRESSURE(1: S, 1)) - MINVAL(PRESSURE(1: S, 1))
        LAYER_PRESSURE_DROP = MAXVAL(PRESSURE) - MINVAL(PRESSURE)
        PRESSURE_RATIO_DROP = GROOVE_PRESSURE_DROP / LAYER_PRESSURE_DROP

        PRINT "(A, / 3(A, F5.3, 3X), /, A)", REPEAT("-", 88), &
            "GROOVE PRESSURE DROP = ", GROOVE_PRESSURE_DROP, &
            "LAYER PRESSURE DROP = ", LAYER_PRESSURE_DROP, &
            "PRESSURE RATIO DROP = ", PRESSURE_RATIO_DROP, REPEAT("-", 88)
    END IF

    PRINT "(A, F8.3, 1X, A)", "EXECUTION TIME: ", STOP_TIME - START_TIME, "SEC"

    DEALLOCATE(RADIUS, ANGLE, PRESSURE, PREVIOUS)

END PROGRAM SOLUTION
