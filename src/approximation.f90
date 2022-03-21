MODULE APPROXIMATION

    PUBLIC :: SINGLE_ACCURACY, DOUBLE_ACCURACY

CONTAINS

    SUBROUTINE SINGLE_ACCURACY(NI, NJ, REGIME, GROOVE_DEPTH, GROOVE_RADIUS, &
            EXTERNAL_RADIUS, ATMOSPHERIC_PRESSURE, RADIUS_STEP, ANGLE_STEP, &
            RADIUS, PRESSURE, PREVIOUS)
        IMPLICIT NONE

        INTEGER :: I, J, S
        INTEGER, INTENT(IN) :: NI, NJ
        DOUBLE PRECISION, INTENT(IN) :: REGIME, GROOVE_DEPTH, GROOVE_RADIUS, &
            EXTERNAL_RADIUS, ATMOSPHERIC_PRESSURE, RADIUS_STEP, ANGLE_STEP
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:, :) :: RADIUS, PRESSURE, PREVIOUS

        S = NINT(GROOVE_RADIUS / EXTERNAL_RADIUS * NI)

        FORALL (J = 1: NJ) PRESSURE(NI, J) = ATMOSPHERIC_PRESSURE
        FORALL (I = 2: NI - 1) PRESSURE(I, NJ) = PRESSURE(I, NJ - 1)
        FORALL (I = S + 1: NI - 1) PRESSURE(I, 1) = PRESSURE(I, 2)

        IF (PREVIOUS(1, 1) >= 0.5) THEN
            PRESSURE(1, 1) = SQRT(2.0 * REGIME * SQRT(PREVIOUS(1, 1) * &
                (1 - PREVIOUS(1, 1))) * 2.0 / 3.0 / GROOVE_DEPTH * RADIUS_STEP + &
                PREVIOUS(2, 1) ** 2)
        ELSE IF (PREVIOUS(1, 1) < 0.5) THEN
            PRESSURE(1, 1) = SQRT(REGIME * 2.0 / 3.0 / GROOVE_DEPTH * RADIUS_STEP + &
                PREVIOUS(2, 1) ** 2)
        END IF

        FORALL (J = 2: NJ) PRESSURE(1, J) = PRESSURE(1, 1)

        DO I = 2, S - 1
            PRESSURE(I, 1) = SQRT((GROOVE_DEPTH * (PREVIOUS(I + 1, 1) ** 2 + &
                PRESSURE(I - 1, 1) ** 2) / RADIUS_STEP ** 2 + 2.0 / RADIUS(I, 1) * &
                PREVIOUS(I, 2) ** 2 / ANGLE_STEP) / (2.0 / RADIUS(I, 1) / ANGLE_STEP + &
                GROOVE_DEPTH * 2.0 / RADIUS_STEP ** 2))
        END DO

        DO J = 2, NJ - 1
            PRESSURE(S, J) = SQRT((GROOVE_DEPTH / RADIUS(S, J) ** 2 * &
                (PREVIOUS(S, J + 1) ** 2 + PRESSURE(S, J - 1) ** 2) / &
                ANGLE_STEP ** 2 + (PRESSURE(S - 1, J) ** 2 + &
                PREVIOUS(S + 1, J) ** 2) / RADIUS_STEP) / (2.0 / RADIUS_STEP + &
                2.0 * GROOVE_DEPTH / RADIUS(S, J) ** 2 / ANGLE_STEP ** 2))
        END DO

        PRESSURE(S, 1) = (2.0 / RADIUS(S, 1) * PREVIOUS(S, 2) / ANGLE_STEP + &
            PRESSURE(S - 1, 1) / RADIUS_STEP) / (1.0 / RADIUS_STEP + &
            2.0 / RADIUS(S, 1) / ANGLE_STEP)

        DO I = 2, NI - 1
            DO J = 2, NJ - 1
                IF (I == S) CYCLE
                PRESSURE(I, J) = SQRT(((RADIUS(I + 1, J) + RADIUS(I, J)) / 2.0 * &
                    PREVIOUS(I + 1, J) ** 2 / RADIUS_STEP ** 2 + (RADIUS(I, J) + &
                    RADIUS(I - 1, J)) / 2.0 * PRESSURE(I - 1, J) ** 2 / &
                    RADIUS_STEP ** 2 + 1.0 / RADIUS(I, J) * (PREVIOUS(I, J + 1) ** 2 + &
                    PRESSURE(I, J - 1) ** 2) / ANGLE_STEP ** 2) / ((RADIUS(I + 1, J) + &
                    RADIUS(I, J)) / 2.0 / RADIUS_STEP ** 2 + (RADIUS(I, J) + &
                    RADIUS(I - 1, J)) / 2.0 / RADIUS_STEP ** 2 + 2.0 / RADIUS(I, J) / &
                    ANGLE_STEP ** 2))
            END DO
        END DO

    END SUBROUTINE SINGLE_ACCURACY

    SUBROUTINE DOUBLE_ACCURACY(NI, NJ, REGIME, GROOVE_DEPTH, GROOVE_RADIUS, &
            EXTERNAL_RADIUS, ATMOSPHERIC_PRESSURE, RADIUS_STEP, ANGLE_STEP, &
            RADIUS, PRESSURE, PREVIOUS)
        IMPLICIT NONE

        INTEGER :: I, J, S
        INTEGER, INTENT(IN) :: NI, NJ
        DOUBLE PRECISION, INTENT(IN) :: REGIME, GROOVE_DEPTH, GROOVE_RADIUS, &
            EXTERNAL_RADIUS, ATMOSPHERIC_PRESSURE, RADIUS_STEP, ANGLE_STEP
        DOUBLE PRECISION, INTENT(OUT), DIMENSION(:, :) :: RADIUS, PRESSURE, PREVIOUS

        S = NINT(GROOVE_RADIUS / EXTERNAL_RADIUS * NI)

        FORALL (J = 1: NJ) PRESSURE(NI, J) = ATMOSPHERIC_PRESSURE
        FORALL (I = 2: NI - 1) PRESSURE(I, NJ) = (4.0 * PRESSURE(I, NJ - 1) - &
            PRESSURE(I, NJ - 2)) / 3.0
        FORALL (I = S + 1: NI - 1) PRESSURE(I, 1) = (4.0 * PRESSURE(I, 2) - &
            PRESSURE(I, 3)) / 3.0

        IF (PREVIOUS(1, 1) >= 0.5) THEN
            PRESSURE(1, 1) = SQRT((2.0 * REGIME * SQRT(PREVIOUS(1, 1) * &
                (1 - PREVIOUS(1, 1))) * 2.0 / 3.0 / GROOVE_DEPTH * 2.0 * &
                RADIUS_STEP + 4.0 * PREVIOUS(2, 1) ** 2 - PREVIOUS(3, 1) ** 2) / 3.0)
        ELSE IF (PREVIOUS(1, 1) < 0.5) THEN
            PRESSURE(1, 1) = SQRT((REGIME * 2.0 / 3.0 / GROOVE_DEPTH * 2.0 * &
                RADIUS_STEP + 4.0 * PREVIOUS(2, 1) ** 2 - PREVIOUS(3, 1) ** 2) / 3.0)
        END IF

        FORALL (J = 2: NJ) PRESSURE(1, J) = PRESSURE(1, 1)

        DO I = 2, S - 1
            PRESSURE(I, 1) = SQRT((GROOVE_DEPTH * (PREVIOUS(I + 1, 1) ** 2 + &
                PRESSURE(I - 1, 1) ** 2) / RADIUS_STEP ** 2 + 1.0 / RADIUS(I , 1) * &
                (4.0 * PREVIOUS(I, 2) ** 2 - PREVIOUS(I, 3) ** 2) / ANGLE_STEP) / &
                (2.0 * GROOVE_DEPTH / RADIUS_STEP ** 2 + 3.0 / RADIUS(I, 1) / &
                ANGLE_STEP))
        END DO

        DO J = 2, NJ - 1
            PRESSURE(S, J) = SQRT(((4.0 * PRESSURE(S - 1, J) ** 2 - &
                PRESSURE(S - 2, J) ** 2) / 2.0 / RADIUS_STEP + &
                (4.0 * PREVIOUS(S + 1, J) ** 2 - PREVIOUS(S + 2, J) ** 2) / &
                2.0 / RADIUS_STEP + GROOVE_DEPTH / RADIUS(S, J) ** 2 * &
                (PREVIOUS(S, J + 1) ** 2 + PRESSURE(S, J - 1) ** 2) / &
                ANGLE_STEP ** 2) / (GROOVE_DEPTH / RADIUS(S, J) ** 2 * &
                2.0 / ANGLE_STEP ** 2 + 3.0 / RADIUS_STEP))
        END DO

        PRESSURE(S, 1) = (1.0 / RADIUS(S, 1) * (4.0 * PREVIOUS(S, 2) - &
            PREVIOUS(S, 3)) / ANGLE_STEP + (4.0 * PRESSURE(S - 1, 1) - &
            PRESSURE(S - 2, 1)) / 2.0 / RADIUS_STEP) / (3.0 / 2.0 / RADIUS_STEP + &
            3.0 / RADIUS(S, 1) / ANGLE_STEP)

        DO I = 2, NI - 1
            DO J = 2, NJ - 1
                IF (I == S) CYCLE
                PRESSURE(I, J) = SQRT(((RADIUS(I + 1, J) + RADIUS(I, J)) / 2.0 * &
                PREVIOUS(I + 1, J) ** 2 / RADIUS_STEP ** 2 + (RADIUS(I, J) + &
                RADIUS(I - 1, J)) / 2.0 * PRESSURE(I - 1, J) ** 2 / &
                RADIUS_STEP ** 2 + 1.0 / RADIUS(I, J) * (PREVIOUS(I, J + 1) ** 2 + &
                PRESSURE(I, J - 1) ** 2) / ANGLE_STEP ** 2) / ((RADIUS(I + 1, J) + &
                RADIUS(I, J)) / 2.0 / RADIUS_STEP ** 2 + (RADIUS(I, J) + &
                RADIUS(I - 1, J)) / 2.0 / RADIUS_STEP ** 2 + 2.0 / RADIUS(I, J) / &
                ANGLE_STEP ** 2))
            END DO
        END DO

    END SUBROUTINE DOUBLE_ACCURACY

END MODULE APPROXIMATION
