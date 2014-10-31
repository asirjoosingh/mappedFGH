      SUBROUTINE PRIMAT (MATRIX,DIM,M,N,VORKO,NACHKO,TAPE,TEXT)
C
C     Zweck: PRIMAT (PRInt MATrix) gibt die in MATRIX stehende Matrix
C            auf die durch TAPE spezifizierte Ausgabeeinheit aus.
C     Autor: Gerald Geudtner (Maerz 1990,
C            Matthias Krack (September 1990)
C
C     MATRIX: Auszugebende Real-Matrix -> MATRIX(DIM,DIM)
C     DIM   : Dimension von MATRIX in der Routine, die PRIMAT aufruft
C     M     : Anzahl der auszugebenen Zeilen von MATRIX (M.LE.DIM)
C     N     : Anzahl der auszugebenen Spalten von MATRIX (N.LE.DIM)
C     VORKO : Anzahl der auszugebenen Vorkommastellen (VORKO.GE.1)
C     NACHKO: Anzahl der auszugebenen Nachkommastellen (NACHKO.GE.0)
C     TAPE  : Nummer der Ausgabeeinheit (TAPE.GE.0)
C             Bei TAPE = 0 wird auf die Standardausgabeeinheit (*)
C             geschrieben.
C     TEXT  : Text fuer die Ueberschrift der Matrix
C
C     Datum der letzten Aenderung: 17.11.1990  (MK)
C
C     ******************************************************************
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER BIS,BREITE,DIM,LINKS,M,N,NACHKO,RECHTS,SPAZAL,TAPE,
     +        VON,VORKO
      CHARACTER TEXT*(*),FORM1*27,FORM2*29,FORM3*28
      REAL*8 MATRIX(DIM,DIM)
C
C     ------------------------------------------------------------------
C
C     *** Definition der variablen Formate ***
C
      FORM1 = '(/, T2, 7X, ??(2X, F??.??))'
C     FORM2 = '(/, T2, 7X, ??(??X, I5, ??X))'
      FORM2 = '(   T2, 7X, ??(??X, I5, ??X))'
      FORM3 = '(T2, I5, 2X, ??(2X, F??.??))'
C
C     *** Ausgabe der Ueberschrift ***
C
      IF (TAPE.EQ.0) THEN
        WRITE (*, 5000) TEXT
 5000   FORMAT (/, T2, A)
      ELSE IF ((TAPE.GT.0).AND.(TAPE.LT.1000)) THEN
        WRITE (TAPE, 5000) TEXT
      ELSE
        WRITE (*, 5010)
 5010   FORMAT (/, T2, 'TAPE = ', I5, /,
     +          /, T2, '*** FEHLER! Unzulaessige Bezifferung der ',
     +                 'Ausgabeeinheit ***')
        GO TO 999
      ENDIF
C
C     *** Ueberpruefung der Parameter M, N und VORKO ***
C
      IF ((M.LT.1).OR.(M.GT.DIM)) THEN
        IF (TAPE.EQ.0) THEN
          WRITE (*, 5020) M
 5020     FORMAT (/, T2, 'M = ', I5, /,
     +            /, T2, '*** FEHLER! Unzulaessiger Wert fuer die ',
     +                   'Zeilenzahl ***')
        ELSE
          WRITE (TAPE, 5020) M
        END IF
        GO TO 999
      END IF
C
      IF ((N.LT.1).OR.(N.GT.DIM)) THEN
        IF (TAPE.EQ.0) THEN
          WRITE (*, 5030) N
 5030     FORMAT (/, T2, 'N = ', I5, /,
     +            /, T2, '*** FEHLER! Unzulaessiger Wert fuer die ',
     +                   'Spaltenzahl ***')
        ELSE
          WRITE (TAPE, 5030) N
        END IF
        GO TO 999
      END IF
C
      IF ((VORKO.LT.1).OR.(VORKO.GT.15)) THEN
        IF (TAPE.EQ.0) THEN
          WRITE (*, 5040) VORKO
 5040     FORMAT (/, T2, 'VORKO = ', I5, /,
     +            /, T2, '*** FEHLER! Unzulaessiger Wert fuer die ',
     +                   'Anzahl der Vorkommastellen ***')
        ELSE
          WRITE (TAPE, 5040) VORKO
        END IF
        GO TO 999
      END IF
C
C     *** Beschreiben der variablen Formate ***
C
      BREITE = VORKO + NACHKO + 4
      SPAZAL = INT(124/BREITE)
      WRITE (FORM2(13:14), 5050) SPAZAL
 5050 FORMAT (I2)
      RECHTS = MAX((NACHKO - 2),1)
      LINKS =  BREITE - RECHTS - 5
      WRITE (FORM2(16:17), 5050) LINKS
      WRITE (FORM2(25:26), 5050) RECHTS
      IF ((NACHKO.GE.0).AND.(NACHKO.LT.16)) THEN
        WRITE (FORM1(13:14), 5050) SPAZAL
        WRITE (FORM3(14:15), 5050) SPAZAL
        WRITE (FORM1(21:22), 5050) BREITE - 2
        WRITE (FORM3(22:23), 5050) BREITE - 2
        WRITE (FORM1(24:25), 5050) NACHKO
        WRITE (FORM3(25:26), 5050) NACHKO
      ELSE
        IF (TAPE.EQ.0) THEN
          WRITE (*, 5060) NACHKO
 5060     FORMAT (/, T2, 'NACHKO = ', I5, /,
     +            /, T2, '*** FEHLER! Unzulaessiger Wert fuer die ',
     +                   ' Anzahl der Nachkommastellen ***')
        ELSE
          WRITE (TAPE, 5060) NACHKO
        END IF
        GO TO 999
      END IF
C
C     *** Ausgabe der Matrix MATRIX im Format M x N ***
C
      DO 20  I = 1, N, SPAZAL
        VON = I
        BIS = MIN((VON + SPAZAL - 1),N)
        IF (TAPE.EQ.0) THEN
          WRITE (*, FORM2) (J, J = VON, BIS)
        ELSE
          WRITE (TAPE, FORM2) (J, J = VON, BIS)
        ENDIF
        DO 10  J = 1, M
          IF (TAPE.EQ.0) THEN
            WRITE (*, FORM3) J, (MATRIX(J,K), K = VON, BIS)
          ELSE
            WRITE (TAPE, FORM3) J, (MATRIX(J,K), K = VON, BIS)
          ENDIF
   10   CONTINUE

CC        IF (TAPE.EQ.0) THEN
CC          WRITE (*, '(/)')
CC        ELSE
CC          WRITE (TAPE, '(/)')
CC        ENDIF

   20 CONTINUE
C
  999 CONTINUE
      END
