      SUBROUTINE RPF501( NCHTDS, CHTDST, TIR_ID, NTYPAR, TYPARR )

c  Copyright (C) 2000-1999
c  By Mechanical Dynamics, Inc. Ann Arbor, Michigan
c
c  All Rights Reserved, This code may not be copied or
c  reproduced in any form, in part or in whole, without the
c  explicit written permission of Mechanical Dynamics, Inc.
c

c DESCRIPTION:
c
c  Reads property file for user tire model based on
c  the Fiala Tire model and initializes the
c  tire parameter array (TYPARR).
c

c ARGUMENT LIST:
c
c    name    type  storage    use           description
c    ======  ====  =========  ===  ====================================
c    NCHTDS   I.S.     -       R   Number of characters in tire 
c                                  property file name.
c    CHTDST   C.A.             R   Tire property file name
c    TIR_ID   I.S.     1       R   Tire GFORCE id
c    NTYPAR   I.S.     1       R   Dimension of TYPARR
c    TYPARR   D.A.   NTYPAR    E   Tire parameter array
c
c   *** Legend:  I integer             S scalar     R referenced
c                D double precision    A array      E evaluated
c                C character
  
C Inputs:

      INTEGER          NCHTDS
      CHARACTER*(*)    CHTDST
      INTEGER          TIR_ID    
      INTEGER          NTYPAR
      
C Outputs:

      DOUBLE PRECISION TYPARR( NTYPAR )
      
C Locals:

C Units conversions:

      CHARACTER*(12)   UNITS(5)
      DOUBLE PRECISION CV2MDL(5)
      DOUBLE PRECISION CV2SI(5)
      DOUBLE PRECISION FCVSI
      DOUBLE PRECISION LCVSI
      DOUBLE PRECISION MCVSI
      DOUBLE PRECISION ACVSI
      DOUBLE PRECISION TCVSI

c Fiala Property File Map

      INCLUDE 'tyrHarty_501.inc'
      
C  RTO variables:

      INTEGER          RETURN_VAL
      DOUBLE PRECISION TMPREAL

C Shape Array RTO Stuff:

      DOUBLE PRECISION TMP1, TMP2
      INTEGER          N_NODES
      INTEGER          ARRPTR
      CHARACTER*80     FORM
      INTEGER          FLEN
      CHARACTER*80     TABLE
      INTEGER          TLEN
      
      LOGICAL          ERRFLG 
      CHARACTER*80     MESSAG
            
c+---------------------------------------------------------------------*
c
c Open the file:

      CALL RTO_OPEN_FILE_F2C ( CHTDST, NCHTDS, RETURN_VAL ) 

      ERRFLG = RETURN_VAL .EQ. 0 
      MESSAG = 'Harty Tyre 501:  No Error opening tire property file.'
      CALL ERRMES ( ERRFLG, MESSAG, TIR_ID, 'STOP' )

      
c Read [UNITS] block from property file:

c  Parameters in the property file may be given in any consistent
c  set of units.  The [UNITS] block identifies those units.
c  During evaluation, however, SI Units are used. So as parameters 
c  are read from the property file they are converted
c  to SI units.

c  SI unit system.
c    LENGTH = meter
c    FORCE  = newton
c    ANGLE  = radians
c    MASS   = kg
c    TIME   = second

c  UNITS(1)-> FORCE  UNITS(2)-> MASS  UNITS(3)-> LENGTH 
c  UNITS(4)-> TIME   UNITS(5)-> ANGLE

      CALL ATRTOU( TIR_ID, UNITS )
      CALL ACUNFN( UNITS, CV2MDL, CV2SI )  
            
      FCVSI = CV2SI(1)			
C   Force Conversion
      MCVSI = CV2SI(2)			
C   Mass Conversion
      LCVSI = CV2SI(3)			
C   Length Conversion
      TCVSI = CV2SI(4)			
C   Time Conversion
      ACVSI = CV2SI(5)			
C   Angle Conversion
      
C**************************  TIRPRP POPULATION  ************************

c Read [MODEL] block:
 
      CALL RTO_READ_REAL_F2C
     . ( 
     .  'MODEL', 5, 'USE_MODE', 8, 
     .   TYPARR( USE_MODE ), RETURN_VAL 
     . )  
         
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG ,
     .  'Harty Tyre 501:  No Use_mode?' 
     .   ,TIR_ID,'STOP')


c Read [DIMENSION] block:

      CALL RTO_READ_REAL_F2C
     . (
     .  'DIMENSION', 9, 'UNLOADED_RADIUS', 15,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG, 
     .  'Harty Tyre 501:  No UNLOADED_RADIUS?' 
     .   ,TIR_ID,'STOP')
     
      TYPARR( UNLOADED_RADIUS ) = TMPREAL  *  LCVSI
      
      CALL RTO_READ_REAL_F2C
     . (
     .  'DIMENSION', 9, 'WIDTH', 5,
     .   TMPREAL, RETURN_VAL 
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG, 'Harty Tyre 501:  No WIDTH?' 
     .               ,TIR_ID,'STOP')
   
      TYPARR( WIDTH )    = TMPREAL  *  LCVSI
      

c Read [PARAMETER] block

      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'VERTICAL_STIFFNESS', 18,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
C-----------------------------------------------------------------------
      CALL ERRMES( ERRFLG, 'Harty Tyre 501:  No VERTICAL_STIFFNESS?' 
     .            ,TIR_ID,'STOP')
      
      TYPARR( VERTICAL_STIFFNESS )   = TMPREAL * (FCVSI / LCVSI)
      
      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'VERTICAL_DAMPING', 16,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG, 'Harty Tyre 501:  No VERTICAL_DAMPING?' 
     .              ,TIR_ID,'STOP')
      
      TYPARR( VERTICAL_DAMPING ) = TMPREAL * (FCVSI * TCVSI / LCVSI)
      

      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'ROLLING_RESISTANCE', 18,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No ROLLING_RESISTANCE?', 
     .   TIR_ID,'STOP')
      
      TYPARR( ROLLING_RESISTANCE )   = TMPREAL
      
      
      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'CMX', 3,
     .   TMPREAL, RETURN_VAL
     . )

      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'rpf501: CMX undefined.',
     .   TIR_ID,'STOP')

      TYPARR( CMX )   = TMPREAL


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'CGAMMA', 6,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No CGAMMA?', 
     .   TIR_ID,'STOP' )
      
      TYPARR( CGAMMA ) = TMPREAL * (FCVSI / ACVSI)
      
      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'UMIN', 4,
     .  TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No UMIN?', 
     .   TIR_ID,'STOP')
  
      TYPARR( UMIN ) = TMPREAL 
      

      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'UMAX', 4,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No UMAX?', 
     .   TIR_ID,'STOP')
      
      TYPARR( UMAX ) = TMPREAL 


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'RELAXATION_LENGTH', 17,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No RELAXATION_LENGTH?', 
     .   TIR_ID,'STOP')
      
      TYPARR( RELAXATION_LENGTH ) = ABS(TMPREAL * LCVSI)


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'ALPHA_CRITICAL', 14,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No ALPHA_CRITICAL?', 
     .   TIR_ID,'STOP')
      
      TYPARR( ALPHA_CRITICAL ) = ABS(TMPREAL * ACVSI)


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'CURVATURE_FACTOR_ANGLE', 22,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No CURVATURE_FACTOR_ANGLE?', 
     .   TIR_ID,'STOP')
      
      TYPARR( curvature_factor_angle ) = ABS(TMPREAL)


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'SCALE_FACTOR_LATERAL', 20,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No SCALE_FACTOR_LATERAL?', 
     .   TIR_ID,'STOP')
      
      TYPARR( SCALE_FACTOR_LATERAL ) = ABS(TMPREAL)


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'RATED_LOAD', 10,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No RATED_LOAD?', 
     .   TIR_ID,'STOP')
      
      TYPARR( rated_load ) = ABS(TMPREAL * MCVSI)


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'SCALE_FACTOR_DIM', 16,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No SCALE_FACTOR_DIM?', 
     .   TIR_ID,'STOP')
      
      TYPARR( scale_factor_dim ) = ABS(TMPREAL)


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'SLIP_RATIO_CRITICAL', 19,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No SLIP_RATIO_CRITICAL?', 
     .   TIR_ID,'STOP')
      
      TYPARR( slip_ratio_critical ) = ABS(TMPREAL)


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'CURVATURE_FACTOR_RATIO', 22,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No CURVATURE_FACTOR_RATIO?', 
     .   TIR_ID,'STOP')
      
      TYPARR( curvature_factor_ratio ) = ABS(TMPREAL)


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'PNEUM_TRAILING_SCALING', 22,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No PNEUM_TRAILING_SCALING?', 
     .   TIR_ID,'STOP')
      
      TYPARR( pneum_trailing_scaling ) = ABS(TMPREAL)


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'PNEUMATIC_LEAD_CAMBER', 21,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No PNEUMATIC_LEAD_CAMBER?', 
     .   TIR_ID,'STOP')
      
      TYPARR( pneumatic_lead_camber ) = ABS(TMPREAL * LCVSI)


      CALL RTO_READ_REAL_F2C
     . (
     .  'PARAMETER', 9, 'LIMIT_CAMBER_ONSET_FRIC', 23,
     .   TMPREAL, RETURN_VAL
     . )
     
      ERRFLG = RETURN_VAL .EQ. 0
      CALL ERRMES( ERRFLG,
     .  'Harty Tyre 501:  No LIMIT_CAMBER_ONSET_FRIC?', 
     .   TIR_ID,'STOP')
      
      TYPARR( limit_camber_onset_fric ) = ABS(TMPREAL * ACVSI)

      n_nodes  = 0
      arrptr   = shape

C READ [SHAPE] BLOCK IF IT EXISTS:

      CALL RTO_START_TABLE_READ_F2C
     . ( 
     .  'SHAPE', 5, FORM, FLEN, RETURN_VAL
     . )
     
      IF ( RETURN_VAL .EQ. 1 ) THEN
      
800     CONTINUE

        CALL RTO_READ_TABLE_LINE_F2C( TABLE, TLEN, RETURN_VAL )
	
        if ( return_val .eq. 1 .and. tlen .gt. 3 ) then

          call act_line_parse (table, tmp1, tmp2, tlen)

          if ( n_nodes .lt. max_shape .and. tlen .eq. 2 ) then
            n_nodes  = n_nodes  + 1        
            typarr( arrptr ) = tmp1
            typarr( arrptr + 1) = tmp2
            arrptr = arrptr + 2
          else

            if ( n_nodes .gt. max_shape) then
	      CALL ERRMES( .true.,
     .          'Harty Tyre 501:  Shape table has more than 10 nodes', 
     .           TIR_ID, 'STOP' )
            endif

            if (tlen .ne. 2) then
              CALL ERRMES( .true.,
     .             'Harty Tyre 501:  Error parsing line of SHAPE table',
     .              TIR_ID, 'STOP' )
            endif

          endif

          goto 800
	  
        endif
 
        typarr( n_shape ) = n_nodes
 
      else

	call usrmes( .true., 
     .   'Harty Tyre 501:  No shape table.  Cylinder will be used' 
     .   ,tir_id, 'WARN')
        	  
      endif
      
C Close tire property file:

      CALL RTO_CLOSE_FILE_F2C ( CHTDST, NCHTDS, RETURN_VAL )
 
      ERRFLG = RETURN_VAL .EQ. 0 
      MESSAG = 'exa_fiaini: Error closing tire property file.'
      CALL ERRMES( ERRFLG, MESSAG, TIR_ID, 'STOP' )

      RETURN
      END
