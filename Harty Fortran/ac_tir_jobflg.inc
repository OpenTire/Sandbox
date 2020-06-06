c*******************************************************************
c
c  Adams/Tire
c  Copyright (C) 1999-Present  MSC.Software Corporation
c  All Rights Reserved.
c
c*******************************************************************
c
c
c Holds the jobflg types for the Standard Tire Interface (STI):
c
c     -1 = No evaluation
c      0 = Normal Mode
c      1 = Inquire:  Subroutine must return the actually used dimensions of
c                    NTYPAR, NDEQVR, NVARS, NIWORK, NWORK
c      2 = First Initialization (new model)
c      3 = Re-Initialization (mid-simulation)
c      4 = Successful Simulaton Step (handled from sen980)
c      5 = Differencing mode
c      6 = End Simulation (handled from acarSDI_fin)

c jobflg = 0   Normal Mode
c jobflg = 1   Subroutine must return the actually used dimensions of
c jobflg = 2   First initialization before starting the simulation
c jobflg = 3   Re-initialization during simulation
c jobflg = 4   Successful step (not used in A/Car)
c jobflg = 5   DFLAG is .TRUE. (A/Car deviation)
c jobflg = 99  Final call of tire model (not used in A/Car)
c
      integer     noeval
      parameter ( noeval = -1 )
c      
      integer     normal
      parameter ( normal = 0 )
c      
      integer     inquire
      parameter ( inquire = 1 )    
c
      integer     init
      parameter ( init  =  2 )
c
      integer     reset
      parameter ( reset  =  3 )
c
      integer     sstep
      parameter ( sstep  =  4 )
c
      integer     diff
      parameter ( diff  =  5 )
c
      integer     endsim
      parameter ( endsim  =  99 )
