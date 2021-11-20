1

   Abaqus 3DEXPERIENCE R2016x                                  Date 20-Nov-2021   Time 12:01:41
   For use by Supplied by Team-SolidSQUAD under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    1     INCREMENT     1     STEP TIME    0.00    


                        S T E P       1     C A L C U L A T I O N   O F   E I G E N V A L U E S

                                            F O R   N A T U R A L   F R E Q U E N C I E S



                                                                                          

     THE LANCZOS EIGENSOLVER IS USED FOR THIS ANALYSIS
     Abaqus WILL COMPUTE UNCOUPLED
      STRUCTURAL AND ACOUSTIC MODES
     NUMBER OF EIGENVALUES                     30
     HIGHEST FREQUENCY OF INTEREST                  1.00000E+18
     MAXIMUM NUMBER OF STEPS WITHIN RUN               35
     BLOCK SIZE FOR LANCZOS PROCEDURE                  7
     THE EIGENVECTORS ARE SCALED SO THAT 
      THE LARGEST DISPLACEMENT ENTRY IN EACH VECTOR
      IS UNITY

     THIS IS A LINEAR PERTURBATION STEP.
     ALL LOADS ARE DEFINED AS CHANGE IN LOAD TO THE REFERENCE STATE

     EXTRAPOLATION WILL NOT BE USED

     CHARACTERISTIC ELEMENT LENGTH     0.500    

     DETAILS REGARDING ACTUAL SOLUTION WAVEFRONT REQUESTED

     DETAILED OUTPUT OF DIAGNOSTICS TO DATABASE REQUESTED

     PRINT OF INCREMENT NUMBER, TIME, ETC., TO THE MESSAGE FILE EVERY     1  INCREMENTS

     COLLECTING MODEL CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     COLLECTING STEP CONSTRAINT INFORMATION FOR OVERCONSTRAINT CHECKS

     CHECK POINT   START OF SOLVER    
  
 COMPUTER PRECISION USED BY ABAQUS      2.22045e-016
 CONVERGENCE CRITERION FOR LANCZOS RUN  1.81899e-012
  
 NUMBER OF EQUATIONS                                                        2028
 NUMBER OF UNRESTRAINED DEGREES OF FREEDOM                                  2028
 HEURISTIC SCALE FOR THE FIRST NONZERO EIGENVALUES (PROBLEM SCALE)  1.14767e-002
 INITIAL SHIFT VALUE                                                0.00000e+000
  
 NUMBER OF CPUS USED BY DIRECT SOLVER     1
  
 NUMBER OF NEGATIVE EIGENVALUES     0
 NUMBER OF POSITIVE EIGENVALUES  2028
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     0
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NUMBER OF SINGULARITIES            3
  
 NUMBER OF ZERO PIVOTS              0
 UPDATED SHIFT VALUE    1.14767e-002
  
 NUMBER OF CPUS USED BY DIRECT SOLVER     1
  
 NUMBER OF NEGATIVE EIGENVALUES     3
 NUMBER OF POSITIVE EIGENVALUES  2025
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     3
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 COMPUTATIONAL INTERVAL FOR THE INITIAL LANCZOS RUN
 
 EIGENVALUES NEEDED TO THE RIGHT OF SHIFT VALUE   30
 EIGENVALUES NEEDED TO THE LEFT OF SHIFT VALUE     3
 RIGHT BOUNDARY OF THE COMPUTATIONAL INTERVAL  3.95574e+037
 LEFT BOUNDARY OF THE COMPUTATIONAL INTERVAL  -3.95574e+037
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   15
 TERMINATION BASED ON THE COST-PER-EIGENVALUE CRITERION
  
 NUMBER OF ACCEPTED EIGENVALUES    8
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    1
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |      3 |        0        3 |     0     0 |     0     0 |      0    30       |
 |-3.95574e+037 1.14767e-002 |-3.95574e+037 1.14767e-002 |-3.95574e+037 3.95574e+037 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   27
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    5
 SHIFT RIGHT        7.34598e+002
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         1.14767e-002
  
 SHIFTING TO THE RIGHT DIRECTION
 NEW SHIFT          7.34598e+002
  
 NUMBER OF NEGATIVE EIGENVALUES   134
 NUMBER OF POSITIVE EIGENVALUES  1894
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS   134
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: EXPAND THE TRUST INTERVAL TO THE RIGHT
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT    126
 RIGHT BOUNDARY OF THE INTERVAL  7.34598e+002
 LEFT BOUNDARY OF THE INTERVAL   3.69150e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |      3 |        0        3 |     0     5 |     0     5 |      0    30       |
 |-3.95574e+037 1.14767e-002 |-3.95574e+037 3.69150e+002 |-3.95574e+037 3.95574e+037 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     0     0 |     0     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 7.34598e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   15
 TERMINATION BASED ON THE COST-PER-EIGENVALUE CRITERION
  
 NUMBER OF ACCEPTED EIGENVALUES   15
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 6.47845e+002 6.49068e+002 6.61119e+002 6.98545e+002
  7.00765e+002 7.10839e+002 7.10839e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  +  |      3 |        0        3 |     0     5 |     0     5 |      0    30       |
 |-3.95574e+037 1.14767e-002 |-3.95574e+037 3.69150e+002 |-3.95574e+037 3.95574e+037 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     0     0 |     0     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 7.34598e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   27
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT        7.34598e+002
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     7
 SHIFT LEFT         1.14767e-002
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          5.08498e+002
  
 NUMBER OF NEGATIVE EIGENVALUES     8
 NUMBER OF POSITIVE EIGENVALUES  2020
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS     8
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT   119
 EIGENVALUES NEEDED TO THE LEFT      0
 RIGHT BOUNDARY OF THE INTERVAL  6.47845e+002
 LEFT BOUNDARY OF THE INTERVAL   5.08498e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0     0 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   21
 TERMINATION BASED ON THE COST-PER-EIGENVALUE CRITERION
  
 NUMBER OF ACCEPTED EIGENVALUES   65
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  5.55943e+002 5.57447e+002 5.57457e+002 5.58383e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.65227e+002 5.72835e+002
  5.73074e+002 5.73121e+002 5.73191e+002 5.73293e+002 5.73515e+002 5.73575e+002
  5.73627e+002 5.82223e+002 5.83745e+002 5.86022e+002 5.88650e+002 5.88650e+002
  5.89297e+002 5.98765e+002 6.10020e+002 6.24455e+002 6.47845e+002 6.49068e+002
  6.61119e+002 6.98545e+002 7.00765e+002 7.10839e+002 7.10839e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    2
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0     0 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   22
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT   25
 SHIFT RIGHT        1.24308e+003
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT         5.08498e+002
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          6.02651e+002
  
 NUMBER OF NEGATIVE EIGENVALUES   123
 NUMBER OF POSITIVE EIGENVALUES  1905
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS   123
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     2
 EIGENVALUES NEEDED TO THE LEFT     67
 RIGHT BOUNDARY OF THE INTERVAL  6.47845e+002
 LEFT BOUNDARY OF THE INTERVAL   5.57457e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0    25 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.57457e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    2
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES   65
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  5.55943e+002 5.57447e+002 5.57457e+002 5.58383e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.65227e+002 5.72835e+002
  5.73074e+002 5.73121e+002 5.73191e+002 5.73293e+002 5.73515e+002 5.73575e+002
  5.73627e+002 5.82223e+002 5.83745e+002 5.86022e+002 5.88650e+002 5.88650e+002
  5.89297e+002 5.98765e+002 6.10020e+002 6.24455e+002 6.47845e+002 6.49068e+002
  6.61119e+002 6.98545e+002 7.00765e+002 7.10839e+002 7.10839e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    3
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0    25 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.57457e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   22
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT        4.40223e+003
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT        -3.69139e+003
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          5.80054e+002
  
 NUMBER OF NEGATIVE EIGENVALUES   115
 NUMBER OF POSITIVE EIGENVALUES  1913
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS   115
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     1
 EIGENVALUES NEEDED TO THE LEFT     66
 RIGHT BOUNDARY OF THE INTERVAL  6.02651e+002
 LEFT BOUNDARY OF THE INTERVAL   5.57457e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0    25 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.57457e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     0     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    2
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES   67
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  5.55943e+002 5.57447e+002 5.57457e+002 5.58383e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.65227e+002 5.71012e+002
  5.72835e+002 5.73074e+002 5.73121e+002 5.73191e+002 5.73293e+002 5.73327e+002
  5.73515e+002 5.73575e+002 5.73627e+002 5.82223e+002 5.83745e+002 5.86022e+002
  5.88650e+002 5.88650e+002 5.89297e+002 5.98765e+002 6.10020e+002 6.24455e+002
  6.47845e+002 6.49068e+002 6.61119e+002 6.98545e+002 7.00765e+002 7.10839e+002
  7.10839e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0    25 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.57457e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     0     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   22
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT        4.40223e+003
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT        -3.69139e+003
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          5.68755e+002
  
 NUMBER OF NEGATIVE EIGENVALUES    95
 NUMBER OF POSITIVE EIGENVALUES  1933
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    95
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT    10
 EIGENVALUES NEEDED TO THE LEFT     54
 RIGHT BOUNDARY OF THE INTERVAL  5.80054e+002
 LEFT BOUNDARY OF THE INTERVAL   5.57457e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0    25 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.57457e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       95       95 |     0     0 |     0     0 |     54    10       |
 | 5.68755e+002 5.68755e+002 | 5.68755e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    1
 SINGULAR BLOCK BETA
  
 NUMBER OF ACCEPTED EIGENVALUES   68
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  5.55943e+002 5.57447e+002 5.57457e+002 5.58383e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.65227e+002 5.71012e+002
  5.72835e+002 5.73062e+002 5.73074e+002 5.73121e+002 5.73191e+002 5.73293e+002
  5.73327e+002 5.73515e+002 5.73575e+002 5.73627e+002 5.82223e+002 5.83745e+002
  5.86022e+002 5.88650e+002 5.88650e+002 5.89297e+002 5.98765e+002 6.10020e+002
  6.24455e+002 6.47845e+002 6.49068e+002 6.61119e+002 6.98545e+002 7.00765e+002
  7.10839e+002 7.10839e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0    25 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.57457e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |       95       95 |     0     0 |     0     0 |     54    10       |
 | 5.68755e+002 5.68755e+002 | 5.68755e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   22
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT        4.40223e+003
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     0
 SHIFT LEFT        -3.69139e+003
  
 SHIFT IS DETERMINED IN ORDER TO CLOSE GAP BETWEEN TRUST INTERVALS
 NEW SHIFT          5.63106e+002
  
 NUMBER OF NEGATIVE EIGENVALUES    94
 NUMBER OF POSITIVE EIGENVALUES  1934
 NUMBER OF NEGATIVE EIGENVALUES NOT ASSOCIATED WITH LAGRANGE MULTIPLIERS    94
 NUMBER OF LAGRANGE MULTIPLIERS     0
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     54
 RIGHT BOUNDARY OF THE INTERVAL  5.63106e+002
 LEFT BOUNDARY OF THE INTERVAL   5.57457e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0    25 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.57457e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |     0     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:   18
 TERMINATION BASED ON THE COST-PER-EIGENVALUE CRITERION
  
 NUMBER OF ACCEPTED EIGENVALUES   83
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  5.55943e+002 5.57447e+002 5.57457e+002 5.58383e+002 5.58383e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.65227e+002 5.71012e+002 5.72835e+002 5.73062e+002 5.73074e+002
  5.73121e+002 5.73191e+002 5.73293e+002 5.73327e+002 5.73515e+002 5.73575e+002
  5.73627e+002 5.82223e+002 5.83745e+002 5.86022e+002 5.88650e+002 5.88650e+002
  5.89297e+002 5.98765e+002 6.10020e+002 6.24455e+002 6.47845e+002 6.49068e+002
  6.61119e+002 6.98545e+002 7.00765e+002 7.10839e+002 7.10839e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0    25 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.57457e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |     0     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   22
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT        5.63106e+002
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT    15
 SHIFT LEFT        -1.71480e+002
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT          5.63106e+002
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     39
 RIGHT BOUNDARY OF THE INTERVAL  5.63106e+002
 LEFT BOUNDARY OF THE INTERVAL   5.08498e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |    15     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    9
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   90
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  5.55943e+002 5.57447e+002 5.57457e+002 5.58383e+002 5.58383e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.65227e+002 5.71012e+002 5.72835e+002 5.73062e+002
  5.73074e+002 5.73121e+002 5.73191e+002 5.73293e+002 5.73327e+002 5.73515e+002
  5.73575e+002 5.73627e+002 5.82223e+002 5.83745e+002 5.86022e+002 5.88650e+002
  5.88650e+002 5.89297e+002 5.98765e+002 6.10020e+002 6.24455e+002 6.47845e+002
  6.49068e+002 6.61119e+002 6.98545e+002 7.00765e+002 7.10839e+002 7.10839e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |    15     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   22
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT        5.63106e+002
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     7
 SHIFT LEFT        -1.71480e+002
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT          5.63106e+002
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     32
 RIGHT BOUNDARY OF THE INTERVAL  5.63106e+002
 LEFT BOUNDARY OF THE INTERVAL   5.08498e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |     7     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    9
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES   97
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  5.55943e+002 5.57447e+002 5.57457e+002 5.58383e+002 5.58383e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.65227e+002 5.71012e+002 5.72835e+002
  5.73062e+002 5.73074e+002 5.73121e+002 5.73191e+002 5.73293e+002 5.73327e+002
  5.73515e+002 5.73575e+002 5.73627e+002 5.82223e+002 5.83745e+002 5.86022e+002
  5.88650e+002 5.88650e+002 5.89297e+002 5.98765e+002 6.10020e+002 6.24455e+002
  6.47845e+002 6.49068e+002 6.61119e+002 6.98545e+002 7.00765e+002 7.10839e+002
  7.10839e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |     7     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   22
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT        5.63106e+002
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     7
 SHIFT LEFT        -1.71480e+002
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT          5.63106e+002
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     25
 RIGHT BOUNDARY OF THE INTERVAL  5.63106e+002
 LEFT BOUNDARY OF THE INTERVAL   5.08498e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |     7     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    9
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES  104
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  5.55943e+002 5.57447e+002 5.57457e+002 5.58383e+002 5.58383e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.65227e+002 5.71012e+002
  5.72835e+002 5.73062e+002 5.73074e+002 5.73121e+002 5.73191e+002 5.73293e+002
  5.73327e+002 5.73515e+002 5.73575e+002 5.73627e+002 5.82223e+002 5.83745e+002
  5.86022e+002 5.88650e+002 5.88650e+002 5.89297e+002 5.98765e+002 6.10020e+002
  6.24455e+002 6.47845e+002 6.49068e+002 6.61119e+002 6.98545e+002 7.00765e+002
  7.10839e+002 7.10839e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |     7     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   22
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT        5.63106e+002
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     7
 SHIFT LEFT        -1.71480e+002
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT          5.63106e+002
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     18
 RIGHT BOUNDARY OF THE INTERVAL  5.63106e+002
 LEFT BOUNDARY OF THE INTERVAL   5.08498e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |     7     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    9
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES  111
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  5.55943e+002 5.57447e+002 5.57457e+002 5.58383e+002 5.58383e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.65227e+002
  5.71012e+002 5.72835e+002 5.73062e+002 5.73074e+002 5.73121e+002 5.73191e+002
  5.73293e+002 5.73327e+002 5.73515e+002 5.73575e+002 5.73627e+002 5.82223e+002
  5.83745e+002 5.86022e+002 5.88650e+002 5.88650e+002 5.89297e+002 5.98765e+002
  6.10020e+002 6.24455e+002 6.47845e+002 6.49068e+002 6.61119e+002 6.98545e+002
  7.00765e+002 7.10839e+002 7.10839e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |     7     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   22
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT        5.63106e+002
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     7
 SHIFT LEFT        -1.71480e+002
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT          5.63106e+002
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT     11
 RIGHT BOUNDARY OF THE INTERVAL  5.63106e+002
 LEFT BOUNDARY OF THE INTERVAL   5.08498e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |     7     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    9
 ALL THE EIGENVALUES ARE INSIDE THE INTERVAL
  
 NUMBER OF ACCEPTED EIGENVALUES  118
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  5.55943e+002 5.57447e+002 5.57457e+002 5.58383e+002 5.58383e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.65227e+002 5.71012e+002 5.72835e+002 5.73062e+002 5.73074e+002 5.73121e+002
  5.73191e+002 5.73293e+002 5.73327e+002 5.73515e+002 5.73575e+002 5.73627e+002
  5.82223e+002 5.83745e+002 5.86022e+002 5.88650e+002 5.88650e+002 5.89297e+002
  5.98765e+002 6.10020e+002 6.24455e+002 6.47845e+002 6.49068e+002 6.61119e+002
  6.98545e+002 7.00765e+002 7.10839e+002 7.10839e+002
  
 FACTORIZATION WITH NEW SHIFT VALUE IS REQUIRED
  
 TRUST INTERVALS AFTER LANCZOS RUN AND COMPRESSION
  
 NUMBER OF TRUST INTERVALS    5
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |     7     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 DATA FOR CHOOSING SHIFT
  
 NUMBER OF EIGENVALUES NEEDED IN THE PRIMARY TRUST INTERVAL   22
  
 NUMBER OF MODES FORMING SENTINEL TO THE RIGHT    0
 SHIFT RIGHT        5.63106e+002
  
 NUMBER OF MODES FORMING SENTINEL TO THE LEFT     7
 SHIFT LEFT        -1.71480e+002
  
 RESTART AT THE SAME SHIFT
 NEW SHIFT          5.63106e+002
  
 NEW LANCZOS RUN IS REQUIRED
  
 TYPE OF RUN: FILLING GAPS
 EIGENVALUES NEEDED TO THE RIGHT     0
 EIGENVALUES NEEDED TO THE LEFT      4
 RIGHT BOUNDARY OF THE INTERVAL  5.63106e+002
 LEFT BOUNDARY OF THE INTERVAL   5.08498e+002
  
 TRUST INTERVALS AFTER FACTORIZATION AND COMPRESSION
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  +  |      8 |        0        8 |     0     0 |     0    50 |      0   119       |
 |-3.95574e+037 5.08498e+002 |-3.95574e+037 5.08498e+002 |-3.95574e+037 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      1 |       94       95 |     0     0 |     7     1 |     54    10       |
 | 5.63106e+002 5.68755e+002 | 5.63106e+002 5.68755e+002 | 5.57457e+002 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 LANCZOS RUN TERMINATED. NUMBER OF STEPS:    9
 ALL THE EIGENVALUES REQUESTED ON THIS LANCZOS RUN ARE FOUND
  
 NUMBER OF ACCEPTED EIGENVALUES  122
  
 ACCEPTED EIGENVALUES
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  5.55943e+002 5.57447e+002 5.57457e+002 5.58383e+002 5.58383e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002
  5.58696e+002 5.58696e+002 5.58696e+002 5.58696e+002 5.65227e+002 5.71012e+002
  5.72835e+002 5.73062e+002 5.73074e+002 5.73121e+002 5.73191e+002 5.73293e+002
  5.73327e+002 5.73515e+002 5.73575e+002 5.73627e+002 5.82223e+002 5.83745e+002
  5.86022e+002 5.88650e+002 5.88650e+002 5.89297e+002 5.98765e+002 6.10020e+002
  6.24455e+002 6.47845e+002 6.49068e+002 6.61119e+002 6.98545e+002 7.00765e+002
  7.10839e+002 7.10839e+002
  
 THE REQUIRED NUMBER OF EIGENVALUES IS FOUND
  
 FINAL TRUST INTERVALS AFTER SUCCESSFUL LANCZOS RUN
  
 NUMBER OF TRUST INTERVALS    4
  
 -------------------------------------------------------------------------------------
 | KEY | INSIDE |      INERTIAS     |  SENTINELS  |   OUTSIDE   |  NEED TO FIND      |
 |      ENDS OF INTERVAL     |     ENDS OF SENTINEL      |   COMPUTATIONAL INTERVAL  |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      134      134 |     7     0 |     7     0 |    126     0       |
 | 7.34598e+002 7.34598e+002 | 6.47845e+002 7.34598e+002 | 3.69150e+002 7.34598e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      123      123 |     0     0 |     0     0 |     67     2       |
 | 6.02651e+002 6.02651e+002 | 6.02651e+002 6.02651e+002 | 5.57457e+002 6.47845e+002 |
 -------------------------------------------------------------------------------------
 |  -  |      0 |      115      115 |     0     0 |     2     0 |     66     1       |
 | 5.80054e+002 5.80054e+002 | 5.80054e+002 5.80054e+002 | 5.57457e+002 6.02651e+002 |
 -------------------------------------------------------------------------------------
 |  +  |     95 |        0       95 |     0     0 |     0     1 |      0    10       |
 |-3.95574e+037 5.68755e+002 |-3.95574e+037 5.68755e+002 |-3.95574e+037 5.80054e+002 |
 -------------------------------------------------------------------------------------
  
 FIRST EIGENVALUE REQUESTED BY THE USER     1
  
 EIGENVALUES REQUESTED BY THE USER
  5.19531e-012 6.07829e-012 6.90491e-012 1.45098e+001 1.55251e+001 1.55251e+001
  1.84571e+002 3.69150e+002 5.38555e+002 5.38565e+002 5.38964e+002 5.38964e+002
  5.40009e+002 5.40013e+002 5.41668e+002 5.41668e+002 5.43812e+002 5.43815e+002
  5.46298e+002 5.46298e+002 5.47188e+002 5.48937e+002 5.48938e+002 5.50065e+002
  5.51544e+002 5.51544e+002 5.53932e+002 5.53933e+002 5.55820e+002 5.55943e+002
  
 EIGENVALUE EXTRACTION USING THE LANCZOS EIGENSOLVER IS COMPLETE
 NUMBER OF LANCZOS RUNS       :  13
 NUMBER OF LANCZOS STEPS      : 128
 NUMBER OF FACTORIZATIONS     :   8
 NUMBER OF BACK SUBSTITUTIONS : 136

     CHECK POINT  END OF SOLVER       

       ELAPSED USER TIME (SEC)      =   1.0000    
       ELAPSED SYSTEM TIME (SEC)    =  0.20000    
       ELAPSED TOTAL CPU TIME (SEC) =   1.2000    
       ELAPSED WALLCLOCK TIME (SEC) =          1
1

   Abaqus 3DEXPERIENCE R2016x                                  Date 20-Nov-2021   Time 12:01:43
   For use by Supplied by Team-SolidSQUAD under license from Dassault Systemes or its subsidiary.

                                                                                 
                                                                                 
 STEP    2     INCREMENT     1     STEP TIME    0.00    
 


                        S T E P       2     M A T R I X   G E N E R A T E
 
  THE FOLLOWING GLOBAL MATRICES WILL BE GENERATED AND WRITTEN TO Job-1-1_X2.sim

    STIFFNESS
    MASS
 
  IF MULTY-POINT CONSTRAINTS ARE SPECIFIED IN THE MODEL, THEY WILL BE APPLIED
 
  THE FOLLOWING GLOBAL MATRICES WILL BE WRITTEN TO THE FOLLOWING TEXT FILES

    STIFFNESS          TO Job-1-1_STIF2.mtx
    MASS               TO Job-1-1_MASS2.mtx
   MATRIX OUTPUT FORMAT - MATRIX INPUT


          THE ANALYSIS HAS BEEN COMPLETED



     ANALYSIS SUMMARY:
     TOTAL OF          2  INCREMENTS
                       0  CUTBACKS IN AUTOMATIC INCREMENTATION
                       0  ITERATIONS INCLUDING CONTACT ITERATIONS IF PRESENT
                     144  PASSES THROUGH THE EQUATION SOLVER OF WHICH 
                       8  INVOLVE MATRIX DECOMPOSITION, INCLUDING
                       0  DECOMPOSITION(S) OF THE MASS MATRIX
                       1  REORDERING OF EQUATIONS TO MINIMIZE WAVEFRONT
                       0  ADDITIONAL RESIDUAL EVALUATIONS FOR LINE SEARCHES
                       0  ADDITIONAL OPERATOR EVALUATIONS FOR LINE SEARCHES
                       0  WARNING MESSAGES DURING USER INPUT PROCESSING
                       0  WARNING MESSAGES DURING ANALYSIS
                       0  ANALYSIS WARNINGS ARE NUMERICAL PROBLEM MESSAGES
                       0  ANALYSIS WARNINGS ARE NEGATIVE EIGENVALUE MESSAGES
                       0  ERROR MESSAGES



     JOB TIME SUMMARY
       USER TIME (SEC)      =   1.3000    
       SYSTEM TIME (SEC)    =  0.30000    
       TOTAL CPU TIME (SEC) =   1.6000    
       WALLCLOCK TIME (SEC) =          2
