# phase-shift-optimization-rspace
Optimizing phase shifts for r-space or coordinate space. Currently working on N2LO for phase shifts.

## Goal as of 03/16/2020

The current goal is to automate Sanjoy's process for calculating phase shifts and low energy parameters.

Currently Sanjoy has to manually edit coefficients, recompile a fortran binary, and then run said fortran binary.

The coefficients are stored in the data file dphqq.d in the following format:

```
c01  1s0     0.126  0.        138.039   0.
fun          3.     11.       1.        0.70
c01  1s0     0.3048 0.        138.039   0.
fun          3.     12.       1.        0.70
c00  1p1    13.319  0.        138.039   0.
fun          3.     11.       1.        0.70
c00  1p1     0.849  0.        138.039   0.
fun          3.     12.       1.        0.70
c10  3s1   -0.5837950.        138.039   0.
fun          3.     11.       1.        0.70
c10  3s1   0.2399   0.        138.039   0.
fun          3.     12.       1.        0.70
s120 3s1   0.00     0.        138.039   0.
fun          3.     14.       1.        0.70
ls0  3s1   0.0      0.        138.039   0.
fun          3.     16.       1.        0.70
c11  3pj   1.28     0.        138.039   0.
fun          3.     11.       1.        0.70
c11  3pj  -0.520    0.        138.039   0.
fun          3.     12.       1.        0.70
ls1  3pj  -0.865    0.        138.039   0.
fun          3.     16.       1.        0.70
s121 3pj   0.00     0.        138.039   0.
fun          3.     14.       1.        0.70
```

where the number 0.126 after c01 1s0 would be the first 1s0 coefficient.

The phase shifts and low energy parameters are calculated and stored in the output file phqq.d in the format:

```
                    d e g r e e s
 elab(mev)    theoretical  experimental  upper   lower               chi**2
 --------------------------------------------------------------------------


  1 s 0       

     0.01        14.51505      0.000    0.000    0.000                 0.00
     0.02        20.03919      0.000    0.000    0.000                 0.00
     0.03        23.99041      0.000    0.000    0.000                 0.00
     1.00        62.00469     62.068   62.098   62.038 n93n            4.45
     5.00        63.55757     63.630   63.710   63.550 n93n            0.82
    10.00        59.86958     59.960   60.070   59.850 n93n            0.68
    25.00        50.75415     50.900   51.090   50.710 n93n            0.59
    50.00        40.27952     40.540   40.820   40.260 n93n            0.87
   100.00        26.40703     26.780   27.160   26.400 n93n            0.96
   150.00        16.80707     16.940   17.350   16.530 n93n            0.11
   200.00         9.51462      8.940    9.330    8.550 n93n            2.17
   250.00         3.73634      1.960    2.330    1.590 n93n           23.05
   300.00        -0.94602     -4.460   -4.030   -4.890 n93n           66.78
 
 low energy parameters    a  =  -23.6688    r  =    2.6848
                                                  total chi**2       100.47
                                                  chi**2/datum        10.05
                                                  for   10 p. of data
```
