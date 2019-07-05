VIB_LEVELS
==========
Tool to solve the Schroedinger equation of a particle in a
one-dimensional potential V(r) which is known at several values of r
and is affected by uncertainties.  The potential is fitted to a
polynomial within a Monte Carlo resampling scheme, and the first four
levels are obtained.

VIB_LEVELS uses the the one-dimensional Schroedinger equation solver
from ANH_QUADRATURE, and some code from POLYFIT.

VIB_LEVELS is a stand-alone Fortran utility which requires a LAPACK
library.  The code can be compiled with:

```
gfortran -o vib_levels vib_levels.f90 -llapack
```

The compile.sh script is provided for reference only.

The code asks for the name of the data file containing the {r,V,dV}
data and the desired expansion order for the polynomial representation
of the potential on standard input.  Optionally, one can provide an
effective mass in Dalton [=a.m.u.; 1 a.u. by default], and a plot file
name root to produce a plot of the data, the fit and the relevant
eigenvalues and eigenfunctions.

Note: for production calculations the hard-coded number of Monte Carlo
samples (nsample) should be increased to 5000 to guarantee an
uncertainty in the uncertainty of less than 1%.

Usage example
=============
The following example uses the cc-pVQZ FCIQMC energies of the carbon
dimer [J. Chem. Phys. 135, 084104 (2011);
https://doi.org/10.1063/1.3624383] as the potential.  The code
produces a plot.dat data file and a plot.gpi gnuplot file intended to
be used with the cairolatex terminal.

```
$ cat curve_fciqmc.dat
1.7008 -75.4360   0.0003
1.8897 -75.6547   0.0003
2.0787 -75.7614   0.0002
2.2677 -75.7987   0.0002
2.3480 -75.80251  0.00008
2.4566 -75.7993   0.0002
2.6456 -75.7798   0.0001
3.0236 -75.7243   0.0002
3.4015 -75.6805   0.0002
3.7795 -75.6454   0.0002
4.1574 -75.6185   0.0002
4.5353 -75.5998   0.0002
4.9133 -75.5881   0.0002
5.2912 -75.5798   0.0003
$ vib_levels 
Enter name of data file containing {r,E,dE} data:
curve_fciqmc.dat

File "curve_fciqmc.dat" contains 14 data lines.

Effective mass in Da [e.g., 0.5 for H2; empty for 1 a.u.]:
6.0

Polynomial expansion order for potential [empty for 4 (= quartic)]:
9

Enter name *root* for plot files [empty to skip plot]:
plot

Expectation values:
  Delta0          :    4.214638069410E-03   9.131764138489E-06
  Delta1          :    1.253310512858E-02   2.621895212204E-05
  Delta2          :    2.069359624749E-02   4.154465185848E-05
  Delta3          :    2.869279856425E-02   5.528378928990E-05
  Virial ratio    :    1.000000000650E+00   1.048993975617E-10
  Orb. exp. order :    2.229000000000E+01   2.637319575906E+00
  r of minimum    :    2.355282891955E+00   7.865957883420E-04
  Minimum energy  :   -7.580264921476E+01   9.635274539802E-05

Parameters reproducing mean (offset by x,y of min in data):
  c_0:   -1.179754474430E-04
  c_1:   -5.793304981712E-03
  c_2:    4.029637121329E-01
  c_3:   -4.769773126791E-01
  c_4:    2.666506504569E-01
  c_5:   -1.044093170722E-01
  c_6:    5.889951538058E-02
  c_7:   -3.225047051166E-02
  c_8:    9.152961852952E-03
  c_9:   -9.753593211581E-04

Made gnuplot plot.gpi file and plain-text plot.dat data file.
```
