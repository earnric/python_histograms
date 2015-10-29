# python_histograms
##Version 1.0
This version generates normalized 2D histograms where the bins are
normalized and weighed by the mass of the star particles. For the
plots where Z<sub>sol</sub> is the x-axis, the weight of the entire star
particle is used. Where the x-axis is Z<sub>sol</sub>/f<sub>pol</sub>
the weight is the mass of the polluted portion of the star particle
(mass * (1-f<sub>pristine</sub>)).
The fully normalized plots divide out the bin area (hist2d) or the bin
width (1d) along with the comoving volume of the simulation (27
Mpc/0.71<sup>3</sup>)

##Initial version 
Python code to generate 2D histograms with 1D histograms along the x &amp; y axes

The starHist2d-MassLogNorm.py program seems to work fine and generates a 2D histogram that 
is log of the weighted mass of the star particles. 

The starHist2d-MassLogNorm-fullNorm.py is supposed to normalize this data by the area of the bins
(d log x / d log y) and the volume of the simulation: 27 Mpc^3 / h^3. This isn't working yet.
