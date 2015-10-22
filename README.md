# python_histograms
Python code to generate 2D histograms with 1D histograms along the x &amp; y axes

The starHist2d-MassLogNorm.py program seems to work fine and generates a 2D histogram that 
is log of the weighted mass of the star particles. 

The starHist2d-MassLogNorm-fullNorm.py is supposed to normalize this data by the area of the bins
(d log x / d log y) and the volume of the simulation: 27 Mpc^3 / h^3. This isn't working yet.
