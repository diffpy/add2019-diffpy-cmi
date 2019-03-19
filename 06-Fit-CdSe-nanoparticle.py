#!/usr/bin/env python
# -*- coding: utf-8 -*-

# We'll need numpy and matplotlib for plotting our results
import numpy as np
import matplotlib.pyplot as plt

# A least squares fitting algorithm from scipy
from scipy.optimize.minpack import leastsq

# DiffPy-CMI modules for building a fitting recipe
from diffpy.Structure import loadStructure
from diffpy.srfit.pdf import PDFContribution
from diffpy.srfit.fitbase import FitRecipe, FitResults

# Files containing our experimental data and structure file
dataFile = "cdse.gr"
structureFile = "cdse.xyz"

# The first thing to construct is a contribution. Since this is a simple
# example, the contribution will simply contain our PDF data and an associated
# structure file. We'll give it the name "cdse"
cdsePDF = PDFContribution("CdSe")

# Load the data and set the r-range over which we'll fit
cdsePDF.loadData(dataFile)
cdsePDF.setCalculationRange(xmin=1, xmax=20, dx=0.01)

# Add the structure from our xyz file to the contribution, since the structure
# model is non-periodic, we need to specify the periodic=False here to get the
# right PDF
cdseStructure = loadStructure(structureFile)
cdsePDF.addStructure("CdSe", cdseStructure, periodic=False)

# The FitRecipe does the work of managing one or more contributions
# that are optimized together.  In addition, FitRecipe configures
# fit variables that are tied to the model parameters and thus
# controls the calculated profiles.
cdseFit = FitRecipe()

# give the PDFContribution to the FitRecipe
cdseFit.addContribution(cdsePDF)

# Here we create variables for the overall scale of the PDF and a delta2
# parameter for correlated motion of neighboring atoms.
cdseFit.addVar(cdsePDF.scale, 1)
cdseFit.addVar(cdsePDF.CdSe.delta2, 5)

# We fix Qdamp based on prior information about our beamline.
cdseFit.addVar(cdsePDF.qdamp, 0.06, fixed=True)

# Since we are calculating PDF from a non-periodic structure, we also need to 
# specify the Qmin to get he correct PDF. The value of Qmin could be the actual
# Qmin in the experiment or the Qmin used in PDF transformation, or some value
# related to the size and the shape of the structure model. Usually a value 
# in (0.5 ~ 1.0) will give reasonable results.
cdsePDF.CdSe.setQmin(1.0)

# The Qmax used in PDF transformation should also be specfied
cdsePDF.CdSe.setQmax(20.0)

# We create the variables of ADP and assign the initial value to them. In this
# example, we use isotropic ADP for all atoms
CdBiso = cdseFit.newVar("Cd_Biso", value=1.0)
SeBiso = cdseFit.newVar("Se_Biso", value=1.0)

# For all atoms in the structure model, we constrain their Biso according to 
# their species  
atoms = cdsePDF.CdSe.phase.getScatterers()
for atom in atoms:
    if atom.element == 'Cd':
        cdseFit.constrain(atom.Biso, CdBiso)
    elif atom.element == 'Se':
        cdseFit.constrain(atom.Biso, SeBiso)

# Now we create a zoomscale factor which stretches the structure model, this is 
# useful when you want to fit the bond length. Note that the relative position
# of atoms are not changed during the refinements
zoomscale = cdseFit.newVar('zoomscale', value=1.0)

# Here is a simple we to assign the zoomscale to the structure. Note that this
# only works for NON-PERIODIC structure 
lattice = cdsePDF.CdSe.phase.getLattice()
cdseFit.constrain(lattice.a, zoomscale)
cdseFit.constrain(lattice.b, zoomscale)
cdseFit.constrain(lattice.c, zoomscale)

# Turn off printout of iteration number.
cdseFit.clearFitHooks()

# We can now execute the fit using scipy's least square optimizer.
print "Refine PDF using scipy's least-squares optimizer:"
print "  variables:", cdseFit.names
print "  initial values:", cdseFit.values
leastsq(cdseFit.residual, cdseFit.values)
print "  final values:", cdseFit.values
print

# Obtain and display the fit results.
cdseResults = FitResults(cdseFit)
print "FIT RESULTS\n"
print cdseResults

# Plot the observed and refined PDF.

# Get the experimental data from the recipe
r = cdseFit.CdSe.profile.x
gobs = cdseFit.CdSe.profile.y

# Get the calculated PDF and compute the difference between the calculated and
# measured PDF
gcalc = cdseFit.CdSe.evaluate()
baseline = 1.1 * gobs.min()
gdiff = gobs - gcalc

# Plot!
plt.figure()
plt.plot(r, gobs, 'bo', label="G(r) data")
plt.plot(r, gcalc, 'r-', label="G(r) fit")
plt.plot(r, gdiff + baseline, 'g-', label="G(r) diff")
plt.plot(r, np.zeros_like(r) + baseline, 'k:')
plt.xlabel(r"$r (\AA)$")
plt.ylabel(r"$G (\AA^{-2})$")
plt.legend()

plt.show()
