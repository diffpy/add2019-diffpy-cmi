import numpy as np
from scipy.optimize import leastsq

from pyobjcryst import loadCrystal
from pyobjcryst.molecule import Molecule
from diffpy.srreal.pdfcalculator import ConstantPeakWidth
from diffpy.srfit.pdf import PDFContribution
from diffpy.srfit.fitbase import FitRecipe, FitResults

nphcrystal = loadCrystal('naphthalene.cif')
nphmol = Molecule(nphcrystal, "naphthalene")
numatoms = nphcrystal.GetNbScatterer()
atoms = [nphcrystal.GetScatterer(i) for i in range(numatoms)]
xyzf = np.array([(a.X, a.Y, a.Z) for a in atoms])
xyzc = np.array([nphcrystal.FractionalToOrthonormalCoords(x, y, z)
                 for x, y, z in xyzf])
xyzcmol = xyzc - xyzc.mean(axis=0)

spC1 = nphcrystal.GetScatteringPower('C1')
for a, (xc, yc, zc)  in zip(atoms, xyzcmol):
    nphmol.AddAtom(xc, yc, zc, spC1, a.GetName())
    nphcrystal.RemoveScatterer(a)
molposition = xyzf.mean(axis=0)
nphcrystal.AddScatterer(nphmol)
nphmol.X, nphmol.Y, nphmol.Z = xyzf.mean(axis=0)

pdfcntb = PDFContribution('pdfcntb')
pdfcntb.loadData('naphthalene.gr')
pdfcntb.qdamp = 0.06
pdfcntb.setCalculationRange(1.1, 25)
pdfcntb.addStructure('nphmol', nphmol, periodic=False)
pdfcntb.addStructure('widecrystal', nphcrystal, periodic=True)
pdfcntb.addStructure('widemolecule', nphcrystal, periodic=False)
pdfcntb.widecrystal._calc.peakwidthmodel = ConstantPeakWidth()
pdfcntb.widemolecule._calc.peakwidthmodel = ConstantPeakWidth()
from naphthalene_functions import fixpeakwidthparameters
fixpeakwidthparameters(pdfcntb)

pdfcntb.setEquation('scale * (nphmol + widecrystal - widemolecule)')

nphfit = FitRecipe()
nphfit.clearFitHooks()
nphfit.addContribution(pdfcntb)

nphfit.addVar(pdfcntb.scale, name='scale')
pcrystal = pdfcntb.widecrystal.phase
# unit cell parameters
nphfit.addVar(pcrystal.a)
nphfit.addVar(pcrystal.b)
nphfit.addVar(pcrystal.c)
# cell-angle beta is in radians in ObjCryst Crystal
# we will refine angle in degrees.
nphfit.newVar('beta', value=np.degrees(pcrystal.beta.value))
nphfit.constrain(pcrystal.beta, 'radians(beta)')
# all carbon species have the same displacement parameter,
# it is sufficient to add constraint for the C1 atom
pmol = pdfcntb.nphmol.phase
nphfit.addVar(pmol.C1.Biso, name='Biso', value=1.0)

# create new variable for intermolecular displacements.
# constrain the fwhm of a Gaussian peak from 2 atoms accordingly.
nphfit.newVar('Binter', value=1.5)
nphfit.constrain(pdfcntb.widecrystal.fwhm, 'sqrt(2 * log(2) * Binter) / pi')
nphfit.constrain(pdfcntb.widemolecule.fwhm, 'sqrt(2 * log(2) * Binter) / pi')

leastsq(nphfit.residual, nphfit.values)
results = FitResults(nphfit)

r = pdfcntb.r.value
gobs = pdfcntb.y.value
gcalc = pdfcntb.evaluate()

from naphthalene_functions import differenceplot
differenceplot(r, gobs, gcalc, baseline=-0.75)
