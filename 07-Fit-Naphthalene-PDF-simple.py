import numpy as np

from pyobjcryst import loadCrystal
from diffpy.srfit.pdf import PDFContribution
from diffpy.srfit.fitbase import Profile, FitRecipe, FitResults

nphcrystal = loadCrystal('naphthalene.cif')

pdfcntb = PDFContribution('pdfcntb')
pdfcntb.loadData('naphthalene.gr')
pdfcntb.qdamp = 0.06
pdfcntb.setCalculationRange(1.1, 25)
pdfcntb.addStructure('nph', nphcrystal)

nphfit = FitRecipe()
nphfit.clearFitHooks()
nphfit.addContribution(pdfcntb)

nphfit.addVar(pdfcntb.scale, name='scale')
nphfit.addVar(pdfcntb.nph.delta2, value=1.0)
nphase = pdfcntb.nph.phase
# unit cell parameters
nphfit.addVar(nphase.a)
nphfit.addVar(nphase.b)
nphfit.addVar(nphase.c)
# cell-angle beta is in radians in ObjCryst Crystal
# we will refine angle in degrees.
nphfit.newVar('beta', value=np.degrees(nphase.beta.value))
nphfit.constrain(nphase.beta, 'radians(beta)')
# all carbon species have the same displacement parameter,
# it is sufficient to add constraint for the C1 atom
nphfit.addVar(nphase.C1.Biso, name='Biso', value=1.0)

from scipy.optimize import leastsq
leastsq(nphfit.residual, nphfit.values)
results = FitResults(nphfit)

r = pdfcntb.r.value
gobs = pdfcntb.y.value
gcalc = pdfcntb.evaluate()

from naphthalene_functions import differenceplot
differenceplot(r, gobs, gcalc, baseline=-0.75)
show()
