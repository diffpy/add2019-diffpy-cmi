def differenceplot(x, yobs, ycalc, baseline=0):
    """Nice plot of data, simulation and difference curve.
    """
    from matplotlib.pyplot import plot, setp, hlines
    ydiff = yobs - ycalc
    rv = plot(x, yobs, 'bo', x, ycalc, 'r-', x, ydiff + baseline, 'g-')
    setp(rv[0], markeredgecolor='blue', markerfacecolor='none')
    hlines(baseline, x.min(), x.max(), linestyles='dashed')
    return rv


def fixpeakwidthparameters(pdfcontribution):
    """Replace parameters for r-dependent peak width with a constant
    peak width 'fwhm'."""
    from diffpy.srfit.fitbase.parameter import ParameterAdapter
    cntb = pdfcontribution
    cntb.unconstrain(cntb.widecrystal.qbroad)
    cntb.unconstrain(cntb.widemolecule.qbroad)
    cntb.removeParameter(cntb.qbroad)
    def addremovepars(pgen):
        pgen.addParameter(ParameterAdapter('fwhm',
            pgen._calc.peakwidthmodel, attr='width'))
        pgen.removeParameter(pgen.delta1)
        pgen.removeParameter(pgen.delta2)
        pgen.removeParameter(pgen.qbroad)
        return
    addremovepars(cntb.widecrystal)
    addremovepars(cntb.widemolecule)
    return
