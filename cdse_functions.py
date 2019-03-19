def differenceplot(recipe, baseline=None, title='', fig=None):
    """Nice plot of data, simulation and difference curve.
    """
    import numpy as np
    from matplotlib.pyplot import setp, subplots
    cntb = next(iter(recipe._contributions.values()))
    x = cntb.profile.x
    yobs = cntb.profile.y
    ycalc = cntb.evaluate()
    ydiff = yobs - ycalc
    if baseline is None:
        baseline = np.min(yobs) - 0.2 * (np.max(yobs) - np.min(yobs))
    _, ax = subplots() if fig is None else subplots(num=fig)
    rv = ax.plot(x, yobs, 'bo', x, ycalc, 'r-', x, ydiff + baseline, 'g-')
    setp(rv[0], markeredgecolor='blue', markerfacecolor='none')
    ax.hlines(baseline, x.min(), x.max(), linestyles='dashed')
    if title:
        ax.set_title(title)
    return rv
