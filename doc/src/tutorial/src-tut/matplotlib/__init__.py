"""
Purpose: make imports from matplotlib.pylab mean scitools.std if
matplotlib is not installed.
"""

import sys

try:
    import matplotlib
    mpl = True
except ImportError:
    mpl = False

if not mpl:
    # Use scitools instead and make it look as matplotlib
    try:
        import scitools.std
        sys.modules['matplotlib.pyplot'] = scitools.std  # essential
    except ImportError:
        print 'Cannot import matplotlib or scitools.std - cannot plot'
        import sys; sys.exit(1)

pyplot = scitools.std  # essential
