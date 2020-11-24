# Loading the python libraries
import scanpy as sc
import scanpy.external as sce
import pickle
import logging
import scanorama
import trvae
from textwrap import wrap
# import trvaep
# from trvaep import pl

# Import user libraries
from gjainPyLib import *
sys.path.append('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/scripts')
# from scrnaseq_module import *

# For X11 display
import matplotlib
# matplotlib.use('TkAgg')
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D      # for 3D projection
from matplotlib.colors import ListedColormap # for sc.palette to colormap
from itertools import combinations           # pairwise combinations

# Third party libraries
sys.path.append('/home/rad/users/gaurav/projects/seqAnalysis/scrnaseq/scripts/scoreCT/src/')
import scorect_api as ct

# Reset random state
np.random.seed(2105)

# For using R inside python
# import rpy2's package module
# Using extensions: To load it each time IPython starts, 
# list it in configuration file: '/home/rad/.ipython/profile_default/ipython_config.py'
import rpy2.robjects as robjects
import rpy2.rinterface_lib.callbacks
from rpy2.robjects import pandas2ri
import anndata2ri

# Ignore R warning messages
rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)

# Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()

# Loading R libraries
from rpy2.robjects.packages import importr
scran        = importr('scran')
RColorBrewer = importr('RColorBrewer')
gam          = importr('gam')
ggplot2      = importr('ggplot2')
plyr         = importr('plyr')
MAST         = importr('MAST')

