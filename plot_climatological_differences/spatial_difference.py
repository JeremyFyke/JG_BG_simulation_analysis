import sys
import gzip
import glob
import os
import re
import numpy as np
from netCDF4 import MFDataset
import matplotlib.pyplot as plt

CaseName="JG_iteration_7"
output_dir='plots_difference'
input_dir='/glade/scratch/jfyke/'+CaseName+'/run/'
variable="SALT"

print "Generating dataset structures"
PreData=MFDataset(input_dir+CaseName+'.pop.h.000*')
PostData=MFDataset(input_dir+CaseName+'.pop.h.014*')

print "Reading data"
Pre=np.mean(PreData.variables[variable][:,0,:,:],axis=0)
Post=np.mean(PostData.variables[variable][:,0,:,:],axis=0)

print "Plotting data"
plt.pcolor(Post-Pre)
plt.colorbar()
plt.show()
