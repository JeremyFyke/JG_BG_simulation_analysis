import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import argparse

plt.close("all")

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument("--bg_case", dest="BG_CaseName", required=True)
parser.add_argument("--bg_years",nargs='+',dest="BG_years", required=True, type=int)
parser.add_argument("--jg_case", dest="JG_CaseName", required=True)
parser.add_argument("--jg_years",nargs='+',dest="JG_years", required=True, type=int)
a = parser.parse_args()

#Month to analyze (all months gets prohibitive)
Month='01'
nbins=100

nGLC=[701,376]
nOCN=(384,320)
nICE=[384,320]
nLND=[192,288]

PrintGLC=1
PrintOCN=1
PrintICE=1
PrintLND=1

def getDims(fname):
    vn=nc.Dataset(fname)
    nc_vars=[var for var in vn.variables]
    nc_dims=[dim for dim in vn.dimensions]
    return nc_vars,nc_dims,vn

def buildArray(ns,YearSpan,Fprefix,Fsuffix):
    nt=len(range(*YearSpan))
    Array=np.zeros((ns[0],ns[1],nt))

    for n,t in enumerate(np.arange(*YearSpan)):
	tlong=format(t,"04")
	fname=Fprefix+tlong+Fsuffix
	#print 'Loading: '+fname
	vn=nc.Dataset(fname)
	Array[:,:,n]=np.squeeze(vn.variables[var])
	vn.close()
    return Array

def PrintFigure(BG_VAR,JG_VAR,Component,VariableName,Outputdir):
    BG_VAR[BG_VAR>1.e10]=0.
    JG_VAR[JG_VAR>1.e10]=0.   
    if np.count_nonzero(BG_VAR) and np.count_nonzero(JG_VAR): #Count # of non-zero values
	plt.close("all")
	iBG=np.nonzero(BG_VAR)
	iJG=np.nonzero(JG_VAR)
	minrange=np.min(BG_VAR[iBG])
	maxrange=np.max(BG_VAR[iBG])
	BG_hist,base=np.histogram(BG_VAR[iBG],
				  bins=nbins,
				  range=(minrange,maxrange),
				  density=True)
	plt.plot(base[:-1],BG_hist,color='red',label='BG')
	JG_hist,base=np.histogram(JG_VAR[iJG],
				  bins=nbins,
				  range=(minrange,maxrange),
				  density=True)
				  
	plt.plot(base[:-1],JG_hist,color='blue',label='JG')
	plt.title(Component+': '+VariableName)
	plt.legend(loc='upper right', shadow=True)		
	if not os.path.exists(Outputdir):
            os.makedirs(Outputdir)			  	  
	outputfile=os.path.join(Outputdir,Component+'_'+var)
	print 'Printing: '+outputfile
	plt.savefig(outputfile)
	
if PrintGLC:
    ###GLC plots###
    #First, get variables and dimensions from an arbitrary history file
    nc_vars,nc_dims,vn=getDims('/glade/scratch/jfyke/'+a.BG_CaseName+'/run/'+a.BG_CaseName+'.cism.h.0001-01-01-00000.nc')
    for var in nc_vars:
	if var not in nc_dims: #exclude plotting dimension arrays
	    tmp=vn.variables[var]
	    if np.array_equal(np.squeeze(np.shape(tmp)),(1,nGLC[0],nGLC[1])): #exclude plotting fields that don't match expected 2D array shape

		Fsuffix='-01-01-00000.nc'

		Fprefix='/glade/scratch/jfyke/'+a.BG_CaseName+'/run/'+a.BG_CaseName+'.cism.h.'
		BG_VAR=buildArray(nGLC,a.BG_years,Fprefix,Fsuffix)

		Fprefix='/glade/scratch/jfyke/'+a.JG_CaseName+'/run/'+a.JG_CaseName+'.cism.h.'
		JG_VAR=buildArray(nGLC,a.JG_years,Fprefix,Fsuffix)

        	VariableName=vn.variables[var].getncattr('long_name')
        	Outputdir=a.JG_CaseName+'_figs'
        	PrintFigure(BG_VAR, JG_VAR, 'Glc', VariableName, Outputdir)

if PrintOCN:
    print 'hello'
    ###OCN plots###    
    #First, get variables and dimensions from an arbitrary history file
    nc_vars,nc_dims,vn=getDims('/glade/scratch/jfyke/'+a.BG_CaseName+'/run/'+a.BG_CaseName+'.pop.h.0001-'+Month+'.nc')
    for var in nc_vars:
        if var not in nc_dims: #exclude plotting dimension arrays
	    tmp=vn.variables[var]
	    if np.array_equal(np.squeeze(np.shape(tmp)),(1,nOCN[0],nOCN[1])): #exclude plotting fields that don't match expected 2D array shape

        	Fsuffix='-'+Month+'.nc'

        	Fprefix='/glade/scratch/jfyke/'+a.BG_CaseName+'/run/'+a.BG_CaseName+'.pop.h.'
		BG_VAR=buildArray(nOCN,a.BG_years,Fprefix,Fsuffix)

		Fprefix='/glade/scratch/jfyke/'+a.JG_CaseName+'/run/'+a.JG_CaseName+'.pop.h.'
		JG_VAR=buildArray(nOCN,a.JG_years,Fprefix,Fsuffix)

		VariableName=vn.variables[var].getncattr('long_name')
		Outputdir=a.JG_CaseName+'_figs'   
		PrintFigure(BG_VAR, JG_VAR, 'Ocn', VariableName, Outputdir)
		
if PrintICE:
    ###ICE plots###
    #First, get variables and dimensions from an arbitrary history file
    nc_vars,nc_dims,vn=getDims('/glade/scratch/jfyke/'+a.BG_CaseName+'/run/'+a.BG_CaseName+'.cice.h.0001-'+Month+'.nc')
    for var in nc_vars:
        if var not in nc_dims: #exclude plotting dimension arrays
	    tmp=vn.variables[var]
	    if np.array_equal(np.squeeze(np.shape(tmp)),(1,nICE[0],nICE[1])): #exclude plotting fields that don't match expected 2D array shape

        	Fsuffix='-'+Month+'.nc'

        	Fprefix='/glade/scratch/jfyke/'+a.BG_CaseName+'/run/'+a.BG_CaseName+'.cice.h.'
		BG_VAR=buildArray(nOCN,a.BG_years,Fprefix,Fsuffix)

		Fprefix='/glade/scratch/jfyke/'+a.JG_CaseName+'/run/'+a.JG_CaseName+'.cice.h.'
		JG_VAR=buildArray(nOCN,a.JG_years,Fprefix,Fsuffix)

		VariableName=vn.variables[var].getncattr('long_name')
		Outputdir=a.JG_CaseName+'_figs'   
		PrintFigure(BG_VAR, JG_VAR, 'Ice', VariableName, Outputdir)

if PrintLND:

    ###LND plots###
    #First, get variables and dimensions from an arbitrary history file
    nc_vars,nc_dims,vn=getDims('/glade/scratch/jfyke/'+a.BG_CaseName+'/run/'+a.BG_CaseName+'.clm2.h0.0001-'+Month+'.nc') 
    for var in nc_vars:
        if var not in nc_dims: #exclude plotting dimension arrays
	    tmp=vn.variables[var]
	    if np.array_equal(np.squeeze(np.shape(tmp)),(1,nLND[0],nLND[1])): #exclude plotting fields that don't match expected 2D array shape

        	Fsuffix='-'+Month+'.nc'

        	Fprefix='/glade/scratch/jfyke/'+a.BG_CaseName+'/run/'+a.BG_CaseName+'.clm2.h0.'
		BG_VAR=buildArray(nLND,a.BG_years,Fprefix,Fsuffix)

		Fprefix='/glade/scratch/jfyke/'+a.JG_CaseName+'/run/'+a.JG_CaseName+'.clm2.h0.'
		JG_VAR=buildArray(nLND,a.JG_years,Fprefix,Fsuffix)

		VariableName=vn.variables[var].getncattr('long_name')
		Outputdir=a.JG_CaseName+'_figs'   
		PrintFigure(BG_VAR, JG_VAR, 'Lnd', VariableName, Outputdir)
