README
This repo contains the scripts to replicate the analyses presented in Grogan et al., 2019, "A new toolbox to distinguish the sources of spatial memory error".

You will need:
	MemToolbox2D. (Grogan, J.P, Fallon, S.J., Zokaei, N., Husain, M., Coulthard, E.J., Manohar, S.G. MemToolbox2D. v1.0.0 (2019). doi:10.5281/zenodo.3355381; Grogan et al., 2019, A new toolbox to distinguish the sources of spatial memory error) 
	
	MemToolbox. (Suchow, J. W., Brady, T. F., Fougnie, D., & Alvarez, G. A. (2013). 
		Modeling visual working  memory with the MemToolbox. Journal of Vision, 
		13(10):9, 1–8.  doi:10.1167/13.10.9. MemToolbox.org)

	The parallel computing toolbox in Matlab
	
	The matlib package (http://www.smanohar.com/matlib.php)
	
	You will also need to make a ./Data and ./Figs folder in the current working folder.
	

Contents:

	col.m	- columnise data
	quantPlot.m - plot data by quantiles (used for some figures)
	
	MemToolbox1DSim.m - simulate and fit 1D models across parameter sweep
	MemToolbox2DSim.m - simulate and fit 2D models across a parameter sweep
	
	MemToolbox2DSimNTrials.m - simulations with different numbers of parameters
	MemToolbox2DSimNDist.m - different numbers of distractors
	MemToolbox2DSimStimConstr.m - Effect of stimulus constraints on parameter recovery

	MemToolbox2DSimConstBias.m - constant translation bias
    MemToolbox2DSimPropBias.m - proportional bias to edge
	MemToolbox2DSimRadBias.m - radial bias to screen centre
	MemToolbox2DSimDistrRadial .m - radial bias to distractor location

	MemToolbox2DSim2AFC.m - 2AFC data simulation and fitting
	MemToolbox2DSimResampling.m - Edge-resamplings vs edge-constraining
	
	MemToolbox2DPosteriors.m  - Plot heatmaps of posterior samples from 1D and 2D fits
	
	MemToolbox2DBehMetrics.m - Compare parameter recovery to behavioural metrics on simulated data
	
	MemToolbox2DFallonPosteriors.m - Plot 1D and 2D posterior samples for parameters from Fallon et al., (2018), Ignoring versus updating in working memory reveal differential roles of attention and feature binding. Cortex. 107, pp. 50-63.
	
	MemToolbox2DPaperFigs.m - Plot all the figures (apart from those using human data) from Grogan et al., 2019 paper