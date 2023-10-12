nxnynz              !dimensions of the grid (always a cube)
----------------------------------------------------------------------
firstiter          !ITMIN    first timestep (if >0, restart from saved)
laststep      	  !ITMAX    last timestep
dumpstats		          !IPRNT1   how often to calculate statistics
dumprestart		          !IPRNT2   how often to write restart files
dumpoutput		          !IWRITE4  how often to write output files
----------------------------------------------------------------------
itmax		        !TMAX     the maximum simulation time
1.0 3           !TRESCALE, NRESCALE  rescale time and # of rescales
0.0             !TSCALAR  time when to start moving scalars in the flow
----------------------------------------------------------------------
1               ! flow type: 0=decay, 1=forced
----------------------------------------------------------------------
reeeeee	          !	RE	 Reynolds number (1/viscosity)
dtttttt		      ! DT	 timestep
----------------------------------------------------------------------
0               ! ISPCV1	 Initial spectrum type (flow field IC)
4               ! mv1	 initial infrared exponent in the spectrum
2.0             ! wm0v1	 initial peak wavenumber in the spectrum
----------------------------------------------------------------------
1	              !	force_type	Forcing type (see README)
2		            ! KFMAX		Forcing range from k=1 to KFMAX
0.5		          ! FAMP 		Imposed energy input
----------------------------------------------------------------------
0	              !	dealias (0:2/3-rule, 1:truncation and phase shift)
----------------------------------------------------------------------
0               ! det_rand  0=generate RN1,RN2, RN3; 1=fix them as below
67748937.       ! RN1
93712139.       ! RN2
10765397.       ! RN3, the random number seeds
----------------------------------------------------------------------
npppppp		            ! np, number of particles, np.gt.0 means laden flows
0 		          ! particle_tracking mechanism (0=trilinear, 1=cint)
0.01	            !	time_p, when to release the particles into the flow
0.0		          ! particle_filter_size, in case we average the velocity
----------------------------------------------------------------------
1		            ! nums, # of passive scalars.
----------------------------------------------------------------------
1      1.0     4.0    2.   0.  ! ## scalars: type, Sc, infrared exp, peak wavenumber, reaction sc (see below)
1      1.0     4.0    4.   0.
1      1.0     4.0    8.   0.
======================================================================
