64              !dimensions of the grid (always a cube)
----------------------------------------------------------------------
000000          !ITMIN    first timestep (if >0, restart from saved)
20000	       	  !ITMAX    last timestep
10		          !IPRNT1   how often to calculate statistics
200		          !IPRNT2   how often to write restart files
200		          !IWRITE4  how often to write output files
----------------------------------------------------------------------
10.0		        !TMAX     the maximum simulation time
1.0 3           !TRESCALE, NRESCALE  rescale time and # of rescales
0.0             !TSCALAR  time when to start moving scalars in the flow
----------------------------------------------------------------------
1               ! flow type: 0=decay, 1=forced
----------------------------------------------------------------------
40.0	          !	RE	 Reynolds number (1/viscosity)
0.0005		      ! DT	 timestep
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
8		            ! np, number of particles, np.gt.0 means laden flows
0 		          ! particle_tracking mechanism 
0.0	            !	time_p, when to release the particles into the flow
0.0		          ! particle_filter_size, in case we average the velocity
----------------------------------------------------------------------
0		            ! nums, # of passive scalars.
----------------------------------------------------------------------
1      1.0     4.0    2. ! ## scalars: type, Sc, infrared exp, peak wavenumber (see below)
1      1.0     4.0    4.
1      1.0     4.0    8.
======================================================================
