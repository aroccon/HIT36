![](vis.gif)


~~~text
__/\\\________/\\\__/\\\\\\\\\\\__/\\\\\\\\\\\\\\\_____/\\\\\\\\\\_____________/\\\\\_
 _\/\\\_______\/\\\_\/////\\\///__\///////\\\/////____/\\\///////\\\________/\\\\////__
  _\/\\\_______\/\\\_____\/\\\___________\/\\\________\///______/\\\______/\\\///_______
   _\/\\\\\\\\\\\\\\\_____\/\\\___________\/\\\_______________/\\\//_____/\\\\\\\\\\\____
    _\/\\\/////////\\\_____\/\\\___________\/\\\______________\////\\\___/\\\\///////\\\__
     _\/\\\_______\/\\\_____\/\\\___________\/\\\_________________\//\\\_\/\\\______\//\\\_
      _\/\\\_______\/\\\_____\/\\\___________\/\\\________/\\\______/\\\__\//\\\______/\\\__
       _\/\\\_______\/\\\__/\\\\\\\\\\\_______\/\\\_______\///\\\\\\\\\/____\///\\\\\\\\\/___
        _\///________\///__\///////////________\///__________\/////////________\/////////_____
~~~

~~~text
Pseudo-spectral code for DNS of Homogenous isotropic turbulence.
Base code is hit3d rearranged and modified.
Main Additions:
  *Paraview files (fields and particles)
  *Automatic compilation and setup (compile.sh)
  *Paralelizzation strategy modified
  *Optimization and code organization
  *Terminal output
  *Organization of the results folder
Main deletions:
  *LES part


Currently working on the following machines:
-OS X (MacBook Pro)
-UNIX (ubuntu 18.04)


Things to be done/updated (ASAP)
-Particles tracking (input, output files, remove cint? only trilinear?)
-Particles output for paraview, reading of the output is fine but files are not generated.
-Compile.sh (automatic setup of input.f90 file)


To run a simulation:
A) Setup the input.f90 file.
B) execute ./compile.sh


Parallelization strategy
Four (three) possible strategies depending on the type of flow considered: single-phase flow or particles-laden flow.
-Single-phase flow:
  *split=0 (no MPI communicator splitting): all task are solving the Eulerian fields and doing the stats (suggested strategy)
  *split=1 (MPI communicator splitting): 2/3 of tasks solve the Eulerian fields and 1/3 do the stats.
Particles-laden flow:
  *split=0 (no MPI communicator splitting): not possible at the moment.
  *split=1 (MPI communicator splitting): 1/2 of tasks solve the Eulerian fields and 1/4 do the stats and 1/4 track the particles.


Output and restart files.
Files containing the Eulerian fields (u_***,v_***,w_***,etc.) and the particle positions (p_***) are stored in set_run/results


Visualization with Paraview
Two possible cases: single-phase and particles-laden flow.
-Single-phase flow:
  *Go to set_run/results/paraview_fields and run go.sh, paraview files are generated in the output folder.
-Particles-laden flow:
  *Go to set_run/results/paraview_fields and run go.sh, paraview files are generated in the output folder.
  *Go to set_run/results/paraview_particles and run go.sh, paraview files are generated in the output folder (not working ATM)


Validation database:
https://torroja.dmt.upm.es/turbdata/agard/chapter3/HOM03/


References:
[1] S.G. Chumakov, "A priori study of subgrid-scale flux of a passive scalar in turbulence", Phys. Rev. E, 78 15563.
[2] S.G. Chumakov, "Scaling properties of subgrid-scale energy dissipation", Phys. Fluids, 19 058104.
[3] L. Machiels, "Predictability of Small-Scale Motion in Isotropic Fluid Turbulence", Phys. Rev. Lett. 79, 3411.
[4] J. Jimenez, A.A. Wray, P.G. Saffman and R.S. Rogallo, "The structure of intense vorticity in isotropic turbulence", J. Fluid Mech., 255, 65-90.
~~~
