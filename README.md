# Optimal transportation networks via Mumford-Shah-type image inpainting 

## General information 

The code implements the method proposed in [1](), Chapter 4.4, for solving the branched transport [2](http://www.uvm.edu/pdodds/research/papers/others/2003/xia2003a.pdf),[3](https://pdfs.semanticscholar.org/d766/7ac83e8dd7c8ce452fe63775a3ddd705efd9.pdf) and urban planning [4](http://www.numdam.org/article/COCV_2005__11_1_88_0.pdf) problem via a convex reformulation as a Mumford–Shah-type image inpainting problem introduced by [5](https://arxiv.org/abs/1601.07402).

The implementation is based on the [QuocMesh software library](https://archive.ins.uni-bonn.de/numod.ins.uni-bonn.de/software/quocmesh/index.html), AG Rumpf, Institute for Numerical Simulation, University of Bonn. A new finite element class and its corresponding operators based on adaptive triangular prism grids have been added. The software is distributed under the conditions of the [Common Development and Distribution License](LICENSE.txt). 

Currently, the software works under Linux only and requires cmake and a gcc compiler version. Some subfunctions admit optional advanced variants requiring the [MOSEK Optimization toolbox](https://www.mosek.com/), which can be included if desired.


## Compile 

To get started, create a clone of the git repository by 

	git clone https://github.com/cdirks/OptimalTransportNetworks

and decompress the archive. The QuocMesh source code is provided in the directory quocmesh, while the binaries will be built to quocBuild. To compile, change to the directory quocBuild, open a terminal and execute

	./go.sh
	make -j<N>

where N is the number of threads, which will configure cmake and write the build files to quocBuild. If everything was installed correctly, there should be no errors occurring.


## Run example 

To run the code, change to the directory quocBuild/projects/OptimalTransportationNetworks and execute

	./OptimalTransportationNetworksOnPrismGrid

which will compute a simple example of branched transport between a single source and two point sinks with equal mass. The results (primal and dual three-dimensional images and delifted solution) will be written to quocBuild/projects/OptimalTransportationNetworks/Results. To change the example or any parameters, open the parameter file Parameters.par in quocmesh/projects/OptimalTransportationNetworks. In the following, we give a brief description of the listed parameters:

	problemtype: Problem type (options are BranchedTransport and UrbanPlanning)
	initialgridlevelxy: Initial xy-level of the grid
	initialgridlevelz: Initial s-level of the grid
	numrefinements: Number of refinement rounds
	maximalgridlevelxy: Maximal xy-level of the grid (corresponding to the smallest possible area of a simplex)
	maximalgridlevelz: Maximal s-level of the grid
	maxiter: Maximal number of iterations in the previous rounds
	maxiterfinalrun: Maximal number of iterations in the last round (in case of 0 refinements, in the first round)
	maxiterproject: Maximal number of Dykstra iterations for the projection onto the set K
	tol: Error tolerance in iteration (corresponding to the value computed in stoppingcriterion)
	tolproject: Error tolerance in Dykstra iterations
	epsilon: Branched transport/urban planning parameters
	a: Urban planning parameter
	taufactor: Factor for rebalancing of primal-dual stepsizes
	refinementtol: Refinement tolerance (elements with an error value greater than refinementtol will be refined)
	example: Transport example (see the list of supported examples in Example.h)
	stoppingcriterion: Iteration stopping criterion (options are PrimalDualGap and PrimalDifference)
	refinementcriterion: Refinement criterion (options are LocalGradient, LocalPrimalDualGap and LocalPrimalDualGapOnly)
	preconditioning: Preconditioning (yes or no)

	
	
[1] Carolin Dirks. Numerical methods for transportation networks. PhD thesis, 2019.

[2] Qinglan Xia. Optimal paths related to transport problems. *Commun. Contemp. Math.*, 5(2):251--279, 2003. 

[3] Francesco Maddalena, Sergio Solimini, and Jean-Michel Morel. A variational model of irrigation patterns. *Interfaces Free Bound.*, 5(4):391--415, 2003.

[4] Alessio Brancolini and Giuseppe Buttazzo. Optimal networks for mass transportation problems. *ESAIM Control Optim. Calc. Var.*, 11(1):88--101 (electronic), 2005.

[5] Alessio Brancolini, Carolin Rossmanith, Benedikt Wirth. Optimal micropatterns in 2D transport networks and their relation to image inpainting. *Archive for Rational Mechanics and Analysis*, 228(1):279--308, 2018.

