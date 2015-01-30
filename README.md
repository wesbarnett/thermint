#Thermodynamic Integration Utility

##About

This utility is designed to be used with GROMACS output from a set of free
energy simulations. This program calculates free energy difference estimates
through thermodynamic integration. Additionally it calculates uncertainty in the
estimate through a bootstrap calculation.

##Prerequisites

You need Boost and GROMACS.

##Installation

    git clone https://github.com/wesbarnett/thermint
    cd thermint 
    mkdir build
    cd build
    cmake ..
    make -j 4
    make install

The GMXPREFIX environmental variable is used to find the GROMACS install
location. Either set it manually or [source the GMXRC file in your GROMACS
installation](http://www.gromacs.org/Documentation/Installation_Instructions#Getting_access_to_GROMACS_after_installation).

##Running

By default the program uses the trapezoid method for the integration, but
Simpson's rule and Gaussian-Legendre quadrature are also available. Beware that
in order to use Gaussian-Legendre quadrature for the integration you must have
run your simulations at the appropriate lambdas. Both the trapezoid rule and
Simpson's method do not require specific lambdas. Equal spacing of simulations
for trapezoid rule and Simpson's method is not required.

When running the integration calculation you must specify the energy files. In
order to do this you must have set separate-dhdl-file as 'no' in your mdp files.
If you had it set as 'yes' (the default), you're out of luck with that set of
simulations (sorry). The easiest way to specify the files is to do the
following: If 'ener' is the common prefix of all of your energy files, do:

    thermint -f ener*.edr
            
Note that dV/dl, which is required for the integration, is read in directly from
the energy files so that the output frequency of these values from your
simulation was determined by the nstenergy variable in your mdp file."

As mentioned above the standard error is estimated using a bootstrap
calculation. The number of iterations in the bootstrap calculation can be
changed with --nboot. By default 1000 iterations are performed to get the
bootstrapped average and then variance. Block bootstrapping is used, such that
each set of dV/dl is split up into several blocks (5 by default) which are
randomly selected with replacement for each free energy calculation. Five blocks
are chosen by default under the assumption that your simulation will be at least
five times the correlation time. The standard error is then estimated to be the
bootstrap standard deviation. The number of blocks can be changed with
--nblocks. Choosing too many blocks will give an uncertainty which is much lower
than the true uncertainty.

When you run the program it will tell you which group (right now coul, vdw,
mass, or restraint) and the lambda values and the file. At the end of the
calculation you'll get a breakdown of the results along with the total free
energy change estimate and uncertainty.
