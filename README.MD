
\chapter{Implementation}

FF\_GenOpt was developed and tested on Debian with Linux kernel 4.9. It is written in Python and works with Python 3.5.3+. It consists of multiple modular components which can be adapted for other fitness functions and use cases. 

\section{Requirements}

The program requires a Python interpreter with \textbf{numpy} and \textbf{matplotlib}. The fitness computation modules requires \textbf{scipy} and a CHARMM installation.

\section{Architecture}

The program consists of multiple files which have to be in the same directory.

\textbf{FF\_GenOpt.py} contains the core module, genetic algorithm, random search algorithm, log management, population management, and parameter management.

\textbf{FF\_GenOpt\_conf.py} contains the configuration file management.

\textbf{FF\_GenOpt\_fitness.py} wraps the fitness function which depends on the file below.

\textbf{FF\_GenOpt\_extern.py} contains modified code from AFMM and Florian Jörg (from the Department of Computational Biological Chemistry at the University of Vienna) for calculating the fitness function.

\section{Configuration}

FF\_GenOpt uses configuration files that can be created and edited in any editor. The configuration files contain the execution command for CHARMM, the path to input files and output files, the optimization algorithm settings, and the parameter names, type, range, and initial values. Similar to Python, comments can be added by using \# in the beginning. A typical configuration file is shown below:

\begin{lstlisting}[breaklines=true]
MDEXEC                      charmm
MDINP                       Molecules/oac/oac.inp
MDOUT                       /dev/shm/opt.log
QMOUT                       Molecules/oac/oac.log
PARAM_FILENAME              FF_GenOpt_results/oac/parameters.str
POPULATION_FILENAME         FF_GenOpt_results/oac/lastPopulation_oac.txt
LOG_FILENAME                FF_GenOpt_results/oac/convergence_oac.log

GENERATIONS                 10000
POPULATION_SIZE             200
MUTATIONS_PER_GENERATION    100
STEP_SIZE                   0.05
CROSSOVER_BLX               0.7
CROSSOVER_SBX               0.25
CROSSOVER_UNIFORM           0.05
MINIMUM_IMPROVEMENT         0.0001
MINIMUM_DIVERSITY           0.0001

#------------------------------------------------
#         name  type        min     max     val
#------------------------------------------------
PARAMETER oacb1 bond 0 600 340
PARAMETER oacb2 bond 0 600 655.5162
PARAMETER oacb3 bond 0 600 316.83
PARAMETER oaca1 angle 0 200 33.0
PARAMETER oaca2 angle 0 200 79.9417
PARAMETER oaca3 angle 0 200 70.0
PARAMETER oaca4 angle 0 200 35.7
PARAMETER oacd1 dihedral 0 30 0.0
PARAMETER oaci1 improper 0 200 71.000
\end{lstlisting}

\texttt{MDINP} points at the CHARMM input file, which contains all commands for reading the force field parameters, computing the normal modes with vibran, and \textit{must point at the same file as} \texttt{PARAM\_FILENAME} (which is rewritten each time the fitness function is computed).

\texttt{MDOUT} contains the output of CHARMM and is only used to compute the fitness function.

\texttt{QMOUT} contains the output of Gaussian with the QM-calculated normal modes. This file should not change during the optimization.

\texttt{POPULATION\_FILENAME} contains the current population with all parameters, fitness values, the vibrational frequencies, and a small CHARMM script for setting the best set of parameters (on top). This file is rewritten after each generation, which enables the user to stop the algorithm at any time and take the best results.

\texttt{LOG\_FILENAME} contains the log file, see section \ref{ch:log}. It can be used to track and analyze the behavior of the algorithm with the current settings.

\texttt{POPULATION\_SIZE} defines the maximum population size. A large population covers a large parameter search space, but the initial convergence is usually slower. High-dimensional problems have larger search spaces and require larger populations than low-dimensional problems. If the population is too small, the population diversity will shrink quickly and computation effort will be wasted on generating new population members. If the population is too large, computation time will be wasted on optimizing a vast number of population members, especially in the beginning.

\texttt{GENERATIONS} is the maximum number of generations to be computed. This number is an upper limit and should not be too small to prevent premature computation stops. Note that computation time and convergence speed per generation can vary significantly and comparing generational progress is not always useful.

\texttt{MUTATIONS\_PER\_GENERATION} is a \textit{suggestion} for the number of crossover and mutation operations per generation. It is also used to control the number of attempts of the random search algorithm. The actual number of mutations and accepted new population members will differ between each generation.

\texttt{STEP\_SIZE} is used exclusively by the random search algorithm and defines the component-wise \textit{maximum} step size when computing a new direction vector. If this value is too small, the search algorithm wastes a lot of computation time with small steps and might get stuck in a small local minimum. If this value is too large then the idea of finding local descending directions by using the negative vector will fail often. For smooth high-dimension problems a smaller step size is preferable (see figure \ref{fig:sqrfn}). However, the force field fitting problem does not have such a smooth surface. A value of 0.05 seems to be a good compromise. In this case at most 5\% of the valid search range for each parameter will be used. It should be noted that this particular parameter is important and should be optimized if the performance is insufficient.

\texttt{CROSSOVER\_BLX}, \texttt{CROSSOVER\_SBX}, and \texttt{CROSSOVER\_UNIFORM} define the probability for choosing either of the crossover operators per mutation operation in the genetic algorithm. The sum of all probabilities \textit{must} be 1.0 or the program will fail. BLX and SBX are also parametrized, but the settings are fixed in code.

\texttt{MINIMUM\_IMPROVEMENT} is a parameter for the random search that defines the minimum relative fitness improvement per step. This parameter is not as important as others, but it prevents wasting computation time on practically irrelevant improvements. 0.01 means that the improvement must be at least 1\%, but the value can be much smaller.

\texttt{MINIMUM\_DIVERSITY} is a parameter that prevents the algorithm from getting stuck in fruitless computations after all population members converged to the same values. Instead of stopping the algorithm here, a new initial population is created but the population member with the best fitness is retained. This value can also be relatively small, e.g. around $10^{-5}$.

Finally, the \texttt{PARAMETER} keyword is used to define the parameter names and parameter type as well as the minimum, maximum, and original values. Genetic algorithms are designed for searching wide parameter ranges, but narrowing the parameter range down will lead to faster results.

\section{Population files}
Population files are generated at the end of each new generation and consist of three parts: All population members with their respective fitness value, the QM vibrational frequencies (multiplied by the QM factor of 0.957) and MM vibrational frequencies in absolute values and as a ratio (for the candidate with the best fitness), and the CHARMM script for setting the parameters of the population member with the best fitness.

Example population file for OAc\textsuperscript{-}:
\begin{lstlisting}
Generation 1244
22.39,[316.2, 47.47, 581.35, 33.73, 75.29, 18.96, 21.96, 0.16, 75.17]
22.65,[316.31, 43.64, 579.66, 33.47, 75.37, 18.79, 22.56, 0.15, 77.39]
...
31.35,[319.67, 72.37, 602.14, 35.75, 68.69, 16.71, 18.39, 0.16, 99.42]
31.41,[315.9, 81.56, 587.97, 35.69, 73.06, 13.74, 16.8, 0.06, 77.63]

   38.0832	:   37.1461
  399.6353	:  399.0209
  578.1499	:  571.3312
  585.1879	:  586.6166
  818.8192	:  830.4762
  922.3170	:  914.9241
  961.9735	:  948.7063
 1225.5265	: 1256.8892
 1301.2723	: 1321.7307
 1368.0702	: 1339.7813
 1379.1796	: 1348.1398
 1668.8337	: 1680.4878
 2848.6474	: 2830.5436
 2910.3859	: 2930.0955
 2933.2282	: 2930.8579

    0.9754
    0.9985
    0.9882
    1.0024
    1.0142
    0.9920
    0.9862
    1.0256
    1.0157
    0.9793
    0.9775
    1.0070
    0.9936
    1.0068
    0.9992

set oacb1 316.20
set oacb2 47.47
set oacb3 581.35
set oaca1 33.73
set oaca2 75.29
set oaca3 18.96
set oaca4 21.96
set oacd1 0.16
set oaci1 75.17
\end{lstlisting}

\section{Log files} \label{ch:log}

When starting an optimization process, the parameter settings are written into the log. The log is in TSV format and can be parsed easily by Excel-like tools or programming languages for analysis and plotting. A new line with 11 columns is appended after each generation. Example log file:

\begin{lstlisting}
Gen	BestF	MeanF	TEvals	GEvals	REvals	t	t/eval
15	61.53	94.68	3092	1883	1209	1558.78	0.50
16	48.35	86.30	3201	1991	1210	1617.76	0.51
17	41.80	83.07	3313	2102	1211	1679.63	0.51
\end{lstlisting}

\begin{lstlisting}
Delta	Improv	Divers
0.00	0.00	582.28
13.19	0.00	2639.27
6.55	0.00	345.59
\end{lstlisting}

The \texttt{Gen} column shows the current generation number. \texttt{BestF} and \texttt{MeanF} represent best fitness and mean fitness of the respective generation. \texttt{GEvals} and \texttt{REvals} show the number of function evaluations triggered by the genetic algorithm and by the rest of the algorithm (e.g. random search). \texttt{TEvals} is the sum of both. \texttt{t} is the total passed computation time in seconds while \texttt{t/eval} is the average time per evaluation. \texttt{Delta} is the improvement since the last generation, \texttt{Improv} is the improvement per evaluation, and \texttt{Divers} is the estimated diversity of the population.

\section{Running FF\_GenOpt}

After creating the configuration file and setting up CHARMM inputs, the program can be started with

\begin{lstlisting}
python3 FF_GenOpt.py configuration_file.cfg
\end{lstlisting}

It can be stopped by killing the process or keeping Ctrl+C pressed. Multiple sessions can be run at the same time remotely via the Linux tool \textbf{tmux}.

\section{Slicing}

The axis-parallel curve visualizations (see figure \ref{fig:oacslice} in appendix A) are created by a simple tool named GridSlicer. It is based on the first version FF\_GenOpt and uses the same configuration files. However, it requires different files to avoid interference with the optimization process and requires additional parameters.

\begin{lstlisting}[breaklines=true]
python3 GridSlicer.py im1.cfg im1/im1.inp /dev/shm/im1.log im1/im1.log im1/parameters.str im1/ 50 21
\end{lstlisting}

The parameters are: configuration file, CHARMM input file, CHARMM output file, log file, CHARMM force field parameters file, image output path, number of focus points (a point that all curves cross), number of computations per curve per dimension. The number of dimensions can be added as an additional parameter to cut down computation time for testing purposes. The output files are in PNG format. 
