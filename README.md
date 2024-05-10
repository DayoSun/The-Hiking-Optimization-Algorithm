# The-Hiking-Optimization-Algorithm

In this project, a novel metaheuristic called ‘The Hiking Optimization Algorithm’ (HOA) is proposed. HOA is inspired by hiking, a popular recreational activity, in recognition of the similarity between the search landscapes of optimization problems and the mountainous terrains traversed by hikers. HOA’s mathematical model is premised on Tobler’s Hiking Function (THF), which determines the walking velocity of hikers (i.e. agents) by considering the elevation of the terrain and the distance covered. THF is employed in determining hikers’ positions in the course of solving an optimization problem. HOA’s performance is demonstrated by benchmarking with 29 well-known test functions (including unimodal, multimodal, fixed-dimension multimodal, and composite functions), three engineering design problems (EDPs), (including I-beam, tension/compression spring, and gear train problems) and two N-P Hard problems (i.e. Traveling Salesman’s and Knapsack Problems). Moreover, HOA’s results are verified by comparison to 14 other metaheuristics, including Teaching Learning Based Optimization (TLBO), Genetic Algorithm (GA), Differential Evolution (DE), Particle Swarm Optimization, Grey Wolf Optimizer (GWO) as well as newly introduced algorithms such as Komodo Mlipir Algorithm (KMA), Quadratic Interpolation Optimization (QIO), and Coronavirus Optimization Algorithm (COVIDOA). In this study, we employ statistical tests such as the Wilcoxon rank sum, Friedman test, and Dunn’s post hoc test for the performance evaluation. HOA’s results are competitive and, in many instances, outperform the aforementioned well-known metaheuristics.


# Main Paper
The Hiking Optimization Algorithm: A novel human-based metaheuristic approach.
https://doi.org/10.1016/j.knosys.2024.111880

#  HOA Algorithm
HOA is located in the HOA_V2.m script. It is based on Tobler's hiking function approach. It has 6 inputs, namely: the objective function, Lower Bounds, Upper Bounds, dimension of the decision variables, number of search agents, and the number of run times. Run the HOA_main.m script 


# Inputs
The major inputs of HOA are the objective function to be optimized and the respective constraints. The objective function and the constraints are stated in the myFitn.m script. In the case of TSP, the objective functions are loaded as separate scripts which may contain the coordinates of the locations/points or a distance matrix.

# Outputs
The key outputs are the values of the decision variables, convergence curves, cost function, or the path cost and path for TSP. The convergence curves and path are in graphical format. However, the convergence curves may be suppressed as the case may be.


# Requirements
The following are the requirements of DSO
- MATLAB 2016a and later versions
- 8GB RAM size
- 256GB Hard disk size
- CPU
The size of the HOA_V2.m script is about 5KB and the HOA_main.m script is about 17KB depending on the objective functions and constraints. The size of all the scripts required to replicate the  results in the paper is about 500kB.

# Design Flow
After setting the parameters of HOA, the design flow is summarised as follows:


# How to use
- State the objective function and constraints in myFitn.m scripts
- Go to the HOA_V2.m script and change the DSO tuning params if you so desire. The params are at default values.
- Go to the HOA_main.m script and set the number of search agents, max. number of iterations, Monte Carlo runs, and the dimension of the problem.
- Run the HOA_main.m script.


# Authors
- Sunday O. Oladejo: School for Data Science and Computational Thinking, University of Stellenbosch, Stellenbosch, South Africa (e-mail: sunday@sun.ac.za)
- Stephen O. Ekwe:  Department of Electrical, Electronic, and Computer Engineering, Cape Peninsula University of Technology, Cape Town, South Africa (e-mail: ekwes@cput.ac.za)
- Seyedali Mirjalili: Centre for Artificial Intelligence Research and Optimisation, Torrens University Australia, Fortitude Valley, Brisbane, QLD 4006, Australia (e-mail: ali.mirjalili@torrens.edu.au)
