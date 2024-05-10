# The-Hiking-Optimization-Algorithm

In this project, a novel metaheuristic called ‘The Hiking Optimization Algorithm’ (HOA) is proposed. HOA is inspired by hiking, a popular recreational activity, in recognition of the similarity between the search landscapes of optimization problems and the mountainous terrains traversed by hikers. HOA’s mathematical model is premised on Tobler’s Hiking Function (THF), which determines the walking velocity of hikers (i.e. agents) by considering the elevation of the terrain and the distance covered. THF is employed in determining hikers’ positions in the course of solving an optimization problem. HOA’s performance is demonstrated by benchmarking with 29 well-known test functions (including unimodal, multimodal, fixed-dimension multimodal, and composite functions), three engineering design problems (EDPs), (including I-beam, tension/compression spring, and gear train problems) and two N-P Hard problems (i.e. Traveling Salesman’s and Knapsack Problems). Moreover, HOA’s results are verified by comparison to 14 other metaheuristics, including Teaching Learning Based Optimization (TLBO), Genetic Algorithm (GA), Differential Evolution (DE), Particle Swarm Optimization, Grey Wolf Optimizer (GWO) as well as newly introduced algorithms such as Komodo Mlipir Algorithm (KMA), Quadratic Interpolation Optimization (QIO), and Coronavirus Optimization Algorithm (COVIDOA). In this study, we employ statistical tests such as the Wilcoxon rank sum, Friedman test, and Dunn’s post hoc test for the performance evaluation. HOA’s results are competitive and, in many instances, outperform the aforementioned well-known metaheuristics.


# Main Paper
The Hiking Optimization Algorithm: A novel human-based metaheuristic approach.
https://doi.org/10.1016/j.knosys.2024.111880

#  HOA Algorithm
HOA is located in the DSO_TPM.m script. It is based on the Two-Process approach. It has 6 inputs, namely: the objective function, Lower Bounds, Upper Bounds, dimension of the decision variables, number of search agents, and the number of run times.
