%% PARTICLE SWARM OPTIMIZATION (PSO)

function [bestfitIter,best_Particle,bestParticle_Position] = PSO(ObjFun,LB,UB,nVar,nPop,MaxIt)
%% Problem Parameters
prob = ObjFun;                  % fitness function
Var = nVar;                     % number of variables
lb = LB;                        % Variable Lower Bound
ub = UB;                        % Variable Upper Bound
cPop = nPop;                    % class population size
Maxiter = MaxIt;                % max no of iteration

%% PSO Parameters
w=0.45;                          % Inertia Weight
c1=1;                         % Personal Learning Coefficient
c2=1;                         % Global Learning Coefficient

%% Starting the Particle Swarm Optimization
f = NaN(cPop,1);                    % vector for fitness fun value of the population
D = Var;                            % length of decision variables

bestfitIter = NaN(Maxiter+1,1);     % vetcor to store best fitness fun value per iteration

P = repmat(lb,cPop,1) + repmat((ub-lb),cPop,1).*rand(cPop,D);   % generation of initial positon of the particles in the population
V = repmat(lb,cPop,1) + repmat((ub-lb),cPop,1).*rand(cPop,D);   % generation of initial velocity of each particle

for q = 1:cPop
    f(q) = prob(P(q,:));        % evaluating the fitness fun value of initial population
end

bestfitIter(1) = min(f);            % vector to store best fitness value of the 0th iteration 

pbest = P;                      % initalize the personal best solution          
f_pbest = f;                    % initalize the fitness of the personal best solution  

[f_gbest,ind] = min(f_pbest);   % initalize the best objective function value 
gbest = P(ind,:);               % dtermine the best solution

%% Iteration Loop (Main Loop)
for t = 1:Maxiter
    for q = 1:cPop
        V(q,:) = w*V(q,:) + c1*rand(1,D).*(pbest(q,:)-P(q,:)) + c2*rand(1,D).*(gbest - P(q,:));     % Determine the new velocity
        
        P(q,:) = P(q,:) + V(q,:);   % update the position
        
        P(q,:) = max(P(q,:),LB);    % bounding the violating to their lower bound
        P(q,:) = min(P(q,:),UB);    % bounding the violating to their upper bound
        
        f(q) = prob(P(q,:));        % determine the fitness of the new solution
        
        if f(q) < f_pbest(q)
           f_pbest(q) = f(q);       % update the fitness function value of the personal best soultion
           pbest(q,:) = P(q,:);     % update the personal best soultion
           
           if f_pbest(q) < f_gbest
           f_gbest = f_pbest(q);        % update the fitness function value of the global best soultion
           gbest = pbest(q,:);          % update the global best soultion
           end
        end
    end
    bestfitIter(t+1) = f_gbest;
end

% Results
best_Particle = f_gbest;
bestParticle_Position = gbest;
