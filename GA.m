%% GENETIC ALGORITHM (GA)

function [bestfitGene,best_Gene,bestGene_Position] = GA(ObjFun,LB,UB,nVar,nPop,MaxIt)
%% Problem Parameters
prob = ObjFun;                  % fitness function
Var = nVar;                     % number of variables
lb = LB;                        % variable Lower Bound
ub = UB;                        % variable Upper Bound
cPop = nPop;                    % class population size
Maxiter = MaxIt;                % max no of iteration


%% GA Parameters
etac = 20;                     % distribution index for crossover
etam = 20;                     % distribution index for mutation
Pc = 0.8;                      % probability of crossover
Pm = 0.2;                      % probability of mutation

%% Starting the Genetic Algorithim
% Initaiization
f = NaN(cPop,1);                 % population vetcor for objective fun value
OffspringObj = NaN(cPop,1);      % population vetcor for fitness fun value

bestfitGene = NaN(Maxiter+1,1);  % vetcor to store best fitness fun value per iteration

% length of decision variables
D = Var;
% generation of initial population (or food source)
P = repmat(lb,cPop,1) + repmat((ub-lb),cPop,1).*rand(cPop,D);

for q = 1:cPop
    f(q) = prob(P(q,:));        % evaluating the fitness fun value
end

bestfitGene(1) = min(f);        % vector to store best fitness value of the 0th iteration

%% Iteration Loop (Main Loop)
for t = 1:Maxiter
    
    %% Tournament Selection
    MatingPool = TournamentSelection(f,cPop);       % Performing tournaments to select best solution
    Parent = P(MatingPool,:);                       % Selecting parents solution
    
    %% Crossover
    offspring = CrossoverSBX(Parent,Pc,etac,lb,ub);
    
    %% Mutation
    offspring = MutationPloy(offspring,Pm,etam,lb,ub);
    
    for j = 1:cPop
        OffspringObj(j) = prob(offspring(j,:));    % Evaluating the fitness of the offspring
    end
    
    CombinedPopulation = [P; offspring];
    [f,ind] = sort([f;OffspringObj]);             % mu + laambda selection makes the first solution the best solution
    
    f = f(1:cPop);          
    P = CombinedPopulation(ind(1:cPop),:);         % extract best fitness population
    
    bestfitGene(t+1) = min(f);
end
% Results
[best_Gene,ind] = min(f);
bestGene_Position = P(ind,:);
end

% function for Tournament Selection
function MatingPool = TournamentSelection(f,cPop)
% Tournament Selction allows each solution to participate exactly twice

MatingPool = NaN(cPop,1);               % vector to store index of parent soultion
indx = randperm(cPop);                  % randomly shuffing the index of population member

for i = 1:(cPop-1)                      % Poolsize in cPop
    Candidate = [indx(i) indx(i+1)];    % selecting one pair of Population member for tournament
    CandidateObj = f(Candidate);
    [~,ind] = min(CandidateObj);        % selecting winner based on minimum fitness value
    MatingPool(i) = Candidate(ind);     % storing the index of the winner
end

% Tournament selection between the last and first member
Candidate = [indx(cPop) indx(1)];
CandidateObj = f(Candidate);
[~,ind] = min(CandidateObj);
MatingPool(cPop) = Candidate(ind);
end

% function for Crossover
function offspring = CrossoverSBX(Parent,Pc,etac,lb,ub)

[NP, Dim] = size(Parent);       % Determine the no. of population and decision variables
indx = randperm(NP);            % Permutating numbers from 1 to cPop
Parent = Parent(indx,:);        % Randomly shuffling parent solutions
offspring = NaN(NP,Dim);        % matrix to store offspring solutions

for i = 1:2:NP                  % selecting parents in pairs for crossover
    r = rand;                   % generate random number to decise if crossover is to occur
    
    if r < Pc                   % checking for crossover probability
        for j = 1:Dim
            r = rand;           % generating random number to detrmiine the beta value
            
            if r<= 0.5
                beta = (2*r)^(1/(etac+1));    % calculate beta value
            else
                beta = (1/(2*(1-r)))^(1/(etac+1));   % calculate beta value
            end
            
            offspring(i,j) = 0.5*(((1+beta)*Parent(i,j))+(1-beta)*Parent(i+1,j));       % generate first offspring solution
            offspring(i+1,j) = 0.5*(((1-beta)*Parent(i,j))+(1+beta)*Parent(i+1,j));     % generate second offspring solution
        end
        
        offspring(i,:) = max(offspring(i,:),lb);     % bounding the violating variable to their lower bound
        offspring(i+1,:) = max(offspring(i+1,:),lb); % bounding the violating variable to their lower bound
        offspring(i,:) = min(offspring(i,:),ub);     % bounding the violating variable to their upper bound
        offspring(i+1,:) = min(offspring(i+1,:),ub); % bounding the violating variable to their upper bound
    else
        offspring(i,:) = Parent(i+1,:);       % copy first parent as first offspring soultion
        offspring(i+1,:) = Parent(i+1,:);     % copy second parent as second offspring soultion
    end
    
end
end

% function for Mutation
function offspring = MutationPloy(offspring,Pm,etam,lb,ub)

[Np,Dim] = size(offspring);

for i = 1:Np
    r=rand;
    if r < Pm                   % checking the mutation probability
        for j = 1:Dim
            r = rand;           % generating random variable to determine each solution
            if r<0.5
                delta = (2*r)^(1/(etam+1)) - 1;                     % calculating delta value
            else
                delta = (1-2*(1-r))^(1/(etam+1)) - 1;               % calculating delta value
            end
            offspring(i,j) = offspring(i,j) + (ub(j)-lb(j))*delta;  % Mutating each variable of offspring solution
        end
        offspring(i,:) = max(offspring(i,:),lb);     % bounding the violating variable to their lower bound
        offspring(i,:) = min(offspring(i,:),ub);     % bounding the violating variable to their upper bound
    end
end
end