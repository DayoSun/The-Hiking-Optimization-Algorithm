%% ARTIFICIAL BEE COLONY OPTIMIZATION (ABC)

function [bestfitIter,best_Bee,bestBee_Position] = ABC(ObjFun,LB,UB,nVar,nPop,MaxIt)
%% Problem Parameters
prob = ObjFun;                  % fitness function
Var = nVar;                     % number of variables
lb = LB;                        % variable Lower Bound
ub = UB;                        % variable Upper Bound
cPop = nPop;                    % class population size
Maxiter = MaxIt;                % max no of iteration


%% ABC Parameters
limit = 5;                      % permissible number of failures

%% Starting the Artificial Bee Colony 
% Initaiization
f = NaN(cPop,1);       % population vetcor for objective fun value
fit = NaN(cPop,1);     % population vetcor for fitness fun value
trial = NaN(cPop,1);   % trial vector

bestfitIter = NaN(Maxiter+1,1);     % vetcor to store best fitness fun value per iteration

% length of decision variables
D = Var;
% generation of initial population (or food source)
P = repmat(lb,cPop,1) + repmat((ub-lb),cPop,1).*rand(cPop,D);

for q = 1:cPop
    f(q) = prob(P(q,:));   % evaluating the objective fun value
    fit(q) = CalFit(f(q));  % evaluating the fitness fun value
end

bestfitIter(1) = min(f);            % vector to store best fitness value of the 0th iteration 

[bestobj,ind] = min(f);   % determine and store the best objective fun value
bestSolve = P(ind,:);     % determine and store the best solution

%% ABC Main Loop
for t = 1:Maxiter
    
    %% Employed bee phase
    for i = 1:cPop
        [trial,P,fit,f] = GenNewSol(prob,lb,ub,cPop,i,P,fit,trial,f,D);
    end
    
    %% onlooker bee phase
    probability = 0.9 * (fit/max(fit))+0.1; % find probability based on fitness fun value
    m = 0; n = 1;
    
    while(m<cPop)
        if(rand<probability(n))
            [trial,P,fit,f] = GenNewSol(prob,lb,ub,cPop,n,P,fit,trial,f,D);
            m = m + 1;    % increase counter per onlooker bee
        end
        n = mod(n,cPop)+1;
    end
    
    % store the best solution
    [bestobj,ind]=min([f;bestobj]);
    CombinedSol = [P;bestSolve];
    bestSolve = CombinedSol(ind,:);
    
    %% Scout Bee Phase
    [val,ind]=max(trial);
    
    if(val>limit)
        trial(ind) = 0;      % reset thr trail value
        P(ind,:) = lb + (ub-lb).*rand(1,D);   % generate a random solution
        f(ind) = prob(P(ind,:));              % determine the objective fun value of new solution
        fit(ind) = CalFit(f(ind));            % determine the fitness fun value of new solution
    end
    
    % store the best solution
    [bestobj,ind]= min([f;bestobj]);
    CombinedSol = [P;bestSolve];
    bestSolve = CombinedSol(ind,:);
    
    bestfitIter(t+1) = bestobj;
end

% store the best solution
[best_Bee,ind]= min([f;bestobj]);
CombinedSol = [P;bestSolve];
bestBee_Position = CombinedSol(ind,:);
end

% function for Calculating fitness
function fit = CalFit(f)

if f>= 0
    fit = 1/(1+f);
else
    fit = 1+abs(f);
end

end

% function for generating new solution
function [trial,P,fit,f] = GenNewSol(prob,LB,UB,nPop,n,P,fit,trial,f,D)

j = randi(D,1);     % randomly select the variable that is to be changed
par = randi(nPop,1);  % randomly select the neighbour

while(par==n)         % ensuring the neigbour is not same as current solution
    par = randi(nPop,1);
end

Xnew = P(n,:);      % variable to generate a new population

Phi = -1+(1-(-1))*rand; % generate a random number between -1 and 1

Xnew(j) = P(n,j) + Phi*(P(n,j) - P(par,j));    % genrating a new solution
Xnew(j) = min(Xnew(j),UB(j));       % bounding the violating variable to their upper bound
Xnew(j) = max(Xnew(j),LB(j));       % bounding the violating variable to their lower bound

objNewSol = prob(Xnew);             % dtermine the new objective fun value
FitnessNewSol = CalFit(objNewSol);  % dtermine the new fitness fun value

if (FitnessNewSol > fit(n))
   P(n,:) = Xnew;                   % new solution enters pool of solutions
   fit(n,:) = FitnessNewSol;        % Update the fitness fun value
   f(n,:) = objNewSol;              % Update the onjective fun value
   trial(n) = 0;                    % reset trial to zero
else
    trial(n) = trial(n) + 1;        % increase the trial counter
end

end
