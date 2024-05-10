%% SANITIZED TEACHNING LEARNING BASED OPTIMIZATION (sTLBO)

function [bestfitIter,bestClass_Value,bestClass_Position,P,f] = sTLBO(ObjFun,LB,UB,nVar,nPop,MaxIt)
%% Problem Parameters
prob = ObjFun;                  % fitness function
Var = nVar;                     % number of variables
lb = LB;                        % Variable Lower Bound
ub = UB;                        % Variable Upper Bound
cPop = nPop;                    % class population size
Maxiter = MaxIt;                % max no of iteration

%% Starting the Teaching learning Optimization
f = NaN(cPop,1);        % vector for fitness fun value of the population
bestfitIter = NaN(Maxiter+1,1);   % vetcor to store best fitness fun value per iteration
D = Var;                % length of decision variables
P = repmat(lb,cPop,1) + repmat((ub-lb),cPop,1).*rand(cPop,D);   % generation of initial population of the class

for q = 1:cPop
    f(q) = prob(P(q,:));   % evaluating the fitness fun value of initial population
end

bestfitIter(1) = min(f);    % stores the best fitness value of the 0th iteration (inital class population)
%% Iteration Loop (Main Loop)
for t = 1:Maxiter
    for i = 1:cPop
        %% Teacher Phase
        Xmean = mean(P);     % determine the mean of the population
        [~,ind] = min(f);    % determine the location of the teacher
        Xbest = P(ind,:);    % copying the solution acting as teacher
        
        TF = randi([1 2],1,1);  % generating either 1 or 2 randomly for teaching factor
        
        Xnew = P(i,:) + rand(1,D).*(Xbest - TF*Xmean);   % generating the new solution
        
        Xnew = min(ub,Xnew);    % bounding the violating to their upper bound
        Xnew = max(lb,Xnew);    % bounding the violating to their lower bound
        
        fnew = prob(Xnew);      % evaluating the fitness of the newly generated solution
        
        if (fnew < f(i))         % greedy selection
            P(i,:) = Xnew;       % include the new solution of population
            f(i) = fnew;         % include the fitness fun value of new solution in population
        end
        
        %% Learner Phase
        part = randi([1 cPop],1,1);     % selection of random parter
        
        % Ensuring that current member is not the partner
        while i == part
            part = randi([1 cPop],1,1);     % selection of random parter
        end
        
        if f(i) < f(part)       % select appropriate euqation to be used in Learner Phase
            Xnew = P(i,:) + rand(1,D).*(P(i,:) - P(part,:));    % generate the new solution
        else
            Xnew = P(i,:) - rand(1,D).*(P(i,:) - P(part,:));    % generate the new solution
        end
        
        Xnew = min(ub,Xnew);     % bounding the violating variable to their upper bound
        Xnew = max(lb,Xnew);     % bounding the violating variable to their lower bound
        
        fnew = prob(Xnew);       % evaluating the fitness of the newly generated solution
        
        if (fnew < f(i))         % greedy selection
            P(i,:) = Xnew;       % include the new solution of population
            f(i) = fnew;         % include the fitness fun value of new solution in population
        end
        
    end
    bestfitIter(t+1) = min(f);
end

% Results
[bestClass_Value,ind] = min(f);
bestClass_Position = P(ind,:);

