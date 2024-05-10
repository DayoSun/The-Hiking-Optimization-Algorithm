%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% THE Hiking Optimizer Algorithm (HOA)                                                   %
% Authors: Sunday Oladejo, Stephen Ekwe, and Mirjalili Seyedali                                                  %
% Version: v2 - TOBLER's HIKING FUNCTION (THF) APPROACH                                     %
% Last Update: 2023-02-09      
% Developed in MATLAB R2023a based on the following paper:
%
%  The Hiking Optimization Algorithm: A novel human-based metaheuristic approach
%  https://doi.org/10.1016/j.knosys.2024.111880
%


%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Best] = HOA_v2(ObjFun,LB,UB,dim,hiker,MaxIter)
%% Problem Parameters
prob = ObjFun;                  % objective function
nVar = dim;                     % dimension of problem
lb = LB;                        % lower bound
ub = UB;                        % upper bound
nPop = hiker;                   % no. of hikers
MaxIt = MaxIter;                % max iteration

%% Pre-allocate
fit = zeros(nPop,1);                    % vectorize variable to store fitness value
Best.iteration = zeros(MaxIt+1,1);      % vectorize variable store fitness value per iteration

%% Start Tobler Hiking Function Optimizer
Pop = repmat(lb,nPop,1) + repmat((ub-lb),nPop,1).*rand(nPop,nVar); % generate initial position of hiker

%% Evaluate fitness
for q = 1:nPop
    fit(q) = prob(Pop(q,:));            % evaluating initial fitness 
end
Best.iteration(1) = min(fit);           % stores initial fitness value at 0th iteration


%% Main Loop
for i = 1:MaxIt
        [~,ind] = min(fit);             % obtain initial fitness
        Xbest = Pop(ind,:);             % obtain golbal best position of initial fitness
%         Lead_Hiker(i) = ind;
%         [~,isd] = max(fit); 
%         Sweaper_Hiker(i) = isd;
    for j = 1:nPop
        
        Xini = (Pop(j,:));              % obtain initial position of jth hiker
    
        theta = randi([0 50],1,1);      % randomize elevation angle of hiker        
        s = tan(theta);                 % compute slope
        SF = randi([1 2],1,1);          % sweep factor generate either 1 or 2 randomly 
        
        Vel = 6.*exp(-3.5.*abs(s+0.05)); % Compute walking velocity based on Tobler's Hiking Function
        newVel =  Vel + rand(1,nVar).*(Xbest - SF.*Xini) ; % determine new position of hiker
        
        
        newPop = Pop(j,:) + newVel;     % Update position of hiker 
        
        newPop = min(ub,newPop);        % bound violating to upper bound
        newPop = max(lb,newPop);        % bound violating to lower bound
        
        fnew = prob(newPop);            % re-evaluate fitness
        if (fnew < fit(j))              % apply greedy selection strategy
            Pop(j,:) = newPop;          % store best position
            fit(j) = fnew;              % store new fitness 
        end
    end
    
     Best.iteration(i+1) = min(fit); % store best fitness per iteration
    disp(['Iteration ' num2str(i+1) ': Best Cost = ' num2str(Best.iteration(i+1))]);
    
    
%     Best.iteration(i+1) = min(fit);     % store best fitness per iteration
%     % Show Iteration Information
%     disp(['Iteration ' num2str(i) ': Best Hike = ' num2str(min(fit))]);    
end
% THFO Solution: Store global best fitness and position
[Best.Hike,idx] = min(fit); 
Best.Position = Pop(idx,:);
end

