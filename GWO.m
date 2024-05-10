 %% GREY WOLF OPTIMIZATION (GWO)

function [bestfitIter,best_Wolf,bestWolf_Position] = GWO(ObjFun,LB,UB,nVar,nPop,MaxIt)
%% Problem Parameters
prob = ObjFun;                  % fitness function
Var = nVar;                     % number of variables
lb = LB;                        % variable Lower Bound
ub = UB;                        % variable Upper Bound
cPop = nPop;                    % class population size
Maxiter = MaxIt;                % max no of iteration

%% GWO Parameters
No_of_Wolf = cPop;
% a = 2;
t = 1; % iteration counter

%% GWO Main Loop
% Initaiization
% Alpha Wolf
alpha_pos = NaN(1,Var); alpha_score = inf;
% Beta Wolf
beta_pos = NaN(1,Var); beta_score = inf;
% Delta Wolf
delta_pos = NaN(1,Var); delta_score = inf;

f = NaN(cPop,1);                    % target vector for fitness fun value of the population
bestfitIter = NaN(Maxiter+1,1);     % vetcor to store best fitness fun value per iteration

% Initialize first population of Grey Wolf
P = Init_Position(No_of_Wolf,Var,ub,lb);

for q = 1:cPop
    f(q) = prob(P(q,:));        % evaluating the fitness fun value of initial population
end

bestfitIter(1) = min(f);            % vector to store best fitness value of the 0th iteration 
%% Calculate and upddate fitness of the three candidate wolves

while t < Maxiter
    for i = 1:No_of_Wolf % for each wolf
        
        % Calculate the fitness function values
        f = prob(P(i,:));
        
        %% update alpha, beta and delta
        % Check Alpha fitness function value
        if f < alpha_score
            alpha_pos = P(i,:);
            alpha_score = f;
        end
        % Check Beta fitness function value
        if f < beta_score
            beta_pos = P(i,:);
            beta_score = f;
        end
        % Check Delta fitness function value
        if f < delta_score
            delta_pos = P(i,:);
            delta_score = f;
        end
    end
    
    
    %% update wolf position
    a = 2-t*(2/Maxiter);
    for i =  1:No_of_Wolf
       for j = 1:Var
          r1 = rand(); r2 = rand();
          A = 2*a*(r1-a);
          C = 2*r2;
          D_alpha = abs(C*alpha_pos(j)-P(i,j));
          X1 = alpha_pos(j) - A* D_alpha;
          
          D_beta = abs(C*beta_pos(j)-P(i,j));
          X2 = beta_pos(j) - A* D_beta;
          
          D_delta = abs(C*delta_pos(j)-P(i,j));
          X3 = delta_pos(j) - A* D_delta;
          
          P(i,j) = (X1 + X2 + X3)/3;
       end
    end
    
    %% Check boundary position and update if necessary
    for i = 1:No_of_Wolf % for each wolf
        No_of_Boundary = size(ub,2);
        
        if No_of_Boundary == 1
           LB_vect = lb.*ones(1,Var);
           UB_vect = ub.*ones(1,Var);
           
        else 
            LB_vect = lb; UB_vect = ub;
        end
        for j = 1:Var
            if P(i,j) > UB_vect(j)
               P(i,j) = UB_vect(j);
            end
            if P(i,j) < LB_vect(j)
               P(i,j) = LB_vect(j);
            end
        end
    end
    
    % iteration
    t = t+1;
    bestfitIter(t+1) = alpha_score;
end
% Results
[best_Wolf,ind] = min(f);
bestWolf_Position = P(ind,:);
end


 
function Position = Init_Position(No_of_Wolf,nVar,UB,LB)

No_of_Boundary = size(UB,2);

if No_of_Boundary == 1
    Position = rand(No_of_Wolf,nVar).*(UB-LB) + LB;
else   % more than one boundary
    for i = 1:No_of_Wolf  % outloop for each Wolf
        for j = 1:nVar  % for each variable
            LB_j = LB(j);
            UB_j = UB(j);
            Position(i,j) = rand()*(UB_j-LB_j) + LB_j;
        end
    end
end


end

