% Real-Coded Simulated Annealing for Discrete and Combinatorial Optimization Problems
function [bestSA_Iter,bestSA_Value,bestSA_Position] = RC_SA(ObjFun,LB,UB,nVar,nPop,MaxIt)

%% Problem Parameters
prob = ObjFun;                  % fitness function
cVar = nVar;                    % number of variables
Var=[1 cVar];                   % Decision Variables Matrix Size
lb = LB;                        % Variable Lower Bound
ub = UB;                        % Variable Upper Bound
cPop = nPop;                    % Population Size
Maxiter = MaxIt;                % max no of iteration

%% SA Parameters
MaxSubIt=20;                    % Maximum Number of Sub-iterations
T0=0.1;                         % Initial Temp.
alpha=0.99;                     % Temp. Reduction Rate
nMove=5;                        % Number of Neighbors per Individual
mu = 0.5;                       % Mutation Rate
sigma = 0.1*(ub-lb);            % Mutation Range (Standard Deviation)

%% Initialization
% Create Empty Structure for Individuals
empty_individual.Position=[];
empty_individual.Cost=[];
% Create Population Array
P=repmat(empty_individual,cPop,1);
% Initialize Best Solution
BestSol.Cost=inf;

% Initialize Population
for i=1:cPop    
    % Initialize Position
    P(i).Position=unifrnd(lb, ub, Var);    
    % Evaluation
    P(i).Cost=prob(P(i).Position);    
    % Update Best Solution
    if P(i).Cost<=BestSol.Cost
        BestSol=P(i);
    end    
end
% Array to Hold Best Cost Values
bestSA_Iter=zeros(Maxiter,1);
% Store Initial Best Cost Ever Found
bestSA_Iter(1)=BestSol.Cost;
% Intialize Temp.
T=T0;

%% SA Main Loop
for it=1:Maxiter    
    for subit=1:MaxSubIt        
        % Create and Evaluate New Solutions
        newP=repmat(empty_individual,cPop,nMove);
        for i=1:cPop
            for j=1:nMove                
                % Create Neighbor
                newP(i,j).Position=Mutate(P(i).Position,mu,sigma(1),lb(1),ub(1));                
                % Evaluation
                newP(i,j).Cost=prob(newP(i,j).Position);                
            end
        end
        newP=newP(:);        
        % Sort Neighbors
        [~, SortOrder]=sort([newP.Cost]);
        newP=newP(SortOrder);        
        for i=1:cPop            
            if newP(i).Cost<=P(i).Cost
                P(i)=newP(i);                
            else
                DELTA=(newP(i).Cost-P(i).Cost)/P(i).Cost;
                lambda=exp(-DELTA/T);
                if rand<=lambda
                    P(i)=newP(i);
                end
            end            
            % Update Best Solution Ever Found
            if P(i).Cost<=BestSol.Cost
                BestSol=P(i);
            end        
        end
    end    
    % Store Best Cost Ever Found
    bestSA_Iter(it+1)=BestSol.Cost;    
    % % Display Iteration Information
    % disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(bestSA_Iter(it+1))]);
    
    % Update Temp.
    T=alpha*T;    
    sigma = 0.98*sigma;    
end
%% Results
% Store Best Answer Ever Found
bestSA_Value=BestSol.Cost;
% Store Best Position Ever Found
bestSA_Position=BestSol.Position;

% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(bestSA_Iter,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost');
% grid on;
end
%% Mutate
function y=Mutate(x,mu,sigma,lb,ub)
    A=(rand(size(x))<=mu);
    J=find(A==1);
    y=x;
    y(J)=x(J)+sigma.*randn(size(J));
    % Clipping
    y=max(y,lb);
    y=min(y,ub);    
end