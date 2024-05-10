% Biogeography-based optimization (BBO) software for minimizing a general function
function [bestSpecie_Iter,bestSpecie_Value,bestSpecie_Position] = BBO(ObjFun,LB,UB,nVar,nPop,MaxIt)

%% Problem Parameters
prob = ObjFun;                  % fitness function
cVar = nVar;                    % number of variables
Var=[1 nVar];                   % Decision Variables Matrix Size
lb = LB;                        % Variable Lower Bound
ub = UB;                        % Variable Upper Bound
cPop = nPop;                    % Number of Habitats (Population Size)
Maxiter = MaxIt;                % max no of iteration

%% BBO Parameters
KeepRate=0.2;                   % Keep Rate
nKeep=round(KeepRate*cPop);     % Number of Kept Habitats
nNew=cPop-nKeep;                % Number of New Habitats
% Migration Rates
mu=linspace(1,0,cPop);          % Emmigration Rates
lambda=1-mu;                    % Immigration Rates
alpha=0.9;
pMutation=0.1;
sigma=0.02*(ub-lb);


%% Starting the Biogeography-based optimization Algorithim
% Empty Habitat
habitat.Position=[];
habitat.Cost=[];
% Create Habitats Array
P=repmat(habitat,cPop,1);
% Initialize Habitats
for i=1:cPop
    P(i).Position=unifrnd(lb,ub,Var);
    P(i).Cost=prob(P(i).Position);
end
% Sort Population
[~, SortOrder]=sort([P.Cost]);
P=P(SortOrder);
% Best Solution Ever Found
BestSol=P(1);
% Array to Hold Best Costs
bestSpecie_Iter=zeros(Maxiter,1);
% Store Initial Best Cost
bestSpecie_Iter(1)=BestSol.Cost;
%% BBO Main Loop
for it=1:Maxiter
    newpop=P;
    for i=1:cPop
        for k=1:cVar
            % Migration
            if rand<=lambda(i)
                % Emmigration Probabilities
                EP=mu;
                EP(i)=0;
                EP=EP/sum(EP);
                
                % Select Source Habitat
                j=RouletteWheelSelection(EP);
                
                % Migration
                newpop(i).Position(k)=P(i).Position(k) ...
                    +alpha*(P(j).Position(k)-P(i).Position(k));
            end
            
            % Mutation
            if rand<=pMutation
                newpop(i).Position(k)=min(newpop(i).Position(k)+sigma*randn);
            end
        end
        
        % Apply Lower and Upper Bound Limits
        newpop(i).Position = max(newpop(i).Position, lb);
        newpop(i).Position = min(newpop(i).Position, ub);
        
        % Evaluation
        newpop(i).Cost=prob(newpop(i).Position);
    end
    
    % Sort New Population
    [~, SortOrder]=sort([newpop.Cost]);
    newpop=newpop(SortOrder);
    
    % Select Next Iteration Population
    P=[P(1:nKeep)
         newpop(1:nNew)];
     
    % Sort Population
    [~, SortOrder]=sort([P.Cost]);
    P=P(SortOrder);
    
    % Update Best Solution Ever Found
    BestSol=P(1);

    % Store iteration
    bestSpecie_Iter(it+1)=BestSol.Cost;
    
    % % Show Iteration Information
    % disp(['Iteration ' num2str(it) ': bestSpecie_Iter = ' num2str(BestCost(it))]);
end

%% Results
% Store Best Cost Ever Found
bestSpecie_Value=BestSol.Cost;

% Store Best Position Ever Found
bestSpecie_Position=BestSol.Position;
    
% figure;
% %plot(BestCost,'LineWidth',2);
% semilogy(BestCost,'LineWidth',2);
% xlabel('Iteration');
% ylabel('Best Cost');
% grid on;

end
% RouletteWheel Selection
function j=RouletteWheelSelection(P)
    r=rand;
    C=cumsum(P);
    j=find(r<=C,1,'first');
end

