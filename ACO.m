% Ant Colony optimization (ACO) Algorithm for Discrete and Combinatorial Optimization Problems
function [bestAnt_Iter,bestAnt_Value,bestAnt_Position] = ACO(ObjFun,LB,UB,cVar,nPop,MaxIt)

%% Problem Parameters
prob = ObjFun;                  % fitness function
nVar = cVar;                    % number of variables
Var=[1 nVar];                   % Decision Variables Matrix Size
lb = LB;                        % Variable Lower Bound
ub = UB;                        % Variable Upper Bound
cPop=nPop;                      % Population Size (Archive Size)
Maxiter = MaxIt;                % max no of iteration

%% ACOR Parameters
nSample=40;                     % Sample Size
q=0.5;                          % Intensification Factor (Selection Pressure)
zeta=1;                         % Deviation-Distance Ratio

%% Starting the Ant Colony Optimization Algorithm
% Create Empty Individual Structure
empty_individual.Position=[];
empty_individual.Cost=[];
% Create Population Matrix
ant=repmat(empty_individual,cPop,1);
% Initialize Population Members
for i=1:cPop    
    % Create Random Solution
    ant(i).Position=unifrnd(lb,ub,Var);    
    % Evaluation
    ant(i).Cost=prob(ant(i).Position);    
end
% Sort Population
[~, SortOrder]=sort([ant.Cost]);
ant=ant(SortOrder);
% Update Best Solution Ever Found
BestSol=ant(1);
% Array to Hold Best Cost Values
bestAnt_Iter=zeros(Maxiter,1);
% Store Initial Best Cost
bestAnt_Iter(1)=BestSol.Cost;
% Solution Weights
w=1/(sqrt(2*pi)*q*cPop)*exp(-0.5*(((1:cPop)-1)/(q*cPop)).^2);
% Selection Probabilities
p=w/sum(w);

%% ACO Main Loop
for it=1:Maxiter    
    % Means
    s=zeros(cPop,nVar);
    for l=1:cPop
        s(l,:)=ant(l).Position;
    end    
    % Standard Deviations
    sigma=zeros(cPop,nVar);
    for l=1:cPop
        D=0;
        for r=1:cPop
            D=D+abs(s(l,:)-s(r,:));
        end
        sigma(l,:)=zeta*D/(cPop-1);
    end    
    % Create New Population Array
    newpop=repmat(empty_individual,nSample,1);
    for t=1:nSample        
        % Initialize Position Matrix
        newpop(t).Position=zeros(Var);        
        % Solution Construction
        for i=1:nVar            
            % Select Gaussian Kernel
            l=RouletteWheelSelection(p);            
            % Generate Gaussian Random Variable
            newpop(t).Position(i)=s(l,i)+sigma(l,i)*randn;            
        end  
        % bounding the violating variable to their upper & lower bound
        newpop(t).Position = min(ub,newpop(t).Position);     
        newpop(t).Position = max(lb,newpop(t).Position);     
        % Evaluation
        newpop(t).Cost=prob(newpop(t).Position);        
    end
    
    % Merge Main Population (Archive) and New Population (Samples)
    ant=[ant;newpop]; %#ok
    
    % Sort Population
    [~, SortOrder]=sort([ant.Cost],'ascend');
    ant=ant(SortOrder);
    
    % Delete Extra Members
    ant=ant(1:cPop);
    
    % Update Best Solution Ever Found
    BestSol=ant(1);
    
    % Store Best Cost
    bestAnt_Iter(it+1)=BestSol.Cost;
    
    % % Show Iteration Information
    % disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(bestAnt_Iter(it))]);    
end

%% Results
% Store Best Cost Ever Found
bestAnt_Value=BestSol.Cost;
% Store Best Position Ever Found
bestAnt_Position=BestSol.Position;

% figure;
% plot(BestCost,'LineWidth',2);
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