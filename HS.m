% Harmony Search Algorithm for Discrete and Combinatorial Optimization Problems
function [bestHS_Iter,bestHS_Value,bestHS_Position] = HS(ObjFun,LB,UB,nVar,nPop,MaxIt)

%% Problem Parameters
prob = ObjFun;                  % fitness function
cVar = nVar;                    % number of variables
Var=[1 cVar];                   % Decision Variables Matrix Size
lb = LB;                        % Variable Lower Bound
ub = UB;                        % Variable Upper Bound
cPop = nPop;                    % Number of Habitats (Population Size)
Maxiter = MaxIt;                % max no of iteration

%% Harmony Search Parametters
nNew=20;                        % Number of New Harmonies
HMCR=0.9;                       % Harmony Memory Consideration Rate
PAR=0.1;                        % Pitch Adjustment Rate
FW=0.02*(ub(1)-lb(1));          % Fret Width (Bandwidth)
FW_damp=0.995;                  % Fret Width Damp Ratio

%% Initialization
% Empty Harmony Structure
empty_harmony.Position=[];
empty_harmony.Cost=[];
% Initialize Harmony Memory
HM=repmat(empty_harmony,cPop,1);
% Create Initial Harmonies
for i=1:cPop
    HM(i).Position=unifrnd(lb,ub,Var);
    HM(i).Cost=prob(HM(i).Position);
end
% Sort Harmony Memory
[~, SortOrder]=sort([HM.Cost]);
HM=HM(SortOrder);
% Update Best Solution Ever Found
BestSol=HM(1);
% Array to Hold Best Cost Values
bestHS_Iter=zeros(Maxiter,1);
%% Store Initial Best Answer
bestHS_Iter(1)=BestSol.Cost;

%% Harmony Search Main Loop
for it=1:Maxiter    
    % Initialize Array for New Harmonies
    NEW=repmat(empty_harmony,nNew,1);    
    % Create New Harmonies
    for k=1:nNew        
        % Create New Harmony Position
        NEW(k).Position=unifrnd(lb,ub,Var);
        for j=1:nVar
            if rand<=HMCR
                % Use Harmony Memory
                i=randi([1 cPop]);
                NEW(k).Position(j)=HM(i).Position(j);
            end            
            % Pitch Adjustment
            if rand<=PAR
                %DELTA=FW*unifrnd(-1,+1);    % Uniform
                DELTA=FW*randn();            % Gaussian (Normal) 
                NEW(k).Position(j)=NEW(k).Position(j)+DELTA;
            end        
        end        
        % Apply Variable Limits
        NEW(k).Position=max(NEW(k).Position,lb);
        NEW(k).Position=min(NEW(k).Position,ub);
        % Evaluation
        NEW(k).Cost=prob(NEW(k).Position);        
    end    
    % Merge Harmony Memory and New Harmonies
    HM=[HM;NEW]; %#ok    
    % Sort Harmony Memory
    [~, SortOrder]=sort([HM.Cost]);
    HM=HM(SortOrder);    
    % Truncate Extra Harmonies
    HM=HM(1:cPop);    
    % Update Best Solution Ever Found
    BestSol=HM(1);    
    % Store Best Cost Ever Found
    bestHS_Iter(it+1)=BestSol.Cost;    
    % % Show Iteration Information
    % disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(bestHS_Iter(it+1))]); 
    
    % Damp Fret Width
    FW=FW*FW_damp;    
end
%% Results
% Store Best Answer Ever Found
bestHS_Value=BestSol.Cost;
% Store Best Position Ever Found
bestHS_Position=BestSol.Position;
end