function [Best_F,Best_P,Conv_curve]=COVIDOA(objFun, LB, UB, nVar, nHiker, MaxItr)


CostFunction = objFun;
MaxIt = MaxItr;
minVal = LB;
maxVal = UB;
nPop = nHiker;
D = nVar;


MR=0.1;
shifttingNo=1;
numOfSubprotiens=2;


VarMin =ones(1,D).*minVal;
VarMax = ones(1,D).*maxVal;
gamma =.5;
beta=.5;
bestsol.Cost = inf;
bestsol.Position=zeros(1,D);
empty_individual.Position=[];
empty_individual.Cost=[];
pop=repmat(empty_individual,nPop,1);

for i = 1:nPop
    
    % Generate Random Solution
    pop(i).Position = unifrnd(VarMin, VarMax, [1 D]);
    
    % Evaluate Solution
    pop(i).Cost = CostFunction(pop(i).Position);
    
    % Compare Solution to Best Solution Ever Found
    if pop(i).Cost < bestsol.Cost
        bestsol.Cost = pop(i).Cost;
        bestsol.Position=pop(i).Position;
    end
end

Conv_curve = nan(MaxIt, 1);
%     bestposition=nan(MaxIt, 1);
%--------- main loop ------------------------------------

for it = 1:MaxIt
    %%--------------------------------------------------------
    %% virus replication phase through frameshifting technique
    %%--------------------------------------------------------
    %%---------- roulette wheel selection----------
    c = [pop.Cost];
    avgc = mean(c);
    if avgc ~= 0
        c = c/avgc;
    end
    probs = exp(-beta*c);
    x=[];
    for k=1:nPop
        parent = pop(RouletteWheelSelection(probs));
        parent.Position = max(parent.Position, VarMin);
        parent.Position = min(parent.Position, VarMax);
        %--------- frameshifting ------------------
        x(k,:)=parent.Position;
        for t=1:numOfSubprotiens
            for i=1:D-shifttingNo
                %  shifting by the value of shifttingNo
                x(k,i)=x(k,i+shifttingNo);
                
            end
            r=unifrnd(minVal, maxVal, [1 shifttingNo]);
            x(k,:)=[x(k,1:D-shifttingNo) r];
            subprotien(t,:)=x(k,:);
        end
        %%---- apply uniform crossover between proteins to generate new solution
        newvirus(k,:)=UniformCrossover(subprotien(1,:),subprotien(2,:),gamma);
        newvirus(k,:)=max(newvirus(k,:), VarMin);
        newvirus(k,:)=min(newvirus(k,:), VarMax);
    end
    % ======== calculating the cost of the new population
    for t=1:nPop
        
        childcost=CostFunction(newvirus(t,:));
        %        numofevaluaions=numofevaluaions+1;
        if childcost < bestsol.Cost
            bestsol.Position = newvirus(t,:);
            bestsol.Cost=childcost;
        end
    end
    
    newPop = repmat(empty_individual, nPop, 1);
    % Mutation
    for l=1:nPop
        for k = 1: D
            R = rand();
            if R < MR
                newvirus(l,k)=minVal+rand*(maxVal-minVal);
            end
            newPop(l).Position=newvirus(l,:);
        end
        newPop(l).Position = max(newPop(l).Position, VarMin);
        newPop(l).Position = min(newPop(l).Position, VarMax);
        newPop(l).Cost=CostFunction(newPop(l).Position);
        if newPop(l).Cost < bestsol.Cost
            bestsol.Position = newPop(l).Position;
            bestsol.Cost=newPop(l).Cost;
        end
        
    end
    
    pop = SortPopulation([pop; newPop]);
    
    % Remove Extra Individuals
    pop = pop(1:nPop);
    Conv_curve(it) = bestsol.Cost;
    Best_F= bestsol.Cost;
    Best_P=bestsol.Position;
    %bestposition(it) = bestsol.Position;
    
    % Display Itertion Information
    %disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(Conv_curve(it))]);
    
end