%=================================================================================
% Komodo Mlipir Algorithm (KMA) version 1.0.0
%
% Developed in MATLAB R2017a based on the following paper:
%
% S. Suyanto, A.A. Ariyanto, A.F. Ariyanto, Komodo Mlipir Algorithm 
% Applied Soft Computing, 2021, 108043, https://doi.org/10.1016/j.asoc.2021.108043
%
% Main function with 2D Visualization
% 
% Updated 23 November 2021 by Prof. Dr. Suyanto, S.T., M.Sc.
% School of Computing, Telkom University
% Jl. Telekomunikasi No 1 Terusan Buah Batu, Bandung 40257, Indonesia
% Email: suyanto@telkomuniversity.ac.id 
% Website: https://suyanto.staff.telkomuniversity.ac.id
%=================================================================================

function [BestIndiv,OptVal,NumEva,fopt,fmean,ProcTime,EvoPopSize] = KMA2D(objFun, LB, UB, nVar, nHiker, MaxItr);



global Nvar FthresholdFX MaxNumEva
global PopSize MinAdaPopSize MaxAdaPopSize
global BigMales BigMalesFX
global Female FemaleFX
global SmallMales SmallMalesFX
global AllHQ AllHQFX
global MlipirRate MutRate MutRadius
global OneElitFX
global Ra
global Rb



Nvar = nVar;        % Dimension can scaled up to thousands for the functions F1-F13, but it is fixed for F14-F23
FthresholdFX= 0;
MaxNumEva = MaxItr;       % Maximum number of evaluations
PopSize = nHiker;   % Population Size (number of Komodo individuals)
Evaluation = objFun;                  % fitness function
Ra = UB;
Rb = LB;

%% Generate the initial population of PopSize Komodo individuals
Pop = PopConsInitialization(PopSize);

%% Calculate all fitness values of the PopSize Komodo individuals
FX = zeros(1,PopSize);
for ii=1:PopSize
    FX(ii) = Evaluation(Pop(ii,:));
end

%% Sort the fitness values of the initial population
[SortedFX,IndFX] = sort(FX, 'ascend');
FX  = SortedFX;
Pop = Pop(IndFX,:);
OneElitFX  = FX(1);   % The best-so-far fitness to check improvement or stagnation 

%% Setting of parameters
MaxAdaPopSize = PopSize * 40;      % Maximum Adaptive Population Size = 5 * 40 = 200
MinAdaPopSize = PopSize * 4;       % Minimum Adaptive PopSize
NumBM         = floor(PopSize/2);  % Number of Big Males
MlipirRate    = (Nvar-1)/Nvar;     % Mlipir Rate is fixed during the first stage of evolution
MutRate       = 0.5;               % Probability of the Female mutation
MutRadius     = 0.5;               % Radius to limit the step of the Female mutation

%--------------------------------------------------------------------------
% First Stage: examining if the benchmark function is simple or complex
%--------------------------------------------------------------------------
IsGlobal    = 0;    % Boolean to check if the global optimum is reached
ImproveRate = 0;    % Improve rate to examine the benchmark function
NumEva      = 0;    % Number of evalutions
Gen         = 0;    % Generation
MaxGenExam1 = 100;  % Maximum generation of the first examination
MaxGenExam2 = 1000; % Maximum generation of the second examination

fopt       = [];    % Best-so-far fitness value in each generation
fmean      = [];    % Mean fitness value in each generation
EvoPopSize = [];    % Population Size used in each generation
GenImprove = 0;     % Generation counter to check the improvement condition

while Gen < MaxGenExam2
    Gen    = Gen + 1;           % Increament of generation
    NumEva = NumEva + PopSize;  % Number of evalutaions
   
    BigMales      = Pop(1:NumBM,:);      % Big Males
    BigMalesFX    = FX(1:NumBM);         % Fitness of Big Males
    Female        = Pop(NumBM+1,:);      % Females
    FemaleFX      = FX(NumBM+1);         % Fitness of Females
    SmallMales    = Pop(NumBM+2:end,:);  % Small Males
    SmallMalesFX  = FX(NumBM+2:end);     % Fitness of Small Males
    
    MoveBigMalesFemaleFirstStage;   % Move the BigMales and Female as well in the first stage
    MoveSmallMalesFirstStage;       % Move (Mlipir) SmallMales in the first stage
    
    Pop = [BigMales ; Female; SmallMales];     % New population after all movements
    FX  = [BigMalesFX FemaleFX SmallMalesFX];  % New fitness after all movements
    
    [SortedFX,IndFX] = sort(FX, 'ascend');     % Sort the new fitness
    FX               = SortedFX;               % Sorted FX
    Pop              = Pop(IndFX,:);           % Sorted Population
    BestIndiv        = Pop(1,:);               % Best (Komodo) Individual
    OptVal           = FX(1);                  % Optimum Value of the given problem

    fopt       = [fopt OptVal];         % Best-so-far fitness value in each generation
    fmean      = [fmean mean(FX)];      % Mean fitness value in each generation
    EvoPopSize = [EvoPopSize PopSize];  % Population Size used in each generation
    
    if OptVal < OneElitFX
        GenImprove  = GenImprove + 1;   % Increament of the GenImprove
        ImproveRate = GenImprove/Gen;   % Calculate the ImproveRate
        OneElitFX   = OptVal;           % Update the OneElitFX
    end
    
    if OptVal <= FthresholdFX  % If the optimum value is equal to or lower than the global optimum
        IsGlobal   = 1;        % The global optimum is reached
        break;                 % Stop the first stage of evolution
    end
   
    if Gen == MaxGenExam1      % If generation is equal to the maximum generation of the first examination
        if ImproveRate < 0.5   % If the Improve Rate is lower than 0.5
            IsGlobal   = 0;    % The global optimum is not reached 
            break;             % Stop the first stage of evolution and run the second stage of evolution
        end
    end
end  % while Gen < MaxGenExam2

%--------------------------------------------------------------------------
% Second Stage: solving the complex benchmark functions using
% self-adaptation of population size
%--------------------------------------------------------------------------
if ~IsGlobal && NumEva <= MaxNumEva
    FirstStagePop = Pop;  % Keep five individuals in the population of the frst stage of evolution
    FirstStageFX  = FX;   % Keep five fitness values in the population of the frst stage of evolution
    SwarmSize     = size(FirstStagePop,1);  % Micro-swarm size: 2 BM, 1 Female, 2 SM
    NumBM         = floor(SwarmSize/2);     % Number of Big Males = floor(5/2) = 2

    IncAdaPopSize = SwarmSize;  % increment of AdaPopSize
    DecAdaPopSize = SwarmSize;  % derement of AdaPopSize
    
    MlipirRate    = 0.5;  % Mlipir Rate is fixed during the second stage of evolution
    MaxGenImprove = 2;    % Maximum Generation to be defined as an improve condition
    MaxGenStagnan = 2;    % Maximum Generation to be defined as a stagnation condition
    
    GenImprove = 0;       % Generation counter to check the improvement condition
    GenStagnan = 0;       % Generation counter to check the stagnation condition
    
    ConsPop   = PopConsInitialization(MaxAdaPopSize-SwarmSize);  % Constrained initialization of 195 individuals
    ConsPopFX = zeros(1,size(ConsPop,1));                        % Initialize fitness values
    for ii=1:size(ConsPop,1)
        ConsPopFX(ii) = Evaluation(ConsPop(ii,:));   % Calculate the fitness value of each individual
    end
    
    Pop           = [FirstStagePop; ConsPop];  % Combine the populations of the first and second stage
    PopSize       = size(Pop,1);               % PopSize is now 5 + 195 = 200 individuals
    FX            = [FirstStageFX ConsPopFX];  % Combine the first and second FX
    OneElitFX     = min(FX);                   % The best-so-far fitness to check improvement or stagnation 
    
    while NumEva < MaxNumEva
        AdaPopSize = size(Pop,1);  % Adaptive PopSize is equal to MaxAdaPopSize = 200 individuals
        AllHQ      = [];           % All High Quality (Big Males)
        AllHQFX    = [];           % All High Quality Fitness
        for ind=1:SwarmSize:AdaPopSize
            MicroSwarm    = Pop(ind:ind+SwarmSize-1,:);
            MicroSwarmFX  = FX(ind:ind+SwarmSize-1);
            [~,IndFX]     = sort(MicroSwarmFX, 'ascend');
            MicroSwarm    = MicroSwarm(IndFX,:);
            MicroSwarmFX  = MicroSwarmFX(IndFX);
            AllHQ         = [AllHQ ; MicroSwarm(1:NumBM,:)];  % All high-quality Big Males
            AllHQFX       = [AllHQFX MicroSwarmFX(1:NumBM)];  % Fitness Values
        end
        for ind=1:SwarmSize:AdaPopSize
            MicroSwarm   = Pop(ind:ind+SwarmSize-1,:);
            MicroSwarmFX = FX(ind:ind+SwarmSize-1);
            [~,IndFX]    = sort(MicroSwarmFX, 'ascend');
            MicroSwarm   = MicroSwarm(IndFX,:);
            MicroSwarmFX = MicroSwarmFX(IndFX);
            BigMales     = MicroSwarm(1:NumBM,:);      % BigMales
            BigMalesFX   = MicroSwarmFX(1:NumBM);      % Fitness of BigMales
            Female       = MicroSwarm(NumBM+1,:);      % Female
            FemaleFX     = MicroSwarmFX(NumBM+1);      % Fitness of Female
            SmallMales   = MicroSwarm(NumBM+2:end,:);  % SmallMales
            SmallMalesFX = MicroSwarmFX(NumBM+2:end);  % Fitness of SmallMales
                        
            MoveBigMalesFemaleSecondStage;  % Move the BigMales and Female as well in the second stage
            MoveSmallMalesSecondStage;      % Move (Mlipir) the SmallMales in the second stage
            
            AllHQ(ind:ind+NumBM-1,:) = BigMales;
            AllHQFX(ind:ind+NumBM-1) = BigMalesFX;
            
            % Resulted new population
            Pop(ind:ind+SwarmSize-1,:) = [BigMales ; Female; SmallMales];
            FX(ind:ind+SwarmSize-1)    = [BigMalesFX FemaleFX SmallMalesFX];
            
            NumEva = NumEva + SwarmSize;   % Number of evalutaions
            
            [~,IndMin] = min(FX);
            OptVal     = FX(IndMin);       % Optimum Value of the given problem
            if OptVal <= FthresholdFX      % Optimum Value reaches the global optimum
                break;
            end
        end
        
        % Random population
        Ind = randperm(length(FX));
        Pop = Pop(Ind,:);
        FX  = FX(Ind);
        
        [~,IndMin] = min(FX);        % Find the minimum (best) value of the fitness FX
        BestIndiv  = Pop(IndMin,:);  % Best Individual (Komodo)
        OptVal     = FX(IndMin);     % The optimum value of the given problem
        fopt       = [fopt OptVal];    % Best-so-far fitness value in each generation
        fmean      = [fmean mean(FX)]; % Mean fitness value in each generation
        
        % if RoundUp(OptVal) <= FthresholdFX
        if OptVal <= FthresholdFX
            break;
        end
        
        %----------------------------------------------------------------------
        % Self-adaptation of population size
        %----------------------------------------------------------------------
        if OptVal < OneElitFX   % consecutive fitness show an improvement
            GenImprove = GenImprove + 1;
            GenStagnan = 0;
            OneElitFX  = OptVal;
        else                    % consecutive fitness values show a stagnation
            GenStagnan = GenStagnan + 1;
            GenImprove = 0;
        end
        
        % If consecutive fitness values show an improvement
        if GenImprove > MaxGenImprove
            AdaPopSize = AdaPopSize - DecAdaPopSize;
            if AdaPopSize < MinAdaPopSize
                AdaPopSize = MinAdaPopSize;
            end            
            [SortedFX,IndFX]  = sort(FX, 'ascend');
            SortedPop         = Pop(IndFX,:);
            Pop               = SortedPop(1:AdaPopSize,:);
            FX                = SortedFX(1:AdaPopSize);
            GenImprove        = 0;
        end
        
        % If consecutive fitness values show a stagnation
        if (GenStagnan > MaxGenStagnan)
            AdaPopSizeOld = size(Pop,1);
            AdaPopSize    = AdaPopSize + IncAdaPopSize;
            NumAddPop     = AdaPopSize - AdaPopSizeOld;
            if AdaPopSize > MaxAdaPopSize
                AdaPopSize = AdaPopSizeOld;
                NumAddPop  = 0;  % No individual added into the population
            end
            if AdaPopSize > AdaPopSizeOld
                NewPop   = zeros(NumAddPop,Nvar);
                NewPopFX = zeros(1,NumAddPop);
                for nn=1:NumAddPop
                    [NewPop(nn,:), NewPopFX(nn)] = AddingPop(BestIndiv);
                end
                Pop(end+1:end+size(NewPop,1),:) = NewPop;
                FX(end+1:end+size(NewPop,1))    = NewPopFX;
                NumEva = NumEva + NumAddPop;   % Number of evalutaions
            else  % NumAddPop == 0; % No adding population
                for nn=1:size(Pop,1)
                    [Pop(nn,:), FX(nn)] = Reposition(Pop(nn,:),FX(nn));
                end
                NumEva = NumEva + size(Pop,1);
            end
            GenStagnan = 0;
        end
       
        RandInd = randperm(size(Pop,1));
        FX  = FX(RandInd);
        Pop = Pop(RandInd,:);
        
        EvoPopSize = [EvoPopSize AdaPopSize];  % Population size in each generation
        Gen = Gen + 1;                         % Increament of generation
    end % while Evolution
end % if ~IsGlobal

ProcTime = toc;


%% Constrained Initialization of Population
function X = PopConsInitialization(PS)
% PS : Population Size
% NL : Number of Locations (at the four corners of the problem landscape)

global Nvar Ra Rb

F1   = [0.01 0.01 0.99 0.99];
F2   = [0.01 0.99 0.01 0.99];
X    = zeros(PS, Nvar);
IndX = 1;
for nn=1:4:PS
    if PS - nn >= 4
        NL = 4;
    else
        NL = PS - nn + 1;
    end
    ss = 1;
    while ss <= NL
        Temp = zeros(1,Nvar);
        for i=1:floor(Nvar/2)
            Temp(:,i) = Rb(i)+(Ra(i)-Rb(i))*(F1(ss)+((rand*2)-1)*0.01);
        end
        for i=floor(Nvar/2)+1:Nvar
            Temp(:,i) = Rb(i)+(Ra(i)-Rb(i))*(F2(ss)+((rand*2)-1)*0.01);
        end        
        X(IndX,:) = Temp;
        IndX = IndX + 1;
        ss = ss + 1;
    end
end


%% Move Big Males and Female in the first stage. The winner mates Female (if the Female wants)
function MoveBigMalesFemaleFirstStage

global BigMales BigMalesFX
global Female FemaleFX
global Nvar

HQ   = BigMales;
HQFX = BigMalesFX;

TempSM   = BigMales;
TempSMFX = BigMalesFX;

for ss=1:size(TempSM,1)
    MaxFolHQ = randi(2);
    VM = zeros(1,Nvar);            % Velocity of a strong male
    RHQ = randperm(size(HQ,1));
    FolHQ = 0;
    for fs=1:length(RHQ)
        ind = RHQ(fs);
        if ind~= ss
            % Semi randomly select sn individual to define an attraction or a distraction
            if HQFX(ind) < TempSMFX(ss) || rand < 0.5
                VM = VM + rand * (HQ(ind,:) - TempSM(ss,:));
            else
                VM = VM + rand * (TempSM(ss,:) - HQ(ind,:));
            end
        end
        FolHQ = FolHQ + 1;
        if FolHQ >= MaxFolHQ
            break
        end
    end
    NewBM = TempSM(ss,:) + VM;         % New Big Males
    NewBM = trimr(NewBM);              % Limit the values into the given dimensional boundaries
    TempSM(ss,:) = NewBM;
    TempSMFX(ss) = Evaluation(NewBM);
end

% Replace the Big Males with the best ones
[BigMales, BigMalesFX] = Replacement(BigMales, BigMalesFX, TempSM, TempSMFX);

WinnerBM  = BigMales(1,:);
WinnerFX  = BigMalesFX(1);

if WinnerFX < FemaleFX || rand < 0.5   % Sexual reproduction
    OffSprings = CrossOver(WinnerBM,Female);
    fx1 = Evaluation(OffSprings(1,:));
    fx2 = Evaluation(OffSprings(2,:));
    
    % Keep the best position of female
    if fx1 < fx2
        if fx1 < FemaleFX
            Female   = OffSprings(1,:);
            FemaleFX = fx1;
        end
    else
        if fx2 < FemaleFX
            Female   = OffSprings(2,:);
            FemaleFX = fx2;
        end
    end
else % Asexual reproduction
    NewFemale = Mutation;
    fx = Evaluation(NewFemale);
    
    % Keep the best position of female
    if fx < FemaleFX
        Female    = NewFemale;
        FemaleFX  = fx;
    end
end


%% Move Big Males and Female in the second stage. The winner mates Female (if the Female wants)
function MoveBigMalesFemaleSecondStage

global BigMales BigMalesFX AllHQ AllHQFX
global Female FemaleFX
global Nvar

if ~isempty(AllHQ)
    GlobalHQ = [BigMales ; AllHQ];
    GlobalHQFX = [BigMalesFX AllHQFX];
else
    GlobalHQ = BigMales;
    GlobalHQFX = BigMalesFX;
end

TempSM   = BigMales;
TempSMFX = BigMalesFX;

for ss=1:size(TempSM,1)
    VM = zeros(1,Nvar);   % Velocity of a strong male
    RHQ = randperm(size(GlobalHQ,1));
    MaxFolHQ = randi(3);
    FolHQ = 0;
    for fs=1:length(RHQ)
        ind = RHQ(fs);
        if ind ~= ss
            % select randomly Individual to define an attraction
            if GlobalHQFX(ind) < TempSMFX(ss) || rand < 0.5
                VM = VM + rand * (GlobalHQ(ind,:) - TempSM(ss,:));
            else % or a distraction
                VM = VM + rand * (TempSM(ss,:) - GlobalHQ(ind,:));
            end
        end
        FolHQ = FolHQ + 1;
        if FolHQ >= MaxFolHQ
            break
        end
    end
    NewBM = TempSM(ss,:) + VM;
    NewBM = trimr(NewBM);
    TempSM(ss,:) = NewBM;
    TempSMFX(ss) = Evaluation(NewBM);
end
[BigMales, BigMalesFX] = Replacement(BigMales, BigMalesFX, TempSM, TempSMFX);

WinnerBM  = BigMales(1,:);
WinnerFX  = BigMalesFX(1);

if WinnerFX < FemaleFX || rand < 0.5   % Sexual reproduction
    OffSprings = CrossOver(WinnerBM,Female);
    fx1 = Evaluation(OffSprings(1,:));
    fx2 = Evaluation(OffSprings(2,:));
    
    % Keep the best position of female
    if fx1 < fx2
        if fx1 < FemaleFX
            Female   = OffSprings(1,:);
            FemaleFX = fx1;
        end
    else
        if fx2 < FemaleFX
            Female   = OffSprings(2,:);
            FemaleFX = fx2;
        end
    end
else % Asexual reproduction
    NewFemale = Mutation;
    fx = Evaluation(NewFemale);
    
    % Keep the best position of female
    if fx < FemaleFX
        Female    = NewFemale;
        FemaleFX  = fx;
    end
end


%% Move (Mlipir) Small Males aside with the maximum Mlipir Rate to do a high-exploitative low-explorative searching
function MoveSmallMalesFirstStage

global MlipirRate  % dimensional rate to "mlipir"
global BigMales
global SmallMales SmallMalesFX % WeakMalesDER
global Nvar

HQ       = BigMales;
TempWM   = SmallMales;
TempWMFX = SmallMalesFX;
MaxFolHQ = 1;

for ww=1:size(SmallMales,1)
    VMlipir = zeros(1,Nvar);          % vector of mlipir velocity
    RHQ = randperm(size(HQ,1));       % Random index of HQ
    FolHQ = 0;                        % Number of followed HQ
    for fs=1:length(RHQ)              % for each female and strong male
        ind = RHQ(fs);                % Index of followed HQ
        A = randperm(Nvar);           % get the attributes of movement
        D = round(MlipirRate*Nvar);   % dimensional size to "mlipir"
        if D >= Nvar
            D = Nvar-1;
        end
        if D < 1
            D = 1;
        end
        M = A(1:D);           % Moving the WM based on the D-attributes
        B = zeros(1,Nvar);    % Initialize Binary pattern
        B(M) = 1;             % Binary pattern
        VMlipir = VMlipir + rand(1,Nvar) .* (HQ(ind,:).*B)-(SmallMales(ww,:).*B);
        FolHQ = FolHQ + 1;
        if FolHQ >= MaxFolHQ
            break
        end
    end
    NewSM = SmallMales(ww,:) + VMlipir;
    NewSM = trimr(NewSM);
    TempWM(ww,:) = NewSM;
    TempWMFX(ww) = Evaluation(NewSM);    
end
SmallMales   = TempWM;
SmallMalesFX = TempWMFX;


%% Move (Mlipir) Small Males aside with MlipirRate = 0.5 to do a low-exploitative high-explorative searching
function MoveSmallMalesSecondStage

global MlipirRate
global BigMales AllHQ
global SmallMales SmallMalesFX
global Nvar

if ~isempty(AllHQ)
    HQ = [BigMales ; AllHQ];
else
    HQ = BigMales;
end

TempWM   = SmallMales;
TempWMFX = SmallMalesFX;

for ww=1:size(SmallMales,1)
    MaxFolHQ = randi(3);              % rand a number 1, 2, or 3
    VMlipir = zeros(1,Nvar);          % vector of mlipir velocity
    RHQ = randperm(size(HQ,1));       % Random index of HQ
    FolHQ = 0;                        % Number of followed HQ
    for fs=1:length(RHQ)              % for each female and strong male
        ind = RHQ(fs);                % Index of followed HQ
        A = randperm(Nvar);           % get the attributes of movement
        D = round(MlipirRate*Nvar);   % dimensional size to "mlipir"
        if D >= Nvar
            D = Nvar-1;
        end
        if D < 1
            D = 1;
        end
        M = A(1:D);           % Moving the WM based on the D-attributes
        B = zeros(1,Nvar);    % Initialize Binary pattern
        B(M) = 1;             % Binary pattern
        VMlipir = VMlipir + rand(1,Nvar) .* (HQ(ind,:).*B)-(SmallMales(ww,:).*B);
        FolHQ = FolHQ + 1;
        if FolHQ >= MaxFolHQ
            break
        end
    end
    NewSM = SmallMales(ww,:) + VMlipir;  
    NewSM = trimr(NewSM);                % Limit the values into the given dimensional boundaries
    TempWM(ww,:) = NewSM;
    TempWMFX(ww) = Evaluation(NewSM);    
end
SmallMales   = TempWM;
SmallMalesFX = TempWMFX;


%% Whole arithmetic crossover
function Offsprings = CrossOver(Parent1,Parent2)

global Nvar

Offsprings = zeros(2,Nvar);    % Initialize Offsprings
for ii=1:Nvar
    rval = rand;               % Generate random value in each dimension
    Offsprings(1,ii) = rval*Parent1(ii) + (1-rval)*Parent2(ii);
    Offsprings(2,ii) = rval*Parent2(ii) + (1-rval)*Parent1(ii);
end
Offsprings(1,:) = trimr(Offsprings(1,:));  % Limit the values into the given dimensional boundaries
Offsprings(2,:) = trimr(Offsprings(2,:));  % Limit the values into the given dimensional boundaries


%% Mutation of the only female
function NewFemale = Mutation

global Female
global Nvar Rb Ra

global MutRate MutRadius

NewFemale = Female;               % Initialize a new Female
MaxStep   = MutRadius*(Ra-Rb);    % Maximum step of the Female mutation
for ii=1:Nvar
    if rand < MutRate             % Check if a random value is lower than the Mutation Rate
        NewFemale(:,ii) = Female(:,ii) + (2*rand-1)*MaxStep(ii);
    end
end
NewFemale = trimr(NewFemale);     % Limit the values into the given dimensional boundaries


%% Adding an individual randomly
function [NewX, NewFX] = AddingPop(X)

global Nvar Rb Ra

NewX  = X + (0.05*levy(1,Nvar,1.5)) .* abs(Ra-Rb);
NewX  = trimr(NewX);       % Limit the values into the given dimensional boundaries
NewFX = Evaluation(NewX);  % Calculate the fitness value


%% Reposition an individual randomly
function [NewX, NewFX] = Reposition(X,FX)

global Nvar Rb Ra
global MutRate MutRadius

TempX   = X;
MaxStep = MutRadius*(Ra-Rb);    % Maximum step of the reposition
for ii=1:Nvar
    if rand < MutRate
        TempX(:,ii) = X(:,ii) + (2*rand-1)*MutRadius*MaxStep(ii);
    end
end

TempX  = trimr(TempX);       % Limit the values into the given dimensional boundaries
TempFX = Evaluation(TempX);  % Calculate the fitness value

if TempFX < FX               % TempFX is better than the original FX
    NewX  = TempX;
    NewFX = TempFX;
else                         % TempFX is worse than or equal to the original FX
    NewX  = X;
    NewFX = FX;
end


%% Replacement: sort the old and new populations and select the best ones
function [Z, FZ] = Replacement(X, FX, Y, FY)
% X  : old population of LX individuals
% FX : old fitness
% Y  : new population
% FY : new fitness
% Z  : survivor of LX individuals
% FZ : survivor fitness

LX   = size(X,1);  % Number of individuals in old population
XY   = [X ; Y];    % Joint individuals of old and new population
FXFY = [FX FY];    % Joint fitness values of old and new population

[SortedVal,SortedInd] = sort(FXFY, 'ascend');   % Sort all fitness values
Z    = XY(SortedInd(1:LX),:);                   % Select the best individuals
FZ   = SortedVal(1:LX);                         % Select the best fitness


%% Limit the values into the given dimensional boundaries
function [Z] = trimr(X)

global Ra Rb Nvar

for ii=1:Nvar
    X(X(:,ii)<Rb(ii),ii) = Rb(ii);
    X(X(:,ii)>Ra(ii),ii) = Ra(ii);
end
Z = X;