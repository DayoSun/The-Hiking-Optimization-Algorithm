%%% Quadratic Interpolation Optimization (QIO) for 23 functions %%%
%--------------------------------------------------------------------------%
% Quadratic Interpolation Optimization (QIO)                               %
% Source codes demo version 1.0                                            %
% The code is based on the following paper:                                %
% W. Zhao, L. Wang, Z. Zhang, S. Mirjalili, N. Khodadadi, Q. Ge, Quadratic % 
% Interpolation Optimization (QIO): A new optimization algorithm based on  % 
% generalized quadratic interpolation and its applications to real-world   % 
% engineering problems, Computer Methods in Applied Mechanics and          %
% Engineering (2023) 116446, https://doi.org/10.1016/j.cma.2023.116446.    %
%--------------------------------------------------------------------------%
function [BestX,BestF,HisBestF]=QIO(BenFunctions,Low,Up,Dim, MaxIt,nPop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FunIndex: Index of function.                       %
% MaxIt: Maximum number of iterations.               %
% nPop: Size of population.                          %
% IndividualsPos: Position of individual population  %
% IndividualsFit: Fitness of individual population   %
% Dim: Dimensionality of prloblem.                   %
% BestX: Best solution found so far.                 %
% BestF: Best fitness corresponding to BestX.        %
% HisBestF: History best fitness over iterations.    %
% Low: Low bound of search space.                    %
% Up: Up bound of search space.                      %
% w1: exploration weight                             %
% w2:exploitation weight                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%[Low,Up,Dim]=FunRange(FunIndex);
IndividualsPos=zeros(nPop,Dim);
IndividualsFit=zeros(nPop,1);
if length(Low)==1
    Low=ones(1,Dim)*Low;
    Up=ones(1,Dim)*Up;
end
for i=1:nPop
    IndividualsPos(i,:)=rand(1,Dim).*(Up-Low)+Low;
    IndividualsFit(i)=BenFunctions(IndividualsPos(i,:));
end

BestF=inf;
BestX=[];

for i=1:nPop
    if IndividualsFit(i)<=BestF
        BestF=IndividualsFit(i);
        BestX=IndividualsPos(i,:);
    end
end

HisBestF=zeros(MaxIt,1);

for It=1:MaxIt
    newIndividualPos=zeros(1,Dim);
    
    for i=1:nPop
        if rand>0.5
            K=[1:i-1 i+1:nPop];
            RandInd=randperm(nPop-1,3);            
            K1=K(RandInd(1));
            K2=K(RandInd(2));
            K3=K(RandInd(3));
            f1=IndividualsFit(i);
            f2=IndividualsFit(K1);
            f3=IndividualsFit(K2);
            for j=1:Dim
                x1=IndividualsPos(i,j);
                x2=IndividualsPos(K1,j);
                x3=IndividualsPos(K2,j);
                % Eq.(25)
                newIndividualPos(j)=GQI(x1,x2,x3,f1,f2,f3,Low(j),Up(j));
            end
            a=cos(pi/2*It/MaxIt);
            b=0.7*a+0.15*a*(cos(5*pi*It/MaxIt)+1);
            % Eq.(27)
            w1=3*b*randn;
            % Exploration, Eq.(26)
            newIndividualPos=newIndividualPos+w1.*(IndividualsPos(K3,:)-...
                newIndividualPos)+round(0.5*(0.05+rand))*(log(rand/(rand)));
        else
            K=[1:i-1 i+1:nPop];
            RandInd=randperm(nPop-1,2);
            K1=K(RandInd(1));
            K2=K(RandInd(2));
            f1=IndividualsFit(K1);
            f2=IndividualsFit(K2);
            f3=BestF;
            for j=1:Dim
                x1=IndividualsPos(K(RandInd(1)),j);
                x2=IndividualsPos(K(RandInd(2)),j);
                x3=BestX(j);
                %Eq.(31)
                newIndividualPos(j)=GQI(x1,x2,x3,f1,f2,f3,Low(j),Up(j));
            end
            %Eq.(32)
            w2=3*(1-(It-1)/MaxIt)*randn;
            rD=randi(Dim);
            %Exploitation, Eq.(30)
            
            newIndividualPos=newIndividualPos+w2*(BestX-round(1+rand)*...
                (Up-Low)/(Up(rD)-Low(rD))*IndividualsPos(i,rD));
            
        end
        newIndividualPos=SpaceBound(newIndividualPos,Up,Low);
        newIndividualFit=BenFunctions(newIndividualPos);
        %Eq.(33)
        if newIndividualFit<IndividualsFit(i)
            IndividualsFit(i)=newIndividualFit;
            IndividualsPos(i,:)=newIndividualPos;
        end
        
    end
    
    for i=1:nPop
        if IndividualsFit(i)<BestF
            BestF=IndividualsFit(i);
            BestX=IndividualsPos(i,:);
        end
    end
    
    HisBestF(It)=BestF;
    
end

