clear All; clc; close All
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% This main script is to test the metaheuristic algorithms with the
% unimodal, multimodal and composite benchmark functions and evaluate their parameters.

% To test the unimodal and multimodal functions ...uncomment lines 30, 31, 64 and
% comment lines   35 and 61
% For composite funcs. uncomment lines 35, 61 & comment lines 30, 31, 64
% Note that to use this script for other objective functions need a minor
% tweak


% BENCHMARK TEST & COMPOSITE FUNCTION
%% Functions

% Unimodal and multimodal functions
Function = {'F1';'F2';'F3';'F4';'F5';'F6';'F7';'F8';'F9';'F10';'F11';'F12';'F13';...
    'F14';'F15';'F16';'F17';'F18';'F19';'F20';'F21';'F22';'F23'};



%Function = {'F34';'F35';'F36';'F37';'F38';'F39'};    


Benchmark = {'Benchmark_500dim.csv'};   % name the accuracy metric  csv  according to the dimensions e.g Benchmark_100dim.csv  for 100 dimension
Compute = {'Compute_500dim.csv'};       % name the computational time metric  csv  according to the dimensions e.g Compute_100dim.csv  for 100 dimension
Rank = {'Rank_500dim.csv'};


Fnn = 1:length(Function);Fnn=Fnn';
%% Tuning Parameters
dim=5;        % dimension of problem
nHiker = 10;   % Number of search agents
run=2;                     % Independent run
MaxItr = 5;  % Maximum numbef of iterations
HH = multiwaitbar(3,[0 0 0 0],{'Please wait','Please wait','Please wait'});

%%

for h = 1:length(dim)
    nVar=dim(h);            % dimension of problem
    for i = 1:length(Function)
        FN= i;
        % Load details of the selected benchmark function
        Function_name=char(Function(i));

        % Composite functions
        %[LB,UB,nVar,fobj,xx,yy,~,xymin]=Get_Modal_Functions(Function_name,nVar,0);
        
        % unimodal and multimodal functions 
        [LB,UB,nVar,fobj,xx,yy,~,xymin]=Test_Functions(Function_name,nVar);
        
        
        lb= LB; ub = UB;  
       
        
        % Check boundary condition
        if length(LB) == 1
            LB = LB.*ones(1,nVar);          % Lower Bounds
            UB = UB.*ones(1,nVar);          % Upper Bounds
        else
            LB; UB;
        end
        % Objective function
        ObjFun = fobj;
        
        multiwaitbar(3,[abs(h/length(dim)), abs(i/length(Function)), abs(i/length(Function))],{'Dimension','Benchmark Function','Run'},HH);
        
        for j = 1:run
            
            %----------------------------Meta-Heauristic Algorithms-------------------------------
            % DSO-based Solution - Proposed Algorithm
             tic; [Best] = HOA_v2(ObjFun,LB,UB,nVar,nHiker,MaxItr);
            THFO_T= toc;
            % sTLBO-based Solution
            tic; [bestfitIter,bestClass_Value,bestClass_Position] = sTLBO(ObjFun,LB,UB,nVar,nHiker,MaxItr);
            TLBO_T= toc;
            % GA-based Solution
            tic; [bestfitGene,best_Gene,bestGene_Position] = GA(ObjFun,LB,UB,nVar,nHiker,MaxItr);
            GA_T= toc;
            % DE-based Solution
            tic; [bestfitPop,best_Pop,bestPop_Position] = DE(ObjFun,LB,UB,nVar,nHiker,MaxItr);
            DE_T= toc;
            % PSO-based Solution
            tic; [bestfit,best_Particle,bestParticle_Position] = PSO(ObjFun,LB,UB,nVar,nHiker,MaxItr);
            PSO_T= toc;
            % ABC-based Solution
            tic; [bestfitBee,best_Bee,bestBee_Position] = ABC(ObjFun,LB,UB,nVar,nHiker,MaxItr);
            ABC_T= toc;
            % GWO-based Solution
            tic; [bestfitWolf,best_Wolf,bestWolf_Position] = GWO(ObjFun,LB,UB,nVar,nHiker,MaxItr);
            GWO_T= toc;
            % SCA-based Solution
            tic; [bestDestination,Destination_position,bestfitDest]=SCA(nHiker,MaxItr,LB,UB,nVar,ObjFun);
            SCA_T= toc;
            % BBO-based Solution
            tic; [bestSpecie_Iter,bestSpecie_Value,bestSpecie_Position] = BBO(ObjFun,LB,UB,nVar,nHiker,MaxItr);
            BBO_T= toc;
            % ACO-based Solution
            tic; [bestAnt_Iter,bestAnt_Value,bestAnt_Position] = ACO(ObjFun,LB,UB,nVar,nHiker,MaxItr);
            ACO_T= toc;
            % RC-SA-based Solution
            tic; [bestSA_Iter,bestSA_Value,bestSA_Position] = RC_SA(ObjFun,LB,UB,nVar,nHiker,MaxItr);
            RC_T= toc;
            % HS-based Solution
            tic; [bestHS_Iter,bestHS_Value,bestHS_Position] = HS(ObjFun,LB,UB,nVar,nHiker,MaxItr);
            HS_T= toc;
            
            % COVIDOA-based Solution
            shifttingNo=1;
            numOfSubprotiens=2;
            MR=0.1;
            tic; %[Best_F,Best_P,Conv_curve]= COVIDOA(ObjFun, LB, UB, nVar, nHiker, MaxItr);
            [Best_F,Best_P,Conv_curve] = COVIDOA(nHiker,MaxItr,lb,ub,nVar,ObjFun,MR,shifttingNo,numOfSubprotiens);
            COVIDOA_T= toc;

            % kOMODO-based Solution
            tic; [BestIndiv,OptVal,NumEva,fopt,fmean,ProcTime,EvoPopSize]= KMA2D(ObjFun, LB, UB, nVar, nHiker, MaxItr, FN);
            KMA2D_T= toc;
            
            % QIO-based solution
            tic;[BestX,BestF,HisBestF]=QIO(ObjFun,lb,ub,nVar,MaxItr,nHiker);
            QIO_T = toc;


            % Obtain cost function
            H1(j,i)=Best.Hike; H2(j,i)=bestClass_Value; H3(j,i)=best_Gene; H4(j,i)=best_Pop;
            H5(j,i)=best_Particle; H6(j,i)=best_Bee; H7(j,i)=best_Wolf; H8(j,i)=bestDestination;
            H9(j,i)=bestSpecie_Value; H10(j,i)=bestAnt_Value; H11(j,i)=bestSA_Value; H12(j,i)=bestHS_Value;
            H13(j,i)=Best_F; H14(j,i)=OptVal; H15(j,i)=BestF;



            % obtain computational time
            C1(j,i)=THFO_T; C2(j,i)=TLBO_T; C3(j,i)=GA_T; C4(j,i)=DE_T; C5(j,i)=PSO_T; C6(j,i)=ABC_T; C7(j,i)=GWO_T; C8(j,i)=SCA_T;
            C9(j,i)=BBO_T; C10(j,i)=ACO_T; C11(j,i)=RC_T; C12(j,i)=HS_T; C13(j,i)=COVIDOA_T;  C14(j,i)=KMA2D_T; C15(j,i)=QIO_T;
            
            % obtain other data
            %P1(j,:)=bestsleepIter'; SAgent(j,:)=mean(Search_Agent); Hitr(:,:,:,j)=H_iter;  %C=permute(FgFit,[1 3 2]); C=reshape(C,[],size(FgFit,2),1);
            %Ah(j,:)=2.*(avgFit).*(1-(j/run));
            %E1(j,:)=mean(Mu_T').*2*(1-(j/run)); %#ok<*UDIM>
            
            % Exploration & expliotation
            %E1(j)=bestSleep_Value;
            %Eo(j)=2*bestSleep_Value*(1-(j/run)); Mu.one(j)=(mean(mean(Search_Agent)))-(j*((mean(mean(Search_Agent)))/run));
            %Mu.two(j)=(min(min(Search_Agent)))-(j*((min(min(Search_Agent)))/run)); Mu.three(j)=(max(max(Search_Agent)))-(j*((max(max(Search_Agent)))/run));
        end
        % obtain data for Wilcoxon rank-sum test with %5 significance
        W1(i)=ranksum((H1(:,i)),[H2(:,i);H3(:,i);H4(:,i);H5(:,i);H6(:,i);H7(:,i);H8(:,i);H9(:,i);H10(:,i);H11(:,i);H12(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W2(i)=ranksum((H2(:,i)),[H1(:,i);H3(:,i);H4(:,i);H5(:,i);H6(:,i);H7(:,i);H8(:,i);H9(:,i);H10(:,i);H11(:,i);H12(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W3(i)=ranksum((H3(:,i)),[H1(:,i);H2(:,i);H4(:,i);H5(:,i);H6(:,i);H7(:,i);H8(:,i);H9(:,i);H10(:,i);H11(:,i);H12(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W4(i)=ranksum((H4(:,i)),[H1(:,i);H2(:,i);H3(:,i);H5(:,i);H6(:,i);H7(:,i);H8(:,i);H9(:,i);H10(:,i);H11(:,i);H12(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W5(i)=ranksum((H5(:,i)),[H1(:,i);H2(:,i);H3(:,i);H4(:,i);H6(:,i);H7(:,i);H8(:,i);H9(:,i);H10(:,i);H11(:,i);H12(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W6(i)=ranksum((H6(:,i)),[H1(:,i);H2(:,i);H3(:,i);H4(:,i);H5(:,i);H7(:,i);H8(:,i);H9(:,i);H10(:,i);H11(:,i);H12(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W7(i)=ranksum((H7(:,i)),[H1(:,i);H2(:,i);H3(:,i);H4(:,i);H5(:,i);H6(:,i);H8(:,i);H9(:,i);H10(:,i);H11(:,i);H12(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W8(i)=ranksum((H8(:,i)),[H1(:,i);H2(:,i);H3(:,i);H4(:,i);H5(:,i);H6(:,i);H7(:,i);H9(:,i);H10(:,i);H11(:,i);H12(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W9(i)=ranksum((H9(:,i)),[H1(:,i);H2(:,i);H3(:,i);H4(:,i);H5(:,i);H6(:,i);H7(:,i);H8(:,i);H10(:,i);H11(:,i);H12(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W10(i)=ranksum((H10(:,i)),[H1(:,i);H2(:,i);H3(:,i);H4(:,i);H5(:,i);H6(:,i);H7(:,i);H8(:,i);H9(:,i);H11(:,i);H12(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W11(i)=ranksum((H11(:,i)),[H1(:,i);H2(:,i);H3(:,i);H4(:,i);H5(:,i);H6(:,i);H7(:,i);H8(:,i);H9(:,i);H10(:,i);H12(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W12(i)=ranksum((H12(:,i)),[H1(:,i);H2(:,i);H3(:,i);H4(:,i);H5(:,i);H6(:,i);H7(:,i);H8(:,i);H9(:,i);H10(:,i);H11(:,i); H13(:,i); H14(:,i); H15(:,i)]);
        W13(i)=ranksum((H13(:,i)),[H1(:,i);H2(:,i);H3(:,i);H4(:,i);H5(:,i);H6(:,i);H7(:,i);H8(:,i);H9(:,i);H10(:,i);H11(:,i); H12(:,i); H14(:,i); H15(:,i)]);
        W14(i)=ranksum((H14(:,i)),[H1(:,i);H2(:,i);H3(:,i);H4(:,i);H5(:,i);H6(:,i);H7(:,i);H8(:,i);H9(:,i);H10(:,i);H11(:,i); H12(:,i); H13(:,i); H15(:,i)]);
        W15(i)=ranksum((H15(:,i)),[H1(:,i);H2(:,i);H3(:,i);H4(:,i);H5(:,i);H6(:,i);H7(:,i);H8(:,i);H9(:,i);H10(:,i);H11(:,i); H12(:,i); H13(:,i); H14(:,i)]);


        %P1=mean(P1); best.sleepIter=P1'; AH=mean(Ah); figure; semilogy(E1(:,1));
        %plotFigs(Function_name,dim(end),i,Hitr,H_miniter,H_maxiter,ini_pop,AH,best,SAgent,xymin,H1,MaxItr,Eo,Mu,E1);
    end
    
    
    %------------------------------------------------ COST FUNCTION --------------------------------------------
    
    % Statistical Analysis (MEAN AND STD)
    %H1 = H1';H2 = H2';H3 = H3';H4 = H4';H5 = H5';H6 = H6';H7 = H7';H8 = H8';H9 = H9';H10 = H10';H11 = H11';H12 = H12';H13 = H13';
    THFO_mean = mean(H1); THFO_STD = std(H1); TLBO_mean = mean(H2); TLBO_STD = std(H2); GA_mean = mean(H3); GA_STD = std(H3);
    DE_mean = mean(H4); DE_STD = std(H4); PSO_mean = mean(H5); PSO_STD = std(H5); ABC_mean = mean(H6); ABC_STD = std(H6);
    GWO_mean = mean(H7); GWO_STD = std(H7); SCA_mean = mean(H8); SCA_STD = std(H8); BBO_mean = mean(H9); BBO_STD = std(H9);
    ACO_mean = mean(H10); ACO_STD = std(H10); RC_mean = mean(H11); RC_STD = std(H11); HS_mean = mean(H12); HS_STD = std(H12);
    COVIDOA_mean = mean(H13); COVIDOA_STD = std(H13); KMA2D_mean = mean(H14); KMA2D_STD = std(H14); QIO_mean = mean(H15); QIO_STD = std(H15);
    
    
    
    % Write AVG & STD of Algorithms to table 1
    M1 = [THFO_mean',THFO_STD']; M2 = [TLBO_mean',TLBO_STD']; M3 = [GA_mean',GA_STD']; M4 = [DE_mean',DE_STD']; %#ok<*SAGROW>
    M5 = [PSO_mean',PSO_STD']; M6 = [ABC_mean',ABC_STD']; M7 = [GWO_mean',GWO_STD']; M8 = [SCA_mean',SCA_STD'];
    M9 = [BBO_mean',BBO_STD']; M10 = [ACO_mean',ACO_STD']; M11 = [RC_mean',RC_STD']; M12 = [HS_mean',HS_STD'];
    M13 = [COVIDOA_mean',COVIDOA_STD']; M14 = [KMA2D_mean',KMA2D_STD']; M15 = [QIO_mean',QIO_STD'];
    
    Stat_table1 = table(Fnn,M1,M2,M3,M4,M5,M6,M7,M8,M9,M10,M11,M12,M13,M14,M15,...
        'VariableNames',{'F','HOA','TLBO','GA','DE','PSO','ABC','GWO','SCA','BBO','ACO','RC_SA','HS', 'COVIDOA','KMA2D','QIO'});
    Stat_table1(1:end,:)
    writetable(Stat_table1,char(Benchmark(h)),'WriteRowNames',true);
    
    %------------------------------------------------ COMPUTATION TIME --------------------------------------------
    
    % Statistical Analysis (MEAN AND STD)
    THFO_com = mean(C1); THFO_cSTD = std(C1); TLBO_com = mean(C2); TLBO_cSTD = std(C2); GA_com = mean(C3); GA_cSTD = std(C3);
    DE_com = mean(C4); DE_cSTD = std(C4); PSO_com = mean(C5); PSO_cSTD = std(C5); ABC_com = mean(C6); ABC_cSTD = std(C6);
    GWO_com = mean(C7); GWO_cSTD = std(C7); SCA_com = mean(C8); SCA_cSTD = std(C8); BBO_com = mean(C9); BBO_cSTD = std(C9);
    ACO_com = mean(C10); ACO_cSTD = std(C10); RC_com = mean(C11); RC_cSTD = std(C11); HS_com = mean(C12); HS_cSTD = std(C12);
    COVIDOA_com = mean(C13); COVIDOA_cSTD = std(C13); KMA2D_com = mean(C14); KMA2D_cSTD = std(C14); QIO_com = mean(C15); QIO_cSTD = std(C15);
    
    
    % Write AVG & STD of Algorithms to table 2
    D1 = [THFO_com',THFO_cSTD']; D2 = [TLBO_com',TLBO_cSTD']; D3 = [GA_com',GA_cSTD']; D4 = [DE_com',DE_cSTD']; %#ok<*SAGROW>
    D5 = [PSO_com',PSO_cSTD']; D6 = [ABC_com',ABC_cSTD']; D7 = [GWO_com',GWO_cSTD']; D8 = [SCA_com',SCA_cSTD'];
    D9 = [BBO_com',BBO_cSTD']; D10 = [ACO_com',ACO_cSTD']; D11 = [RC_com',RC_cSTD']; D12 = [HS_com',HS_cSTD'];
    D13 = [COVIDOA_com',COVIDOA_cSTD']; D14 = [KMA2D_com',KMA2D_cSTD']; D15 = [QIO_com',QIO_cSTD'];
    
    
    Stat_table2 = table(Fnn,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10,D11,D12,D13,D14,D15,...
        'VariableNames',{'F','HOA','TLBO','GA','DE','PSO','ABC','GWO','SCA','BBO','ACO','RC_SA','HS','COVIDOA','KMA2D', 'QIO'});
    Stat_table2(1:end,:)
    writetable(Stat_table2,char(Compute(h)),'WriteRowNames',true);
    
    %------------------------------------------------ Wilcoxon RANK-SUM --------------------------------------------
    
    % Write AVG & STD of Algorithms to table 2
    Stat_table3 = table(Fnn,W1',W2',W3',W4',W5',W6',W7',W8',W9',W10',W11',W12',W13',W14', W15',...
        'VariableNames',{'F','HOA','TLBO','GA','DE','PSO','ABC','GWO','SCA','BBO','ACO','RC_SA','HS', 'COVIDOA','KMA2D','QIO'});
    Stat_table3(1:end,:)
    writetable(Stat_table3,char(Rank(h)),'WriteRowNames',true);
    
    %------------------------------------------------ OTHER PLOTs --------------------------------------------
    
    
    % Compare normalized min and max values of benchmark test functions.....
% % % %     bar_mean1=[min(H1(1)),min(H2(1)),min(H3(1)),min(H4(1)),min(H5(1)),min(H6(1)),min(H7(1)),min(H8(1)),min(H9(1)),min(H10(1)),min(H11(1)),min(H12(1))];
% % % %     bar_mean2=[min(H1(4)),min(H2(4)),min(H3(4)),min(H4(4)),min(H5(4)),min(H6(4)),min(H7(4)),min(H8(4)),min(H9(4)),min(H10(4)),min(H11(4)),min(H12(4))];
% % % %     bar_mean3=[min(H1(9)),min(H2(9)),min(H3(9)),min(H4(9)),min(H5(9)),min(H6(9)),min(H7(9)),min(H8(9)),min(H9(9)),min(H10(9)),min(H11(9)),min(H12(9))];
% % % %     bar_mean4=[min(H1(10)),min(H2(10)),min(H3(10)),min(H4(10)),min(H5(10)),min(H6(10)),min(H7(10)),min(H8(10)),min(H9(10)),min(H10(10)),min(H11(10)),min(H12(10))];
% % % %     bar_mean5=[min(H1(14)),min(H2(14)),min(H3(14)),min(H4(14)),min(H5(14)),min(H6(14)),min(H7(14)),min(H8(14)),min(H9(14)),min(H10(14)),min(H11(14)),min(H12(14))];
% % % %     bar_mean6=[min(H1(16)),min(H2(16)),min(H3(16)),min(H4(16)),min(H5(16)),min(H6(16)),min(H7(16)),min(H8(16)),min(H9(16)),min(H10(16)),min(H11(16)),min(H12(16))];
% % % %     % plots
% % % %     figure('Position',[284   214   660   290])
% % % %     subplot(3, 6, 1:3,'align'); bar(bar_mean1); title('F1')
% % % %     ylabel('f(min)'); set(gca,'xticklabel',{'HOA','TLBO','GA','DE','PSO','ABC','GWO','SCA','BBO','ACO','RC-SA','HS'});
% % % %     subplot(3, 6, 4:6,'align'); bar(bar_mean2); title('F4')
% % % %     ylabel('f(min)'); set(gca,'xticklabel',{'HOA','TLBO','GA','DE','PSO','ABC','GWO','SCA','BBO','ACO','RC-SA','HS'});
% % % %     subplot(3, 6, 7:9,'align'); bar(bar_mean3); title('F9')
% % % %     ylabel('f(min)'); set(gca,'xticklabel',{'HOA','TLBO','GA','DE','PSO','ABC','GWO','SCA','BBO','ACO','RC-SA','HS'});
% % % %     subplot(3, 6, 10:12,'align'); bar(bar_mean4); title('F10')
% % % %     ylabel('f(min)'); set(gca,'xticklabel',{'HOA','TLBO','GA','DE','PSO','ABC','GWO','SCA','BBO','ACO','RC-SA','HS'});
% % % %     subplot(3, 6, 13:15,'align'); bar(bar_mean5); title('F14')
% % % %     ylabel('f(min)'); set(gca,'xticklabel',{'HOA','TLBO','GA','DE','PSO','ABC','GWO','SCA','BBO','ACO','RC-SA','HS'});
% % % %     subplot(3, 6, 16:18,'align'); bar(bar_mean6); title('F16')
% % % %     ylabel('f(min)'); set(gca,'xticklabel',{'HOA','TLBO','GA','DE','PSO','ABC','GWO','SCA','BBO','ACO','RC-SA','HS'});
% % % %     
% % % %     % Save figure to folder
% % % %     % home acct
% % % %     %FolderName = 'C:\Users\sid\Dropbox\MATLAB FILES\2021 MATWORK\HOA\Figures';
% % % %     % sch acct
% % % %     FolderName = 'C:\Dayo\THF Optimization\Figures';
% % % %     [~, file]  = fileparts(Function_name);  % Remove extension
% % % %     saveas(gca, fullfile(FolderName, [file, '.fig']));  % Append .fig
% % % %     %close(gcf);
    
    % clear place holderS
    H1=[];H2=[];H3=[];H4=[];H5=[];H6=[];H7=[];H8=[];H9=[];H10=[];H11=[];H12=[]; H13 = []; H14 = []; H15 = [];
    C1=[];C2=[];C3=[];C4=[];C5=[];C6=[];C7=[];C8=[];C9=[];C10=[];C11=[];C12=[]; C13 = []; C14 = []; C15 = [];
    W1=[];W2=[];W3=[];W4=[];W5=[];W6=[];W7=[];W8=[];W9=[];W10=[];W11=[];W12=[]; W13 = []; W14 = []; W15 = [];
    
end
exist = 2;
toc
delete(HH.figure)
clear('HH')