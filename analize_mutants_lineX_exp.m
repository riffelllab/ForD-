%% Code to analyze data from MUTANTS LINE_X experiments
% - Line 1 and Line 4 are Opsin 1 allele mutants. 
%   They are the same type of Opsin 1 mutation but their activity in the wind tunnel is different. M1 (Opsin 1 mutants) activity in the wind tunnel is also different from line 1 and line 4.
% - Line 2 is Opsin 1 & 2 double mutants.  
% - Line 3 is Opsin 1 & 2 double mutants. 
%   Both lines (2 and 3)  have different gRNA mutations. Can it be that our m0 experiments were done with the mutants from line 3? 
% - Line 5 is Wild Type.
%   I will group these data wit our current WT mosquitoes
% - Line 6 is Opsin 2 allele mutants.

% WARNING 1: The Wild Type experiments (wt) were done using the Riffell Lab
% aedes aegypti lines (instead of UC Santa Barbara)
% WARNING 2: Experiments done after MArch 23rd 2020 are 3 hours duration
% instead of 5 (1hr AIR - 1hr CO2 - 1hr PostCO2)

%%function datasetList= load_all_dataset()
close all;
clear variables;
format long;
%opengl('save', 'software');


% Information related to the data to analyze
workspace= 'D:\Diego\MatlabWorkspace\AnalysisFLYDRA_v_08272020';

inputPath= 'D:\Diego\MatlabWorkspace\MosquitoVisionProject\INPUT_Data\';
outputPath= 'D:\Diego\MatlabWorkspace\MosquitoVisionProject\OUTPUT_Graphs\';
subFolder= 'FLYDRA_Trajectories_mutants_lineX\';
%subFolder= 'FLYDRA_Trajectories_mutants_lineX\AdditionalExp\';

% Load the name of the files containing the experiments datasets stored in
% subFolder
inputPath= strcat(inputPath,subFolder);
cd(inputPath)
filesList=dir('*mainbrain.h5');

% Move to the MATLAB workspace
cd(workspace);
%Create the outputFolder folder
outputFolder = load_output_folder(outputPath, subFolder);

% Variables to use in this script
loadFullDataset=true;
fps=60.0; %90.0;
% Lower time threshold (in seconds) to plot a trajectory or not
flightTimeLimit= 1.5;
%Select the Control Color value
baseColor= 'white';  %select between: black, white, R-UHE
baseColorIndexList=zeros(length(filesList),1);
testedColorList= cell(length(filesList),1);
% Test section X,Y,Z limits (start_x= -lim_x and end_x= +lim_x --> Total x distance is lim_x*2)
% These lim values are the current size of the Test Section (DO NOT CHANGE IT)
lim_x= 0.9144;
lim_y= 0.3048;
lim_z= 0.6096;

fileNameList={};
typeGenderList={};
%% ===========================  
% Load and clean the raw data from Flydra
for expIndex= 1:length(filesList)
        
        % Load file name to work with 
        fileName= filesList(expIndex).name;
        filePath= strcat(inputPath,fileName);
        fileNameList=vertcat(fileNameList, {fileName(1:15)});
        %disp (filePath);
        % Load all the information from the FLYDRA .h5 file    
        [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z]= load_data_from_file(filePath, loadFullDataset);
        %tempData= [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z];
        
        % Load the experiment settings (cues positions & exp timestamps)
        [cuesSetup,ts_startAIR,ts_startCO2,ts_endCO2,ts_endAIR, mType, mGender]=  load_exp_settings(fileName);            
        % --- Clean the data loaded from the H5 file ---
        % Erase the X,Y and Z values that are out of the test section limits
        [attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z]= erase_pos_outside_WT_v2(attr_id, attr_time, attr_frame, attr_x, attr_y, attr_z, lim_x, lim_y, lim_z);
        % Erase the entries recorded before and after the mass flow
        % controller script is launched and the entries with timeStamp == 0.0.
        expIndexes= find(attr_time(:) >= ts_startAIR & attr_time(:) <= ts_endAIR & attr_time(:) ~=0.0);
        % Erase all insectIDs that have been detected in only 1 frame (or
        % below 0.x seconds (is it worth it to keep it?)
        
        % --- ---

        %Fill the dataset with the data from current experiment
        dataset(expIndex).fileName=    fileName;
        dataset(expIndex).type=        mType;
        dataset(expIndex).gender=       mGender;
        dataset(expIndex).expCues=    cuesSetup;
        
        dataset(expIndex).attr_id=     attr_id(expIndexes);
        dataset(expIndex).attr_time=   attr_time(expIndexes);
        dataset(expIndex).attr_frame=  attr_frame(expIndexes);
        dataset(expIndex).attr_x=      attr_x(expIndexes);
        dataset(expIndex).attr_y=      attr_y(expIndexes);
        dataset(expIndex).attr_z=      attr_z(expIndexes);

        % Estimate the timestamps values for odor stimulus ON and OFF
        % attr_time contains timestamp epochs (in seconds)
        indexPrevCO2= find(dataset(expIndex).attr_time(:) < ts_startCO2);
        indexWithCO2= find(dataset(expIndex).attr_time(:) >= ts_startCO2 & dataset(expIndex).attr_time(:) < ts_endCO2-1);
        indexPostCO2= find(dataset(expIndex).attr_time(:) >= ts_endCO2);
                
        dataset(expIndex).stim(indexPrevCO2,1)= {'AIR'};
        dataset(expIndex).stim(indexWithCO2,1)= {'CO2'};
        dataset(expIndex).stim(indexPostCO2,1)= {'postCO2'};
        
        % add the type and gender to list (used in the generate report fct
        typeGenderList= vertcat(typeGenderList, [{dataset(expIndex).type}, {dataset(expIndex).gender}]); 
        
        if length(cuesSetup(:,1)) == 3
            %Load the position where the baseCue is and the color tested in experiment
            baseCueIndex=find(strcmp(cuesSetup([2,3],1), baseColor));
            baseColorIndexList(expIndex)= baseCueIndex;
            testedColorList(expIndex)= cuesSetup(find(~strcmp(cuesSetup([2,3],1), baseColor))+1,1);
        end
end

%Clear temporary variables from workspace
clear fileName loadFullDataset indexPrevCO2 indexWithCO2 indexPostCO2
clear attr_id attr_time attr_frame attr_x attr_y attr_z
clear cuesSetup ts_startAIR ts_startCO2 ts_endCO2 ts_endAIR mType mGender
clear filesList
% ===========================
% ===========================


%% =================== 
% Estimate relative activity by adding all trajectories duration and then
% dividing it by the total number of trajectories. (different estimations 
% for prior/during/post CO2) 
% WARNING: WE are dividing each of the total number of trajectories prior,
% during and post by the total number of trajectories prior CO2
activityEstList= zeros(length(dataset), 6);
relFlightActDuration= zeros(length(dataset), 3);
xNames= cell(1,length(dataset));
for expIndex=1:length(dataset)
    %Find the indexes for AIR, CO2 and postCO2
    indexPrevCO2 = find(strcmp(dataset(expIndex).stim(:), 'AIR'));
    indexWithCO2 = find(strcmp(dataset(expIndex).stim(:),'CO2'));
    indexPostCO2 = find(strcmp(dataset(expIndex).stim(:),'postCO2'));
    % columns 1/3/5 has the # of trajectories detected per PrevCO2/WithCO2/PostCO2
    % and 2/4/6 contains the sum of these trajectories duration (in seconds)
    [activityEstList(expIndex,1), activityEstList(expIndex,2)]= estimate_relative_flt_activity_v2([dataset(expIndex).attr_id(indexPrevCO2) dataset(expIndex).attr_time(indexPrevCO2)], flightTimeLimit);
    [activityEstList(expIndex,3), activityEstList(expIndex,4)]= estimate_relative_flt_activity_v2([dataset(expIndex).attr_id(indexWithCO2) dataset(expIndex).attr_time(indexWithCO2)], flightTimeLimit);
    [activityEstList(expIndex,5), activityEstList(expIndex,6)]= estimate_relative_flt_activity_v2([dataset(expIndex).attr_id(indexPostCO2) dataset(expIndex).attr_time(indexPostCO2)], flightTimeLimit);
    
    %relFlightActDuration(expIndex,:)= [activityEstList(expIndex,2)/activityEstList(expIndex,2), activityEstList(expIndex,4)/activityEstList(expIndex,2), activityEstList(expIndex,6)/activityEstList(expIndex,2);];
    relFlightActDuration(expIndex,:)= [activityEstList(expIndex,2)/activityEstList(expIndex,2), activityEstList(expIndex,4)/activityEstList(expIndex,2), activityEstList(expIndex,6)/activityEstList(expIndex,4);];
 
    xNames{expIndex}= strcat(dataset(expIndex).fileName(1:8),' - (',dataset(expIndex).type, ' - ', dataset(expIndex).gender, ')');
end

figure()
bar(relFlightActDuration(1:end,:));
l= cell(1,3);
l{1}='Prior CO2'; l{2}='With CO2'; l{3}= 'Post CO2';
legend(l, 'Location', 'northeast');
if length(dataset) > 10
    set(gca,'XtickLabel', xNames(2:2:end));
else
    set(gca,'XtickLabel', xNames);
end
if contains(version, 'R2019')
    xtickangle(90)    
end
%title('Estimated flight activity compared to prior-CO2 activity');
title('Estimated flight activity compared to previous activity');
% ===========================
% ===========================


%% ===========================    
% Count the amount of trajectories detected while stim is prior, during and
% post CO2
trajSummary= zeros(length(dataset), 9);
for expIndex= 1:length(dataset)
    %Find the indexes fro AIR, CO2 and postCO2
    indexPrevCO2 = find(strcmp(dataset(expIndex).stim(:), 'AIR'));
    indexWithCO2 = find(strcmp(dataset(expIndex).stim(:),'CO2'));
    indexPostCO2 = find(strcmp(dataset(expIndex).stim(:),'postCO2'));
    
    % Count the total number of trajectories longer than flightTimeLimit
    % seconds. AS each row value in dataset.attr_id(:) is related to 1
    % frame we only need to check in how many frames that ID has been
    % detected
    disp(strcat(' --- FILE: ',dataset(expIndex).fileName, ' --- '));
    disp(' --- Data prev CO2 --');  
    [trajSummary(expIndex,1), trajSummary(expIndex,2), trajSummary(expIndex,3)]= count_trajectories([dataset(expIndex).attr_id(indexPrevCO2) dataset(expIndex).attr_time(indexPrevCO2)], flightTimeLimit);
    disp(' --- Data during CO2 --');
   [trajSummary(expIndex,4), trajSummary(expIndex,5), trajSummary(expIndex,6)]= count_trajectories([dataset(expIndex).attr_id(indexWithCO2), dataset(expIndex).attr_time(indexWithCO2)], flightTimeLimit);
    disp(' --- Data post CO2 --');
    [trajSummary(expIndex,7), trajSummary(expIndex,8), trajSummary(expIndex,9)]= count_trajectories([dataset(expIndex).attr_id(indexPostCO2), dataset(expIndex).attr_time(indexPostCO2)], flightTimeLimit);
    disp(' =====================================');
end
clear indexPrevCO2 indexWithCO2 indexPostCO2
% ===========================
% ===========================


%% Generate your own speed calculator

for odorChecked=[{'AIR'},{'CO2'},{'postCO2'}]
    for expIndex=1:length(dataset)
        odorIndex= find(strcmp(dataset(expIndex).stim, odorChecked));
        data=[dataset(expIndex).attr_id(odorIndex), dataset(expIndex).attr_time(odorIndex), dataset(expIndex).attr_x(odorIndex), dataset(expIndex).attr_y(odorIndex), dataset(expIndex).attr_z(odorIndex)];
        %select all the IDs from the vector
        uniqueID=unique(data(:,1));
        %Transpose from a column matrix to a row matrix
        uniqueID=uniqueID';
        % Initialice local variables
        totalTraj=0; 
        flightTimeList= nan(length(uniqueID),1);
        totalFlightTime=0;
        speedList=[];
        objIDsList=[];
        % For each insect ID, check if 
        for index= 1:length(uniqueID)
            objID= uniqueID(index);
            %Load the frames where appears the current objID
            objTime= find(data(:,1) == objID);
            % Estimate the duration of each insect ID trajectory
            duration= get_trajectory_duration(data(objTime(:),2));

            %Consider only if flight duration is bigger than flightTimeLimit seconds
            if duration >= flightTimeLimit
                % Calculate the lenght traveled by trajectory ID in each axis
                distX= sum(abs(diff(data(objTime(:),3))));
                distY= sum(abs(diff(data(objTime(:),4))));
                distZ= sum(abs(diff(data(objTime(:),5))));

                % Estimate the mean velocity
                distT= sqrt(distX^2 + distY^2 + distZ^2);
                speed= distT/duration;

                speedList= [speedList, speed]; 
                objIDsList= [objIDsList, objID];
            end
        end
        
        disp(strcat('-- ', int2str(expIndex), ' ---- ', odorChecked))
        disp(strcat('max speed:', sprintf('%.6f',max(speedList))))
        
        % Store the mean speed of each trajectory ID
        if strcmp(odorChecked, 'AIR')
            meanSpeedList(expIndex, 1)= mean(speedList);
            meanSpeedExp(expIndex).IDsAIR= objIDsList;
            meanSpeedExp(expIndex).speedsAIR= speedList;
        elseif strcmp(odorChecked, 'CO2')
            meanSpeedList(expIndex, 2)= mean(speedList);
            meanSpeedExp(expIndex).IDsCO2= objIDsList;
            meanSpeedExp(expIndex).speedsCO2= speedList;
        elseif strcmp(odorChecked, 'postCO2')
            meanSpeedList(expIndex, 3)= mean(speedList);
            meanSpeedExp(expIndex).IDsPostCO2= objIDsList;
            meanSpeedExp(expIndex).speedsPostCO2= speedList;
        end       

    end
end 

plot(meanSpeedList)
lg=[{'AIR'}, {'CO2'}, {'postCO2'}]
legend(lg, 'Location', 'northeast');
%% ==========================
%% ==========================


%% ===========================
% Plot the heatmaps for each of the experiments. There are 3 figures for
% each experiment heatmap prior/during/post CO2being released

%Parameters related to the heatmaps
topValueNormalized= 0.0001; %0.0003;
imgTitle= strcat(dataset(expIndex).fileName(1:15),'_',insType,'_',dataset(expIndex).gender);
plotCues= false;     %Control Flag to plot the odor/visual cues
    

for expIndex= 1:length(dataset)
    %Find the indexes for AIR, CO2 and postCO2
    indexPrevCO2 = find(strcmp(dataset(expIndex).stim(:), 'AIR'));
    indexWithCO2 = find(strcmp(dataset(expIndex).stim(:),'CO2'));
    indexPostCO2 = find(strcmp(dataset(expIndex).stim(:),'postCO2'));
    
    % Plot the heatmap in for the trajectories inside the test section
    if dataset(expIndex).type(1) =='m'
        insType= 'opsin_';
        if dataset(expIndex).type(2) =='0'
            insType= strcat(insType,'1&2');    
        else
            insType= strcat(insType, dataset(expIndex).type(2));
        end
    elseif dataset(expIndex).type(1) == 'l'
        insType= strcat('line_',dataset(expIndex).type(2));
    else
        insType= 'wild_type';
    end
    plot_XY_XZ_heatmaps_v5(plotCues,dataset(expIndex).attr_x(indexPrevCO2), dataset(expIndex).attr_y(indexPrevCO2), dataset(expIndex).attr_z(indexPrevCO2), dataset(expIndex).expCues, topValueNormalized, strcat(num2str(expIndex),'_1_',imgTitle,'_heatmap_prevCO2'));    
    plot_XY_XZ_heatmaps_v5(plotCues,dataset(expIndex).attr_x(indexWithCO2), dataset(expIndex).attr_y(indexWithCO2), dataset(expIndex).attr_z(indexWithCO2), dataset(expIndex).expCues, topValueNormalized, strcat(num2str(expIndex),'_2_',imgTitle,'_heatmap_withCO2'));
    plot_XY_XZ_heatmaps_v5(plotCues,dataset(expIndex).attr_x(indexPostCO2), dataset(expIndex).attr_y(indexPostCO2), dataset(expIndex).attr_z(indexPostCO2), dataset(expIndex).expCues, topValueNormalized, strcat(num2str(expIndex),'_3_',imgTitle,'_heatmap_postCO2'));   
end    
clear indexPrevCO2 indexWithCO2 indexPostCO2 imgTitle insType plotCues
% ===========================
% ===========================


%% ===========================
% Count insects (of particular type and gender) near visualClue when the CO2 is being released
% An insect ID will be counted as many times as it passes over the volume

% => (A.) This loop runs for each type of mosquito (wt,mX, lX), gender (m, f) or odor used (AIR, CO2, postCO2) <= 
% => To generate several datafiles, run it several times changing the values checkedType, checkedGender and checkedOdor <=  
%mosqTypes= ['all', 'wt','l1','l2','l3','l4','l5','m0', 'm1','m2'];
%checkedType= mosqTypes(1);
%checkedGender= 'f';
checkedOdor= 'AIR';
dataFolder= 'analysisData_oldVersion\';
xNames= cell(1,length(dataset));
countListWithCO2= zeros(length(dataset),length(dataset(1).expCues(2:end,1)));
for expIndex= 1:length(dataset)
    % disp(strcat(' - Working with expIndex:',num2str(expIndex)));
    % Create the new file Excel file name
    disp(strcat(' * Working with group: ',{' '}, dataset(expIndex).fileName(1:15)));
    fileName= strcat(dataset(expIndex).fileName(1:15), '_countsInsideCueVol_',checkedOdor,'.xlsx');
    
    % Returns vector with the blobal counts per each cue position [counts in pos1, counts in pos2, ..., counts in posN]
    countListWithCO2(expIndex,:)= count_insect_in_volume_v6(strcat(outputPath, outputFolder, dataFolder), fileName, dataset(expIndex), flightTimeLimit, checkedOdor);
    
    % Pick the experiment dates to plot later as the X axis legend
    %xNames{expIndex}= strcat(dataset(expIndex).fileName(1:8),' - (',dataset(expIndex).type, ' - ', dataset(expIndex).gender, ')');
end
clear checkedType checkedGender checkedOdor dataFolder
% => END (A.) LOOP <=

% ===========================
% ===========================


%% =================== PI VALUES =====================
% Generates the PreferenceIndex value per each experiment date
piList= generate_PI(countListWithCO2, baseColorIndexList);
%Change the NaN by 0
piList(isnan(piList))=0;
figure()
bar(piList)
set(gca,'XtickLabel', xNames);
if contains(version, 'R2019')
    xtickangle(45)    
end
title('Preference Index towards the tested color');

% For PI values grouped by mutants type, check plot_PI_for_mutants_exp.m
% ===========================
% ===========================



%% =============== GENERATE SUMMARY REPORT ===============
% generate report with all the information for this type of experiments
sheetName= strcat('MUTANTS_with_base_color_',baseColor);
fileName= strcat('Add_Exp_Data_Mosquito_Mutants_Project_Report_countsInVol_v6','.xls');
generate_report_for_mutants(outputPath, fileName, sheetName, fileNameList, typeGenderList, trajSummary, countListWithCO2, baseColorIndexList, testedColorList, piList)
%
% ===========================
% ===========================


%% ------------------------------------------------------------------------
% METHODS TO WORK WITH ALL THE DATA PREVIOUSLY GENERATED IN SUMMARIES/FILES
% -------------------------------------------------------------------------


%% =============== LOAD SUMMARY REPORT ===============
% Load the report created above with all the information for all experiments already pre-analyzed
sheetName= strcat('MUTANTS_with_base_color_',baseColor);
fileName= strcat('Add_Exp_Data_Mosquito_Mutants_Project_Report_countsInVol_v6','.xls');
[smryTrajectories, smryCounts, smryPIs, expList]= load_summary_report_mutants(outputPath, fileName, sheetName);
% ===========================
% ===========================


%% =================== WORKING WITH INDIVIDUAL FILES PER EXPERIMENT
% Load instants where insect are inside a cue volume and group them by time after the odor was released
%load the set of files (AIR,CO2 or postCO2) that you want to use
if strcmp(sheetName,'WT_black_vs_white_to_compareWithMutants')
    outputFolder= 'exp_white_black\';
end

% => (B.) This loop runs for each type of odor (AIR, CO2, postCO2) do you want to load the data files <= 
% => To load datafiles related to different odor stim, run it several times changing the checkedOdor value <=  
% Select the type of odor / part of experiment to load (AIR, CO2, postCO2)
odorChecked='AIR';
if (exist('outputPath', 'var') &&  exist('outputFolder', 'var'))
    % Load files
    filesPath= strcat(outputPath, outputFolder, 'analysisData_oldVersion\');
    cd(filesPath);
    filesList=dir(strcat('*countsInsideCueVol_',odorChecked,'.xlsx'));
    cd(workspace);
    % Define number of groups to generate
    numGrps=120; %=120 =60 sec per group || =4 =1800 sec per grp
    %find the counts near visual cue over time and per experiment
    [filesNameList, p1, p2, p3, p4, t1, t2, IDsInP1, IDsInP2]= load_insect_data_per_time_groups(numGrps, filesPath, filesList);
    % save this information in a .CSV file in the project folder
    save_counts_over_time_in_file(p1, p2, p3, p4)
else
    disp('  * Warning: file path and/or folder not specified in the variabes workspace ');
end
%If loading AIR or postCO2 data, assign it to its correct variable
switch(odorChecked)
    case 'CO2'
        p1IDsCO2= IDsInP1;
        p2IDsCO2= IDsInP2;
        p1CO2= p1;
        p2CO2= p2;
        t1CO2= t1;
        t2CO2= t2;
        
        % If working the the additional black vs white data done with wt from Riffell Lab (for those experiment done after March 23 2020)
        % For the 3hours exp, groups the subset of 30seconds in couples to make 1 minute each subset
        if strcmp(sheetName,'WT_black_vs_white_to_compareWithMutants')
            p1CO2((end-4):end,:)= horzcat(p1CO2((end-4):end,1:2:end)+p1CO2((end-4):end,2:2:end), zeros(length(p1CO2((end-4):end,1)),60));
            p2CO2((end-4):end,:)= horzcat(p2CO2((end-4):end,1:2:end)+p2CO2((end-4):end,2:2:end), zeros(length(p2CO2((end-4):end,1)),60));
            t1CO2((end-4):end,:)= horzcat(t1CO2((end-4):end,1:2:end)+t1CO2((end-4):end,2:2:end), zeros(length(t1CO2((end-4):end,1)),60));
            t2CO2((end-4):end,:)= horzcat(t2CO2((end-4):end,1:2:end)+t2CO2((end-4):end,2:2:end), zeros(length(t2CO2((end-4):end,1)),60));
        end
    case'AIR'
        p1IDsAIR= IDsInP1;
        p2IDsAIR= IDsInP2;
        p1AIR= p1;
        p2AIR= p2;
        t1AIR= t1;
        t2AIR= t2;
        % The aclimatation part is ALWAYS 1 hour duration (for all exp types),
        % groups the subset of 30seconds in couples to make 1 minute each subset
        p1AIR= horzcat(p1AIR(:,1:2:end)+p1AIR(:,2:2:end));
        p2AIR= horzcat(p2AIR(:,1:2:end)+p2AIR(:,2:2:end));
        t1AIR= horzcat(t1AIR(:,1:2:end)+t1AIR(:,2:2:end));
        t2AIR= horzcat(t2AIR(:,1:2:end)+t2AIR(:,2:2:end));
    case'postCO2'
        p1IDsPostCO2= IDsInP1;
        p2IDsPostCO2= IDsInP2;
        p1PostCO2= p1;
        p2PostCO2= p2;
        t1PostCO2= t1;
        t2PostCO2= t2;
        % If working the the additional black vs white data done with wt from Riffell Lab (for those experiment done after March 23 2020)
        % For the 3hours exp, groups the subset of 30seconds in couples to make 1 minute each subset
        if strcmp(sheetName,'WT_black_vs_white_to_compareWithMutants')
            p1PostCO2((end-4):end,:)= horzcat(p1PostCO2((end-4):end,1:2:end)+p1PostCO2((end-4):end,2:2:end), zeros(length(p1PostCO2((end-4):end,1)),60));
            p2PostCO2((end-4):end,:)= horzcat(p2PostCO2((end-4):end,1:2:end)+p2PostCO2((end-4):end,2:2:end), zeros(length(p2PostCO2((end-4):end,1)),60));
            t1PostCO2((end-4):end,:)= horzcat(t1PostCO2((end-4):end,1:2:end)+t1PostCO2((end-4):end,2:2:end), zeros(length(t1PostCO2((end-4):end,1)),60));
            t2PostCO2((end-4):end,:)= horzcat(t2PostCO2((end-4):end,1:2:end)+t2PostCO2((end-4):end,2:2:end), zeros(length(t2PostCO2((end-4):end,1)),60));
        end
end
%If you are only using 2 Visual Cues
clear p1 p2 p3 p4 t1 t2 numGrps odorChecked IDsInP1 IDsInP2
% => END (B.) LOOP <=
% ==================================================
% ==================================================    


%% ====================== (I) % of trajectories NEAR BOTH CUES (per mutation type)
% If working the RAW data : TrajSummary contains information about the trajectories. Column values are: 
% [Number of Traj AIR, Avg duration traj AIR, Max. duration traj AIR, ...
%   Number of Traj CO2, Avg duration traj CO2, Max. duration traj CO2, ...   
%       Number of Traj PostCO2, Avg duration traj PostCO2, Max. duration traj PostCO2]
% If working with SUMMARY REPORT data, use the smryTrajectorys table
uniqueTypes= unique(expList.mosqType);
typeIndex=1;
figure()
for typeValue=uniqueTypes'
    clear temp tempNorm clear trajNearCues
    i= find(strcmp(expList.mosqType, typeValue));
    temp= p1IDsAIR(i)+p2IDsAIR(i);
    tempNorm= bsxfun(@rdivide,(temp*100),smryTrajectories.numOfTrajAIR(i));
    trajNearCues(:,1)= tempNorm;
    temp= p1IDsCO2(i)+p2IDsCO2(i);
    tempNorm= bsxfun(@rdivide,(temp*100),smryTrajectories.numOfTrajCO2(i));
    trajNearCues(:, 2)= tempNorm;
    %temp= p1IDsPostCO2(i)+p2IDsPostCO2(i);
    %tempNorm= bsxfun(@rdivide,(temp*100),smryTrajectories.numOfTrajPostCO2(i));
    %trajNearCues(:, 3)= tempNorm;

    %Generate the boxplot for this type of mutation
    subplot(3,3,typeIndex)
    if (length(i) == 1)
        plot(trajNearCues, '*r', 'MarkerSize', 5);
    else
        boxplot([trajNearCues(:,1), trajNearCues(:,2),trajNearCues(:,3)])
    end
    title(strcat('% of trajectories detected near BOTH cues (',typeValue,')'));
    xNames=[{'PriorCO2'}, {'WithCO2'},{'PostCO2'}];
    set(gca,'XtickLabel', xNames);
    ylabel('% of Trajectories')
    typeIndex= typeIndex+1;    
end
% ==================================================
% ==================================================    


%% WARNING: In following analysis the postCO2 is not analyzed as we are not
% sure how much the CO2 remains can affect insect behavior


% ====================== (II) % of trajectories NEAR TESTED CUE (per mutation type)
% If working the RAW data : TrajSummary contains information about the trajectories. Column values are: 
% [Number of Traj AIR, Avg duration traj AIR, Max. duration traj AIR, ...
%   Number of Traj CO2, Avg duration traj CO2, Max. duration traj CO2, ...   
%       Number of Traj PostCO2, Avg duration traj PostCO2, Max. duration traj PostCO2]
% If working with SUMMARY REPORT data, use the smryTrajectorys table
uniqueTypes= unique(expList.mosqType);
typeIndex=1;
figure()
%tempValListAIR= zeros(length(dataset),1);
%tempValListCO2= zeros(length(dataset),1);
%tempValListPostCO2= zeros(length(dataset),1);
for typeValue=uniqueTypes'
    clear temp tempNorm clear trajNearTestedCue
    i= find(strcmp(expList.mosqType, typeValue));
    tempAIR=[];
    tempCO2=[];
%    tempPostCO2=[];
    for v=i'
        % group the data in function of the base color cue position
        if (smryCounts.baseColor_Pos(v) == 1)
            tempAIR= [tempAIR; p2IDsAIR(v)];
            tempCO2= [tempCO2; p2IDsCO2(v)];
%            tempPostCO2= [tempPostCO2; p2IDsPostCO2(v)];
        else
            tempAIR= [tempAIR; p1IDsAIR(v)];
            tempCO2= [tempCO2; p1IDsCO2(v)];
%            tempPostCO2= [tempPostCO2; p1IDsPostCO2(v)];
        end
    end
    % Calculate the % of trajectories interested in the tested color cue
    tempNorm= bsxfun(@rdivide,(tempAIR*100),smryTrajectories.numOfTrajAIR(i));
    trajNearTestedCue(:,1)= tempNorm;
    %tempValListAIR(i)= tempNorm;
    tempNorm= bsxfun(@rdivide,(tempCO2*100),smryTrajectories.numOfTrajCO2(i));
    trajNearTestedCue(:, 2)= tempNorm;
    %tempValListCO2(i)= tempNorm;
%    tempNorm= bsxfun(@rdivide,(tempPostCO2*100),smryTrajectories.numOfTrajPostCO2(i));
%    trajNearCues(:, 3)= tempNorm;
%    tempValListPostCO2(i)= tempNorm;
    
    %Generate a SCATTER PLOT for this type of mutation
    subplot(3,3,typeIndex)
    scatter(ones(size(trajNearTestedCue,1),1), trajNearTestedCue(:,1), 'jitter', 'on', 'MarkerFaceColor', [0, 0.5, 0]);
    hold on
    scatter(ones(size(trajNearTestedCue,1),1), trajNearTestedCue(:,2), 'jitter', 'on', 'MarkerFaceColor', [0.5, 0, 0]);
%    scatter(ones(size(trajNearCues,1),1), trajNearCues(:,3), 'jitter', 'on', 'MarkerFaceColor', [0, 0, 0.5]);

    % plot also the means values 
    y= [mean(trajNearTestedCue(:,1)), mean(trajNearTestedCue(:,1))];
    plot([0.9, 1.1], y, '--g', 'MarkerSize',15,'MarkerEdgeColor',[0 0.5 0]);
    y= [mean(trajNearTestedCue(:,2)), mean(trajNearTestedCue(:,2))];
    plot([0.9, 1.1], y, '--r', 'MarkerSize',15,'MarkerEdgeColor',[0.5 0 0]);
%    y= [mean(trajNearCues(:,3)), mean(trajNearCues(:,3))];
%    plot([0.9, 1.1], y, '--b', 'MarkerSize',15,'MarkerEdgeColor',[0 0 0.5]);
    hold off
    lg= [{' With AIR'}, {'With CO2'}];
    legend(lg);
    title(strcat('% of trajectories detected near black cue (',typeValue,')'));
    xlim([0,2]);
    ylim([0, 25]);
    ylabel('% of Trajectories')
    typeIndex= typeIndex+1;    
end
clear tempAIR tempCO2 tempPostCO2 tempNorm
clear lg y i typeValue tipeIndex
% ==================================================
% ==================================================    



%% ====================== (III) % of time NEAR TESTED CUE per mutation type)
% the matrix tXAIR, txCO2 and txPostCO2 have the amount of time spent by
% all trajectories for a given time slot of 60 seconds.

%Obtain the total amount of time per experiment (AIR and CO2 only)
totalT1AIR= sum(t1AIR,2);
totalT2AIR= sum(t2AIR,2);
totalT1CO2= sum(t1CO2,2);
totalT2CO2= sum(t2CO2,2);
%totalT1PostCO2= sum(t1PostCO2,2);
%totalT2PostCO2= sum(t2PostCO2,2);

%Initialize parameters
totalFlightTimeAIR= zeros(size(expList.date,1),1);
totalFlightTimeCO2= zeros(size(expList.date,1),1);
%totalFlightTimePostCO2= zeros(size(expList.date,1),1);

% => (C.) This loop runs for each type of odor (AIR, CO2, postCO2) you want to analyze <= 
% => To work with the data related to different odor stim, run it several times changing the checkedOdor value (AIR, CO2, postCO2)<=  
odorChecked= 'CO2';
for expIndex= 1:length(dataset)  
    %Find the indexes for AIR, CO2 and postCO2
    indexOdor = find(strcmp(dataset(expIndex).stim(:), odorChecked)); 
    tempDataset= [dataset(expIndex).attr_id(indexOdor) dataset(expIndex).attr_time(indexOdor)];
    %select all the IDs from the dataset
    uniqueID=unique(tempDataset(:,1));
    %Transpose from a column matrix to a row matrix
    uniqueID=uniqueID';
    tempFlightTime=0;
    %tempCount=0;
    % For each trajectory ID check if it is longer than flightTimeLimit
    % and (if yes) add it to the total for the experiment
    for index= 1:length(uniqueID)
        objID= uniqueID(index);
        %Load the frames where appears the current objID
        objTime= find(tempDataset(:,1) == objID);
        % Estimate the duration of each insect ID trajectory
        duration= get_trajectory_duration(tempDataset(objTime,2));
        %Consider only if flight duration is bigger than flightTimeLimit seconds
        if (duration >= flightTimeLimit)
            tempFlightTime= tempFlightTime+duration;
            %tempCount= tempCount+1;
        end
    end
    if strcmp(odorChecked, 'CO2')
        totalFlightTimeCO2(expIndex)= tempFlightTime;
    elseif strcmp(odorChecked, 'postCO2')
        totalFlightTimePostCO2(expIndex)= tempFlightTime;
    else
        %AIR
        totalFlightTimeAIR(expIndex)= tempFlightTime;
    end
    %disp(strcat('total counts--',odorChecked,'---', num2str(tempCount)));
end
% => END (C.) LOOP <=

%Group by Line_X type, normalize and plot
uniqueTypes= unique(expList.mosqType);
typeIndex=1;
figure()
%tempValListAIR=zeros(length(dataset),1);
%tempValListCO2=zeros(length(dataset),1);
for typeValue=uniqueTypes'
    clear temp tempNorm clear timeNearBlackCue
    i= find(strcmp(expList.mosqType, typeValue));
    tempAIR=[];
    tempCO2=[];
    tempPostCO2= [];
    for v=i'
        if (smryCounts.baseColor_Pos(v) == 1)
            tempAIR= [tempAIR; totalT2AIR(v)];
            tempCO2= [tempCO2; totalT2CO2(v)];
%            tempPostCO2= [tempPostCO2; totalT2PostCO2(v)];
        else
            tempAIR= [tempAIR; totalT1AIR(v)];
            tempCO2= [tempCO2; totalT1CO2(v)];
%            tempPostCO2= [tempPostCO2; totalT1PostCO2(v)];
            
        end
    end
    tempNorm= bsxfun(@rdivide,(tempAIR*100),totalFlightTimeAIR(i));
    timeNearBlackCue(:,1)= tempNorm;
    %tempValListAIR(i)= tempNorm;
    tempNorm= bsxfun(@rdivide,(tempCO2*100),totalFlightTimeCO2(i));
    timeNearBlackCue(:, 2)= tempNorm;
   % tempValListCO2(i)= tempNorm;
%     tempNorm= bsxfun(@rdivide,(tempPostCO2*100),totalFlightTimePostCO2(i));
%     timeNearBlackCue(:, 3)= tempNorm;
    
    %Generate a SCATTER PLOT for this type of mutation
    %subplot(3,3,typeIndex)
    scatter(ones(size(timeNearBlackCue,1),1), timeNearBlackCue(:,1), 'jitter', 'on', 'MarkerFaceColor', [0, 0.5, 0]);
    hold on
    scatter(ones(size(timeNearBlackCue,1),1), timeNearBlackCue(:,2), 'jitter', 'on', 'MarkerFaceColor', [0.5, 0, 0]);
%    scatter(ones(size(timeNearBlackCue,1),1), timeNearBlackCue(:,3), 'jitter', 'on', 'MarkerFaceColor', [0, 0, 0.5]);

    % plot also the means values 
    y= [mean(timeNearBlackCue(:,1)), mean(timeNearBlackCue(:,1))];
    plot([0.9, 1.1], y, '--g', 'MarkerSize',15,'MarkerEdgeColor',[0 0.5 0]);
    y= [mean(timeNearBlackCue(:,2)), mean(timeNearBlackCue(:,2))];
    plot([0.9, 1.1], y, '--r', 'MarkerSize',15,'MarkerEdgeColor',[0.5 0 0]);
%     y= [mean(timeNearBlackCue(:,3)), mean(timeNearBlackCue(:,3))];
%     plot([0.9, 1.1], y, '--b', 'MarkerSize',15,'MarkerEdgeColor',[0 0 0.5]);

    hold off
    lg= [{' With AIR'}, {'With CO2'}];
    %lg= [{' With AIR'}, {'With CO2'}, {'With AIR (postCO2)'}];
    legend(lg);
    title(strcat('% of trajectory time detected near black cue (',typeValue,')'));
    xlim([0,2]);
    ylim([0, 25]);
    %xNames=[{'PriorCO2'}, {'WithCO2'},{'PostCO2'}];
    %set(gca,'XtickLabel', xNames);
    ylabel('% of Trajectories')
    typeIndex= typeIndex+1;
    disp(strcat('AIR -- ', typeValue));
    disp(strcat('--indexes: ',num2str(i)));
    disp(strcat('% time trajs spent in cue:',num2str(timeNearBlackCue(:,1))));
    disp('CO2')
    disp(strcat('--indexes: ',num2str(i)));
    disp(strcat('% time trajs spent in cue:',num2str(timeNearBlackCue(:,2))));
%     disp('Post CO2')
%     disp(strcat('--indexes: ',num2str(i)));
%     disp(strcat('% time trajs spent in cue:',num2str(timeNearBlackCue(:,3))));
%     
end

clear tempAIR tempCO2 tempNorm
clear lg y i typeValue tipeIndex
clear totalT1AIR totalT2AIR totalT1CO2 totalT2CO2
% ==================================================
% ==================================================    


%% ====== % OF TIME NEAR CUE (INDIVIDUALLY FOR EACH TRAJ NEAR CUE) ====
% Generate % of time spent near cue per time of trajectory (for each
% trajecotry individually)

%(1.) Estimate time spent near visual cues for each one of trhe
%trajectories
odorChecked='CO2';
% Load files
filesPath= strcat(outputPath, outputFolder,'analysisData_oldVersion\');
%filesPath= strcat(outputPath, outputFolder,'analysisData_newVersion\');
cd(filesPath);
filesList=dir(strcat('*countsInsideCueVol_',odorChecked,'.xlsx'));
cd(workspace);
for fileIndex= 1:length(filesList)       
    fileName= filesList(fileIndex).name;
    disp(strcat(' - Working with file: ', {' '}, fileName));
    filesName(fileIndex)= {fileName(1:8)};
    % data loaded structure: [vCuePosition	objID	timeStamp NaN initialTimeStamp finalTimeStamp]
    dataFromExcel = table2array(readtable(strcat(filesPath, fileName)));
    expData(fileIndex).name= fileName;
    % Initialize the correct fields in function of the odor used
    if strcmp(odorChecked, 'AIR')
        expData(fileIndex).listIDsInCueAIR=[];
        expData(fileIndex).listTimeInCueAIR=[];
        expData(fileIndex).listTrajTimeAIR=[];
        expData(fileIndex).percentNearCueAIR=[];
    elseif strcmp(odorChecked, 'CO2')
        expData(fileIndex).listIDsInCueCO2=[];
        expData(fileIndex).listTimeInCueCO2=[];
        expData(fileIndex).listTrajTimeCO2=[];
        expData(fileIndex).percentNearCueCO2=[];
    end
    if isempty(dataFromExcel)
        disp(strcat('   -> No counts detected in file: ',{' '}, fileName));
    else
        %Initialize the object that will contain the information for the current experiment            
        tmpListTime=[];
        tmpIndex= 1;
        % To load timestamp data (time spent by insect near thew cue
        % Sorted the the data by cue visited ID timestamp
        sortedData= sortrows(dataFromExcel(:,1:3));
        if length(unique(sortedData(:,1)))== 1
            %If we have counts in only one position, check if is the Tested cue
            if baseColorIndexList(expIndex) == unique(sortedData(:,1))
                continue
            end
        else
            % Split the data in function to the position values (1, 2)
            split= find(sortedData(:,1)==2, 1);
            % Find the first appearence of the position 2 
            if baseColorIndexList(fileIndex)== 1
                %Base color was placed in position 1 (+Y axis) in current experiment
                sortedData= sortedData(split:end,:);
            else
                %Base color was placed in position 2 (-Y axis) in current experiment
                sortedData= sortedData(1:(split-1),:);
            end   
        end
        %Estimates the time spend for each ID in the pos where was TESTED cue
        tmpIndex=1;
        for objID= unique(sortedData(:,2))'
            idIndexes= find(sortedData(:,2) ==objID);
            % Calculate the Delta time between timestamps values for the id 
            tsDiff= diff(sortedData(idIndexes,3));
            % find which  sequential counts didn't happen "~=consecutively"
            out= find(tsDiff > 0.05);
            t=0;
            k=1;
            if any(out)
                for i= out
                    % the insect left the volume at a given moment
                    t= t + sum(tsDiff(k:(i-1)));
                    k= i+1;
                end
                t= t + sum(tsDiff(k:end));
            else
                %the insects has been inside the volume all the time
                t= sum(tsDiff);
            end
            tmpListTime(1,tmpIndex)= t;
            tmpIndex= tmpIndex+1;
        end
        %Add the IDs and the their time spent visiting the cue to the structure
        if strcmp(odorChecked, 'AIR')
            expData(fileIndex).listIDsInCueAIR=unique(sortedData(:,2))';
            expData(fileIndex).listTimeInCueAIR=tmpListTime;
  
        elseif strcmp(odorChecked, 'CO2')
            expData(fileIndex).listIDsInCueCO2=unique(sortedData(:,2))';
            expData(fileIndex).listTimeInCueCO2=tmpListTime;

        end
    end
end
clear tempIDs indexes dataFromExcel tmpIndex tmpListTime t k idIndexes tsDiff        
 

%(2. A) Once we have generated the time spent near the cues (1.), estimate 
% the amount of time spent by each of the trajectories ids FOR TRAJ WITH
% DURATION LONGER THAN 1.5 SEC ONLY
% WARNING: Data for AIR and CO2 must be loaded first
odorChecked='CO2';
for expIndex= 1:length(dataset)
    % Initialize temp var
    tmpListTrajTime=[];
    %Find the indexes for AIR, CO2 and postCO2 to create a working subdataset
    indexOdor = find(strcmp(dataset(expIndex).stim(:), odorChecked)); 
    tempDataset= [dataset(expIndex).attr_id(indexOdor) dataset(expIndex).attr_time(indexOdor)];
    % For each trajectory ID check if it is longer than flightTimeLimit
    % and (if yes) add it to the total for the experiment
    if strcmp(odorChecked, 'AIR')
        tmpListIDs= expData(expIndex).listIDsInCueAIR;
    elseif strcmp(odorChecked, 'CO2')
       tmpListIDs= expData(expIndex).listIDsInCueCO2;
    end
    for tmpIndex= 1:length(tmpListIDs)
        objID= tmpListIDs(tmpIndex);
        %Load the frames where appears the current objID
        objTime= find(tempDataset(:,1) == objID);
        % Estimate the duration of each insect ID trajectory
        duration= get_trajectory_duration(tempDataset(objTime,2));
        %Consider only if flight duration is bigger than flightTimeLimit seconds
        if (duration <= flightTimeLimit)
              duration=NaN;
        end
        tmpListTrajTime(1,tmpIndex)= duration;
    end
    if strcmp(odorChecked, 'AIR')
        IDsToKeep= ~isnan(tmpListTrajTime);
        expData(expIndex).listIDsInCueAIR= expData(expIndex).listIDsInCueAIR(IDsToKeep);
        expData(expIndex).listTimeInCueAIR= expData(expIndex).listTimeInCueAIR(IDsToKeep);
        expData(expIndex).listTrajTimeAIR= tmpListTrajTime(IDsToKeep);
    elseif strcmp(odorChecked, 'CO2')
        IDsToKeep= ~isnan(tmpListTrajTime);
        expData(expIndex).listIDsInCueCO2= expData(expIndex).listIDsInCueCO2(IDsToKeep);
        expData(expIndex).listTimeInCueCO2= expData(expIndex).listTimeInCueCO2(IDsToKeep);
        expData(expIndex).listTrajTimeCO2= tmpListTrajTime(IDsToKeep);
    end
end
clear tmpIndex tmpListTrajTime


% (2.B) == SCATTER PLOT of the data generated before for AIR and CO2
% WARNING: Data for AIR and CO2 must be loaded first
mosqTypes={'m0';'m2';'wt';'m0'}; %type of mosquito for each experiment date
uniqueTypes={'wt';'m0';'m2'};
%mosqTypes={'m0';'m0';'m1';'m2';'l1';'l4';'l3';'l5';'l2';'l6';'l1';'l6';'l2';'l5';'wt';'l3';'l4';'l5';'l3';'l4';'wt'};
%uniqueTypes={'wt';'m0';'m2'};%;'m1';'l1';'l2';'l3';'l4';'l5';'l6'};
figure()
for typeIndex=1:size(uniqueTypes,1)
    % Group experiments data in function of type of mosquitoes
    typeValue= uniqueTypes(typeIndex);
    i= find(strcmp(mosqTypes, typeValue));
    timeTrajNearCueAIR=[];
    timeTrajNearCueCO2=[];
    for expIndex=i'
        % FOR AIR
        tmpNorm= bsxfun(@rdivide, (expData(expIndex).listTimeInCueAIR*100),expData(expIndex).listTrajTimeAIR);
        timeTrajNearCueAIR= [timeTrajNearCueAIR, tmpNorm];
        % Add percentage values to the strucutre
        expData(expIndex).percentNearCueAIR= tmpNorm;
        %FOR CO2
        tmpNorm= bsxfun(@rdivide, (expData(expIndex).listTimeInCueCO2*100),expData(expIndex).listTrajTimeCO2);
        timeTrajNearCueCO2= [timeTrajNearCueCO2, tmpNorm];
        % Add percentage values to the strucutre
        expData(expIndex).percentNearCueCO2= tmpNorm;
    end
    
    %Save values in table
    timeInCuePerTypeNorm(typeIndex).type= typeValue; 
    timeInCuePerTypeNorm(typeIndex).timeTrajNearCueAIR= timeTrajNearCueAIR;
    timeInCuePerTypeNorm(typeIndex).timeTrajNearCueCO2= timeTrajNearCueCO2;
    %Generate a SCATTER PLOT for this type of mutation
    subplot(3,3,typeIndex)
    scatter(ones(size(timeTrajNearCueAIR,2),1), timeTrajNearCueAIR, 10, 'jitter', 'on', 'MarkerFaceColor', [0, 0.5, 0], 'MarkerFaceAlpha', 0.5);
    hold on
    scatter(ones(size(timeTrajNearCueCO2,2),1)*2, timeTrajNearCueCO2, 10, 'jitter', 'on', 'MarkerFaceColor', [0.5, 0, 0], 'MarkerFaceAlpha', 0.5);
    % plot also the means values 
    y= [mean(timeTrajNearCueAIR), mean(timeTrajNearCueAIR)];
    plot([0.9, 1.1], y, 'g', 'LineWidth',5,'MarkerEdgeColor',[0 0.5 0]);
    y= [mean(timeTrajNearCueCO2), mean(timeTrajNearCueCO2)];
    plot([1.9, 2.1], y, 'r', 'LineWidth',5,'MarkerEdgeColor',[0.5 0 0]);
    hold off
    lg= [{' With AIR'}, {'With CO2'}];
    %lg= [{' With AIR'}, {'With CO2'}, {'With AIR (postCO2)'}];
    legend(lg);
    title(strcat('% time near black cue per trajectory (',typeValue,')'));
    xlim([0,3]);
    ylim([0, 100]);
    ylabel('% of time')
    %typeIndex= typeIndex+1;
end    
% ====================================
%=====================================


%%  ==== FIND AND STORE ALL THE TRAJECTORIES DURATION 
%       DETECTED FOR A ODOR STIM. ====
odorChecked='CO2';
for expIndex= 1:length(dataset)
    % Initialize temp var
    tmpIndex=1;
    tempDataset= [];
    tmpListTrajTime= [];
    tmpListTimeAIR= [];
    tmpListIDsCO2= [];
    %Find the indexes for AIR, CO2 and postCO2 to create a working subdataset
    indexOdor = find(strcmp(dataset(expIndex).stim(:), odorChecked)); 
    tempDataset= [dataset(expIndex).attr_id(indexOdor) dataset(expIndex).attr_time(indexOdor)];
    % Select all the IDs from the dataset
     tmpListIDs= unique(tempDataset(:,1));
     countTraj=1;
     for index= 1:length(tmpListIDs)
        objID= tmpListIDs(index);
        %Load the frames where appears the current objID
        objTime= find(tempDataset(:,1) == objID);
        % Estimate the duration of each insect ID trajectory
        duration= get_trajectory_duration(tempDataset(objTime,2));
        %Consider only if flight duration is bigger than flightTimeLimit seconds
        if (duration >= flightTimeLimit)
            if strcmp(odorChecked, 'AIR')
                expDataTrajDuration(expIndex).listTrajTimeAIR(1,countTraj)= duration;
                countTraj=countTraj + 1;
           elseif strcmp(odorChecked, 'CO2')
                expDataTrajDuration(expIndex).listTrajTimeCO2(1,countTraj)= duration;
                countTraj=countTraj + 1;
            end
        end
     end
end

%SAVE THE DURATION IN AN EXCEL FILE
fileNameAIR= 'addData_exp_Duration_All_Traj_Detected_AIR.xlsx';
fileNameCO2= 'addData_exp_Duration_All_Traj_Detected_CO2.xlsx';
outputFileAIR= strcat(outputPath,outputFolder,fileNameAIR);
outputFileCO2= strcat(outputPath,outputFolder,fileNameCO2);

for expIndex= 1:length(dataset)
    tAIR=  table(expDataTrajDuration(expIndex).listTrajTimeAIR', 'VariableNames', {'traj_duration_AIR'});
    tCO2 = table(expDataTrajDuration(expIndex).listTrajTimeCO2', 'VariableNames', {'traj_duration_CO2'});

    sheetName= dataset(expIndex).fileName(1:8);
    %writetable([expDataTrajDuration(1).listTrajTimeAIR'; expDataTrajDuration(1).listTrajTimeCO2'], outputFile, 'sheet',expSheet);
    writetable(tAIR, outputFileAIR, 'sheet', sheetName);
    writetable(tCO2,  outputFileCO2, 'sheet', sheetName);
end
% ===========================
% ===========================
