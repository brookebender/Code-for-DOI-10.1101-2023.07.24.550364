clear all; close all; clc;
%% user defined variables
%Get data from "Analyzed files" folder
dataDir = '/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';
%Save data to "Analyzed files" folder
saveDir = '/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';
experimentName = 'BB_P1_2_8_23';
%% change directory and import data
cd(dataDir)
files = dir('*-groupStats*');                                               %find all processed group data in directory (BE and Co)
clear dataDir
%% put relevant individual data into group structure
for f = 1:size(files,1)                                                     %cycle through files (f)
    load(files(f).name);                                                    %load file (f)

    if groupname ==1                                                   
        STIM.bin = bin;
        STIM.individual = groupEvent;
        STIM.traces = groupSummary;
        STIM.AUC = groupAUC;
        elseif groupname == 2                                                     
        FREarly.bin = bin;
        FREarly.individual = groupEvent;
        FREarly.traces = groupSummary;
        FREarly.AUC = groupAUC;
        elseif groupname == 3                                                     
        SOEarly.bin = bin;
        SOEarly.individual = groupEvent;
        SOEarly.traces = groupSummary;
        SOEarly.AUC = groupAUC;
        elseif groupname == 4                                                     
        FRMid.bin = bin;
        FRMid.individual = groupEvent;
        FRMid.traces = groupSummary;
        FRMid.AUC = groupAUC;
        elseif groupname == 5                                                   
        SOMid.bin = bin;
        SOMid.individual = groupEvent;
        SOMid.traces = groupSummary;
        SOMid.AUC = groupAUC;
        elseif groupname == 6                                                  
        FRLate.bin = bin;
        FRLate.individual = groupEvent;
        FRLate.traces = groupSummary;
        FRLate.AUC = groupAUC;
        elseif groupname == 7                                                   
        SOLate.bin = bin;
        SOLate.individual = groupEvent;
        SOLate.traces = groupSummary;
        SOLate.AUC = groupAUC;
        elseif groupname == 8                                                   
        FRextcueext.bin = bin;
        FRextcueext.individual = groupEvent;
        FRextcueext.traces = groupSummary;
        FRextcueext.AUC = groupAUC;
        elseif groupname == 9                                                   
        SOextcueext.bin = bin;
        SOextcueext.individual = groupEvent;
        SOextcueext.traces = groupSummary;
        SOextcueext.AUC = groupAUC;
        elseif groupname == 10                                                   
        extcueext.bin = bin;
        extcueext.individual = groupEvent;
        extcueext.traces = groupSummary;
        extcueext.AUC = groupAUC;
        elseif groupname == 11 
        FRreactcueext.bin = bin;
        FRreactcueext.individual = groupEvent;
        FRreactcueext.traces = groupSummary;
        FRreactcueext.AUC = groupAUC;
        elseif groupname == 12 
        SOreactcueext.bin = bin;
        SOreactcueext.individual = groupEvent;
        SOreactcueext.traces = groupSummary;
        SOreactcueext.AUC = groupAUC;
        elseif groupname == 13 
        reactcueext.bin = bin;
        reactcueext.individual = groupEvent;
        reactcueext.traces = groupSummary;
        reactcueext.AUC = groupAUC;
        elseif groupname == 14
        FRreinst.bin = bin;
        FRreinst.individual = groupEvent;
        FRreinst.traces = groupSummary;
        FRreinst.AUC = groupAUC;
        elseif groupname == 15
        FRreinstReact.bin = bin;
        FRreinstReact.individual = groupEvent;
        FRreinstReact.traces = groupSummary;
        FRreinstReact.AUC = groupAUC;
        elseif groupname == 16
        SOreinst.bin = bin;
        SOreinst.individual = groupEvent;
        SOreinst.traces = groupSummary;
        SOreinst.AUC = groupAUC;
        elseif groupname == 17
        SOreinstReact.bin = bin;
        SOreinstReact.individual = groupEvent;
        SOreinstReact.traces = groupSummary;
        SOreinstReact.AUC = groupAUC;
        elseif groupname == 18
        cocaine.bin = bin;
        cocaine.individual = groupEvent;
        cocaine.traces = groupSummary;
        cocaine.AUC = groupAUC;
        elseif groupname == 19
        saline.bin = bin;
        saline.individual = groupEvent;
        saline.traces = groupSummary;
        saline.AUC = groupAUC;
        elseif groupname == 20
        SOEarlyExcl.bin = bin;
        SOEarlyExcl.individual = groupEvent;
        SOEarlyExcl.traces = groupSummary;
        SOEarlyExcl.AUC = groupAUC;
        elseif groupname == 21
        FRLateExcl.bin = bin;
        FRLateExcl.individual = groupEvent;
        FRLateExcl.traces = groupSummary;
        FRLateExcl.AUC = groupAUC;
    else
        disp('Error; hit ctrl c')                                           %if the group is something else, there is an error
        pause
    end
    clear bin groupSummary groupEvent group groupAUC groupData
end
clear f
%% change directory and save
cd(saveDir)
clearvars -except extcueext STIM FREarly SOEarly FRMid SOMid FRLate SOLate FRextcueext SOextcueext reactcueext FRreinst SOreinst FRreinstReact SOreinstReact cocaine saline SOEarlyExcl FRLateExcl experimentName window
save(experimentName)