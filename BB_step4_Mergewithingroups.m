clear all; close all; clc;
%% user defined variables
%Get data from "Analyzed files" folder
dataDir = '/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';
%Save data to "Analyzed files" folder
saveDir = '/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';

%% Cycle through all groups and merge together 
%Create a string of all the groups you want to look at by filename end
groupnames= ["*STIM"; "*FREarly"; "*SOEarly"; "*FRMid"; "*SOMid"; "*FRLate"; "*SOLate"; "*FRextcueext"; "*SOextcueext"; "*extcueext"; "*FRreactcueext"; "*SOreactcueext"; "*reactcueext"; "*FRreinst"; "*FRreinstReact"; "*SOreinst"; "*SOreinstReact"; "*cocaine"; "*saline"];
for g=1:size(groupnames,1)                                                  %cycle through each group to be merged
groupstring=groupnames (g, 1);                                              %find current group name
groupchar=convertStringsToChars(groupstring);                               %Convert to character for saving purposes
groupname=g;                                                                %Number value for group
fileEnd=groupnames(g,1)+"-events*";                                         %Find all events files that end in groupnames
%% change directory and import data
cd(dataDir)
files = dir(fileEnd); %find all processed invididual data in directory
%% put relevant individual data into group structure
for f = 1:size(files,1)                                                     %for file (f) through the number of '-events' files
    load(files(f).name, 'animalID')                                         %load the animal ID
%     animalID = (['m',animalID]);                                          %make sure there's a letter in the animalID so it can be stored in a structure properly
    animalID = animalID(find(~isspace(animalID)));
    groupData.(animalID) = load(files(f).name,...                           %put event, eventMean, eventSEM, and PF Feeding structures into for each animal ID into the groupData structure
        'eventMean','AUCmean', 'inclusion','groupname');                      
    clear animalID
end
load(files(f).name,'window')                                                %load the window information (need to carry through to final step)
clear f files
%% change directory and save
cd(saveDir)
clearvars -except groupData window inclusion groupchar groupname groupnames dataDir saveDir
save([groupchar, '-merged'])
clear groupData
end