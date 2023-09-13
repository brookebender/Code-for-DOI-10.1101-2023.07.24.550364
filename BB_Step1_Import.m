%% first import raw trace data
%File where photo data is located *REPLACE PATH/FILE NAME*
filedir="/Users/brookebender/Documents/MATLAB/Photometry/BB_P2/Results/BB_P2_M3_SO_Veh_PhotoData.csv";
%Import photometry data into proper format
[Frame, timestamp, RawF_465_F1, RawF_465_F2, RawF_410_F1, RawF_410_F2, RawF_560_F1, RawF_560_F2] = importfile2(filedir, [2, Inf]);
%Set save file to Imported data file
saveDir = '/Users/brookebender/Documents/MATLAB/Photometry/BB_P2/Imported data';
cd(saveDir);
%Name animalID *REPLACE RAT_DATE_PROGRAM* (This needs to be different for each recording session)
animalID = 'M3_SO_Veh'; 
%Save as "animalID-Imported"
save([animalID, '-Imported']);
%Clear photometry data file
clear filedir
%% Import event files
%File where event data is located *REPLACE PATH/FILE NAME*
filedir="/Users/brookebender/Documents/MATLAB/Photometry/BB_P2/Results/BB_P2_M3_SO_Veh_EventData.csv";
%Import event data into proper format
[Event, StartTimes, EndTimes] = importevents(filedir, [2, Inf]);
%Set save file to Imported data file
cd(saveDir);
%Save as "animalID-Events"
save([animalID, '-Events']);
clear