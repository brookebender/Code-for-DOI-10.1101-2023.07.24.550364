clear all; close all; clc;
%% User Modified Variables
%Get data from "Analyzed files" folder
dataDir = '/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';
%Save data to "Analyzed files" folder
saveDir = '/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';

%pulling in data after completing Step 4. 
group = 'STIM';
fileName = '*STIM-merged';


AUCwin=[0, 1]; %Enter the window size in seconds relative to event onset to calculate the AUC (event onset = time 0, FR 3 sec after event would be [0 3])
FRameRate = 30;                                                           %FRame rate of FP camera; it will be 1/2 the set FRame rate if in interleaved mode, or the set FRame rate if in constant
FRameDuration = 1/FRameRate;
%% change directory and import data
cd(dataDir)
load(fileName)
clear groupchar groupnames

%% AUC calculations for each trace surrounding events based on AUCwin

%The AUC structure created with this should have an individual AUC value
%for each session x color x event x fiber x trial 

fnID=fieldnames(groupData);                                             
fnWindow = AUCwin;
twocolor={'signal','red'};
fnRegion = {'F1','F2'};

for a = 1:size(fnID)
    for q=1:length(twocolor)
        fnEventType = fieldnames(groupData.(fnID{a}).eventMean.(twocolor{q}));
        for e = 1:size(fnEventType,1)                                               %cycle through event type (e)
            winDuration = ceil(abs(fnWindow(1,1)) +...             %calculate the max number of samples in a given window of interest based on the window size and FRameDuration
                (fnWindow(1,2))/FRameDuration);
            for r = 1:length(fnRegion);                              %cycle through regions (r)
                    tempEvent = groupData.(fnID{a}).eventMean.(twocolor{q}).(fnEventType{e}).post.(fnRegion{r});
                    tempData = transpose(tempEvent);        %put indexed trace data into a temporary vector
                    if size(tempData,1)>=winDuration
                        tempData=tempData(1:winDuration);                             %occasionally getting a PETH that has an extra FRame on the end, FR adjusting length. 
                        AUCstats.mean.(fnID{a}).(twocolor{q}).(fnEventType{e}).(fnRegion{r})=...
                             trapz(tempData);
                        peakstats.mean.(fnID{a}).(twocolor{q}).(fnEventType{e}).(fnRegion{r})=...
                             max(tempData);
                    else
                end
                
            end
            
        end
    end
end

clear winDuration FRameDuration temp* idx a e r y ySize w p k len* tempmot

%% make arrays for each group



%Make list of all sessions included
sessionslist=fieldnames(AUCstats.mean);
for b=1:size(sessionslist);
    eventlist=fieldnames(AUCstats.mean.(sessionslist{b}).red);
    for c=1:size(eventlist);
        tempevent=eventlist{c};
        templabel=sessionslist{b};
        AUC_red.F1.(tempevent)(b,1)=AUCstats.mean.(templabel).red.(tempevent).F1;
        AUC_red.F2.(tempevent)(b,1)=AUCstats.mean.(templabel).red.(tempevent).F2;
        AUC_signal.F1.(tempevent)(b,1)=AUCstats.mean.(templabel).signal.(tempevent).F1;
        AUC_signal.F2.(tempevent)(b,1)=AUCstats.mean.(templabel).signal.(tempevent).F2;
    end
end

%Make list of all sessions included
sessionslist=fieldnames(peakstats.mean);
for b=1:size(sessionslist);
    eventlist=fieldnames(peakstats.mean.(sessionslist{b}).red);
    for c=1:size(eventlist);
        tempevent=eventlist{c};
        templabel=sessionslist{b};
        peak_red.F1.(tempevent)(b,1)=peakstats.mean.(templabel).red.(tempevent).F1;
        peak_red.F2.(tempevent)(b,1)=peakstats.mean.(templabel).red.(tempevent).F2;
        peak_signal.F1.(tempevent)(b,1)=peakstats.mean.(templabel).signal.(tempevent).F1;
        peak_signal.F2.(tempevent)(b,1)=peakstats.mean.(templabel).signal.(tempevent).F2;
    end
end

%% Save
cd(saveDir)
save([group, '-AUCandpeaks'])