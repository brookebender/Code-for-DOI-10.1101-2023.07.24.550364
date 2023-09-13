clear all; close all; clc;
%% user modified variables
% dataDir = 'C:\Users\Sierra\Documents\MATLAB\Photometry Data laptop\Analyzed Files\PFCGABASNFR_analyzed';
dataDir = '/Users/brookebender/Documents/MATLAB/Photometry/Imported data';
saveDir='/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';
% saveDir = 'C:\Users\Sierra\Documents\MATLAB\Photometry Data laptop\Analyzed Files\PFCGABASNFR_analyzed';
EventFile = 'M6_5_31_STIMDay2-events.mat'; 
TracesFile = 'M6_5_31_STIMDay2-traces.mat';
%^REPLACE FILE NAMES FOR EACH FILE


%% change directory and import data
cd(dataDir)
load(EventFile)
clear dataDir EventFile


Event=cellstr(Event); %Convert Event names into a cell array
EventIntervals=[StartTimes EndTimes]; %Taking the times and TTL high/low output and putting them together

if any(strcmp(Event,'House Light'))
idxGP1 = find(strcmp(Event, 'House Light')); %find each occurance of each event and place in new array... there's a faster way to do this but i'll figure it out later
Input1=[idxGP1(1,:) idxGP1(end)]; % taking the first through last occurance of events and making a new array.
TTLGP1=EventIntervals((Input1(1,1):Input1(1,2)),:); %New array with just the TTL timestamps that correspond to each GP port
TS.HouseLighttimes=TTLGP1(1:end,1);
else
end

%BB edits
TS.Tonetimes=TTLGP1(:,1)+30;                                     %tones occur 30 seconds after houselight turns on
TS.CueLighttimes=TTLGP1(:,1)+60;                                 %cues occur 60 seconds after houselight turns on
TS.Levertimes=TTLGP1(:,1)+90;                                    %lever insertions occur 90 seconds after houselight turns on
TS.CueandTonetimes=TTLGP1(:,1)+120;                              %cues+tones occur 120 seconds after houselight turns on
startofsession=TTLGP1(1,1)-30;                                  %session starts 30 seconds before first houselight on time
TS.Infusiontimes=startofsession+720;
for a=[2:4];
    TS.Infusiontimes(a,1)=TS.Infusiontimes(1,1)+60*(a-1);
end 

%BB edits
fields = fieldnames(TS);
TS = rmfield(TS, fields(structfun(@isempty, TS))); %remove any empty arrays in TS

% idxGP4 = find(strcmp(Event, ' GPIO-4'));
% Input4=[idxGP4(1,:) idxGP4(end)];
% TTLGP4=EventIntervals((Input4(1,1):Input4(1,2)),:);

clear Input1 Event idxGP1

dataDir = '/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';
load(TracesFile)
%% initiate variables (user modifiable)
window.whole = [-3, 30];                                                     %set the window in seconds surrounding whole, pre-event, and post-event periods
window.pre = [-3, 0];
window.post = [0, 30];


frameRate = 30;                                                           %frame rate of FP camera; it will be 1/2 the set frame rate if in interleaved mode, or the set frame rate if in constant
frameDuration = 1/frameRate;                                                %get frame duration by getting reciprocal of frame rate
clear frameRate
%% get all calcium data surrounding behavioral events and stimuli

fnEventType = fieldnames(TS);                                               %get fieldnames of TS
fnWindow = fieldnames(window);
fnRegion = fieldnames(traces.signal);                                              %get fieldnames of signal
twocolor={'signal','red'};

% exist motionTS var
% if ans==1

motioncheck=exist('motionTS', 'var');
if motioncheck==1
motionTSd=motionTS./30;                                                     %Putting the artifact timestamps into seconds

     %Find event TS that happen during artifacts and remove from analysis
     for d=1:size(fnEventType,1)
         for m=1:size(motionTSd,1)
             tempmot=TS.(fnEventType{d});
             remove=find(tempmot>=motionTSd(m,1)& tempmot<=motionTSd(m,2));
             tempmot(remove)=[];
             TS.(fnEventType{d})=tempmot;
         end
     end
else
end
     if isempty(inclusion)                                                   %getting rid of pesky prblem later on, remove this if you actually need to exclude traces.
         inclusion=1;
     else
     end
     
for g=1:size(fnEventType,1)                 %If a TS becomes empty after deleting events that occur in a motion artifact, this fixes it
    if isempty(TS.(fnEventType{g}))
        field=fnEventType(g);
        TS=rmfield(TS,field);
    else
    end
end

fnEventType = fieldnames(TS);  
     
     for p=1:length(twocolor)
         for e = 1:size(fnEventType,1)                                               %cycle through event type (e)
             for w = 1:size(fnWindow,1) %cycle through window type (w)
                 winDuration = ceil((abs(window.(fnWindow{w})(1,1)) +...             %calculate the max number of samples in a given window of interest based on the window size and frameDuration
                     window.(fnWindow{w})(1,2))/frameDuration);
                 for r = 1:size(fnRegion,1)                                          %cycle through regions (r)
                     for y = 1:size(TS.(fnEventType{e}),1)                           %cycle through every (y) timestamp for an event category (e)
                         tempTS = timestamp;
                         idx = tempTS >= (TS.(fnEventType{e})(y,1) +...              %index the timestamps in the window of a given event
                             window.(fnWindow{w})(1,1)) & tempTS <=...
                             (TS.(fnEventType{e})(y,1) + window.(fnWindow{w})(1,2));
                         event.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r})...       %prefill the matrix with NaNs (solves ragged trial problem)
                             (y,1:winDuration) = NaN;
                         tempData = traces.(twocolor{p}).(fnRegion{r}).zscorefourier(idx);                       %put indexed trace data into a temporary vector
%                          tempData=tempData(1:winDuration);                             %occasionally getting a PETH that has an extra frame on the end.
                         ySize = size(tempData,1);                                   %calculate the number of columns in the working/incoming trial
                         event.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r})...       %fill the NaN'ed row with the indexed values through the size of that trial
                             (y,1:ySize) = tempData;
                         AUC.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r})(y,1) =...  %calculate area under the curve
                             trapz(tempData);
                     end

                     len(r,1) = length(event.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...       %store the length of each region
                         .(fnRegion{r}));
                 end
                 lenMax = max(len);                                                  %find the maximum length of all the regions
                 for k = 1:size(fnRegion,1)                                          %for region r, determine if the size is smaller than the maximum length
                     if size(event.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...
                             .(fnRegion{k}),2) < lenMax
                         event.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{k})...       %if it is smaller, put a NaN in the maximum length column
                             (:,lenMax) = NaN;
                     else
                     end
                 end
             end
         end
     end

clear winDuration frameDuration temp* idx e r y ySize w p k len*

%% calculate calcium trace means for behavioral events and stimuli
for p=1:length(twocolor)
    for e = 1:size(fnEventType,1)                                               %cycle through event type (e)
        for w = 1:size(fnWindow,1)                                              %cycle through window type (w)
            for r = 1:size(fnRegion,1)
                 nu=size(event.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r}),1);                                              
                if nu==1
                  eventMean.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r})= event.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r});
                else
                eventMean.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r}) =...     %calculate the mean of each event type (e) for each region (r), excluding the NaNs
                    nanmean(event.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r}));
                idx=~isnan(event.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r})); %index the number of values that are not NaNs
                n = sum(idx,1);                                                             %get the number of non-NaN values for proper n
                stdDev = nanstd(event.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...         %calculate the standard deviation of each event type (e) for each region (r), excluding the NaNs
                    .(fnRegion{r}));
                eventSEM.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r}) =...      %calculate the SEM using the NaN stand deviation and non-NaN n
                    stdDev./sqrt(n);
                end
            end 
        end
                    %figure
                    %plot(eventMean.(twocolor{p}).(fnEventType{e}).whole.F1)
    end
end

clear e w r idx n stdDev p
%% calculate AUC mean and SEM
for p=1:length(twocolor)
for e = 1:size(fnEventType,1)                                               %cycle through event types (e)
    for w = 1:size(fnWindow,1)                                              %cycle through window type (w)
        for r = 1:size(fnRegion,1)                                          %cycle through region (r)
            AUCmean.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r}) =...       %calculate the mean of each AUC
                nanmean(AUC.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r}));
            idx = ~isnan(AUC.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r})); %index the number of values that are not NaNs
            n = sum(idx,1);                                                 %get the number of non-NaN values for proper n
            stdDev = nanstd(AUC.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...           %calculate the standard deviation of each AUC, excluding the NaNs
                .(fnRegion{r}));              
            AUCSEM.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r}) =...        %calculate the SEM using the NaN stand deviation and non-NaN n
                stdDev./sqrt(n);          
        end
    end
end
end

clear e w r p idx n stdDev fnEventType fnWindow fnRegion

%% enter group
userInput = input('Enter 1 for FR, 2 for SO: ');
disp(userInput)
if userInput == 1
    training = 1; %FR
    groupname='FR'; 
elseif userInput == 2
    training = 2; %SO
    groupname='SO';
else
    disp('Error; redo group entry ')
    pause
end
clear userInput

%% enter sex
userInput = input('Enter 1 for male, 2 for female: ');
disp(userInput)
if userInput == 1
    sex = 1; %male
    sexname='male'; 
elseif userInput == 2
    sex = 2; %female
    sexname='female';
else
    disp('Error; redo sex entry ')
    pause
end
clear userInput

%% enter training
    trainingname='STIM';


%% clear variables and save data
cd(saveDir)
clearvars -except animalID event* window group inclusion AUC* TS traces sexname groupname trainingname
save([animalID, sexname, groupname, trainingname, '-events'])