clear all; close all; clc;
%% user modified variables
%Get data from "Imported data" folder
dataDir = '/Users/brookebender/Documents/MATLAB/Photometry/Imported data';
%Set save folder to "Analyzed files" folder
saveDir='/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';
% *REPLACE RAT_DATE_PROGRAM below*
EventFile = 'F6_3_3_FR5Day2-events.mat'; 
TracesFile = 'F6_3_3_FR5Day2-traces.mat';

%% change directory and import data
cd(dataDir)
%Load event file from "Imported data" folder
load(EventFile)
clear dataDir EventFile

%% Extract event times
Event=cellstr(Event); %Convert Event names into a cell array
EventIntervals=[StartTimes EndTimes]; %Taking the times and TTL high/low output and putting them together

if any(strcmp(Event,'House Light'))
idxGP1 = find(strcmp(Event, 'House Light')); %find each occurance of each event and place in new array
Input1=[idxGP1(1,:) idxGP1(end)]; % taking the first through last occurance of events and making a new array
TTLGP1=EventIntervals((Input1(1,1):Input1(1,2)),:); %New array with just the TTL timestamps that correspond to each 
TS.HouseLightofftimes=TTLGP1(1:end,2); %Puts houselight off times into an array
TS.HouseLightontimes=TTLGP1(1:end,1); %Puts houselight on times into an array
else
end

if any(strcmp(Event,'Active Lever'))
idxGP2 = find(strcmp(Event, 'Active Lever'));
Input2=[idxGP2(1,:) idxGP2(end)];
TTLGP2=EventIntervals((Input2(1,1):Input2(1,2)),:);
TS.ActiveLevertimes=TTLGP2(1:end,1); %array has timestamps for all active lever presses
else
end

if any(strcmp(Event,'Inactive Lever'))
idxGP3 = find(strcmp(Event, 'Inactive Lever'));
Input3=[idxGP3(1,:) idxGP3(end)];
TTLGP3=EventIntervals((Input3(1,1):Input3(1,2)),:);
TS.Inacttimes=TTLGP3(1:end,1); %array has timestamps for all inactive lever presses
else
end

if exist('TTLGP2', 'var');      %Prevent the next step from occuring if there are no active lever presses
    actcounter=1;
else
    actcounter=0;
end
if actcounter>0
timeuntilnextpress=diff(TTLGP2(:,1));        %Make an array where each value is the time after each lever press
timesincelastpress=[timeuntilnextpress(1:0,:); TTLGP2(1,1)-TTLGP1(1,1); timeuntilnextpress(0+1:end,:)];  %Make an array where each value is the time since the last press, with the first value time between when the session started and the first press occurred
timeuntilnextpress(end+1,:)=(TTLGP1(end,2));                                                 %Make the last value when the house light last turns off
isolatedpresses=find((timesincelastpress(:,1)>=3)&(timeuntilnextpress(:,1)>=5));           %Find lever presses with no press in 3 seconds prior or 5 seconds after
TS.IsolatedActiveLeverPresses=TTLGP2(isolatedpresses,1); %array has timestamps for active lever presses with no press in 3 seconds prior or 5 seconds after
else 
end

if exist('TTLGP3', 'var');      %Prevent the next step from occuring if there are no inactive lever presses
    inactcounter=1;
else
    inactcounter=0;
end

if inactcounter>0
timeuntilnextpress2=diff(TTLGP3(:,1));        %Make an array where each value is the time after each lever press
timesincelastpress2=[timeuntilnextpress2(1:0,:); TTLGP3(1,1)-TTLGP1(1,1); timeuntilnextpress2(0+1:end,:)];  %Make the first value time since session start (which is the first value)
timeuntilnextpress2(end+1,:)=(TTLGP1(end,2));                                                 %Make the last value when the house light last turns off
isolatedpresses2=find((timesincelastpress2(:,1)>=3)&(timeuntilnextpress2(:,1)>=5));           %Find lever presses with no press in 3 seconds prior or 5 seconds after
TS.IsolatedInactiveLeverPresses=TTLGP3(isolatedpresses2,1); %array has timestamps for inactive lever presses with no press in 3 seconds prior or 5 seconds after
else
end

if any(strcmp(Event,'Cue Light'))
idxGP4 = find(strcmp(Event, 'Cue Light'));
Input4=[idxGP4(1,:) idxGP4(end)];
TTLGP4=EventIntervals((Input4(1,1):Input4(1,2)),:);
TS.CueLighttimes=TTLGP4(1:end,1); %array has timestamps for all cue light occurances
else
end

if exist('TTLGP4', 'var');      %Prevent the next step from occuring if there are no cue lights
    cuecounter=1;
else
    cuecounter=0;
end

if cuecounter>0
timeuntilnextcue=diff(TTLGP4(:,1));        %Make an array where each value is the time after each cue
timesincelastcue=[timeuntilnextcue(1:0,:); TTLGP4(1,1)-TTLGP1(1,1); timeuntilnextcue(0+1:end,:)];  %Make the first value time since session start (which is the first value)
timeuntilnextcue(end+1,:)=(TTLGP1(end,2));                                                 %Make the last value when the house light last turns off
isolatedcues=find((timesincelastcue(:,1)>=3)&(timeuntilnextcue(:,1)>=5));           %Find cue with no cue in 3 seconds prior or 5 seconds after
TS.IsolatedCues=TTLGP4(isolatedcues,1); %array has timestamps for cues with no cue in 3 seconds prior or 5 seconds after

lengthofcue=TTLGP4(:,2)-TTLGP4(:,1);        %make a matrix that subtracts cue start times from cue off times to find length of cue
cueslong=find(lengthofcue(:,1)>18);         %make a matrix with the positions of cues longer than 18 seconds (~20 second cues)
cuesshort=find(lengthofcue(:,1)<2);        %make a matrix with the positions of cues shorter than 2 seconds (~1 second cues)
TS.LongCueLighttimes=TTLGP4(cueslong,1);    %array has timestamps for all long cue lights
TS.ShortCueLighttimes=TTLGP4(cuesshort,1);  %array has timestamps for all short cue lights
else
end

if isfield (TS, 'IsolatedCues')
    if isfield (TS, 'ShortCueLighttimes')
IC=TS.IsolatedCues;    %Matrix of all isolated cues
SCLT=TS.ShortCueLighttimes;    %Matrix of all shortcuetimes
posmatch1=ismember(IC,SCLT);  %Rows with 1 indicate matching value
posIC=find(posmatch1==1);  %Indicates position of isolated cues that are also short
TS.IsolatedShortCues=IC(posIC, 1);   %array has timestamps for short cues with no cue in 3 seconds prior or 5 seconds after
else
end
else
end

if isfield (TS, 'IsolatedCues')
    if isfield (TS, 'LongCueLighttimes')
IC=TS.IsolatedCues;    %Matrix of all isolated cues
LCLT=TS.LongCueLighttimes;    %Matrix of all shortcuetimes
posmatch1=ismember(IC,LCLT);  %Rows with 1 indicate matching value
posIC=find(posmatch1==1);  %Indicates position of isolated cues that are also short
TS.IsolatedLongCues=IC(posIC, 1);   %array has timestamps for short cues with no cue in 3 seconds prior or 5 seconds after
else
end
else
end

if actcounter>0      %Prevent next step from happening if there are no active lever presses
    if cuecounter>0
B = repmat(TS.ActiveLevertimes,[1 length(TS.CueLighttimes)]); %Find timestamps in in TS.ActiveLevertimes that match up with timestamps in TS.CueLighttimes
[minValue,closestIndex] = min(abs(B-TS.CueLighttimes'),[],1);
closestValue = TS.ActiveLevertimes(closestIndex');
TS.UnrfActiveLevertimes=setdiff(TS.ActiveLevertimes,closestValue); %Array has timestamps for all unreinforced active lever presses
else
  TS.UnrfActiveLevertimes=TS.ActiveLevertimes;
    end
else 
end

if isfield (TS, 'UnrfActiveLevertimes')
    if isfield (TS, 'IsolatedActiveLeverPresses')
UAL=TS.UnrfActiveLevertimes;    %Matrix of all unreinforced active lever presses
IAL=TS.IsolatedActiveLeverPresses;    %Matrix of all isolated active lever presses
posmatch=ismember(UAL,IAL);  %Rows with 1 indicate matching value
posUAL=find(posmatch==1);  %Indicates position of unreinforced active lever presses that are also isolated
TS.IsolatedUnrfActiveLeverPresses=UAL(posUAL, 1);   %array has timestamps for unreinforced active lever presses with no press in 3 seconds prior or 5 seconds after
else
end
else
end

fields = fieldnames(TS);
TS = rmfield(TS, fields(structfun(@isempty, TS))); %remove any empty arrays in TS

clear Input1 Input2 Input3 Input4 Event idxGP1 idxGP2 idxGP3 idxGP4 actcounter B closestIndex closestValue cue* endtimes eventintervals inactcounter isolated* lengthofcue posUAL posmatch time* UAL IAL

dataDir = '/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';
cd(dataDir)
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
userInput = input('Enter 1 for early, 2 for mid, 3 for late: ');
disp(userInput)
if userInput == 1
    training = 1; %early
    trainingname='early'; 
elseif userInput == 2
    training = 2; %mid
    trainingname='mid';
elseif userInput==3
    training = 3; %late
    trainingname='late';
else
    disp('Error; redo training entry ')
    pause
end
clear userInput


%% clear variables and save data
%Set save folder to "Analyzed files" folder
saveDir='/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';
cd(saveDir)
clearvars -except animalID event* window group inclusion AUC* TS traces groupname sexname trainingname
save([animalID, sexname, groupname, trainingname, '-events'])






