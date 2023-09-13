clear all; close all; clc;
%% user defined variables
%Get data from "Analyzed files" folder
dataDir = '/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';
%Save data to "Analyzed files" folder
saveDir = '/Users/brookebender/Documents/MATLAB/Photometry/Analyzed files';

%fileNames = ["*STIM-merged.mat"; "*FREarly-merged.mat"; "*SOEarly-merged.mat"; "*FRMid-merged.mat"; "*SOMid-merged.mat"; "*FRLate-merged.mat"; "*SOLate-merged.mat"; "*FRextcueext-merged.mat"; "*SOextcueext-merged.mat"; "*extcueext-merged.mat"; "*FRreactcueext-merged.mat"; "*SOreactcueext-merged.mat"; "*reactcueext-merged.mat"; "*FRreinst-merged.mat"; "*FRreinstReact-merged.mat"; "*SOreinst-merged.mat"; "*SOreinstReact-merged.mat"; "*cocaine-merged.mat"; "*saline-merged.mat"];
%for g=1:size(fileNames,1)                                                  %cycle through each merged file
%% change directory and import data
%cd(dataDir)
%fileName=fileNames(g,1);
%load(fileName);
%clear fileName

fileName = '*SOReinst-merged.mat';
%% change directory and import data
cd(dataDir)
load(fileName)
clear dataDir fileName
%% create trial name matrix
binDur = 0.2;                                                              %bin duration in seconds; user's choice, but works best if it's a factor of the frame rate
frameRate = 30;                                                           %frame rate for each signal; numerator is camera frame rate, denominator is number of interleaved channels
binFrameNum = binDur*frameRate;                                             %generate the number of frames for each bin by multiplying the bin duration and frame rate
clear binDur frameRate
%% get largest size for each event and event mean for each animal
fnID = fieldnames(groupData);

for i = 1:size(fnID,1)                                                      %for rat ID 
  twocolor={'signal','red'}; 
  fnEventType = fieldnames(groupData.(fnID{i}).eventMean.signal);
    for p=1:size(twocolor,2)              %get fieldnames for event type
for e = 1:size(fnEventType,1)                                           %for event type e
 
    fnWindow = fieldnames(groupData.(fnID{i}).eventMean.(twocolor{p})...              %get fieldnames for window type
            .(fnEventType{e}));
        for w = 1:size(fnWindow,1)                                          %for window type w
            fnRegion = fieldnames(groupData.(fnID{i}).eventMean.(twocolor{p})...          %get fieldnames for region
                .(fnEventType{e}).(fnWindow{w}));
            len.(fnEventType{e})(i,w) = size(groupData.(fnID{i})...         %put the length of one region (all will the same) into a length matrix for that event and window type; this will contain the lengths for each mouse by the end
                .eventMean.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{1}),2);
        end
end
    end
end
fnLen = fieldnames(len);
for l = 1:size(fnLen,1)                                                     %for every event type, find the maximum length across sessions
    for x = 1:size(len.(fnLen{l}),2)                                        %cycle throug each window type (x, in columns)
        lenMax.(fnLen{l})(1,x) = max(len.(fnLen{l})(:,x));              
    end
end
clear len fnLen i e l w p   
    
   
%% put individual trace means into group matrix, excluding no signal regions
for i = 1:size(fnID,1)    
    twocolor={'signal','red'}; 
    for p=1:size(twocolor,2)%cycle through animal IDs (i)
        fnEventType = fieldnames(groupData.(fnID{i}).eventMean.(twocolor{p}));
    for e = 1:size(fnEventType,1)                                           %cycle through event types (e)
        for w = 1:size(fnWindow,1)                                          %cycle through window types (w)
            for r = 1:size(fnRegion,1)                                      %cycle through region (r)
                yLen = lenMax.(fnEventType{e})(1,w);                        %get the maximum length across sessions for that event type
                groupEvent.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r})...  %prefill groupEvent with NaNs to the appropriate maximum length
                    (i,1:yLen) = NaN; 
                if groupData.(fnID{i}).inclusion(r,1) == 1  %if the inclusion value for the region is one, the signal is good; if not (0), exclude that region
                    tempData = groupData.(fnID{i}).eventMean.(twocolor{p})...             %put working data into temporary vector
                        .(fnEventType{e}).(fnWindow{w}).(fnRegion{r});
                    groupEvent.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...
                        .(fnRegion{r})(i,1:size(tempData,2)) = tempData; %put each animal's (i) event Mean trace data for event type (e), window type (w), and region (r) into the i row of a group matrix
                else 
                end
%                 if groupEvent.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r})...  %prefill groupEvent with NaNs to the appropriate maximum length
%                     (i,1)==0
%                 end
            end
        end
    end
    end
end
clear i e w r fnData lenMax tempData yLen p

           twocolor= fieldnames(groupEvent);  %Fixes problem where rows with event that didn't exist are filled with 0 instead of NaN
           for p=1:size(twocolor,1)
               fnEventType = fieldnames(groupEvent.(twocolor{p}));
               for e = 1:size(fnEventType,1)
                   for w = 1:size(fnWindow,1)
                       for r = 1:size(fnRegion,1)
                           tempData=groupEvent.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r});
                           for s=1:size(tempData,1)
                               if sum(tempData(s,:))==0
                               tempData(s,:)=NaN;
                               groupEvent.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r})=tempData;
                               else   
                               end
                           end
                       end
                   end
               end
           end
           
clear i e w r fnData lenMax tempData yLen p

%% put individual AUC means into group matrix, excluding no signal regions
for i = 1:size(fnID,1)                                                      %cycle through animal IDs (i)
    twocolor={'signal','red'};
    for p=1:size(twocolor,2)
         fnEventType = fieldnames(groupData.(fnID{i}).eventMean.(twocolor{p}));
    for e = 1:size(fnEventType,1)                                           %cycle through event types (e)
        for w = 1:size(fnWindow,1)                                          %cycle through window types (w)
            for r = 1:size(fnRegion,1)                                      %cycle through region (r)
                groupAUC.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...                  %prefill the AUC with NaN
                    .(fnRegion{r})(i,1) = NaN;
                if groupData.(fnID{i}).inclusion(r,1) == 1                  %if the inclusion value for the region is one, the signal is good; if not (0), exclude that region
                    groupAUC.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...              %put each animal's (i) AUC for event type (e), window type (w), and region (r) into the i row of a group matrix
                        .(fnRegion{r})(i,1) = groupData.(fnID{i})...
                    .AUCmean.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r});   
                else
                end
            end
        end
    end
    end
end
clear i e w r p

           twocolor= fieldnames(groupAUC);   %Fixes problem where rows with event that didn't exist are filled with 0 instead of NaN
           for p=1:size(twocolor,1)
               fnEventType = fieldnames(groupAUC.(twocolor{p}));
               for e = 1:size(fnEventType,1)
                   for w = 1:size(fnWindow,1)
                       for r = 1:size(fnRegion,1)
                           tempData=groupAUC.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r});
                           for s=1:size(tempData,1)
                               if sum(tempData(s,:))==0
                               tempData(s,:)=NaN;
                               groupAUC.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r})=tempData;
                               else   
                               end
                           end
                       end
                   end
               end
           end
           
           clear i e w r p
%% bin event traces
for e = 1:size(fnEventType,1)                                               %cycle through event type (e)
    twocolor={'signal','red'};
    for p=1:size(twocolor,2)
    for w = 1:size(fnWindow,1)                                              %cycle through window types (w)
        for r = 1:size(fnRegion,1)                                          %cycle through region (r)
            tempData = groupEvent.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...         %put data into temporary vector    
                .(fnRegion{r});              
            binNum = size(tempData,2)/binFrameNum;                          %get the number of bins by dividing the temporary vector by the binFrameNum
            for i = 1:size(tempData,1)                                      %cycle through number of columns (i)
                for x = 1:binNum
                    x1 = x*binFrameNum - (binFrameNum-1);                   %cycle through the bin numbers (x) and make x1 the start point and x2 the endpoint for each bin
                    x2 = x*binFrameNum;
                    bin.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r})...     %make bin (x) the nanmean of the data between x1 and x2
                        (i,x) = nanmean(tempData(i,x1:x2));
                end
            end
        end
    end
    end
end
clear tempData e w r i x x1 x2 binNum binFrameNum
%% get mean and SEM for all group trace data

for k=1:size(fnID,1)                                                %makes fneventtype include all events even if they don't appear for every rat
    tempevents=fieldnames(groupData.(fnID{k}).eventMean.signal);
    tempsize=size(tempevents,1);
    eventsize(k,1)=tempsize;
    allevents=max(eventsize);
    if tempsize == allevents
    fnEventType=tempevents;
    else
    end
end

for e = 1:size(fnEventType,1)                                               %cycle through event type (e)
    for p=1:size(twocolor,2)
    for w = 1:size(fnWindow,1)                                              %cycle through window types (w)
        for r =1:size(fnRegion,1)                                           %cycle through region (r)
            tempData =  groupEvent.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...        %put working data into temporary vector
                .(fnRegion{r});             
            groupSummary.mean.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...             %put the mean (excluding NaNs) of event type (e), region (r) into a summary structure
                .(fnRegion{r}) = nanmean(tempData);
            idx = ~isnan(tempData);                                         %index all the non-Nan numbers
            n = sum(idx,1);                                                 %sum all the indexed numbers to get the non-NaN n
            stdDev = nanstd(tempData);                                      %get the standard deviation of tempData, excluding NaNs
            groupSummary.SEM.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...              %calculate the SEM and put it into a summary structure
                .(fnRegion{r})= stdDev./sqrt(n); 
        end
    end    
    end
end
clear e w r tempData idx n stdDev p
%% get mean and SEM for AUC
for e = 1:size(fnEventType,1)                                               %cycle through event type (e)
    for p=1:size(twocolor,2)
    for w = 1:size(fnWindow,1)                                              %cycle through window types (w)
        for r = 1:size(fnRegion,1)                                          %cycle through region (r)
            tempData = groupAUC.(twocolor{p}).(fnEventType{e}).(fnWindow{w})...           %put working data into temporary vector
                .(fnRegion{r});
            groupAUC.mean.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r}) =... %put the mean (excluding NaNs) of event type (e), region (r) into a summary structure
                nanmean(tempData);
            idx = ~isnan(tempData);                                         %index all the non-Nan numbers
            n = sum(idx,1);                                                 %sum all the indexed numbers to get the non-NaN n
            stdDev = nanstd(tempData);                                      %get the standard deviation of tempData, excluding NaNs
            groupAUC.SEM.(twocolor{p}).(fnEventType{e}).(fnWindow{w}).(fnRegion{r}) =...  %calculate the SEM and put it into a summary structure
                stdDev./sqrt(n);
        end
    end
    end
end
clear e w r n fnEventType stdDev fnRegion fnWindow idx p
%% change directory and save
cd(saveDir)
clearvars -except group* window bin dataDir saveDir fileNames
save([groupchar, '-groupStats'])
%end

