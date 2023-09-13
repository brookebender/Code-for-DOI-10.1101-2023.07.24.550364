%%Updated Ahmari lab code from Zoe to work with Plexon data output files.
%Remember to replace Rat_Date_Program in line 12 - must run for each session
close all 
clear 

%% user modified variables
%Get data from "Imported data" folder
dataDir = '/Users/brookebender/Documents/MATLAB/Photometry/BB_P2/Imported data';
%Set save folder to "Analyzed files" folder
saveDir = '/Users/brookebender/Documents/MATLAB/Photometry/BB_P2/Analyzed files';
% *REPLACE RAT_DATE_PROGRAM below*
fileName = 'F1_FR3_CNO-Imported.mat';                                               


slidingTime = 120;                                                          %time in seconds for sliding window normalization
motionAvg = 50; %number of values to include for motion artifact average interpolation value
%% change directory and import data
cd(dataDir)
load(fileName)
clearvars -except RawF_465_F1 RawF_465_F2 RawF_410_F1 RawF_410_F2 RawF_560_F1 RawF_560_F2 Frame timestamp motionAvg slidingTime animalID saveDir

% find Nan values that can sometimes occur during the last timepoint
TF(1,:)=isnan(RawF_410_F1);
TF(2,:)=isnan(RawF_465_F1);
TF(3,:)=isnan(RawF_560_F1);
TF(4,:)=isnan(RawF_410_F2);
TF(5,:)=isnan(RawF_465_F2);
TF(6,:)=isnan(RawF_560_F2);

% If the final value of at least one channel is a Nan value, erase the last timepoint
emptyvar=isempty(TF);
if emptyvar==0
RawF_410_F1(end,:)=[];
RawF_465_F1(end,:)=[];
RawF_560_F1(end,:)=[];
RawF_410_F2(end,:)=[];
RawF_465_F2(end,:)=[];
RawF_560_F2(end,:)=[];
Frame(end,:)=[];
timestamp(end,:)=[];
else
end
%% Create alignTraces structure 
alignTraces.signal.F1=RawF_465_F1;  %"signal" is the green channel (dLight in my case), F1=fiber 1=DLS, and F2=fiber2=DMS
alignTraces.isosbestic.F1=RawF_410_F1;
alignTraces.red.F1=RawF_560_F1;
alignTraces.signal.F2=RawF_465_F2;
alignTraces.isosbestic.F2=RawF_410_F2;
alignTraces.red.F2=RawF_560_F2;
fnRegion = fieldnames(alignTraces.signal); %fnRegion refers to Fiber 1 or Fiber 2 (DLS or DMS)
fnSignal=fieldnames(alignTraces);
%% indicate which regions have signal
% Plot the entire signals for F1 (Figure 1) and F2 (Figure 2)
figure
plot(alignTraces.isosbestic.F1)
hold on
plot(alignTraces.signal.F1)
plot(alignTraces.red.F1)
legend('isosbestic','signal','red')

figure
plot(alignTraces.isosbestic.F2)
hold on
plot(alignTraces.signal.F2)
plot(alignTraces.red.F2)
legend('isosbestic','signal','red')

% Generate "inclusion" variable where position (1,1) indicates whether there is signal in F1 and position (2,1) indicates whethere there is signal in F2
userInput = input('1 for signal in F1, 0 for no signal in F1: ');%indicate if the region has signal
inclusion = userInput;
userInput = input('1 for signal in F2, 0 for no signal in F2: ');%indicate if the region has signal
inclusion(2,1) = userInput;
    
%store index of traces you should include (1) or exclude (0) for future use
close all
fourier=0;
expdetrend=0;
clear r userInput expInput fourierInput
%% ensure all timestamp and signal vectors are the same length
for r = 1:size(fnRegion, 1)
    num(1,r) = min([length(alignTraces.signal.(fnRegion{r}))...
        ,length(alignTraces.isosbestic.(fnRegion{r})),length(alignTraces.red.(fnRegion{r}))]);                   %determine which trace (green, red, or isosbestic) has the fewest frames/shortest length
end
len = min(num);
for r = 1:size(fnRegion,1)                                                  %make all trace vectors equal to that length and put into new separated structure
    traces.signal.(fnRegion{r}).raw =...
        alignTraces.signal.(fnRegion{r})(1:len);
    traces.isosbestic.(fnRegion{r}).raw =...
        alignTraces.isosbestic.(fnRegion{r})(1:len);
    traces.red.(fnRegion{r}).raw=alignTraces.red.(fnRegion{r})(1:len);
    
end
timestampsTS=timestamp(1:len);
clear alignTraces alignTracesTS r len num

%% build a low pass filter
lpFilt = designfilt('lowpassiir','FilterOrder',6, ...
    'PassbandFrequency',3, ...
    'SampleRate',30); %zoes was 13.333, but the plexon is 30. You can double check by dividing the length of your Frame vector (which is the total# of frames) by 30(fps), should = the length of your session in seconds.
%% forward and reverse filter the traces (avoids phase shift)
for r = 1:size(fnRegion,1)                                                  %filter trace for each region (r)
    traces.signal.(fnRegion{r}).filt = filtfilt(lpFilt,traces.signal.(fnRegion{r}).raw);
    traces.isosbestic.(fnRegion{r}).filt=filtfilt(lpFilt,traces.isosbestic.(fnRegion{r}).raw);
    traces.red.(fnRegion{r}).filt = filtfilt(lpFilt,traces.red.(fnRegion{r}).raw);
    figure
    plot(traces.signal.(fnRegion{r}).filt)
    hold on
    plot(traces.isosbestic.(fnRegion{r}).filt)
    plot(traces.red.(fnRegion{r}).filt)
    legend('signal','isosbestic','red')
    pause
end
clear r lpFilt
close all

%% identify motion artifacts for removal
for r = 1:size(fnRegion,1)                                                  %cycle through region r
    tempSignal = traces.signal.(fnRegion{r}).filt;
    tempSignal=tempSignal';
    tempIso = traces.isosbestic.(fnRegion{r}).filt;
    tempIso=tempIso';
    tempRed = traces.red.(fnRegion{r}).filt;
    tempRed=tempRed';
    figure                                                                  
    plot(tempIso-mean(tempIso))
    hold on
    plot(tempSignal-mean(tempSignal))  
    plot(tempRed-mean(tempRed))
    userInput = input('Are there motion artifacts? 1 for yes, 2 for no: '); %indicate if there are motion artifact
    if userInput == 1                                                       %if there are motion artifacts
        motionNum = input...
            ('Number of motion artifacts: ');                               %indicate the number of motion artifacts
        for i = 1:motionNum
            motionSize = input('Press 1 for shift, 2 for artifact: ');      %a shift indicates a jump in signal and will be removed and evened out. Artifacts are only removed.
            motion.type(r,i) = motionSize;
            disp('Zoom to start of motion artifact, then press enter ');    %zoom to start of motion artifact i and then press enter
            pause
            [startX,~] = ginput(1);                                         %select start with crosshairs
            disp('Zoom to end of motion artifact, then press enter ');      %zoom to end of motion artifact i and then press enter
            pause
            [stopX,~] = ginput(1);                                          %select end with crosshairs
            close
            if startX<=1
                startTS=1;
            else
            startTS = floor(startX);                                        %round down to nearest integer for start column number
            end
            stopTS = ceil(stopX);                                           %round up to nearest integer for stop column number
            motionTS(i,:)= [startTS stopTS];
            if stopTS > size(tempSignal,2)
                stopTS = size(tempSignal,2);
            end
            if startTS==1
                stopGroup = tempSignal...                                      %get motionAvg number of signals after stop column
                    (1,stopTS:(stopTS+motionAvg));
            stopIso=tempIso(1,stopTS:(stopTS+motionAvg));
                stopRed=tempRed(1,stopTS:(stopTS+motionAvg));
             replaceMean = mean(stopGroup);   
             replaceMeanIso=mean(stopIso);   
             replaceMeanRed=mean(stopRed);
             tempSignal(1,startTS:stopTS) = replaceMean;
             tempIso(1,startTS:stopTS) = mean(stopIso);
             tempRed(1,startTS:stopTS) = mean(stopRed);
            else
            startGroup = tempSignal...                                      %get motionAvg number of signals before start column
                (((startTS-motionAvg):startTS));
            startIso=tempIso(((startTS-motionAvg):startTS));    
            startRed=tempRed(((startTS-motionAvg):startTS)); 

            if stopTS+motionAvg >= size(tempSignal,2)
                if motionSize == 1
                    replaceMean = mean(startGroup);
                    tempSignal(1,startTS:end) = mean(startGroup);
                    replaceMeanIso=mean(startIso);
                    tempIso(1,startTS:end)=mean(startIso);
                      replaceMeanRed=mean(startRed);
                    tempRed(1,startTS:end)=mean(startRed);
                elseif motionSize == 2
                    tempSignal(1,startTS:end) =...                           %replace motion artifact values (between start and stop) with the average of the replacement value vector
                        mean(startGroup);
                    tempIso(1,startTS:end)=mean(startIso);
                    tempRed(1,startTS:end)=mean(startRed);
                end
            else
                stopGroup = tempSignal...                                      %get motionAvg number of signals after stop column
                    (1,stopTS:(stopTS+motionAvg));
                stopIso=tempIso(1,stopTS:(stopTS+motionAvg));
                stopRed=tempRed(1,stopTS:(stopTS+motionAvg));
                if motionSize == 1
                    replaceMean = mean(startGroup) - mean(stopGroup);
                    replaceMeanIso = mean(startIso)- mean(stopIso);
                    replaceMeanRed = mean(startRed)- mean(stopRed);
                    
                    tempSignal(1,stopTS:end) = tempSignal(1,stopTS:end) +...
                        replaceMean;
                    tempIso(1,stopTS:end) = tempIso(1,stopTS:end) + replaceMeanIso;
                    tempRed(1,stopTS:end) = tempRed(1,stopTS:end) + replaceMeanRed;
                    
                    replaceGroup = [startGroup (stopGroup + replaceMean)];
                    replaceGroupIso = [startIso (stopIso + replaceMeanIso)];
                    replaceGroupRed = [startRed (stopRed + replaceMeanRed)];
                    
                    tempSignal(1,startTS:stopTS) = mean(replaceGroup);
                    tempIso(1,startTS:stopTS) = mean(replaceGroupIso);
                    tempRed(1,startTS:stopTS) = mean(replaceGroupRed);
                    
                elseif motionSize == 2
                    replaceGroup = [startGroup stopGroup];                      %make new vector with both start and stop values
                    replaceGroupIso=[startIso stopIso];
                     replaceGroupRed=[startRed stopRed];
                    tempSignal(1,startTS:stopTS) =...                           %replace motion artifact values (between start and stop) with the average of the replacement value vector
                        mean(replaceGroup);
                   tempIso(1,startTS:stopTS)=mean(replaceGroupIso);
                   tempRed(1,startTS:stopTS)=mean(replaceGroupRed);
                end
            end
            end
            traces.signal.(fnRegion{r}).motion = tempSignal';
            traces.isosbestic.(fnRegion{r}).motion = tempIso';
            traces.red.(fnRegion{r}).motion = tempRed';
            motion.start(r,i) = startTS;
            motion.stop(r,i) = stopTS;
            clear tempSignal tempSub tempRed
            
            figure
            tempSignal = traces.signal.(fnRegion{r}).motion';
            tempIso= traces.isosbestic.(fnRegion{r}).motion';
            tempRed= traces.red.(fnRegion{r}).motion';
            plot(tempSignal-mean(tempSignal))                                     %plot new trace with corrected motion artifact
            hold on
            plot(tempIso-mean(tempIso))
            plot(tempRed-mean(tempRed))
            pause
        end
    else
        motion.type(r,1) = NaN;
        motion.start(r,1) = NaN;
        motion.stop(r,1) = NaN;
        traces.signal.(fnRegion{r}).motion = tempSignal';
        traces.isosbestic.(fnRegion{r}).motion=tempIso';
        traces.red.(fnRegion{r}).motion=tempRed';
        close all
    end
end
fnMotion = fieldnames(motion);
for m = 1:size(fnMotion,1)
    idx = motion.(fnMotion{m}) == 0;
    motion.(fnMotion{m})(idx) = NaN;
end
close all
clear  r i userInput start* stop* replaceGroup replaceGroupIso replaceGroupRed replaceMean replaceMeanIso replaceMeanRed motionNum temp* m fnMotion idx motionSize 

%% Fit isosbestic to signal channel and do df/f calculation
%fit isosbestic to signal using a least squares polynomial fit of degree 1.
%then use line equation: fitted isosbestic=a*isosbestic+b
for r=1:size(fnRegion,1)
        isotrace=traces.isosbestic.(fnRegion{r}).motion;
        signaltrace=traces.signal.(fnRegion{r}).motion;

  p=polyfit((isotrace),(signaltrace),1);
  a=p(1); b=p(2);
isofitgreen=a.*isotrace+b;

figure
plot(signaltrace)
hold on
plot(isofitgreen)
legend('signal','fit')

traces.isosbestic.(fnRegion{r}).fittedgreen=isofitgreen;

%% df/f calculation: signal - fitted control / fitted control

df=signaltrace-isofitgreen;
dff=df./isofitgreen;

figure
plot(dff)   

traces.signal.(fnRegion{r}).dff=dff;
end

%%Fit isosbestic to red channel an do df/f calculation
for r=1:size(fnRegion,1)
        isotrace=traces.isosbestic.(fnRegion{r}).motion;
        redtrace=traces.red.(fnRegion{r}).motion;
     p=polyfit((isotrace),(redtrace),1);
     a=p(1); b=p(2);
      isofitred=a.*isotrace+b; 
figure
plot(redtrace)
hold on
plot(isofitred)
legend('signal','fit')

traces.isosbestic.(fnRegion{r}).fittedred=isofitred;

%% df/f calculation: red - fitted control / fitted control
df=redtrace-isofitred;
dff=df./isofitred;

figure
plot(dff)   

traces.red.(fnRegion{r}).dff=dff;
end
pause
close all
%% z-score entire trace to itself
for r = 1:size(fnRegion,1)                                                  %z-score trace for each region (r)
    traces.signal.(fnRegion{r}).zscore = zscore(traces.signal.(fnRegion{r}).dff);
    figure
    plot(traces.signal.(fnRegion{r}).zscore);
%     pl=traces.signal.(fnRegion{r}).zscore(1:1000);
%     plot(pl)
end
for r = 1:size(fnRegion,1)                                                  %z-score trace for each region (r)
    traces.red.(fnRegion{r}).zscore = zscore(traces.red.(fnRegion{r}).dff);
    figure
    plot(traces.red.(fnRegion{r}).zscore);
%     pl=traces.red.(fnRegion{r}).zscore(1:1000);
%     plot(pl)
end
clear r

%% fourier 5 detrending if necessary
twocolor={'signal','red'};
for p=1:length(twocolor)
    for r = 1:size(fnRegion,1)
tempData = traces.(twocolor{p}).(fnRegion{r}).zscore;    %put data into temporary structure
[xData, yData] = prepareCurveData( [], tempData );                  %get x and y components of trace

ft = fittype( 'fourier5' );                                         %set up fit type fourier 5
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0 0 0 0 0 0 0 0 0 0 0 8.05175281242979e-05];

[fitresult, ~] = fit( xData, yData, ft, opts );                     %fit model to data

figure( 'Name', 'untitled fit 1' );                                 %plot model and data
h = plot( fitresult, xData, yData );
legend( h, 'tempData', 'untitled fit 1', 'Location', 'NorthEast' );
ylabel tempData
grid on

f = fitresult;                                                      %reassign fitresult structure to f

fitline = f.a0 + f.a1*cos(xData*f.w) + f.b1*sin(xData*f.w)...       %create fitline based on model coefficients
    + f.a2*cos(2*xData*f.w) + f.b2*sin(2*xData*f.w) +...
    f.a3*cos(3*xData*f.w) + f.b3*sin(3*xData*f.w) +...
    f.a4*cos(4*xData*f.w) + f.b4*sin(4*xData*f.w) +...
    f.a5*cos(5*xData*f.w) + f.b5*sin(5*xData*f.w);

zfourier = (tempData-fitline)';        %subtract fitline from trace
plot(zfourier)
traces.(twocolor{p}).(fnRegion{r}).zscorefourier=zfourier';
    end
end
close all
%% clear variables and save data
%Set save folder to "Analyzed files" folder
saveDir = '/Users/brookebender/Documents/MATLAB/Photometry/BB_P2/Analyzed files';
cd(saveDir);
clearvars -except animalID traces* signal inclusion motion motionTS timestamp timestampTS Frame
save([animalID, '-traces']);
