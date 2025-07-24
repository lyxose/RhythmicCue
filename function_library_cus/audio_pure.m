% old verision; unsupport different tSOA
clear;
clc;
sca;
KbName('UnifyKeyNames'); 

if ~exist('./Data', 'dir')
    mkdir('./Data');
end
if ~exist('./Data/Interrupted','dir')
    mkdir('./Data/Interrupted')
end
saveRaw = true;  % whether to save all variables, if not, only the tabel of results will be saved.
DTstr = datestr(datetime, 'yyyymmddTHHMM');

addpath('function_library_cus');

[keyIsDown, ~, keyCode] = KbCheck;
if sum(keyCode) ~= 0    % make sure only one key is pressed
    error('Check keyboard hardware!!! A certain key is pressed! %s',KbName(keyCode))
end
% get SubjInfo
[groupID, subjID, session, location, subjName, subjGender, subjAge, gstgAmp, seqID] = InformationBox;

% block sequence should be randomized across subjects
seqTypes = {[1 2 3 4],[2 1 3 4],[3 2 1 4],[4 2 3 1]
            [1 2 4 3],[2 1 4 3],[3 2 4 1],[4 2 1 3]
            [1 3 2 4],[2 3 1 4],[3 1 2 4],[4 3 2 1]
            [1 3 4 2],[2 3 4 1],[3 1 4 2],[4 3 1 2]
            [1 4 2 3],[2 4 1 3],[3 4 2 1],[4 1 2 3]
            [1 4 3 2],[2 4 3 1],[3 4 1 2],[4 1 3 2]};

% Set parameters
InitializePsychSound;
disp(struct2table(PsychPortAudio('GetDevices')));
deviceID = input('Choose correct audio device and input its DeviceIndex: ');
sampRate = input('Input its DefaultSampleRate: ');

typeSeq  = seqTypes{seqID};              % sequence of schedule conditions, should be balanced across subjects  
triNum   = 90;                     % trial number of each schedule condition, should be an integer multiple of length(tSOA)
pretNum  = 90;                     % trial number of threshold stage
SOA      = 0.5;                    % standard SOA in Periodic Predictable condition
rdSOA    = [0.4 0.9];              % the range of the random SOA in AP3 and AU condition
stiD     = 0.05;                   % duration of each beep
ramp     = 0.004;                  % Fade in and fade out
ITIs     = [0.8 1.4];              % inter trial interval range (randomly selected in each trial)
maxRT    = 2;                      % skip to next trial in 2s (facilitating post-target EEG analysis)
tSOA     = [SOA/2 SOA 3*SOA/2];    % the interval between last cue and the target (based on tSOAp distribution) 
tSOAp    = [1/3   1/3   1/3  ];    % the probability distribution of each tSOA condition
cFreq    = 1750;                   % pitch of cue
tFreq    = [2000 1500];            % pitch of target
keys = {'UpArrow', 'DownArrow'};   % response keyName
keyCodes = KbName(keys);           % response keyCodes
noiseAmp = 0.1;                    % 90% amplitude is remained for titrating target loudness
cueAmp   = 0.4;                    % fixed amplitude which should be clear enough for all person
if gstgAmp == 0
    gstgAmp  = 0.08;                   % guessed amplitude of target (as the start point of staircase)
end
ampStep  = gstgAmp/10;             % staircase step

% stimulus schedule 
% * first value indicates the SOA between 1st and 2nd cue
% * if there are 7 interval value, the cue will show for 8 times
% * if a certain interval value is NaN, this interval will be randomly
%   drawn from a uniform distribution defined by rdSOA
% 1. periodic predictable
PP  = [1 1 1 1 1 1 1].*SOA;
% 2. aperiodic predictable (accelerating)
AP1 = [1/0.65 1/0.70 1/0.75 1/0.8 1/0.85 1/0.90 1/0.95];
% 3. aperiodic predictable (random pairs, 2N-1,)
AP2 = [1 nan 1 nan 1 nan 1 nan 1 nan 1 nan 1 nan];
% 4. aperiodic unpredictable (random)
AU  = [nan nan nan nan nan nan nan];

cueTypes = {'PP','AP1','AP2','AU'};
cSOAs = {PP,AP1,AP2,AU};
% trial-results table
[~,block] = expRun.generateTrialList('ID',nan,'cueType', ...
    nan,'t0',nan,'ITI',nan,'cSOA',nan,'tSOA',tSOA,'tFreq', ...
    1:length(tFreq),'tgAmp',nan,'tgTime',nan,'RT',nan, ...
    'judge',0,'soaSeed',nan,'noiseSeed',nan);
% threshold stage
repTimes = pretNum / height(block);  
preblock = repmat(block, repTimes, 1);  
results = preblock(randperm(height(preblock)), :);  
results.cueType = repmat({'AU'},pretNum,1);
results.cSOA = repmat({AU},pretNum,1);
% formal task
repTimes = triNum / height(block);  
block = repmat(block, repTimes, 1);  
for i = typeSeq(1:4)
    newblock = block(randperm(height(block)), :);
    newblock.cueType = repmat({cueTypes{i}},triNum,1);
    newblock.cSOA = repmat({cSOAs{i}},triNum,1);
    results = [results; newblock];
end
results.ITI = rand(height(results),1)*diff(ITIs) + ITIs(1);
results.soaSeed = randi(2^32-1,height(results),1);
results.noiseSeed = randi(2^32-1,height(results),1);
results.ID = transpose(-pretNum+1:triNum*4);

try
%% staircase titrating task
pahandle = PsychPortAudio('Open', deviceID, 1, 3, sampRate, 2);
scr = max(Screen('Screens'));
[wpnt,winRect] = Screen('OpenWindow',scr,127);

% titrate white noise volume
PsychPortAudio('Volume',pahandle,0.5);
while 1
    WN = noiseAmp.*(2.*rand(1,sampRate)-1);
    PsychPortAudio('FillBuffer', pahandle, [WN; WN]);
    PsychPortAudio('Start', pahandle, 1, 0, 1);
    inp = input('test the next volume(0~1), or input 0 to use current volume:  ');
    PsychPortAudio('Stop',pahandle);
    if inp == 0
        break
    else
        PsychPortAudio('Volume',pahandle,inp);
    end
end

% play example stimulus, then check
streams{1} = genStream(0,0,cFreq,rand(1),tFreq(1),maxRT,stiD,sampRate,noiseAmp,cueAmp,0.25,ramp);
streams{2} = genStream(0,0,cFreq,rand(1),tFreq(2),maxRT,stiD,sampRate,noiseAmp,cueAmp,0.25,ramp);
while 1
    inp = input('Which target to play? (1/2, 0 to skip)');
    if inp > 0 
        PsychPortAudio('FillBuffer', pahandle, [streams{inp}; streams{inp}]);
        PsychPortAudio('Start', pahandle, 1, 0, 1);
        PsychPortAudio('Stop', pahandle,1);
    else
        break
    end
end
% generate embedded stream
% quest or 'one up three down' staircase to get mixture ratio
tgAmp = gstgAmp;
corrCount = 0;  % counting correct times for staircase procedure
for i = 1:pretNum
    cSOA = results.cSOA{i};
    nanIndices = isnan(cSOA);  % find nan
    if any(nanIndices)
        numNaNs = sum(nanIndices);  % count nan
        originalRNG = rng;  % Preserve global random number generator state 
        rng(results.soaSeed(i));  
        cSOA(nanIndices) = rand(1, numNaNs) * diff(rdSOA) + rdSOA(1);  % random SOA in rdSOA range
        rng(originalRNG);  % Maintain global random number consistency
    end
    results.cSOA{i}=cSOA;
    tSOA = results.tSOA(i);
    noiseSeed = results.noiseSeed(i);
    ITI = results.ITI(i);
    stream = genStream(ITI, cSOA, cFreq, tSOA, tFreq(results.tFreq(i)), maxRT, stiD, sampRate, noiseAmp, cueAmp, tgAmp, ramp, noiseSeed);
    PsychPortAudio('FillBuffer', pahandle, [stream; stream]);
    t0 = PsychPortAudio('Start', pahandle, 1, 0, 1);  % no-repeat, start rightnow, get the actual timing
    tt = GetSecs;
    results.t0(i) = t0;
    tgTime = t0+ITI+sum(cSOA)+tSOA;
    results.tgTime(i) = tgTime;
    WaitSecs('UntilTime',tgTime);
    timeout = tgTime + maxRT;  % 
    RT = nan;
    while GetSecs < timeout
        [keyIsDown, keyT, keyCode] = KbCheck;
        if sum(keyCode) == 1    % make sure only one key is pressed
            if keyCode(KbName(keys{results.tFreq(i)})) % correct key is pressed
                results.judge(i) = 1;
            elseif keyCode(KbName(keys{3-results.tFreq(i)})) % wrong key is pressed
                results.judge(i) = 0;
            end
            if any(keyCode(keyCodes))
                RT = keyT-tgTime;
                results.RT(i) = RT;
                break
            end
        end
        WaitSecs(0.005);
    end
    [startTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop', pahandle,1);
%     GetSecs-tt-(ITI+sum(cSOA)+tSOA+maxRT)
    % staircase: one-up three-down
    if results.judge(i)==0
        tgAmp = tgAmp+ampStep;
        corrCount = 0;
    else
        corrCount = corrCount + 1;
        if corrCount==3
            tgAmp = tgAmp-ampStep;
            corrCount = 0;
        end
    end     
    results.tgAmp(i)=tgAmp;
    fprintf('#%.0f  %.4fs, "%s", judge-%.0f, temporalErr-%.4fs\n',results.ID(i),RT,KbName(keyCode),results.judge(i),estStopTime-timeout)
end
figure;
plot(results.tgAmp(results.ID<0));
% update table
SubjInfo = readtable('./Data/SubjInfo.csv');
rowIdx = find(SubjInfo.subjID == subjID & SubjInfo.groupID == groupID,1);
SubjInfo(rowIdx,'threshold') = {tgAmp};
writetable(SubjInfo,'./Data/SubjInfo.csv');

input('Continue to Formal Task?')


%% main experiment
for i = pretNum + 1:4*triNum
    if mod(i-pretNum, triNum) == 1 
        input('Continue to Next Block?') % free rest time
    end
    cSOA = results.cSOA{i};
    nanIndices = isnan(cSOA);  % find nan
    if any(nanIndices)
        numNaNs = sum(nanIndices);  % count nan
        originalRNG = rng;  % Preserve global random number generator state 
        rng(results.soaSeed(i));  
        cSOA(nanIndices) = rand(1, numNaNs) * diff(rdSOA) + rdSOA(1);  % random SOA in rdSOA range
        rng(originalRNG);  % Maintain global random number consistency
    end
    results.cSOA{i}=cSOA;
    tSOA = results.tSOA(i);
    noiseSeed = results.noiseSeed(i);
    ITI = results.ITI(i);
    stream = genStream(ITI, cSOA, cFreq, tSOA, tFreq(results.tFreq(i)), maxRT, stiD, sampRate, noiseAmp, cueAmp, tgAmp, ramp, noiseSeed);
    PsychPortAudio('FillBuffer', pahandle, [stream; stream]);
    tt = GetSecs;
    t0 = PsychPortAudio('Start', pahandle, 1, 0, 1);  % no-repeat, start rightnow, get the actual timing
    results.t0(i) = t0;
    tgTime = t0+ITI+sum(cSOA)+tSOA;
    results.tgTime(i) = tgTime;
    WaitSecs('UntilTime',tgTime);
    timeout = tgTime + maxRT;  % 
    RT = nan;
    while GetSecs < timeout
        [keyIsDown, keyT, keyCode] = KbCheck;
        if sum(keyCode) == 1    % make sure only one key is pressed
            if keyCode(KbName(keys{results.tFreq(i)})) % correct key is pressed
                results.judge(i) = 1;
            elseif keyCode(KbName(keys{3-results.tFreq(i)})) % wrong key is pressed
                results.judge(i) = 0;
            end
            if any(keyCode(keyCodes))
                RT = keyT-tgTime;
                results.RT(i) = RT;
                break
            end
        end
        WaitSecs(0.005);
    end
    [startTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop', pahandle,1);

    results.tgAmp(i)=tgAmp;
    fprintf('#%.0f  %.4fs, "%s", judge-%.0f, temporalErr-%.4fs\n',results.ID(i),RT,KbName(keyCode),results.judge(i),estStopTime-timeout)
end

%%
PsychPortAudio('Close',pahandle);
save(sprintf('./Data/Result_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr),"results")
if saveRaw
    save(sprintf('./Data/EXP_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
end
catch me
    save(sprintf('./Data/Interrupted/EXPINT_Sub%.0f_Ses%.0f_%s_%s_%s',subjID, session, location, subjName, DTstr))
    PsychPortAudio('Close');
end

