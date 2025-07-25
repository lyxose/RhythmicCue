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
instFolder = './Instructions';

% get SubjInfo
[groupID, subjID, subjName, subjGender, subjAge, gstgAmp, seqID] = InformationBox('A');
headDist = input('Distance between eyes and screen (cm):');

% block sequence should be randomized across subjects
seqTypes = {[1 2 3 4],[2 1 3 4],[3 2 1 4],[4 2 3 1]
            [1 2 4 3],[2 1 4 3],[3 2 4 1],[4 2 1 3]
            [1 3 2 4],[2 3 1 4],[3 1 2 4],[4 3 2 1]
            [1 3 4 2],[2 3 4 1],[3 1 4 2],[4 3 1 2]
            [1 4 2 3],[2 4 1 3],[3 4 2 1],[4 1 2 3]
            [1 4 3 2],[2 4 3 1],[3 4 1 2],[4 1 3 2]};
cueTypes = {'PP','AP1','AP2','AU'};

% stimulus schedule 
% * first value indicates the SOA between 1st and 2nd cue
% * if there are 7 interval value, the cue will show for 8 times
% * if a certain interval value is NaN, this interval will be randomly
%   drawn from a uniform distribution defined by rdSOA
SOA      = 0.5;                    % standard SOA in Periodic Predictable condition
rdSOA    = [0.4 0.9];              % the range of the random SOA in AP3 and AU condition
tSOAp    = [1   1   1];    % the probability distribution of each tSOA condition (still unused)
% 1. periodic predictable
cSOAs.PP  = [1 1 1 1 1 1 1].*SOA;
tSOAs.PP = [1/2 2/2 3/2].*SOA;
% 2. aperiodic predictable (control the last beat)
cSOAs.AP1 = [1.56 1.48 1.40 1.32 1.24 1.16 1.08].*SOA;
tSOAs.AP1 = [1/2 2/2 3/2].*SOA;
% 3. aperiodic predictable (control average time)
cSOAs.AP2 = [1.24 1.16 1.08 1.00 0.92 0.84 0.76].*SOA;
tSOAs.AP2 = 0.68.*[1/2 2/2 3/2].*SOA;
% 4. aperiodic unpredictable (random)
cSOAs.AU  = [nan nan nan nan nan nan nan];
tSOAs.AU  = [1/2 2/2 3/2].*SOA;

% Set parameters
InitializePsychSound;
DeviceTable = struct2table(PsychPortAudio('GetDevices'));
PriorList = DeviceTable(DeviceTable.HostAudioAPIId>2 & DeviceTable.NrOutputChannels>0,{'HostAudioAPIName','DeviceIndex','DeviceName','DefaultSampleRate'});
disp(PriorList);
if height(PriorList)==0
    warning('No good timing device can be found!!')
    disp(DeviceTable(:,{'HostAudioAPIName','DeviceIndex','DeviceName','DefaultSampleRate'}));
end
deviceID = input(['Choose correct audio device and input its DeviceIndex ' ...
                  '\n(priority: ASIO > WASAPI > WDM-KS > DS > MNE): ']);
sampRate = DeviceTable.DefaultSampleRate(DeviceTable.DeviceIndex==deviceID);
typeSeq  = seqTypes{seqID};        % sequence of schedule conditions, should be balanced across subjects  
triNum   = 90;                     % trial number of each schedule condition, should be an integer multiple of length(tSOA)
pretNum  = 60;                      % trial number of threshold stage
stiD     = 0.05;                   % duration of each beep
ramp     = 0.004;                  % Fade in and fade out
ITIs     = [0.8 1.4];              % inter trial interval range (randomly selected in each trial)
maxRT    = 2;                      % skip to next trial in 2s (facilitating post-target EEG analysis)
cFreq    = 1750;                   % pitch of cue
tFreq    = [2000 1500];            % pitch of target
keys = {'UpArrow', 'DownArrow'};   % response keyName
keyCodes = KbName(keys);           % response keyCodes
noiseAmp = 0.1;                    % 90% amplitude is remained for titrating target loudness
cueAmp   = 0.4;                    % fixed amplitude which should be clear enough for all person
fixSize  = 1.2;                    % diameter of fixation dot, in degree 
if gstgAmp == 0
    gstgAmp  = 0.08;                   % guessed amplitude of target (as the start point of staircase)
end
ampStep  = gstgAmp/10;             % staircase step

% trial-results table
[~,block] = expRun.generateTrialList('ID',nan,'cueType', ...
    nan,'t0',nan,'ITI',nan,'cSOA',nan,'tSOA',tSOAs.AU,'tFreq', ...
    1:length(tFreq),'tgAmp',nan,'tgTime',nan,'RT',nan, ...
    'judge',0,'soaSeed',nan,'noiseSeed',nan);
% threshold stage
repTimes = pretNum / height(block);  
preblock = repmat(block, repTimes, 1);  
results = preblock(randperm(height(preblock)), :);  
results.cueType = repmat({'AU'},pretNum,1);
results.cSOA = repmat({cSOAs.AU},pretNum,1);
% formal task
repTimes = triNum / height(block);  
block = repmat(block, repTimes, 1);  
for type = cueTypes
    newblock = block(randperm(height(block)), :);
    newblock.cueType = repmat(type,triNum,1);
    newblock.cSOA = repmat({cSOAs.(type{1})},triNum,1);
    newblock.tSOA = repmat(transpose(tSOAs.(type{1})),triNum/length(tSOAs.(type{1})),1);
    results = [results; newblock];
end
results.ITI = rand(height(results),1)*diff(ITIs) + ITIs(1);
results.soaSeed = randi(2^32-1,height(results),1);
results.noiseSeed = randi(2^32-1,height(results),1);
results.ID = transpose(-pretNum+1:triNum*4);
% get AU cSOA by permutate AP2
% permutation way selected from specific subset
% Exclude sequences with three consecutive increasing/decreasing values
% Exclude sequences with three interval in arithmetic progression
% Generate all permutations of 1-7
P = perms(1:7);
% Caculate the diff. metrix of P
first = diff(P,[],2);
% Mark local arithmetic sequences
arithm = ~diff(first,[],2);
% Consecutive increases/decreases
second = diff(sign(first),[],2);
% Three consecutive monotonic changes
% Note: Consecutive 2/-2 never occur, so consecutive zeros must be zero differences
third = diff(second,[],2);
% All eligible permutations (reject rows with any arithmetic sequence OR any three consecutive monotonic changes):
idx = ~(any(~third,2) | any(arithm,2));
ALTs = P(idx,:);
% For each trial, randomly select one permutation sequence
AUidx = strcmp(results.cueType,'AU');
auSeqs = ALTs(randi(size(ALTs,1),sum(AUidx),1),:);
for i = transpose(find(AUidx))
    results.cSOA(i) = {cSOAs.AP2(ALTs(randi(size(ALTs,1)),:))};
end

try
%% staircase titrating task
pahandle = PsychPortAudio('Open', deviceID, 1, 3, sampRate, 2);

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
streams{1} = genStream(0,cSOAs.AP2(ALTs(randi(size(ALTs,1)),:)),cFreq,rand(1),tFreq(1),maxRT,stiD,sampRate,noiseAmp,cueAmp,0.25,ramp);
streams{2} = genStream(0,cSOAs.AP2(ALTs(randi(size(ALTs,1)),:)),cFreq,rand(1),tFreq(2),maxRT,stiD,sampRate,noiseAmp,cueAmp,0.25,ramp);
while 1
    inp = input('Which target to play? (1/2, 0 to skip)');
    if inp > 0
        PsychPortAudio('FillBuffer', pahandle, [streams{inp}; streams{inp}]);
        PsychPortAudio('Start', pahandle, 1, 0, 1);
        PsychPortAudio('Stop', pahandle,1);
        timeout = GetSecs+2;
        while GetSecs < timeout
            [keyIsDown, keyT, keyCode] = KbCheck;
            if sum(keyCode) == 1    % make sure only one key is pressed
                if keyCode(KbName(keys{inp})) % correct key is pressed
                    disp('Correct!')
                    break
                elseif keyCode(KbName(keys{3-inp})) % wrong key is pressed
                    disp('Wrong!')
                    break
                end
            end
            WaitSecs(0.005);
        end
        if GetSecs > timeout
            disp('Time OUT!')
        end
    else
        break
    end
end
% generate embedded stream
% quest or 'one up three down' staircase to get mixture ratio
tgAmp = gstgAmp;
corrCount = 0;  % counting correct times for staircase procedure
scr = max(Screen('Screens'));
[w,winRect] = Screen('OpenWindow',scr,127);
scWidth = Screen('DisplaySize',scr)/10; % in cm
ut = UT(scWidth,winRect(3),headDist,false);
dotpRad = ut.deg2pix(fixSize)/2; % radius of fixation dot in pixcel
dotRect = [-dotpRad,-dotpRad,dotpRad,dotpRad]+[winRect(3),winRect(4),winRect(3),winRect(4)]./2;

[keyIsDown, ~, keyCode] = KbCheck;
if sum(keyCode) ~= 0    % make sure only one key is pressed
    error('Check keyboard hardware!!! A certain key is pressed! %s',KbName(keyCode))
end

for i = 1:pretNum
    Screen('FillOval',w,0,dotRect);
    Screen('Flip',w);
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
    % tt = GetSecs;
    results.t0(i) = t0;
    tgTime = t0+ITI+sum(cSOA)+tSOA;
    results.tgTime(i) = tgTime;
    Screen('Flip',w,t0+mean(ITIs));   % central dot indicates that the trial is started, which disappeared after the min ITI
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
        checkend;
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
    fprintf('#%.0f  %.4fs, "%s", judge-%.0f, tgAmp-%.3f, temporalErr-%.4fs\n',results.ID(i), RT, KbName(keyCode), results.judge(i), tgAmp, estStopTime-timeout)
end
figure;
plot(results.tgAmp(results.ID<0));
% update table
SubjInfo = readtable('./Data/SubjInfo.csv');
rowIdx = find(SubjInfo.subjID == subjID & SubjInfo.groupID == groupID,1);
SubjInfo(rowIdx,'thresholdA') = {tgAmp};
writetable(SubjInfo,'./Data/SubjInfo.csv');

% input('Continue to Formal Task?')


%% main experiment
for i = pretNum + (1:4*triNum)
    if mod(i-pretNum, triNum) == 1 % each block rest 10s+
        showInstruc_audioRest(w,'Rest',instFolder,'space','backspace',10);
    elseif mod(i-pretNum, 10) == 1 % each 10 trial rest 1s+
        showInstruc_audioRest(w,'Rest',instFolder,'space','backspace',1);
    end
    Screen('FillOval',w,0,dotRect);
    Screen('Flip',w);
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
    % tt = GetSecs;
    t0 = PsychPortAudio('Start', pahandle, 1, 0, 1);  % no-repeat, start rightnow, get the actual timing
    results.t0(i) = t0;
    tgTime = t0+ITI+sum(cSOA)+tSOA;
    Screen('Flip',w,t0+mean(ITIs));   % central dot indicates that the trial is started, which disappeared after the min ITI
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
        checkend;
        WaitSecs(0.005);
    end
    [startTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop', pahandle,1);

    results.tgAmp(i)=tgAmp;
    fprintf('#%.0f  %.4fs, "%s", judge-%.0f, temporalErr-%.4fs\n',results.ID(i),RT,KbName(keyCode),results.judge(i),estStopTime-timeout)
end

%%
sca;
PsychPortAudio('Close',pahandle);
save(sprintf('./Data/A_Result_G%.0f_Sub%.0f_%s_%s',groupID, subjID, subjName, DTstr),'results')
if saveRaw
    save(sprintf('./Data/A_EXP_G%.0f_Sub%.0f_%s_%s',groupID, subjID, subjName, DTstr))
end
catch me
    sca;
    save(sprintf('./Data/Interrupted/A_EXPINT_G%.0f_Sub%.0f_%s_%s',groupID, subjID, subjName, DTstr))
    PsychPortAudio('Close');
end



