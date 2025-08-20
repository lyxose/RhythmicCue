clear;
clc;
sca;
KbName('UnifyKeyNames'); 
Screen('Preference', 'SkipSyncTests', 1);
% PsychDebugWindowConfiguration;
% DisableKeysForKbCheck(KbName('f22'));
checkThreshStage = true; % whether to quit the PTB screen to check the threshold stage performance before the formal task

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

WaitSecs(0.2);
[keyIsDown, ~, keyCode] = KbCheck;
if sum(keyCode) ~= 0    % make sure only one key is pressed
    error('Check keyboard hardware!!! A certain key is pressed: %s',KbName(keyCode))
end

% get SubjInfo
[groupID, subjID, subjName, subjGender, subjAge, gstgAmp, seqID, formator] = InformationBox('A');
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
cSOAs.PP  = [1    1    1    1    1    1    1   ].*SOA;
tSOAs.PP = [1/2 2/2 3/2].*SOA;
% 2. aperiodic predictable (control the last beat)
cSOAs.AP1 = [2.05 1.90 1.75 1.60 1.45 1.30 1.15].*SOA;
tSOAs.AP1 = [1/2 2/2 3/2].*SOA;
% 3. aperiodic predictable (control average time)
cSOAs.AP2 = [1.45 1.30 1.15 1.00 0.85 0.70 0.55].*SOA;
tSOAs.AP2 = 0.4.*[1/2 2/2 3/2].*SOA;
% 4. aperiodic unpredictable (random)
cSOAs.AU  = [nan nan nan nan nan nan nan];
tSOAs.AU  = [1/2 2/2 3/2].*SOA;
% Repeat tSOA according to the proportion specified by tSOAp
% 1. Generate a tSOA sequence repeated according to the specified proportion
tSOAs.PP  = repelem(tSOAs.PP,  tSOAp);
tSOAs.AP1 = repelem(tSOAs.AP1, tSOAp);
tSOAs.AP2 = repelem(tSOAs.AP2, tSOAp);
tSOAs.AU  = repelem(tSOAs.AU,  tSOAp);

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
triNum   = 60;                     % trial number of each schedule condition, should be an integer multiple of length(tSOA)
catTriR  = 0;                      % catch trial rate of whole experiment
checkPer = 15;                     % check performance each "checkPer" trials
pretNum  = 90;                     % pretrial number of threshold stage
stiD     = 0.05;                   % duration of each beep
ramp     = 0.005;                  % Fade in and fade out
ITIs     = [0.8 1.4];              % inter trial interval range (randomly selected in each trial)
maxRT    = 2;                      % skip to next trial in 2s (facilitating post-target EEG analysis)
cFreq    = 1750;                   % pitch of cue
tFreq    = [2000 1500];            % pitch of target
keys = {'UpArrow', 'DownArrow'};   % response keyName
keyCodes = KbName(keys);           % response keyCodes
noiseAmp = 0.1;                    % 90% amplitude is remained for titrating target loudness
cueAmp   = 0.4;                    % fixed amplitude which should be clear enough for all person
fixSize  = 0.162*2;                % diameter of fixation dot, in degree 
if gstgAmp == 0
    gstgAmp  = 0.05;               % guessed amplitude of target (as the start point of staircase)
end
ampStep  = gstgAmp/10;             % staircase step of target amplitude
dynaStep = 0.8;                    % dynamicly decrease after each reverse (set to 1 to keep stepsize consistent)
textSize = 30;

% check parameters setting
if triNum < 1 || mod(triNum*(1-catTriR), sum(tSOAp))*length(tFreq)~=0 
    [~,catTriN] = rat(catTriR);
    error('triNum*(1-catTriR)/2 must be divisible by both sum(tSOAp) * length(tFreq)! try %.0f', 2*lcm(catTriN,sum(tSOAp)*length(tFreq)));
end
if pretNum < 1 || mod(pretNum*(1-catTriR), sum(tSOAp))*length(tFreq)*length(cueTypes) ~=0 ...
               || mod(pretNum*catTriR, length(tSOAs)) ~=0
    [~,catTriN] = rat(catTriR);
    error('pretNum*(1-catTriR) must be divisible by sum(tSOAp) * length(tFreq) * length(cueTypes), and pretNum*catTriR by length(cueTypes)! try %.0f', lcm(catTriN,sum(tSOAp)*length(tFreq)*length(cueTypes)));
end

%% trial-results table
preCatNum = pretNum * catTriR; % catch trial number in pretrial stage
preTruNum = pretNum-preCatNum;

[~,block] = expRun.generateTrialList('ID',nan,'cueType', ...
    nan,'t0',nan,'ITI',nan,'cSOA',{[]},'tSOA',1:length(tSOAs.AU),'tFreq', ...
    1:length(tFreq),'tgAmp',nan,'tgTime',nan,'RT',nan, 'Key',{''},...
    'judge',0,'soaSeed',nan,'noiseSeed',nan);
[~,catchs] = expRun.generateTrialList('ID',nan,'cueType', ...
    nan,'t0',nan,'ITI',nan,'cSOA',{[]},'tSOA',1:length(tSOAs.AU),'tFreq', ...
    0,'tgAmp',0,'tgTime',nan,'RT',nan, 'Key',{''},...
    'judge',nan,'soaSeed',nan,'noiseSeed',nan);
% threshold stage
repTimes = preTruNum / height(block);  
repCTimes = preCatNum / height(catchs);
preblock = [repmat(block, repTimes, 1); repmat(catchs, repCTimes, 1)];
results = preblock(randperm(height(preblock)), :);  
results.cueType = repmat({'AU'},pretNum,1);
results.cSOA = repmat({cSOAs.AU},pretNum,1);
for i = 1:pretNum
    results.tSOA(i) = tSOAs.AU(results.tSOA(i));
end

% formal task

[~,block] = expRun.generateTrialList('ID',nan,'cueType', ...
    nan,'t0',nan,'ITI',nan,'cSOA',{[]},'tSOA',tSOAs.AU,'tFreq', ...
    1:length(tFreq),'tgAmp',nan,'tgTime',nan,'RT',nan, 'Key',{''},...
    'judge',nan,'soaSeed',nan,'noiseSeed',nan);
[~,catchs] = expRun.generateTrialList('ID',nan,'cueType', ...
    nan,'t0',nan,'ITI',nan,'cSOA',{[]},'tSOA',1:length(tSOAs.AU),'tFreq', ...
    0,'tgAmp',0,'tgTime',nan,'RT',nan, 'Key',{''},...
    'judge',nan,'soaSeed',nan,'noiseSeed',nan);
triNum = triNum / 2; % split to 2 part 
CatNum = triNum * catTriR;
TruNum = triNum - CatNum;
repTimes = TruNum / height(block);  
block = repmat(block, repTimes, 1);   
% concetnate catch trials
catBlock = repmat(catchs, CatNum / height(catchs), 1);
catBlock.tFreq(:) = 0;
catBlock.tgAmp(:) = 0;
block = [block; catBlock];
for type = [cueTypes, cueTypes]
    newblock = block;
    newblock.cueType = repmat(type,triNum,1);
    newblock.cSOA = repmat({cSOAs.(type{1})},triNum,1);
    newblock.tSOA = repmat(transpose(tSOAs.(type{1})),triNum/length(tSOAs.(type{1})),1);
    results = [results; newblock(randperm(height(block)), :)];
end
triNum = triNum * 2;
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
%% instruction and threshold stage by quest
pahandle = PsychPortAudio('Open', deviceID, 1, 3, sampRate, 2);

% titrate white noise volume
PsychPortAudio('Volume',pahandle,0.004);% 0.003 for 604-5 ; 0.004 for 604-4 with TANGMAI earphone
while 1
    WN = noiseAmp.*(2.*rand(1,2.*sampRate)-1);
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
while 1
    inp = input('Which target to play? (Enter 1/2 to choose, 0 to skip): ');
    if inp > 0
        Tseq= cSOAs.AP2(ALTs(randi(size(ALTs,1)),:));
        ITI = min(ITIs)+rand*diff(ITIs);
        tSOA = tSOAs.AU(randi(length(tSOAs.AU)));
        stream = genStream(min(ITIs),ITI,Tseq,cFreq,tSOA,tFreq(inp),maxRT,stiD,sampRate,noiseAmp,cueAmp,mean([cueAmp,gstgAmp]),ramp);
        PsychPortAudio('FillBuffer', pahandle, [stream; stream]);
        t0 = PsychPortAudio('Start', pahandle, 1, 0, 1);
        tgTime = t0 + ITI + sum(Tseq)+tSOA;
        timeout = tgTime+maxRT;
        WaitSecs('UntilTime',tgTime);
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
        PsychPortAudio('Stop', pahandle,1);
    else
        break
    end
end
% generate embedded stream
% 'one up two down' staircase to get mixture ratio
corrCount = 0;  % counting correct times for staircase procedure
scr = max(Screen('Screens')); % 1; for 604-4
[w,winRect] = Screen('OpenWindow',scr,127);
scWidth = Screen('DisplaySize',scr)/10; % in cm
Screen('TextSize',w,textSize);
ut = UT(scWidth,winRect(3),headDist,false);
dotpRad = ut.deg2pix(fixSize)/2; % radius of fixation dot in pixcel
dotRect = [-dotpRad,-dotpRad,dotpRad,dotpRad]+[winRect(3),winRect(4),winRect(3),winRect(4)]./2;

lastChange = 0; 
changeIdx = 0;
HideCursor(scr);
if scr==0
    checkScreen = 1; % whether to show the performance screen 
else
    checkScreen = -1;
end


% Parameters for QUEST
pThreshold = 0.5;   % Target performance level for the threshold
beta = 3.5;         % Steepness of the Weibull function
delta = 0.01;       % Lapse rate (proportion of errors due to inattention)
gamma = 0.05;        % Guess rate (for 2AFC, it's 50%)
range = 5;          % Range of plausible log10(amplitude) values to test
grain = 0.01;       % Granularity of the tested values

% Initial guess for the threshold amplitude, converted to log10 scale
tGuess = log10(gstgAmp);
% Standard deviation of the initial guess
tGuessSd = 2;

q = QuestCreate(tGuess, tGuessSd, pThreshold, beta, delta, gamma, grain, range);
q.normalizePdf = 1; % Recommended for better performance
q = QuestUpdate(q, log10(0.001), 0);             % initialization with no prior data
q = QuestUpdate(q, log10(cueAmp), 1); % initialization with known over-threshold cue amplitude

for i = 1:pretNum
    if mod(i, checkPer) == 1 && i>1% each 10 trial rest 1s+
        if checkScreen == 1
            oper = showTrialStats(w, i, checkPer, results, 0, 'Return', 'BackSpace', instFolder, 'Check');
            checkScreen = checkScreen*oper; % if oper == -1, then change to rest screen
        else
            oper = showInstruc_Rest(w, 'Rest', instFolder, 'space', 'BackSpace', 1);
            checkScreen = checkScreen*oper; % if oper == -1, then change to check screen
        end
    end
    Screen('FillOval',w,0,dotRect);
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
    if results.tFreq(i)~=0 % not catch trial
        thistFreq = tFreq(results.tFreq(i));
        tgAmp = 10^QuestQuantile(q); % get the next amplitude to test
    else
        thistFreq = 0;
        tgAmp = 0;
    end
    stream = genStream(min(ITIs), ITI, cSOA, cFreq, tSOA, thistFreq, maxRT, stiD, sampRate, noiseAmp, cueAmp, tgAmp, ramp, noiseSeed);
    PsychPortAudio('FillBuffer', pahandle, [stream; stream]);
    t0 = PsychPortAudio('Start', pahandle, 1, 0, 1);  % no-repeat, start rightnow, get the actual timing
    Screen('Flip',w);

    % tt = GetSecs;
    results.t0(i) = t0;
    tgTime = t0+ITI+sum(cSOA)+tSOA;
    results.tgTime(i) = tgTime;
    Screen('Flip',w,t0+min(ITIs));   % central dot indicates that the trial is started, which disappeared after the min ITI
    WaitSecs('UntilTime',tgTime);
    timeout = tgTime + maxRT;  % 
    RT = nan;
    while GetSecs < timeout
        [keyIsDown, keyT, keyCode] = KbCheck;
        if any(keyCode(keyCodes)) && sum(keyCode) == 1 % one of specified keys is pressed
            if results.tFreq(i)~=0 % if not catch trial, check response
                if keyCode(KbName(keys{results.tFreq(i)})) % correct key is pressed
                    results.judge(i) = 1;
                elseif keyCode(KbName(keys{3-results.tFreq(i)})) % wrong key is pressed
                    results.judge(i) = 0;
                end
                RT = keyT-tgTime;
                results.RT(i) = RT;
                results.Key{i} = KbName(keyCode);
                break
            else % catch trial, just wait for timeout
                results.judge(i) = 0;
                RT = keyT-tgTime;
                results.RT(i) = RT;
                results.Key{i} = KbName(keyCode);
                break
            end
        end
        checkend;
        WaitSecs(0.005);
    end
    if ~(any(keyCode(keyCodes)) && sum(keyCode) == 1) % timeout without correct response
        results.Key{i} = '_NaN_';
        if results.tFreq(i)~=0 % if not catch trial, judge timeout as wrong
            results.judge(i) = 0;
        else % catch trial, judge timeout as correct
            results.judge(i) = 1; % no response is correct for catch trial
        end
    end
    [startTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop', pahandle,1);
    fprintf('Pre-%s  #%.0f  %.4fs, "%s", judge-%.0f, tgAmp-%.3f, temporalErr-%.4fs\n',results.cueType{i}, results.ID(i), RT, KbName(find(keyCode,1)), results.judge(i), tgAmp, estStopTime-timeout)
    results.tgAmp(i)=tgAmp;
    % quest update
    if results.tFreq(i)~=0 % if not catch trial, update quest
        q = QuestUpdate(q, log10(tgAmp), results.judge(i));
    end
    
end
%% threshold stage results
% visualization pretrials 
figure;
plot(results.tgAmp(results.tgAmp>0));
% update table
SubjInfo = readtable('./Data/SubjInfo.csv','Format',formator);
rowIdx = find(SubjInfo.subjID == subjID & SubjInfo.groupID == groupID,1);
SubjInfo(rowIdx,'thresholdA') = {tgAmp};
writetable(SubjInfo,'./Data/SubjInfo.csv');

if checkThreshStage
    sca;
    ShowCursor(scr);
    y ='y';
    n ='n';
    while true
        inp = input("Ask the subject's confidence and check the tgAmp curve. \nEnter 'n' to abort or enter 'y' to continue.");
        if strcmp(inp,'y')
            break
        elseif strcmp(inp,'n')
            error('Experiment aborted by user!')
        end
    end
    [w,winRect] = Screen('OpenWindow',scr,127.5);
    Screen('TextSize',w,textSize);
    HideCursor(scr);
end

tgAmp =10^  QuestQuantile(q); % get the final threshold amplitude
%% main experiment
if scr==0
    checkScreen = 1; % whether to show the performance screen 
else
    checkScreen = -1;
end
for i = pretNum + (1:4*triNum)
    if mod(i-pretNum, triNum) == 1 && (i-pretNum)>1% each block rest 10s+
        oper = showInstruc_Rest(w,'Rest',instFolder,'space','backspace',10);
        checkScreen = checkScreen*oper; % if oper == -1, then change check/rest screen setting
    elseif mod(i-pretNum, checkPer) == 1 && (i-pretNum)>1% each 15 trial rest 2s+
        if checkScreen == 1
            oper = showTrialStats(w, i, checkPer, results, pretNum, 'Return', 'BackSpace', instFolder, 'Check');
            checkScreen = checkScreen*oper; % if oper == -1, then change to rest screen
        else
            oper = showInstruc_Rest(w, 'Rest', instFolder, 'space', 'BackSpace', 1);
            checkScreen = checkScreen*oper; % if oper == -1, then change to check screen
        end
    end
    Screen('FillOval',w,0,dotRect);
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
    if results.tFreq(i)~=0
        thistFreq = tFreq(results.tFreq(i));
    else
        thistFreq = 0;
    end
    stream = genStream(min(ITIs), ITI, cSOA, cFreq, tSOA, thistFreq, maxRT, stiD, sampRate, noiseAmp, cueAmp, tgAmp, ramp, noiseSeed);
    PsychPortAudio('FillBuffer', pahandle, [stream; stream]);
    t0 = PsychPortAudio('Start', pahandle, 1, 0, 1);  % no-repeat, start rightnow, get the actual timing
    Screen('Flip',w);

    tgTime = t0+ITI+sum(cSOA)+tSOA;
    timeout = tgTime + maxRT;  % 
    results.t0(i) = t0;
    Screen('Flip',w,t0+min(ITIs));   % central dot indicates that the trial is started, which disappeared after the min ITI
    results.tgTime(i) = tgTime;
    WaitSecs('UntilTime',tgTime);
    RT = nan;
    while GetSecs < timeout
        [keyIsDown, keyT, keyCode] = KbCheck;
        if any(keyCode(keyCodes)) && sum(keyCode) == 1 % one of specified keys is pressed
            if results.tFreq(i)~=0 % if not catch trial, check response
                if keyCode(KbName(keys{results.tFreq(i)})) % correct key is pressed
                    results.judge(i) = 1;
                elseif keyCode(KbName(keys{3-results.tFreq(i)})) % wrong key is pressed
                    results.judge(i) = 0;
                end
                RT = keyT-tgTime;
                results.RT(i) = RT;
                results.Key{i} = KbName(keyCode);
                break
            else % catch trial, just wait for timeout
                results.judge(i) = 0;
                RT = keyT-tgTime;
                results.RT(i) = RT;
                results.Key{i} = KbName(keyCode);
                break
            end
        end
        checkend;
        WaitSecs(0.005);
    end
    if ~(any(keyCode(keyCodes)) && sum(keyCode) == 1) % timeout without correct response
        results.Key{i} = '_NaN_';
        if results.tFreq(i)~=0 % if not catch trial, judge timeout as wrong
            results.judge(i) = 0;
        else % catch trial, judge timeout as correct
            results.judge(i) = 1; % no response is correct for catch trial
        end
    end
    [startTime, endPositionSecs, xruns, estStopTime] = PsychPortAudio('Stop', pahandle,1);
    results.tgAmp(i)=tgAmp;
    fprintf('%s  #%.0f  %.4fs, "%s", judge-%.0f, temporalErr-%.4fs\n', results.cueType{i}, results.ID(i),RT,KbName(find(keyCode,1)),results.judge(i),estStopTime-timeout)
end

%%
sca;
PsychPortAudio('Close',pahandle);
writetable(results,sprintf('./Data/A_Result_G%.0f_Sub%.0f_%s_%s.csv',groupID, subjID, subjName, DTstr))
if saveRaw
    % save all variables
    save(sprintf('./Data/A_EXP_G%.0f_Sub%.0f_%s_%s',groupID, subjID, subjName, DTstr))
    % save matlab cmd output
    diary(sprintf('./Data/A_EXP_G%.0f_Sub%.0f_%s_%s.txt',groupID, subjID, subjName, DTstr));
end
catch me
    sca;
    PsychPortAudio('Close');
    ShowCursor(scr);
    disp(me);
    disp(me.stack);
    save(sprintf('./Data/Interrupted/A_EXPINT_G%.0f_Sub%.0f_%s_%s',groupID, subjID, subjName, DTstr))
end



