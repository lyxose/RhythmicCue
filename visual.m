
clear;
clc;
sca;
KbName('UnifyKeyNames'); 
% PsychDebugWindowConfiguration;
% DisableKeysForKbCheck(KbName('f22'));

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
    error('Check keyboard hardware!!! A certain key is pressed! %s',KbName(keyCode))
end

% get SubjInfo
[groupID, subjID, subjName, subjGender, subjAge, gstgAmp, seqID] = InformationBox('V');
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
triNum   = 60;                     % trial number of each schedule condition, should be an integer multiple of length(tSOA)
pretNum  = 180;                    % maximum pretrial number of threshold stage
revTimes = 10;                     % stop threshold stage when reaching this reverse times
stiD     = 0.1;                    % duration of each flash
ramp     = 0.004;                  % Fade in and fade out
ITIs     = [0.8 1.4];              % inter trial interval range (randomly selected in each trial)
maxRT    = 2;                      % skip to next trial in 2s (facilitating post-target EEG analysis)
cFreq    = 1750;                   % pitch of cue
tTilt    = [45 135];               % target orientation conditions
keys = {'UpArrow', 'DownArrow'};   % response keyName
keyCodes = KbName(keys);           % response keyCodes
noiseAmp = 0.2;                    % background contrast
cueAmp   = 0.3;                    % fixed contrast which should be clear enough for all person
GaborWidth = 1.2; % target width in degree
GaborSF = 2;        % cycle per degree
bgWidth = 15;       % background width, in degree
if gstgAmp == 0
    gstgAmp  = 0.008;               % guessed contrast of target (as the start point of staircase)
end
ampStep  = gstgAmp/10;             % staircase step of target amplitude
dynaStep = 0.8;                    % dynamicly decrease after each reverse (set to 1 to keep stepsize consistent)

% trial-results table
[~,block] = expRun.generateTrialList('ID',nan,'cueType', ...
    nan,'t0',nan,'ITI',nan,'cSOA',nan,'tSOA',tSOAs.AU,'tTilt', ...
    1:length(tTilt),'tgAmp',nan,'tgTime',nan,'RT',nan, ...
    'judge',0,'soaSeed',nan,'noiseSeed',nan);
% threshold stage
repTimes = pretNum / height(block);  
preblock = repmat(block, repTimes, 1);  
results = preblock(randperm(height(preblock)), :);  
results.cueType = repmat({'AU'},pretNum,1);
results.cSOA = repmat({cSOAs.AU},pretNum,1);
% formal task
triNum = triNum / 2; % split to 2 part 
repTimes = triNum / height(block);  
block = repmat(block, repTimes, 1);  
for type = [cueTypes,cueTypes]
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
%% staircase titrating task
pahandle = PsychPortAudio('Open', deviceID, 1, 3, sampRate, 2);

PsychPortAudio('Volume',pahandle,0.025);% 604-4 with TANGMAI earphone

scr = 1;%max(Screen('Screens'));
[wpnt,winRect] = Screen('OpenWindow',scr,127);
scWidth = Screen('DisplaySize',scr)/10; % in cm
ut = UT(scWidth,winRect(3),headDist,false);
dotpRad = ut.deg2pix(GaborWidth)/2; % radius of fixation dot in pixcel
dotRect = [-dotpRad,-dotpRad,dotpRad,dotpRad]+[winRect(3),winRect(4),winRect(3),winRect(4)]./2;
tgCenter = round(winRect(3:4)/2);

% generate stimuli
cueTexture = genStimTex(wpnt, winRect, ut, cueAmp, tgCenter, GaborSF, GaborWidth, 90);
tgTexture = zeros(1,length(tTilt));
for i = 1:length(tTilt)
    tgTexture(i) = genStimTex(wpnt, winRect, ut, cueAmp, tgCenter, GaborSF, GaborWidth, tTilt(i));
end

% play example stimulus, then check
while 1
    disp('Which target to play? (Press 1/2 to choose, 0 to skip): ');
    while 1
        [keyIsDown, keyT, keyCode] = KbCheck;
            if keyCode(KbName('1!')) ||  keyCode(KbName('1'))% key '1' is pressed
                inp = 1;
                break
            elseif keyCode(KbName('2@')) || keyCode(KbName('2')) %  key '2' is pressed
                inp = 2;
                break
            elseif keyCode(KbName('0)')) || keyCode(KbName('0')) %  key '0' is pressed
                inp = 0;
                break
            end
        checkend;
        WaitSecs(0.005);
    end
    if inp > 0
        beep = MakeBeep(cFreq, min(ITIs));
        PsychPortAudio('FillBuffer', pahandle, [beep;beep]);
        t0 = PsychPortAudio('Start', pahandle, 1, 0, 1);  % no-repeat, start rightnow, get the actual timing
        Screen('FillOval',wpnt,0,dotRect);
        Screen('Flip',wpnt);
        Screen('Flip',wpnt,t0+min(ITIs));   % central dot indicates that the trial is started, which disappeared after the min ITI

%         tempSeq = GetSecs + rand*diff(ITIs)+min(ITIs) + [0,cumsum(cSOAs.AP2(ALTs(randi(size(ALTs,1)),:)))];
        tempSeq = GetSecs + rand*diff(ITIs)+min(ITIs) + [0,cumsum(cSOAs.AP2)];
        fbTs = zeros(1,length(tempSeq));
        for i = 1:length(tempSeq)
            Screen('Drawtexture', wpnt, cueTexture)
            fbTs(i) = Screen('Flip', wpnt, tempSeq(i));
            Screen('Flip', wpnt, tempSeq(i)+stiD);
        end
        Screen('Drawtexture', wpnt, tgTexture(inp))
        tgTime = Screen('Flip', wpnt, fbTs(i) + tSOAs.AU(randi(length(tSOAs.AU))));
        Screen('Flip', wpnt, tgTime+stiD);

        disp('Temporal Error = ')
        disp(fbTs-tempSeq);
        
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
            checkend;
            WaitSecs(0.005);
        end
        if GetSecs > timeout
            disp('Time OUT!')
        end
        WaitSecs('UntilTime', timeout);
    else
        break
    end
end
% 'one up two down' staircase to get mixture ratio
tgAmp = gstgAmp;
corrCount = 0;  % counting correct times for staircase procedure

lastChange = 0; 
for i = 1:pretNum
    if mod(i, 10) == 1 && i>1% each 10 trial rest 1s+
        showInstruc_audioRest(wpnt,'Rest',instFolder,'space','backspace',1);
    end
    % get trial parameters    
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

    % close old texture then draw new textures to save computer memory 
    for j = 1:length(tTilt)
        Screen('Close',tgTexture(j))
        tgTexture(j) = genStimTex(wpnt, winRect, ut, tgAmp, tgCenter, GaborSF, GaborWidth, tTilt(j));
    end

    % Show fixation with pure beep
    beep = MakeBeep(cFreq, min(ITIs));
    PsychPortAudio('FillBuffer', pahandle, [beep;beep]);
    t0 = PsychPortAudio('Start', pahandle, 1, 0, 1);  % no-repeat, start rightnow, get the actual timing
    Screen('FillOval',wpnt,0,dotRect);
    Screen('Flip',wpnt);
    Screen('Flip',wpnt,t0+min(ITIs));   % central dot indicates that the trial is started, which disappeared after the min ITI

    % display the sequence of stimulus
    tempSeq = GetSecs + ITI + [0,cumsum(cSOA)];
    fbTs = zeros(1,length(tempSeq));

    % draw cue
    for j = 1:length(tempSeq)
        Screen('Drawtexture',wpnt, cueTexture)
        fbTs(j) = Screen('Flip', wpnt, tempSeq(j));
        Screen('Flip', wpnt, tempSeq(j)+stiD);
    end
    % draw texture
    Screen('Drawtexture', wpnt, tgTexture(results.tTilt(i)))
    tgTime = Screen('Flip', wpnt, fbTs(j) + tSOA);
    Screen('Flip', wpnt, tgTime+stiD);
    results.tgTime(i) = tgTime;
    timeout = tgTime + maxRT;  % 

    RT = nan;
    while GetSecs < timeout
        [keyIsDown, keyT, keyCode] = KbCheck;
        if sum(keyCode) == 1    % make sure only one key is pressed
            if keyCode(KbName(keys{results.tTilt(i)})) % correct key is pressed
                results.judge(i) = 1;
            elseif keyCode(KbName(keys{3-results.tTilt(i)})) % wrong key is pressed
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
    
    WaitSecs('UntilTime', timeout);
    
    % staircase: one-up two-down
    if results.judge(i)==0
        tgAmp = tgAmp+ampStep;
        corrCount = 0;
        if lastChange(end)~=-1
            lastChange = [lastChange,-1];
            ampStep = dynaStep*ampStep;
        end
    else
        corrCount = corrCount + 1;
        if corrCount==2
            tgAmp = tgAmp-ampStep;
            corrCount = 0;
            if lastChange(end)~=1
                lastChange = [lastChange,1];
                ampStep = dynaStep*ampStep;
            end
        end
    end     
    results.tgAmp(i)=tgAmp;
    fprintf('PreAU  #%.0f  %.4fs, "%s", judge-%.0f, tgAmp-%.4f, sumTemporalErr-%.4fs\n', results.ID(i), RT, KbName(keyCode), results.judge(i), tgAmp, sum(abs(tempSeq-fbTs)))
    % Check if reversed for enough times and keeped this amp once, then break 
    if length(lastChange) >= 2+revTimes && corrCount~=0 % 2+ because of the 0-head and the first step should not be considered
        break
    end
end

% visualization pretrials 
figure;
plot(results.tgAmp(1:i));

% IF reach maxium pretrial number
if i == pretNum
    sca;
    warning('Not enough reverse times(%d)!!', revTimes)
    input('Press Ctrl+C to abort or Enter to continue.')
    [wpnt,winRect] = Screen('OpenWindow',scr,127);
end
% update table
SubjInfo = readtable('./Data/SubjInfo.csv');
rowIdx = find(SubjInfo.subjID == subjID & SubjInfo.groupID == groupID,1);
SubjInfo(rowIdx,'thresholdA') = {tgAmp};
writetable(SubjInfo,'./Data/SubjInfo.csv');

% input('Continue to Formal Task?')
% generate target frames for all formal trials
tgTexture = zeros(1,length(tTilt));
for j = 1:length(tTilt)
    tgTexture(j) = genStimTex(wpnt, winRect, ut, tgAmp, tgCenter, GaborSF, GaborWidth, tTilt(j));
end

%% main experiment

for i = pretNum + (1:4*triNum)
    if mod(i, 10) == 1 && i>1% each 10 trial rest 1s+
        showInstruc_audioRest(wpnt,'Rest',instFolder,'space','backspace',1);
    end
    % get trial parameters    
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

    % Show fixation with pure beep
    beep = MakeBeep(cFreq, min(ITIs));
    PsychPortAudio('FillBuffer', pahandle, [beep;beep]);
    t0 = PsychPortAudio('Start', pahandle, 1, 0, 1);  % no-repeat, start rightnow, get the actual timing
    Screen('FillOval',wpnt,0,dotRect);
    Screen('Flip',wpnt);
    Screen('Flip',wpnt,t0+min(ITIs));   % central dot indicates that the trial is started, which disappeared after the min ITI
    
    % display the sequence of stimulus
    tempSeq = GetSecs + ITI + [0,cumsum(cSOA)];
    fbTs = zeros(1,length(tempSeq));

    % draw cue
    for j = 1:length(tempSeq)
        Screen('Drawtexture',wpnt,cueTexture)
        fbTs(j) = Screen('Flip', wpnt, tempSeq(j));
        Screen('Flip', wpnt, tempSeq(j)+stiD);
    end
    % draw target
    Screen('Drawtexture', wpnt, tgTexture(results.tTilt(i)))
    tgTime = Screen('Flip', wpnt, fbTs(j) + tSOA);
    Screen('Flip', wpnt, tgTime+stiD);
    results.tgTime(i) = tgTime;
    timeout = tgTime + maxRT;  % 

    RT = nan;
    while GetSecs < timeout
        [keyIsDown, keyT, keyCode] = KbCheck;
        if sum(keyCode) == 1    % make sure only one key is pressed
            if keyCode(KbName(keys{results.tTilt(i)})) % correct key is pressed
                results.judge(i) = 1;
            elseif keyCode(KbName(keys{3-results.tTilt(i)})) % wrong key is pressed
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
    WaitSecs('UntilTime', timeout);

%     disp('Temporal Error = ')
%     disp(fbTs-tempSeq);
    
    results.tgAmp(i)=tgAmp;
    fprintf('%s  #%.0f  %.4fs, "%s", judge-%.0f, tgAmp-%.4f, sumTemporalErr-%.4fs\n', results.cueType(i) ,results.ID(i), RT, KbName(keyCode), results.judge(i), tgAmp, sum(abs(tempSeq-fbTs)))
end

%%
sca;
PsychPortAudio('Close',pahandle);
writetable(results,sprintf('./Data/V_Result_G%.0f_Sub%.0f_%s_%s.csv',groupID, subjID, subjName, DTstr))
if saveRaw
    save(sprintf('./Data/V_EXP_G%.0f_Sub%.0f_%s_%s',groupID, subjID, subjName, DTstr))
end
catch me
    sca;
    disp(me);
    disp(me.stack);
    save(sprintf('./Data/Interrupted/V_EXPINT_G%.0f_Sub%.0f_%s_%s',groupID, subjID, subjName, DTstr))
    PsychPortAudio('Close');
end



