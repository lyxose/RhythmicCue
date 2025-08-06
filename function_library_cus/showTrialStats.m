function oper = showTrialStats(wptr, i, showLastN, results, skipFirstNum, passKey, errKey, instFolder, instName)
    % Show recent trial performance statistics 
    % 
    % This function displays the recent N trials' accuracy, overall accuracy, and the key sequence
    % 
    % Parameters:
    % wptr        : Window pointer (returned by Screen('OpenWindow'))
    % i           : Current trial index
    % showLastN   : Number of recent trials to show detail and statistics for
    % results     : Table containing trial results with fields 'judge' and 'Key'
    % skipFirstNum: Number of pretrials before the formal experiment, which is used to skip the first few trials
    % passKey     : Key name for oper = 1 (e.g., 'UpArrow')
    % errKey      : Key name for oper = -1 (e.g., 'DownArrow')
    % instFolder  : Folder path containing the instruction image files
    % instName    : Partial name to match instruction image files

    % 1. draw background image
    % Read the list of image file names from the specified folder
    afiles = dir(fullfile(instFolder, 'EmptyMatchedStruct'));
    for ftype={'*.png', '*.jpg', '*.jpeg'}
        files = dir(fullfile(instFolder, ftype{1})); % Assuming JPEG images
        afiles = [afiles; files];
    end
    fileNames = {afiles.name};
    
    % Match the file path based on the given name
    matchedFiles = fileNames(contains(fileNames, instName));
    
    % Check if there is at least one match
    if isempty(matchedFiles)
        error('No file matched found for the instruction image.');
    end
    
    % Get the first matched image file path and draw it
    imgPointer = 1;
    imgPath = fullfile(instFolder, matchedFiles{imgPointer});

    % 2. 计算正确率
    idxLastN = (i-showLastN):(i-1);
    idxLastN = idxLastN(idxLastN > skipFirstNum); 
    if isempty(idxLastN)
        accLastN = NaN;
    else
        accLastN = mean(results.judge(idxLastN),'omitnan');
    end
    idxAll = (skipFirstNum+1):(i-1);
    if isempty(idxAll)
        accAll = NaN;
    else
        accAll = mean(results.judge(idxAll),'omitnan');
    end

    % 3. 获取最近LastN次按键序列
    keySeq = '';
    for k = idxLastN
        thisKey = results.Key{k};
        if length(thisKey) >= 1
            thisKey = thisKey(1); % 只取第一个字符
        else
            thisKey = ' '; % 如果没有按键，使用空格表示
        end
        keySeq = [keySeq thisKey];
    end

    % 4. 绘制文本
    infoStr = sprintf('ACC of Last %.0f: %.2f\n    total ACC: %.2f\n    Act: %s', showLastN, accLastN, accAll, keySeq);
    [tex,~] = drawCentImg(wptr, imgPath, 'fit');
    Screen('DrawText', wptr, infoStr, 30, 30, [0 0 0]);
%     Screen('DrawTexture',wptr, tex,[],centRect);
    Screen('Flip', wptr);
    % 5. 检查按键
    oper = 0;
    while oper == 0
        [keyIsDown, ~, keyCode] = KbCheck;
        if keyIsDown
            if ismember(passKey, KbName(keyCode))
                oper = 1;
            elseif ismember(errKey, KbName(keyCode))
                oper = -1;
            end
        end
        checkend;
        WaitSecs(0.02);
    end    
    Screen('Close', tex);
end
