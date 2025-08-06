% 需要兼容听觉和视觉任务：
% 固定2个session，自动记录session顺序：session1 session2 填入A/V，取消输入
% threshold在表格中也需要包含两列，可通过判断任务标记字符，选择采用的threshold列
% 不要修改note内容，增加一列note
% staircase不要手动固定，而是根据此前实验结果平均值进行估计
% 修改保存时引用的session为字符格式
% 需要生成两种seqType，且保证不相同

% block sequence should be randomized across subjects
% seqTypes = {[1 2 3 4],[2 1 3 4],[3 2 1 4],[4 2 3 1]
%             [1 2 4 3],[2 1 4 3],[3 2 4 1],[4 2 1 3]
%             [1 3 2 4],[2 3 1 4],[3 1 2 4],[4 3 2 1]
%             [1 3 4 2],[2 3 4 1],[3 1 4 2],[4 3 1 2]
%             [1 4 2 3],[2 4 1 3],[3 4 2 1],[4 1 2 3]
%             [1 4 3 2],[2 4 3 1],[3 4 1 2],[4 1 3 2]};
function [groupID, subjID, subjName, subjGender, subjAge, threshold, seqTypeID, formator] = InformationBox(modal)
    if ~exist('./Data', 'dir')
        mkdir('./Data');
    end    
    infoFilePath = './Data/SubjInfo.csv';
    threhKey = ['threshold',modal];
    dateTime = datestr(datetime,'yyyymmddTHHMM');
    % describe all var. in SubjInfo table                    
    Tag =         {'Header',     'InputPrompt',                                                     'DefaultAnswer',  'A',   'V',    'Type'};
    defaultCell = {'groupID',    [],                                                                1,                 0,     0,     '%d'    
                   'subjID',     [],                                                                1,                 0,     0,     '%d'    
                   'location',   'Enter location, e.g. IP | IBP | BFU | BNU | PKU ... (or Oth):',   'IP',              1,     1,     '%s'    
                   'subjName',   'Enter Subject Name (QuanPin):',                                   'MingZipinying',   1,     1,     '%s'    
                   'subjGender', 'Enter Subject Gender  [Man: M; Woman: F]:'                        'M',               1,     1,     '%s'    
                   'subjAge',    'Enter Subject Age:',                                              22,                1,     1,     '%d'    
                   'thresholdA', 'Enter tgAmp Threshold (unknown: 0):',                             0,                 1,     0,     '%f'    
                   'thresholdV', 'Enter tgAmp Threshold (unknown: 0):',                             0,                 0,     1,     '%f'    
                   'seqTypeA',   'Enter Block Seq. Type (1~24):',                                   0,                 1,     0,     '%d'    
                   'seqTypeV',   'Enter Block Seq. Type (1~24):',                                   0,                 0,     1,     '%d'    
                   'dateTA',     'Enter Date in "yyyymmddTHHMM" form:'                              [],                1,     0,     '%s'    
                   'dateTV',     'Enter Date in "yyyymmddTHHMM" form:'                              [],                0,     1,     '%s'    
                   'session1',   [],                                                                [],                0,     0,     '%s'    
                   'session2',   [],                                                                [],                0,     0,     '%s'    
                   'Note1',      [],                                                                [],                0,     0,     '%s'    
                   'Note2',      [],                                                                [],                0,     0,     '%s'};       
    defaultTable  = cell2table(defaultCell,"VariableNames",Tag);    
    defaultTable.Properties.RowNames = defaultTable.Header;
    defaultTable.Header = [];
    header = transpose(defaultTable.Properties.RowNames);
    formator = [defaultTable.Type{:}];
    exampleInfo = transpose(defaultTable.DefaultAnswer);
    
    modalMark = transpose(defaultTable.(modal)==1);     % select variables will be set in specific modal
    seqType = ['seqType',modal];               % to get the correct colomn for each modal
    oseqType = ['seqType',setdiff('VA',modal)];% get the other modal label to get the corresponding data
    dateT = sprintf('dateT%s',modal);                   % to get the correct colomn for each modal

    groupID = input('Enter Group ID (int):');

    newSubj = true;
    reRun = false;


    %% Get default answer
    
    % threshold set to 0 for unknown situation
    threshold = 0;

    % first subject, creat table
    if ~exist(infoFilePath, 'file')
        subjID  = input('Enter Subject ID (int, max: 99):');
        seqTypeID = 1;
        rowIdx = 1;

    % following subjects, match existing information
    else
        SubjInfo = readtable(infoFilePath,'Format',formator);
        % Get group members
        groupSubset = SubjInfo.groupID == groupID;
        if any(groupSubset)
            subjID = input(sprintf('Enter Subject ID (next: %.0f):', max(SubjInfo.subjID(groupSubset))+1)); %     
        else
            subjID  = input('Enter Subject ID (int, max: 99):');
        end
        % Find corresponding row
        rowIdx = find(SubjInfo.subjID == subjID & SubjInfo.groupID == groupID,1);
        % IF existing
        if ~isempty(rowIdx)
            exampleInfo = table2cell(SubjInfo(rowIdx,:));
            % IF second session
            if ~strcmp(SubjInfo.session1(rowIdx),modal)
                newSubj = false;
                % choose a different seqType
                seqTypeID = balanCond([SubjInfo{groupSubset,seqType};SubjInfo{rowIdx,oseqType}],SubjInfo{~groupSubset,seqType},1:24);
            else
                % IF re-Run in same condition, do not change any setting (except date)
                reRun = true;
                seqTypeID = exampleInfo{strcmp(header,seqType)};
            end
        % IF new one
        else
            rowIdx = height(SubjInfo)+1;
            seqTypeID = balanCond(SubjInfo{groupSubset,seqType},SubjInfo{~groupSubset,seqType},1:24);
        end
        % get average threshold as the guessed staricase start point
        priorThreIdx = groupSubset & SubjInfo.(threhKey)~=0;
        if any(priorThreIdx) && ~reRun
            threshold = mean(SubjInfo{priorThreIdx, threhKey});
            exampleInfo{strcmp(header,threhKey)} = threshold;
        end

    end    
    % Record date and time of start running
    exampleInfo{strcmp(header,dateT)}=dateTime;
    exampleInfo{strcmp(header,seqType)}=seqTypeID;
    
    %% Show information box
    prompt = defaultTable.InputPrompt(modalMark);
    name=sprintf('Info_%s_G%d_Subj%d', modal, groupID, subjID);
    numlines=1;
    defaultanswer = exampleInfo(modalMark);
    defaultanswer = cellfun(@(x) num2str(x), defaultanswer, 'UniformOutput', false);
    answer   = inputdlg(prompt,name,numlines,defaultanswer);

    %% Write SubjInfo.csv table
    exampleInfo{strcmp(header,'groupID')} = groupID;
    exampleInfo{strcmp(header,'subjID')} = subjID;
    % IF this is 1st session or its re-Run
    if newSubj
        exampleInfo{strcmp(header,'session1')} =modal;
    % This is the 2nd session
    else
        exampleInfo{strcmp(header,'session2')} =modal;
    end
    exampleInfo(modalMark) = transpose(answer);
    if ~exist(infoFilePath, 'file')
        SubjInfo = exampleInfo;
    else
        SubjInfo= table2cell(SubjInfo);
        SubjInfo(rowIdx,:) = exampleInfo;
    end
    
    % get vars for output
    subjName = exampleInfo{strcmp(header,'subjName')};
    subjGender = exampleInfo{strcmp(header,'subjGender')};
    subjAge = exampleInfo{strcmp(header,'subjAge')};
    threshold = str2double(exampleInfo{strcmp(header,threhKey)});

    % Rewrite the table file
    writecell([header;SubjInfo],infoFilePath);
return