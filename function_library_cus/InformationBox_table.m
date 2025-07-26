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
function [groupID, subjID, session, location, subjName, subjGender, subjAge, threshold, seqType] = InformationBox(modal)
    if ~exist('./Data', 'dir')
        mkdir('./Data');
    end    
    infoFilePath = './Data/SubjInfo.csv';
    dateTime = datestr(datetime,'yyyymmddTHHMM');
    Tag =         {'Header',     'InputPrompt',                                 'DefaultAnswer',  'A',   'V'};

    defaultCell = {'groupID',    [],                                            1,                 0,     0 
                   'subjID',     [],                                            1,                 0,     0 
                   'location',   'Enter location, IP or PKU or ... (or Oth):',  'IP',              1,     1 
                   'subjName',   'Enter Subject Name (QuanPin):',               'MingZipinying',   1,     1 
                   'subjGender', 'Enter Subject Gender  [Man: M; Woman: F]:'    'M',               1,     1 
                   'subjAge',    'Enter Subject Age:',                          22,                1,     1 
                   'thresholdA', 'Enter tgAmp Threshold (unknown: 0):',         0,                 1,     1 
                   'thresholdV', 'Enter tgAmp Threshold (unknown: 0):',         0,                 1,     1 
                   'seqTypeA',   'Enter Block Seq. Type (1~24):',               [],                1,     0 
                   'seqTypeV',   'Enter Block Seq. Type (1~24):',               [],                0,     1 
                   'dateTA',     'Enter Date in "yyyymmddTHHMM" form:'          [],                1,     0 
                   'dateTV',     'Enter Date in "yyyymmddTHHMM" form:'          [],                0,     1 
                   'session1',   [],                                            [],                0,     0 
                   'session2',   [],                                            [],                0,     0 
                   'Note1',      [],                                            [],                0,     0 
                   'Note2',      [],                                            [],                0,     0  };     
    defaultTable  = cell2table(defaultCell,"VariableNames",Tag);    
    defaultTable.Properties.RowNames = defaultTable.Header;
    defaultTable.Header = [];
    header = transpose(defaultTable.Properties.RowNames);
    exampleInfo = transpose(defaultTable.DefaultAnswer);
    
    modalMask = transpose(defaultTable.(modal)==1);     % select variables will be set in specific modal
    seqType = sprintf('seqType%s',modal);               % to get the correct colomn for each modal
    oseqType = sprintf('seqType%s',setdiff('VA',modal));% get the other modal label to get the corresponding data
    dateT = sprintf('dateT%s',modal);                   % to get the correct colomn for each modal

    groupID = input('Enter Group ID (int):');
    


    %% Get default answer

    % first subject, creat table
    if ~exist(infoFilePath, 'file')
        subjID  = input('Enter Subject ID (int, max: 99):');
        SubjInfo = cell2table(exampleInfo, 'VariableNames', header,'Format',formator);
        SubjInfo.(seqType) = 1;     % MUST NOT set in defaultCell in prevent from miss-marking unused seqType
        SubjInfo.(dateT) = dateTime;
        rowIdx = 1;
%         defaultanswer = table2cell(SubjInfo(rowIdx, modalMask));

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
            if ~strcmp(SubjInfo.session1(rowIdx),modal)
                % choose a different seqType
                SubjInfo{rowIdx,seqType} = balanCond([SubjInfo{groupSubset,seqType};SubjInfo{rowIdx,oseqType}],SubjInfo{~groupSubset,seqType},1:24);
            end
            % IF re-Run in same condition, do not change setting
        % IF new one
        else
            rowIdx = height(SubjInfo)+1;
            SubjInfo = table2cell(SubjInfo);
            SubjInfo(rowIdx,:) = exampleInfo;
            % Select least frequent seqType (1-24) in groupSubset; 
            SubjInfo{rowIdx,seqType} = balanCond(SubjInfo{groupSubset,seqType},SubjInfo{~groupSubset,seqType},1:24);

        end
    end    
%         SubjInfo= table2cell(SubjInfo);
     

    
    %% Show information box
    prompt = defaultTable.InputPrompt(modalMask);
    name=sprintf('Info_%s_G%d_Subj%d', modal, subjID);
    numlines=1;
    defaultanswer = table2cell(SubjInfo(rowIdx, modalMask));
    defaultanswer = cellfun(@(x) num2str(x), defaultanswer, 'UniformOutput', false);
    answer   = inputdlg(prompt,name,numlines,defaultanswer);

    %% Write SubjInfo.csv table
    SubjInfo.groupID = groupID;
    SubjInfo.subjID = subjID;  
    % IF this is 1st session or its re-Run
    if isempty(SubjInfo.session1{rowIdx}) || SubjInfo.session1{rowIdx}==modal
        SubjInfo.session1{rowIdx}=modal;
    % This is the 2nd session
    else
        SubjInfo.session2{rowIdx}=modal;
    end
    SubjInfo= table2cell(SubjInfo);
    SubjInfo{rowIdx,modalMask} = transpose(answer);
    

%     session  = str2double(answer{1});  
%     location = answer{2};   % IP for Institute of Psychology, CAS for other institutd in CAS, Oth for other subjects
%     subjName  = answer{3};
%     subjGender= answer{4};
%     subjAge   = str2double(answer{5});
%     threshold = str2double(answer{6});  % measured contrast threshold 
%     seqType = str2double(answer{7});
%     THRdate = str2double(answer{8});
%     % update table
%     SubjInfo(rowIdx,:) = {groupID, subjID, session, {location}, {subjName}, {subjGender}, subjAge, threshold, seqType, THRdate, {''}};
    writecell([header;SubjInfo],infoFilePath);

return