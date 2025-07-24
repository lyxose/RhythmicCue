function [groupID, subjID, session, location, subjName, subjGender, subjAge, threshold, seqType] = InformationBox
    if ~exist('./Data', 'dir')
        mkdir('./Data');
    end    
    infoFilePath = './Data/SubjInfo.csv';
    dateNow = str2double(datestr(datetime,'yyyymmdd'));
    exampleInfo = {1, 1, 1, 'IP', 'MingZipinying', 'M', 22, 0, 1, dateNow,''};
    prompt={'Enter Exp. Session:',...
            'Enter location, IP or PKU or ... (or Oth):',...
            'Enter Subject Name (QuanPin):',...
            'Enter Subject Gender  [Man: M; Woman: F]:',...
            'Enter Subject Age:',...
            'Enter tgAmp Threshold (unknown: 0):',...
            'Enter Block Seq. Type (1~24):',...
            'Enter Date in "yyyymmdd" form:'};
    name='Experimental Information';
    numlines=1;
    
    % block sequence should be randomized across subjects
    % seqTypes = {[1 2 3 4],[2 1 3 4],[3 2 1 4],[4 2 3 1]
    %             [1 2 4 3],[2 1 4 3],[3 2 4 1],[4 2 1 3]
    %             [1 3 2 4],[2 3 1 4],[3 1 2 4],[4 3 2 1]
    %             [1 3 4 2],[2 3 4 1],[3 1 4 2],[4 3 1 2]
    %             [1 4 2 3],[2 4 1 3],[3 4 2 1],[4 1 2 3]
    %             [1 4 3 2],[2 4 3 1],[3 4 1 2],[4 1 3 2]};
    groupID = input('Enter Group ID (int):');
    % check existance
    if ~exist(infoFilePath, 'file')
        header = {'groupID', 'subjID', 'session', 'location', 'subjName', 'subjGender', 'subjAge', 'threshold', 'seqType', 'THRdate', 'Note'};
        SubjInfo = cell2table(exampleInfo, 'VariableNames', header);
        subjID  = input('Enter Subject ID (int, max: 99):');
        rowIdx = 1;
        defaultanswer = exampleInfo(3:end-1);
    else
        SubjInfo = readtable(infoFilePath);
        SubjInfo.Note = string(SubjInfo.Note);
        groupSubset = SubjInfo.groupID == groupID;
        if any(groupSubset)
            subjID = input(sprintf('Enter Subject ID (next: %.0f):', max(SubjInfo.subjID(groupSubset))+1)); %     
        else
            subjID  = input('Enter Subject ID (int, max: 99):');
        end
        rowIdx = find(SubjInfo.subjID == subjID & SubjInfo.groupID == groupID,1);
        if ~isempty(rowIdx)
            defaultanswer = table2cell(SubjInfo(rowIdx, 3:end)); 
            defaultanswer{1} = defaultanswer{1}+1; % default as the next session number to prevent data overwrite or comfusing
            defaultanswer{7} = mod(defaultanswer{7}+1,24);  % change a block sequence type
        else
            rowIdx = height(SubjInfo)+1;
            defaultanswer = exampleInfo(3:end-1);
            % Select least frequent seqType (1-24) in groupSubset; 
            usedTypes = histcounts(SubjInfo.seqType(groupSubset),1:25);
            typeCandi = find(usedTypes == min(usedTypes));  % use find(xx,1) can make the type be selected in numerical order
            rest_seqType = typeCandi(randi(numel(typeCandi))); % random choice if ties
            defaultanswer{7} = rest_seqType;
        end    
    end

    defaultanswer = cellfun(@(x) num2str(x), defaultanswer, 'UniformOutput', false);
    answer   = inputdlg(prompt,name,numlines,defaultanswer);
    session  = str2double(answer{1});  
    location = answer{2};   % IP for Institute of Psychology, CAS for other institutd in CAS, Oth for other subjects
    subjName  = answer{3};
    subjGender= answer{4};
    subjAge   = str2double(answer{5});
    threshold = str2double(answer{6});  % measured contrast threshold 
    seqType = str2double(answer{7});
    THRdate = str2double(answer{8});
    % update table
    SubjInfo(rowIdx,:) = {groupID, subjID, session, {location}, {subjName}, {subjGender}, subjAge, threshold, seqType, THRdate, {''}};
    writetable(SubjInfo,infoFilePath);
return