% make certain condition balanced across subjects by 3 constrain
% 1. not repeat in the same subject
% 2. not repeat in the group
% 3. repeat the other group to pair

function suggLabel = balanCond(usedCondLabels, pairCondLabels, labelRange)
    if isstring(labelRange)
        usedCounts = arrayfun(@(x) sum(strcmp(usedCondLabels,x)),labelRange);
        pairCounts = arrayfun(@(x) sum(strcmp(pairCondLabels,x)),labelRange);
    elseif isnumeric(labelRange)
        usedCounts = arrayfun(@(x) sum(usedCondLabels==x),labelRange);
        pairCounts = arrayfun(@(x) sum(pairCondLabels==x),labelRange);
    else
        error('Unsupport form of labelRange')
    end
    diffCounts = usedCounts-pairCounts;
    miniCounts = min(diffCounts);
    ALTs = labelRange(diffCounts==miniCounts);

    if ~isempty(ALTs)
        suggLabel = ALTs(randi(numel(ALTs)));
    else
        error('No aviliable candidate! check input of balanCond')
    end
end

