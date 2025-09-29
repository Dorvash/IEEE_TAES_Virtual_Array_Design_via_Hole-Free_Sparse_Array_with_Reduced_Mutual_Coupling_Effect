function [indices, objfunVal] = E_SMART(DOF,varargin)

tableResult     = true;

for i = 1:length(varargin)
    if strcmpi(varargin{i}, 'result') && i < length(varargin)
        if strcmpi(varargin{i+1}, 'off')
            tableResult = false;
            break
        end
    end
end

array           = zeros(1,DOF+1);

if DOF >= 5
    array(1) = 1;
    array(2) = 1;
    array(end) = 1;
else
    errordlg('Choose DOF greater or equal to 5','Error')
    return
end

indices                 = find(array == 1) - 1;
differences             = [0; sort(unique(abs(diff(nchoosek(indices, 2), 1, 2))))];
DOFvec                  = (0:DOF)';
objfunVal               = zeros(1e3,1);
cnt                     = 1;
objfunVal(cnt)          = sum(DOFvec);

while (numel(differences) ~= DOF+1) || any((differences) ~= DOFvec)

    cnt = cnt + 1;
    notCov = setdiff(DOFvec,differences);
    nextDiff = max(notCov);
    posCan = intersect([indices + nextDiff, indices - nextDiff], DOFvec);

    newDiffs = abs(posCan - indices);
    newDiffs_MC = newDiffs(all((newDiffs~=1)'),:);
    posCan_MC = posCan(all((newDiffs~=1)'),:);

    uniqueDiffs = arrayfun(@(i)fliplr(setdiff(newDiffs_MC(i,:), differences)), (1:length(posCan_MC)).', 'UniformOutput', false);
    uniqueDiffsLen = cellfun(@length, uniqueDiffs);
    maxLocs = find(max(uniqueDiffsLen) == uniqueDiffsLen);
    

    if length(maxLocs) >= 2
        diffsVal = sum(cell2mat([uniqueDiffs(maxLocs)]) .* (DOF-1).^(length(uniqueDiffs{maxLocs(1)})-1:-1:0),2);
        [~, idx] = max(diffsVal);
        maxLocs = maxLocs(idx);
    end
    
    differences = union(differences,uniqueDiffs{maxLocs});
    array(posCan_MC(maxLocs) + 1) = 1;
    indices = find(array == 1) - 1;

    objfunVal(cnt) = sum(notCov);

end

objfunVal(cnt+2:end)    = [];

if tableResult
    fprintf('================== Results ==================\n');
    fprintf('DOF = %d, Total Number of Sensors = %d\n', DOF, length(indices));
    fprintf('=============================================\n');
end

end