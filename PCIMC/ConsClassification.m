function [EdgePop, EdgeCons] = ConsClassification(Problem, ArchiveCon, Population, N, EdgeCon, ArchivePopNum, numCon)
    
    Archive = [];
    for i = EdgeCon
        Archive = [Archive, ArchiveCon{1, i}];
    end
    U = [Archive, Population];
    [FrontNo, ~] = NDSort(U.objs, inf);
    PopulationNo = FrontNo(length(Archive)+1:end);
    targetNo = FrontNo(1:end-N);
    targetNo = reshape(targetNo, ArchivePopNum, numCon)';
    Range = zeros(numCon, 2);
    refMin = min(PopulationNo);
    refMax = max(PopulationNo);
    for i = 1:numCon
        Range(i, 1) = min(targetNo(i, :));
        Range(i, 2) = max(targetNo(i, :));
    end
    Range = [Range; [refMin, refMax]];

    type4 = [];
    type3 = [];
    type_1_2 = [];
    for i = 1:numCon
        currentMin = Range(i, 1);
        currentMax = Range(i, 2);
        if currentMin > refMax && currentMax > refMax
            type4 = [type4; i];
        elseif (currentMin <= refMin) && (currentMax <= refMax)
            type_1_2 = [type_1_2; i];
        else
            type3 = [type3; i];
        end
    end
    
    if type4
        EdgeCons = pickCon(Range(type4, :), type4);
    elseif type3
        EdgeCons = pickCon(Range(type3, :), type3);
    else
        EdgeCons = 1:numCon;
    end

    EdgePop = [];
    for i = EdgeCons
        EdgePop = [EdgePop ArchiveCon{1, i}];
    end

end

function EdgeCons = pickCon(Range, index)
    EdgeCons = [];
    Max = max(Range(:, 2));
    MaxRange = Range(Range(:, 2) == Max, :);
    Min = min(MaxRange(:, 1));
    for i = 1:size(Range, 1)
        if Range(i, 1) >= Min
            EdgeCons = [EdgeCons, index(i)];
        end
    end
end