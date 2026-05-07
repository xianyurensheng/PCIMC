function Population = UpdateArchive(Population, N, curCon, numCon)
    
    if isempty(Population)
        return;
    end
    % [~,b]=unique(Population.objs, 'rows');
    % Population = Population(1, b);
    PopObj = Population.objs;
    PopCon = Population.cons;
    if length(curCon) < numCon
        PopCon = PopCon(:, curCon);
    end
    [FrontNo, ~] = NDSort([PopObj, sum(max(0,PopCon),2)], inf);
    FeasibleIndex = sum(max(0,PopCon), 2) == 0;
    TempFrontNo = FrontNo(FeasibleIndex);
    Feasible = Population(FeasibleIndex);
    if ~isempty(Feasible)
        Population = Feasible(TempFrontNo == 1);
        FrontNo = TempFrontNo(TempFrontNo == 1);
    end
    Population = FeasibleUpdate(Population, FrontNo, N); 
    

end

function Population = FeasibleUpdate(Population, FrontNo, N)
    if length(Population) > N
        FrontNoSort = sort(FrontNo);
        MaxFNo = FrontNoSort(N);
        Next = FrontNo < MaxFNo;
        Population1 = Population(Next);
        Last = FrontNo == MaxFNo;
        
        Population2 = Population(Last);
        if sum(Last)~= N - sum(Next)
        Population2 = Truncation(Population2, Population, N-sum(Next));
        end
    
        Population = [Population1, Population2]; 
    end
end

function Population = Truncation(Population,PopAll,N)

    %% Truncation
    Zmin       = min(PopAll.objs,[],1);
    PopObjTemp = Population.objs;
    PopObj = (PopObjTemp -repmat(Zmin,length(Population),1))./(repmat(max(PopAll.objs),length(Population),1)-repmat(Zmin,length(Population),1));

    
    Distance = pdist2(PopObj,PopObj);
    Distance(logical(eye(length(Distance)))) = inf;
    Del = false(1,size(PopObj,1));
    while sum(Del) < size(PopObj,1)-N
        Remain   = find(~Del);
        Temp     = sort(Distance(Remain,Remain),2);
        [~,Rank] = sortrows(Temp);
        Del(Remain(Rank(1))) = true;
    end
    
      Population = Population(Del==0);
end
