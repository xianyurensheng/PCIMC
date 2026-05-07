classdef PCIMC < ALGORITHM
% <multi> <real/binary/permutation><constrained>
% xxx algorithm

    methods
        function main(Algorithm, Problem)
            gen = 2;
            state = 0;
            change_threshold = 1e-3;

            Population = Problem.Initialization();
            numCon = size(Population(1, 1).con, 2);
            Fitness = CalFitness(Population.objs);
            Objvalues(1) = sum(sum(Population.objs, 1));

            ArchiveEdge = [];
            ArchiveCon = [];
            ArchivePopNum = max(floor(Problem.N / numCon), floor(sqrt(Problem.N)));
            for i = 1:numCon
                ConFitness = CalFitness(Population.objs, Population.cons, i);
                ConIndex = TournamentSelection(2, ArchivePopNum, ConFitness);
                ArchiveCon{1, i} = Population(ConIndex);
                ArchiveCon{2, i} = ConFitness(ConIndex);
            end
            
            ns = 0;
            finishCon = 0;
            InnerCon = [];
            EdgeCon = 1:numCon;

            epsilon0 = inf;

            alpha = 0.1;

            if ns == 0
                ns = ns + 1;
            end
            
            while Algorithm.NotTerminated(Population)

                Offspring = OperatorDE_pbest_1_main(Population, Problem.N, Problem, Fitness, alpha);

                if finishCon == 0
                    for i = EdgeCon
                        ArchMatingPool = TournamentSelection(2, 2*ArchivePopNum, ArchiveCon{2,i});
                        ArchiveOffspring = OperatorDE(Problem, ArchiveCon{1,i}, ArchiveCon{1,i}(ArchMatingPool(1:end/2)), ArchiveCon{1,i}(ArchMatingPool(end/2+1:end)));
                        OffspringCon = Offspring.cons;
                        feasibleOff = Offspring(OffspringCon(:, i) <= 0);
                        [ArchiveCon{1,i}, ArchiveCon{2,i}] = EnvironmentalSelection([ArchiveCon{1, i}, ArchiveOffspring, feasibleOff], ArchivePopNum, i, numCon);
                    end
                    [Population, Fitness] = EnvironmentalSelection([Population, Offspring], Problem.N, InnerCon, numCon);
                elseif finishCon >= numCon
                    ArchiveEdgeLength = length(ArchiveEdge);
                    P1 = randperm(ArchiveEdgeLength);
                    P2 = randperm(ArchiveEdgeLength);
                    ArchiveOffspring = OperatorDE(Problem, ArchiveEdge, ArchiveEdge(P1), ArchiveEdge(P2));
                    ArchiveEdge = UpdateArchive([ArchiveEdge, Offspring, ArchiveOffspring], Problem.N, EdgeCon, numCon);
                    cp        = (-log(epsilon0)-6)/log(1-0.5);
                    epsilon   = epsilon0*(1-Problem.FE/Problem.maxFE)^cp;
                    [Population, Fitness] = Improve_E_EnvironmentalSelection([Population, Offspring, ArchiveOffspring], Problem.N, epsilon);
                end


                stage = finishCon == numCon;
                if finishCon < numCon
                    Objvalues(gen) = sum(sum(abs(Population.objs), 1));
                    Archive = [];
                    for k = 1:numCon
                        Archive = [Archive; ArchiveCon{1,k}];
                    end
                    state = is_stable(Objvalues, gen, Population, Archive, Problem.N, ArchivePopNum, change_threshold, Problem.M, stage);
                elseif finishCon == numCon
                    Objvalues(gen) = sum(sum(abs(ArchiveEdge.objs), 1));
                    state = is_stable(Objvalues, gen, ArchiveEdge, Archive, Problem.N, ArchivePopNum, change_threshold, Problem.M, stage);
                end

                if state > 0 && finishCon == 0   
                    [ArchiveEdge, EdgeCon] = ConsClassification(Problem, ArchiveCon, Population, Problem.N, EdgeCon, ArchivePopNum, numCon);
                    InnerCon = 1:numCon;
                    finishCon = numCon;
                    cons = Population.cons;
                    cons(cons <= 0) = 0;
                    conss = sum(cons(:, EdgeCon),2);
                    epsilon0 = max(conss);
                    if epsilon0 == 0
                        epsilon0 = 1;
                    end
                    Objvalues(gen-1) = 0;
                elseif state > 0 && finishCon == numCon
                    EdgeCon = 1:numCon;
                    finishCon = numCon + 1;
                end
                
                if finishCon == 0 && Problem.FE >= Problem.maxFE * 0.7
                    cons  = Population.cons;
                    cons(cons <= 0) = 0;
                    conss = sum(cons,2);
                    epsilon0 = max(conss);
                    if epsilon0 == 0
                        epsilon0 = 1;
                    end
                    for i = numCon
                        ArchiveEdge = [ArchiveEdge ArchiveCon{1,i}];
                    end
                    [Population, Fitness] = EnvironmentalSelection([Population, ArchiveEdge], Problem.N, EdgeCon, numCon);
                    finishCon = numCon + 1;
                end

                gen = gen+1;
                if Problem.FE >=  Problem.maxFE
                    c = numCon+1;
                end


            end
        end
    end
end

function result = is_stable(Objvalues, gen, Population, Archive, N, ArchivePopNum, change_threshold, M, type)
    result = 0;
    [FrontNo, ~] = NDSort(Population.objs, N);
    NC = size(find(FrontNo==1), 2);
    if NC == N
        max_change = abs(Objvalues(gen)-Objvalues(gen-1));
        % change_threshold = change_threshold * abs(((Objvalues(gen) / N))/(M))*10^(M-2);
        change_threshold = change_threshold * abs((Objvalues(gen) )/( N* M))*10^(M-2);
        if max_change <= change_threshold
            result = result + 1;
            if type == 0
                archiveNum = size(Archive, 1);
                for k = 1:archiveNum
                    [FrontNo, ~] = NDSort(Archive(k,:).objs, N);
                    NC = size(find(FrontNo==1), 2);
                    if NC == ArchivePopNum
                        result = result + 1;
                    end
                end
                result = result - archiveNum;
            end
        end
    end

end