function [ Offspring ] = OperatorDE_pbest_1_main(Population, popsize, Problem, Fitness, p)

    permutation = randperm(Problem.N);
    r0          = permutation;
    [r1,r2,~]   = gnR1R2R3(Problem.N, r0);
    
    array = permutation(1:popsize);
    
    [~, indBest] = sort(Fitness, 'ascend');
    pNP          = max(round(p * Problem.N), 2);    % choose at least two best solutions  
    randindex    = ceil(rand(1, popsize) * pNP);    % select from [1, 2, 3, ..., pNP]
    randindex    = max(1, randindex);               % to avoid the problem that rand = 0 and thus ceil(rand) = 0
    pbest        = Population(indBest(randindex));  % randomly choose one of the top 100p% solutions
    
    Offspring        = OperatorDE_pbest_1(Problem,Population(array).decs,pbest.decs,Population(r1(1:popsize)).decs,Population(r2(1:popsize)).decs);
    Offspring     = Problem.Evaluation(Offspring);
end