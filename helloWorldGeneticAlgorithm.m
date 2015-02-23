%% Standard Genetic Algorithm -- Solves For A User Input String
% Works Cited: "Efficiently Vectorized Code for Population Based
% Optimization Algorithms" by Oliver Rice & Rickard Nyman

% Author: Grant Williams (2015)
% Questions? Email me at: grant.williams.okstate.edu


clear;close all;clc;
tic

%% Select Target String
target  = 'Hello, world';
% *Can Be Any String With Any Values and Any Length!*

%% Parameters                    
popSize = 1000;                                 % Population Size
genome  = length(target);                       % Genome Size
mutRate = .01;                                  % Mutation Rate
S       = 4;                                    % Tournament Size
best    = Inf;                                  % Initialize Best (arbitrarily large)
MaxVal  = max(double(target));                  % Max Integer Value Needed
ideal   = double(target);                       % Convert Target to Integers
selection = 1;                                  % 0: Roulette
                                                % 1: Tournament
                                                % 2: 50% Truncation

crossover = 1;                                  % 0: Uniform crossover
                                                % 1: 1 point crossover
                                                % 2: 2 point crossover
%% Initialize Population
Pop = round(rand(popSize,genome)*(MaxVal-1)+1); % Creates Population With Corrected Genome Length

for Gen = 1:1e6                                 % A Very Large Number Was Chosen, But Shouldn't Be Needed
    
    %% Fitness
    
    % The fitness function starts by subtracting each element of each
    % member of the population from each element of the target string. The
    % function then takes the absolute value of the differences and sums
    % each row and stores the function as a mx1 matrix.
    
    F = sum(abs(bsxfun(@minus,Pop,ideal)),2);       
    
     
    
    % Finding Best Members for Score Keeping Reasons
    [current,currentGenome] = min(F);
    
    % Stores New Best Values and Prints New Best Scores
    if current < best
        best = current;
        bestGenome = Pop(currentGenome,:);
        
        fprintf('Gen: %d  |  Fitness: %d  |  ',Gen, best);
        disp(char(bestGenome));
    elseif best == 0
        break
    end

    %% Selection
    if selection == 0
    % Roulette
    [~,W] = min(ones(popSize,1)*(rand(1,2*popSize))>((cumsum(F)*ones(1,2*popSize)/sum(F))),[],1);
    
    elseif selection == 1
    % Tournament
    T = round(rand(2*popSize,S)*(popSize-1)+1);                     % Tournaments
    [~,idx] = min(F(T),[],2);                                       % Index to Determine Winners         
    W = T(sub2ind(size(T),(1:2*popSize)',idx));                     % Winners
    
    elseif selection == 2
    % 50% Truncation
    [~,V] = sort(F,'descend');                                      % Sort Fitness in Ascending Order
    V = V(popSize/2+1:end);                                         % Winner Pool
    W = V(round(rand(2*popSize,1)*(popSize/2-1)+1))';               % Winners    
    end
    
    %% Crossover
    
    if crossover == 0
    % Uniform Crossover
    idx = logical(round(rand(size(Pop))));                          % Index of Genome from Winner 2
    Pop2 = Pop(W(1:2:end),:);                                       % Set Pop2 = Pop Winners 1
    P2A = Pop(W(2:2:end),:);                                        % Assemble Pop2 Winners 2
    Pop2(idx) = P2A(idx);                                           % Combine Winners 1 and 2
    
    elseif crossover == 1
    % 1-Point Crossover
    Pop2 = Pop(W(1:2:end),:);                                       % New Population From Pop 1 Winners
    P2A = Pop(W(2:2:end),:);                                        % Assemble the New Population
    Ref = ones(popSize,1)*(1:genome);                               % The Reference Matrix
    idx = (round(rand(popSize,1)*(genome-1)+1)*ones(1,genome))>Ref; % Logical Indexing
    Pop2(idx) = P2A(idx);                                           % Recombine Both Parts of Winners
    
    elseif crossover == 2
    % 2-Point Crossover
    Pop2 = Pop(W(1:2:end),:);                                       % New Pop is Winners of old Pop
    P2A  = Pop(W(2:2:end),:);                                       % Assemble Pop2 Winners 2
    Ref  = ones(popSize,1)*(1:genome);                              % Ones Matrix
    CP   = sort(round(rand(popSize,2)*(genome-1)+1),2);             % Crossover Points
    idx = CP(:,1)*ones(1,genome)<Ref&CP(:,2)*ones(1,genome)>Ref;    % Index
    Pop2(idx)=P2A(idx);                                             % Recombine Winners
    end
    %% Mutation 
    idx = rand(size(Pop2))<mutRate;                                 % Index of Mutations
    Pop2(idx) = round(rand([1,sum(sum(idx))])*(MaxVal-1)+1);        % Mutated Value
    
    %% Reset Poplulations
    Pop = Pop2;
   
end

toc
