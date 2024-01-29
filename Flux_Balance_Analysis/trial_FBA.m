%Load the file containing context specific GEM models
addpath(genpath('/home/shweta/Documents/MATLAB/RAVEN-v2.8.0'))

load("All_Gem_Models.mat")


for i = 1:150

    disp(['Model: ' num2str(i) ' of ' num2str(150)])

    %To add a reaction
    rxnsToAdd.rxns = {'ATP_Hydrolysis'};
    rxnsToAdd.equations = {'ATP[c] + H2O[c] => ADP[c] + H+[c] + Pi[c]'};
    All_GEM_Models{i}=addRxns(All_GEM_Models{i},rxnsToAdd,3)
    
    constructEquations(All_GEM_Models{i}, 'ATP_Hydrolysis')

    All_GEM_Models{i} = setParam(All_GEM_Models{i}, 'obj', 'ATP_Hydrolysis', 1);

    % Perform FBA
    solution = solveLP(All_GEM_Models{i});

    % Display the flux distribution
    disp(solution.x)
    sol=solution.x

    %numeric matrix to txt file with ID names
    FBAresult = strcat(All_GEM_Models{i}.id, '_FBA', '.txt');
    dlmwrite(FBAresult, sol, '\t');

    %Fetching the Reaction IDs in each model
    RxnID=All_GEM_Models{i}.rxns
    filename = strcat(All_GEM_Models{i}.id, '_RxnID', '.txt');
    writecell(RxnID,filename)
end