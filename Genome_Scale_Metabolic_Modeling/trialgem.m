%Loading Human-GEM into MATLAB
addpath(genpath('/home/shweta/Documents/MATLAB/Human-GEM/'))
addpath(genpath('/home/shweta/Documents/MATLAB/RAVEN-2.7.11'))
load('Human-GEM.mat');

%load the model as a structure named ihuman
ihuman

%prepare the reference model for use with ftINIT
prepDataHumanGEM = prepHumanModelForftINIT(ihuman, false, '/home/shweta/Documents/MATLAB/Human-GEM/data/metabolicTasks/metabolicTasks_Essential.txt', '/home/shweta/Documents/MATLAB/Human-GEM/model/reactions.tsv')
save('prepDataHumanGEM.mat', 'prepDataHumanGEM')

%Prepare the transcriptomic data
TPM_data = readtable('/home/shweta/Documents/MATLAB/TPM_NEW.txt');
[~, n] = size(TPM_data);
numSamp = n-1;
TPM_data

%Extract information from the table into a structure
data_struct.genes = TPM_data{:, 1};
data_struct.tissues = TPM_data.Properties.VariableNames(2:n);
data_struct.levels = TPM_data{:, 2:n};
data_struct.threshold = 1;
data_struct

%Run ftINIT for generating the context specific models
load('prepDataHumanGEM.mat');
All_GEM_Models = cell(numSamp, 1);
for i = 1:numSamp
    disp(['Model: ' num2str(i) ' of ' num2str(2)])
    All_GEM_Models{i} = ftINIT(prepDataHumanGEM, data_struct.tissues{i}, [], [], data_struct, {}, getHumanGEMINITSteps('1+0'), false, true);
    Name=data_struct.tissues{i};
    All_GEM_Models{i}.id = Name;
end
save('All_Gem_Models.mat', 'All_GEM_Models', '-v7.3')
load('All_Gem_Models.mat')