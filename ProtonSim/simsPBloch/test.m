disp('--------------------Running MRI Sim--------------------');
tic
sched = findResource('scheduler','type','jobmanager','name','MRISIM-JM','LookupURL','cluster1');

pjob = createParallelJob(sched);
set(pjob, 'MaximumNumberOfWorkers', 3);

%% for G5
% set(pjob,'FileDependencies',{'/Users/nilesh/Monte Carlo/DCT-Oct2006/G5 simsP patients final'});    %

%% for Zaphod
% set(pjob,'FileDependencies',{'G:\Nilesh\Monte Carlo\DCT-Oct2006\G5 sims\G5 simsP patients final'});

%% Generalized, FileDependenciesString defined in BatchScriptSim
set(pjob,'FileDependencies',{'testCall.m' 'testF.m'});

simTask = createTask(pjob, @testCall, 1, {});

submit(pjob);
waitForState(pjob);

outputP = getAllOutputArguments(pjob)
toc

destroy(pjob);
