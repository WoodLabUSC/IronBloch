function sched = getJobManagerInfo(varagin)

[~,hostname] = system('hostname');
hostname = cellstr(hostname);


if strfind(hostname{1},'hpc-login')
  sched = findResource('scheduler','type','torque');
  sched.ClusterMatlabRoot = '/usr/usc/matlab/default';
  sched.SubmitArguments = '-l walltime=23:59:00';
elseif strfind(hostname{1},'tinkertoy')
  sched = findResource('jobmanager','LookupURL','localhost');
else
  sched = findResource('scheduler','type','jobmanager','LookupURL','localhost');
end

