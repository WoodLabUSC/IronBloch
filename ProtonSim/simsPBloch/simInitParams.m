%% Eamon Doyle, CHLA/USC
%% Written:  6/5/2013
%% Set up parameters for simulation
%% simParams contains everything needed to define the physiological structure of what is being simulated
%% jobParams contains variables that set up file paths and other administrative things that are functionally required but vary with hardware that won't change the output of the simulation overall
%% In future revisions, this will also accept and parse a configuration file


function [jobParams simParams] = simInitParams(varargin)

p = inputParser;
p.addOptional('B0',3);
p.addOptional('jobParams',struct);
p.addOptional('simParams',struct);
p.parse(varargin{:});

jobParams = p.Results.jobParams;
simParams = p.Results.simParams;


[discard,hostname] = system('hostname');
clear discard;
if strfind(hostname,'hpc-login')
	jobParams.resultsDir = '/auto/rcf-proj/jw3/eamondoy/1Tfields';
	jobParams.FileDependenciesString = '/auto/rcf-proj/jw3/eamondoy/research_code/Iron/Cluster/ProtonSim/simsPBloch';
	jobParams.resultsSaveLoc = 'auto/rcf-proj/jw3/eamondoy/results';
elseif strfind(hostname,'terrapin')
	jobParams.resultsDir = '~/Documents/1Tfields';  %this is actually only where the mag fields are stored
	jobParams.FileDependenciesString = '/home/eamon/Documents/research_code/Iron/Cluster/ProtonSim/simsPBloch';
	jobParams.resultsSaveLoc = '~/Documents/';
else
	error('No file paths provided for this host');
end

jobParams.sched = getJobManagerInfo();


%%%%%%%%%%%%%%%%%%%%% Nominal parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Field will be generated and saved with these parameters
%% Alterations can be done during proton simulation since they are only
%% multiplying factors.

simParams.wetToDryWtRatio = 4.1;  %% WDR, from Zuyderhoudt et al.
%% previous simulations were performed at 3.5, unless specified.

%% To get the iron concentration, keeping approximately the same multiplying 
%% factor as the patients for similar R2 range, (FE*delX/volFrac).
simParams.delX = 9.74e-8;       %% 4:1 hemosiderin/ferritin mix


%% All fields will be generated for 1T, appropriate field multipliers will
%% be used during MRI sim below.
simParams.B0 = p.Results.B0;     

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~isfield(simParams,'numProtons')
    simParams.numProtons = 1000;
end


%%%%%%%%%%%%%%%%%%% Parameters that can be altered for geometry generation

simParams.cellBiasFlag = 4;   %% 0 for uniform,
%% 1 for inter-cellular gaussian bias,
%% 2 for inter-cellular patient specific anisotropy,
%% 3 for patient specific inter-cellular anisotropy
%% along with inter-particle distance measures.
%% 4, same as 3 but with sinusoids.
%% 5, same as 4 but here only sinusoids are filled, when spheres cannot be
%% accommodated, they are spilled into the hepatocyes. (simulating Type 4
%% hemochromatosis).

if simParams.cellBiasFlag == 0
	simParams.cellSigma = [];     %% specified but not used !
elseif simParams.cellBiasFlag == 1
	simParams.cellSigma = 10;     %% std dev of cell-iron pdf
elseif simParams.cellBiasFlag == 2
	simParams.cellSigma = [];     %% specified but not used !
elseif simParams.cellBiasFlag == 3
	simParams.cellSigma = [];     %% specified but not used !
elseif simParams.cellBiasFlag == 4
	simParams.cellSigma = [];     %% specified but not used !    
elseif simParams.cellBiasFlag == 5
	simParams.cellSigma = [];     %% specified but not used !        
end


%% Dont change for now...
simParams.NNfactor = 1;   %% nearest neighbor conversion factor, distance = NNfactor* (from NN histogram)
simParams.NNfactorS = 1;    %% sinusoids

%% Distribution Variability Type
simParams.ironDist.distributionVariabilityType = 1;
%% 3, 4, 5, 6 are available but are not going to be used.

%% Variability is std. dev or RMS error
switch simParams.ironDist.distributionVariabilityType

	case 1
		%% 1 : Mean distribution parameters
		simParams.ironDist.Size_Spread_FE_Variability = 0;
		simParams.ironDist.Anisotropy_Spread_FE_Variability = 0;
		simParams.ironDist.Anisotropy_Shape_FE_Variability = 0;
		simParams.ironDist.NN_Shape_FE_Variability = 0;

		
	case 2
		%% 2 : Random
		simParams.ironDist.Size_Spread_FE_Variability = 0.013;      %% stdev
		simParams.ironDist.Anisotropy_Spread_FE_Variability = 0.064;
		simParams.ironDist.Anisotropy_Shape_FE_Variability = 0.412;
		simParams.ironDist.NN_Shape_FE_Variability = 1.221;
		

	case 3
		%% 3 : Object Size Variability, Spread
		simParams.ironDist.Size_Spread_FE_Variability = 0.013;      %% multiples of stdev = 0.012896
		simParams.ironDist.Anisotropy_Spread_FE_Variability = 0;
		simParams.ironDist.Anisotropy_Shape_FE_Variability = 0;
		simParams.ironDist.NN_Shape_FE_Variability = 0;

simParams
	case 4
		%% 4 : Nearest Neighbor (NN) Variability, Shape
		simParams.ironDist.NN_Shape_FE_Variability = 1.221;
		simParams.ironDist.Size_Spread_FE_Variability = 0;
		simParams.ironDist.Anisotropy_Spread_FE_Variability = 0;
		simParams.ironDist.Anisotropy_Shape_FE_Variability = 0;


	case 5
		%% 5 : Anisotropy Variability, Spread
		simParams.ironDist.Anisotropy_Spread_FE_Variability = 0.064;
		simParams.ironDist.Size_Spread_FE_Variability = 0;
		simParams.ironDist.Anisotropy_Shape_FE_Variability = 0;
		simParams.ironDist.NN_Shape_FE_Variability = 0;


	case 6
		%% 6 : Anisotropy Variability, Shape
		simParams.ironDist.Anisotropy_Shape_FE_Variability = 0.412;
		simParams.ironDist.Size_Spread_FE_Variability = 0;
		simParams.ironDist.Anisotropy_Spread_FE_Variability = 0;
		simParams.ironDist.NN_Shape_FE_Variability = 0;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% patientID = {'rw';'gh';'he';'sc';'gz';'vm';...
%              'vl';'bb';'mf';'mb';'lq';'yp';...
%              'bs';'bj';'ts';'na';'bl';'cl';'sp';'es'};
% FE = [1.3 1.4 4.4 4.6 5.9 7.8 ...
%       12.7 13.4 14.8 16.2 16.6 19.2 ... 
%       23.1 25.5 29.0 29.6 30.0 32.9 35.4 57.8];  

simParams.sim_box_side = 80;     % um
simParams.hepatocyte_dim = 20;    % in um

simParams.fieldGridStep = 0.5;    % um for every grid point, 0.5 is nyquist for single sphere, however
simParams.spill = 2;      % um, to avoid NAN's during interpolation of boundary points

