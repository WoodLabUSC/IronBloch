function [patientInfo,patientDir] = PreparePatientInfo(resultsDir,cellBiasFlag,cellSigma,FE,distributionVariabilityType,Size_Spread_FE_Variability,Anisotropy_Spread_FE_Variability,Anisotropy_Shape_FE_Variability,NN_Shape_FE_Variability)
%% Preparing patientInfo.mat for given FE
%% lysosome distribution is computed using gamma function
%% volume fraction is estimated by linear relation with FE
%% 1 October 2007
%% Nilesh Ghugre, CHLA/USC
%% based on final lysosome size histograms
%% inter-cellular iron anisotropy
%% inter-particle distance histrogram

R = linspace(0.05,1.6,100);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Vol Frac Relationship is linear

volFrac = 0.0036813 + 0.001274*FE;

%% size distribution
%% gamma and beta have a linear trend with FE


mybeta = 0.0916947 + 0.000733*FE + Size_Spread_FE_Variability;   % scaling or spread parameter


mygamma = 2.7240626 + 1.9327553*mybeta; % shape parameter

mu   = 0.05; % initial location parameter, decided by min(R)
amp_scale = 1;    % amplitude scaling

% figure;plot(FE,mygamma,'go');
% figure;plot(FE,mybeta,'go');
% figure;plot(FE,volFrac,'go');

for i = 1:length(FE)
    
    x0 = [mygamma mybeta(i) mu amp_scale];
    GammaDist(i,:) = gammafit(x0,R);
    % GammaDist(i,:) = GammaDist(i,:)./max(GammaDist(i,:));
    GammaCDF(i,:) = cumsum( GammaDist(i,:)./sum(GammaDist(i,:)) );
    
    patientInfo(i).r = R;
    patientInfo(i).freq = GammaDist(i,:);
    patientInfo(i).dist_est_cumsum = GammaCDF(i,:);
    patientInfo(i).id = i;
    patientInfo(i).FE = FE(i);
    patientInfo(i).volFrac = volFrac(i)*100;
    
    % patientInfo(i).name = patientID{i};
        
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inter-cellular iron anisotropy
%% Taken from G:\Nilesh\Iron patients\LM-iron studies\cell
%% analysis\cellDriver2.m

% model: change with FE
shapeEst = 1.024473 + 0.0178712*FE + Anisotropy_Shape_FE_Variability;
scaleEst = exp(-3.270029 + 0.3448803*log(FE)) + Anisotropy_Spread_FE_Variability;
ampEst = 0.0754587 + 0.0039904*FE;

t1 = linspace(0,1,100);

for i = 1:length(FE)
    
    x_fit = [shapeEst(i) scaleEst(i) 0 ampEst(i)];
    Mfit = gammafit2(x_fit,t1);
    Mfit(isinf(Mfit)) = 0;
    
    patientInfo(i).normFE = t1;
    patientInfo(i).cellIronCDF = cumsum(Mfit)./max(cumsum(Mfit));
    patientInfo(i).cellIronPDF = Mfit;
    
    % figure;plot(patientInfo(i).normFE,patientInfo(i).cellIronCDF);
    
    clear Mfit x_fit;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inter-particle distance histograms
% model: change with FE
shapeEst = exp(-0.130851 + 0.4763741*log(FE)) + NN_Shape_FE_Variability;
scaleEst = exp(0.7404703 - 1.4030494*log(shapeEst));
ampEst = exp(0.9947519 - 0.2470084*log(FE));

t1 = linspace(0,15,100);

for i = 1:length(FE)
    
    x_fit = [shapeEst(i) scaleEst(i) 0 ampEst(i)];
    Mfit = gammafit2(x_fit,t1);
    Mfit(isinf(Mfit)) = 0;
    
    patientInfo(i).nnDist = t1;
    patientInfo(i).nnCDF = cumsum(Mfit)./max(cumsum(Mfit));
    patientInfo(i).nnPDF = Mfit;
    
    % figure;plot(patientInfo(i).nnDist,patientInfo(i).nnCDF);
    
    clear Mfit x_fit;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sinusoidal Iron Fraction.

for i = 1:length(FE)
    if FE(i) < 100
        % sinusoidIronFraction(i) = exp(4.3593308 - 0.3660853*log(FE(i))) ./ 100;
        
        %% using hepatocyte since hepatocyte iron can probably be scored
        %% more reliably than sinusoids.
        hepatocyteIronFraction(i) = exp(3.5669345 + 0.1857491*log(FE))./100;
        
        sinusoidIronFraction(i) = 1-hepatocyteIronFraction(i);
        
        if sinusoidIronFraction(i) > 1
            sinusoidIronFraction(i) = 1;
        end
    else
        sinusoidIronFraction(i) = 0.3;
    end
    patientInfo(i).sinusoidIronFraction = sinusoidIronFraction(i);
    
end


%%% Save initial parameters

for i = 1:length(FE)
   
        patientInfo(i).Anisotropy_Shape_FE_Variability = Anisotropy_Shape_FE_Variability;
        patientInfo(i).Size_Spread_FE_Variability = Size_Spread_FE_Variability;
        patientInfo(i).Anisotropy_Spread_FE_Variability = Anisotropy_Spread_FE_Variability;
        patientInfo(i).NN_Shape_FE_Variability = NN_Shape_FE_Variability;
    
end



%%% Preparing and creating appropriate directories to save data

%% 0 for uniform,
%% 1 for inter-cellular gaussian bias,
%% 2 for inter-cellular patient specific anisotropy,
%% 3 for patient specific inter-cellular anisotropy
%% along with inter-particle distance measures.
%% 4, same as 3 but with sinusoids.

if cellBiasFlag == 0
    patientDir = [resultsDir '/FE-' num2str(FE) '/UniformDistribution'];
elseif cellBiasFlag == 1
    patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution1/' num2str(cellSigma)];
elseif cellBiasFlag == 2
    patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution2'];
elseif cellBiasFlag == 3    
    switch distributionVariabilityType        
        case 1
            patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution3/Mean'];                  
        case 2
            patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution3/Size/Spread/' num2str(Size_Spread_FE_Variability)];
        case 3
            patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution3/NN/Shape/' num2str(NN_Shape_FE_Variability)];
        case 4  
            patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution3/Anisotropy/Spread/' num2str(Anisotropy_Spread_FE_Variability)];
        case 5    
            patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution3/Anisotropy/Shape/' num2str(Anisotropy_Shape_FE_Variability)];
    end    
elseif cellBiasFlag == 4
    switch distributionVariabilityType        
        case 1
            patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution4/Mean'];                  
        case 2
            patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution4/Size/Spread/' num2str(Size_Spread_FE_Variability)];
        case 3
            patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution4/NN/Shape/' num2str(NN_Shape_FE_Variability)];
        case 4  
            patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution4/Anisotropy/Spread/' num2str(Anisotropy_Spread_FE_Variability)];
        case 5    
            patientDir = [resultsDir '/FE-' num2str(FE) '/BiasedDistribution4/Anisotropy/Shape/' num2str(Anisotropy_Shape_FE_Variability)];
    end    
end

mkdir(patientDir);





