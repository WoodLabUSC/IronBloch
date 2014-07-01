%% Simulation on a cluster (parallel computation)
%% 23 Oct 2006
%% Nilesh Ghugre, CHLA/USC

%% random initialization of rand
rand('twister',sum(100*clock));

t = 0:step:interval;

%% logarithmic echo spacing
%% these are actually the times when 180 pulse is applied and not when echo
%% is formed, echo is formed at 2*TE.
TE = logspace(log(0.1)/log(10),log(30)/log(10),15);
% D = 3.125;  % micron^2/msec, this is the value for H20 in 1.5% agarose, taken from
% Ingolf Sack et al., Magnetic resonance elastography and
% diffusion-weighted iamging of the sol/gel phase transition in agarose,
% JMR 166 (2004) 252-261
% This has been used for comparison with synthetic compound behavior

%% Normal computation
% tic;[signal, signalSE] = simulate(sphereInfo,patientInfo,delBzGridP,patientIndx,fieldGridStep,spill,numProtons,step,interval,TE,D);toc
if 0
disp('--------------------Running MRI Sim--------------------');
tic
sched = findResource('scheduler','type','jobmanager','name','default_jobmanager','LookupURL','localhost');

pjob = createParallelJob(sched);
set(pjob, 'MaximumNumberOfWorkers', numWorkersSim);

%% for G5
% set(pjob,'FileDependencies',{'/Users/nilesh/Monte Carlo/DCT-Oct2006/G5 simsP patients final'});    %

%% for Zaphod
% set(pjob,'FileDependencies',{'G:\Nilesh\Monte Carlo\DCT-Oct2006\G5 sims\G5 simsP patients final'});

%% Generalized, FileDependenciesString defined in BatchScriptSim
set(pjob,'FileDependencies',{FileDependenciesString});

simTask = createTask(pjob, @simulateP, 2, {sphereInfo,sim_box_side,patientInfo,delBzGridP,patientIndx,fieldGridStep,spill,numProtons,step,interval,TE,D,hepatocyte_dim,cellBoundaryFlag,cellBiasFlag,sinusoidBoundaryFlag});
tic
submit(pjob);
waitForState(pjob,'finished');
toc
outputP = getAllOutputArguments(pjob);
toc

signal = outputP{1,1};
signalSE =  outputP{1,2};

destroy(pjob);
clear outputP simTask pjob sched;
end %this loop for testing.  destroy immediately so you don't forget, eamon

disp('running the new one');
[signalOut, signalSEOut]  = simulateP(sphereInfo,sim_box_side,patientInfo,delBzGridP,patientIndx,fieldGridStep,spill,numProtons,step,interval,TE,D,hepatocyte_dim,cellBoundaryFlag,cellBiasFlag,sinusoidBoundaryFlag);
%disp('running the old one');
%[signalOut, signalSEOut]  = simulatePoldrework(sphereInfo,sim_box_side,patientInfo,delBzGridP,patientIndx,fieldGridStep,spill,numProtons,step,interval,TE,D,hepatocyte_dim,cellBoundaryFlag,cellBiasFlag,sinusoidBoundaryFlag);

disp('---------------------End of Sim--------------------');

 
if simResultsFlag == 1

    % step  = 0.0005; %msec        very important factor !!!
    % interval = 60; %msec
    % t = 0:step:interval;
    % TE = logspace(log(0.1)/log(10),log(30)/log(10),15);
    
    %%%%%%%%%%%%%%%
    %%%% RESULTS
    %% T2star
    totSignal  = (1/numProtons) * sum(signal,1);           % complex fid
    % figure;plot(t,abs(totSignal));

    T2s_est = [1 10 20 50];

    for u=1:length(T2s_est)
        [S0(u),T2s(u),C(u),Res_1(u)] = fitexp_mc(t,abs(totSignal),T2s_est(u));
        % s_fit = expc([S0 T2s(u) C],t);
        % figure;plot(t,abs(totSignal));
        % hold on;
        % plot(t,s_fit,'r');
        % hold off;
    end

    [P,Q] = min(Res_1);     % choose fit with least residual

    simR2s(patientIndx) = 1000/(T2s(Q));

    %% T2 SE

    %% TE defined before main for-loop
    % signalSE = zeros(numProtons,length(t));
    % totPhaseSE = zeros(numProtons,length(t));    %% total phase accrual for T2 SE

    %% extracting echoes

    for m = 1:length(TE)

        TE_indx = round(TE(m)/step);

        totSignalSE = (1/numProtons) * signalSE(m,:);

        if(m==length(TE))       %% checking if end of observation time frame is reached
            [signalSEecho(m),signalSEindx(m)] = max(real(totSignalSE(2*TE_indx-100:2*TE_indx)));
        else
            [signalSEecho(m),signalSEindx(m)] = max(real(totSignalSE(2*TE_indx-100:2*TE_indx+100)));
        end

    end

    % figure;plot(t,abs(totSignalSE));

    tEcho = TE * 2;

    % figure;plot(tEcho,signalSEecho);

    T2_est = [1 10 20 50 80 120];


    for u = 1:length(T2_est)
        [S0,T2(u),C,Res_2(u)] = fitexp_mc(tEcho,signalSEecho,T2_est(u));
%         s_fit = expc([S0 T2(u) C],tEcho);
%         figure;plot(tEcho,signalSEecho,'o-');
%         hold on;
%         plot(tEcho,s_fit,'r');
%         hold off;
    end

    [P,Q] = min(Res_2);

    simR2(patientIndx) = 1000/T2(Q);

end



if CPMGsimResultsFlag == 1

    echo_time = [0.1 0.2 0.3 0.4 0.5 0.6 0.8 1 2 3 4 5 6 10 20];
    step  = 0.0005; %msec        very important factor !!!
    interval = 60; %msec
    t = 0:step:interval;


    for  echo_loop = 1:length(echo_time)

        clear tau;
        echo_spacing = echo_time(echo_loop);    % renamed, echo_time is actually echo_spacing
        tau = echo_spacing:echo_spacing*2:58;
        tau_indx = round(tau./step);

        totSignalCPMG = (1/numProtons) * signalCPMG(echo_loop,:);

        %%% Extracting peak amplitudes at echos

        q = [1 round((tau + echo_spacing)./step)];
        % figure;plot(t(q),abs(signalCPMG(echo_loop,q)));

        % Monoexp
        T2_est = [1 10 20 50 80 120];

        for p = 1:length(T2_est)
            [S0,T2,C,Res] = fitexp_mc(t(q),abs(signalCPMG(echo_loop,q)),T2_est(p));
            T2_mono(p) = T2;
            Res_mono(p) = Res;
            S0_mono(p) = S0;
            C_mono(p) = C;

        end

        [Y,I] = min(Res_mono.^2);

        T2_cpmg(echo_loop) = T2_mono(I);
        S0_cpmg(echo_loop) = S0_mono(I);
        C_cpmg(echo_loop)  = C_mono(I);

        ME(echo_loop).S0 = S0_mono(I);
        ME(echo_loop).T2 = T2_mono(I);
        ME(echo_loop).C = C_mono(I);

%                 % Calculate fitted function values
%                 s_fit = expc([S0_mono(I) T2_mono(I) C_mono(I)],t);
%         
%                 figure;
%                 plot(t(q),abs(signalCPMG(echo_loop,q)),'bo');
%                 hold on;
%                 % plot(t,abs(S));
%                 plot(t,s_fit,'r');
%                 hold off;

        % Biexp
        T2_low_est  = 1;%[1 2 10 15];
        T2_high_est = T2_cpmg(echo_loop);

        for p = 1:length(T2_low_est)

            T2_est = [T2_low_est(p) T2_high_est];
            [S0,T2,C,Res] = fitbiexp_mc(t(q),abs(signalCPMG(echo_loop,q)),T2_est);
            T2_bi(p,:) = T2;
            Res_bi(p) = Res;
            S0_bi(p,:) = S0;
            C_bi(p) = C;
        end

        [Y,I] = min(Res_bi.^2);


        BE(echo_loop).S0 = S0_bi(I,:);
        BE(echo_loop).T2 = T2_bi(I,:);
        BE(echo_loop).C  = C_bi(I);

        %                     % Calculate fitted function values
        %                     s_fit = biexpc([S0_bi(I,1) T2_bi(I,1) S0_bi(I,2) T2_bi(I,2) C_bi(I)],t);
        %
        %                     figure;
        %                     plot(t(q),abs(signalCPMG(echo_loop,q)),'go');
        %                     hold on;
        %                     % plot(t,abs(S));
        %                     plot(t,s_fit,'r');
        %                     hold off;

        % non-exponential model

        mono_est = [S0_cpmg(echo_loop) T2_cpmg(echo_loop) C_cpmg(echo_loop)];

        [S0,T2,a,c] = fit_nonexp(t(q),abs(signalCPMG(echo_loop,q)),mono_est);
        NE(echo_loop).S0 = S0;
        NE(echo_loop).T2 = T2;
        NE(echo_loop).a = a;
        NE(echo_loop).c = c;

        %                         data.tau = echo_spacing;
        %                         s_fit = nonexp([S0 T2 a c],t,data);
        %                         figure;
        %                         plot(t(q),abs(signalCPMG(echo_loop,q)),'o');
        %                         hold on;
        %                         % plot(t,abs(S));
        %                         plot(t,s_fit,'-r');
        %                         hold off;
% 


    end     % echo_loop


    %%% NE Results

    for p = 1:size(NE,2)-1
        temp_a = NE(p).a;
        this_a(p) = temp_a;
    end

    echo_time1 = echo_time(1:length(this_a));

    % figure;plot(echo_time1,log(this_a),'-o');

    figure;plot(echo_time1,(this_a),'-o');

    FE = patientInfo(patientIndx).FE;
    this_L = (0.22./this_a)*(FE/3.5)*1.5;
    %% FE is dry wt. , the formula uses wet wt., 1.5 is B0
    figure;plot(echo_time1,this_L,'-o');
    title('NE');
    ylabel('L (um)');
    xlabel('Echo time (ms)');

    
    %% Monoexp T2

    echo_time1 = echo_time(1:size(T2_cpmg,2));

    figure; plot( echo_time, 1000./T2_cpmg(:),'o:' );
    title('T2 cpmg monoexp');
    xlabel('tau (msec)');
    ylabel('R2');

    %% Biexp Amplitudes and T2

    for p = 1:size(BE,2)-1

        Amp(p,:)  = BE(p).S0;
        BiT2(p,:) = BE(p).T2;

    end

    echo_time1 = echo_time(1:size(Amp,1));

    figure; plot( echo_time1', Amp(:,1),'o:r' );
    hold on;
    plot( echo_time1', Amp(:,2),'o:' );
    title('Amp: cpmg Biexp');
    xlabel('tau (msec)');
    ylabel('Amp');

    figure; plot( echo_time1, BiT2(:,1),'o:r' );
    hold on;
    plot( echo_time1, BiT2(:,2),'o:' );
    title('T2: cpmg Biexp');
    xlabel('tau (msec)');
    ylabel('T2');


end     % if

