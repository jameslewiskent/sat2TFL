function [ITs,PP_FAs,IT_FAs,Segment_Sizes] = epg_Sandwich(TRs,T1s,T2,ratios,Velocity,Angle,Diff_co,Train_Size,Segment_Factor)
% EPG function to simulate image trains for Option 1 B1 mapping scheme
% Controls looping structure over different flip angles, T1s etc.
%
%   INFO:
%       Second adapation to Chung et al.'s B1 mapping method. Further
%       segmented image train and proton density image 'sandwiched' to
%       before the pre-pulse.
%
%	INPUT:
%       TRs - Array of Repetition times to simulate (s)
%       T1s - Array of T1 relaxation times to simulate (s)
%       T2 - T2 relaxation time (s)
%       ratios - Equal to 0 to produce a 'parameter' map of PP_FA vs IT_FA OR Array of fixed
%       ratios of PP_FA to IT_FA e.g. [5 10 20 30].
%       Velocity - Velocity of coherent flow (m/s)
%       Angle - Angle of coherent flow (rad)
%       Diff_co - Diffusion co-efficient (m^2/s) (isotropic)
%       Train_Size - Number of phase encode steps, i.e. how many total RF pulses
%       in an image train, which will be zero-filled up to Matrix_size
%
%
%	OUTPUT:
%       ITs - 1x2 Cell containing Image Train 1 and 2
%       PP_FAs - Array of Pre-pulse Flip Angles used to simulate (degrees)
%       IT_FAs - Array of Image Train Flip Angles used to simulate (degrees)
%       Segment_Sizes - Number of segments acquired and phase encodes in each segment(centrically reordered)
%
% J. Kent. 2021. Using B.Hargreaves EPG Functions.
nFreqSamples = 1;
RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
PP_FAs = linspace(1,180,180); % Pre-pulse Flip Angles in Degrees was 1 180
IT_TE = 1.78e-3; % Image Train Echo Time (s)
IT_TR = 4.2e-3; % Image Train Repitition Time (s)
Dummy_Scans = 2; % Add 'dummy' scans which help achieve steady-state
phi1 = zeros(1,Dummy_Scans+Segment_Factor); % Pre-Pulse RF phase (rad)
RF_Spoiling_Increment_Checks = 20;
phi2 = cumsum((0:(ceil(Train_Size/Segment_Factor)*2*Dummy_Scans +2*Train_Size+Segment_Factor+Dummy_Scans+RF_Spoiling_Increment_Checks))*50).*(pi./180); % Imaging phases (rad) (for RF spoiling- need enough for optional dummy scans too)
phi2 = phi2(RF_Spoiling_Increment_Checks:end);


prep_spoils = 5.* ones(1,Dummy_Scans+Segment_Factor+1); % 2*2.^(0:Segment_Factor) % Number of unit gradients to move through after pre-pulse (spoiling- need enough for optional dummy scans too)
train_spoils = 5.* ones(1,ceil(Train_Size/Segment_Factor));
noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
man_spoil = 1; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
HR_TR = 0; % Vary the TR as if it were ECG gated to a heartrate, 0 - off, 1 - on. Essentially, not strict TR but range e.g. a 1 second TR may be simulated as 1.2-1.4 s within the simulation. Vary randomly normally distributed.

Ejection_Frac = 0; % Ejection fraction on or off- reset magnetisation
magtrack_flag = zeros(1,5); % Set to 1 to prevent magnetisation tracking in sequence, otherwise will graph results for 6 different PP_FA's
RF_Phase = 0; % Used for mTx synthetic body simulations

% ---------------------------------------------------------------------- %
% -------------            END of USER input              -------------- %
% ---------------------------------------------------------------------- %

% Console output and checks etc.

% Calculate number of phase encodes in each segment and check
Train_Size_Tot = Train_Size;
Segment_Sizes = zeros(1,Segment_Factor);
for Seg_n = 1:Segment_Factor
    Segment_Sizes(Seg_n)= ceil(Train_Size_Tot./(Segment_Factor+1-Seg_n));
    Train_Size_Tot = Train_Size - sum(Segment_Sizes,'all');
end
if sum(Segment_Sizes,'all') ~= Train_Size
    error('ERROR: epg_func: Number of phase encodes in segments does not equal size of train requested.')
elseif Segment_Factor == 1
    disp(['Number of phase encodes in image train: ',num2str(Train_Size)]);
else
    disp(['Number of phase encodes in image train: ',num2str(Train_Size),', per segment: ',num2str(Segment_Sizes)]);
end

%
% if Segment_Factor == 4
%     train_spoils = 2.^(0:Train_Size/Segment_Factor);
% else
%     train_spoils = [1:1:1*(Train_Size/Segment_Factor)];
% end

% Notificaitons to user
if man_spoil == 1
    disp('WARNING: Assuming perfect spoiling')
end

if HR_TR == 1
    disp('Heartrate variable TR option is on, TR period non-constant')
end

if RF_Spoiling == 0
    disp('RF spoiling is turned off')
    phi2 = zeros(size(phi2));
end


% End of checks and command output, following is control of loops through T1/TR/Ratio values etc.

    % ----------------------- Parameter Data Sim ----------------------- %

if length(ratios) == 1 && ratios(1) == 0
    magtrack_flag = ones(1,size(magtrack_flag,2)); % Don't magtrack
    IT_FAs = linspace(0.1,20,200); % Image Train Flip Angles in Degrees
    IT1 = zeros(length(PP_FAs),length(IT_FAs),length(TRs),length(T1s),Train_Size);
    IT2 = zeros(length(PP_FAs),length(IT_FAs),length(TRs),length(T1s),Train_Size);

    for PP_FAs_n = 1:length(PP_FAs)
        PP_FA = PP_FAs(PP_FAs_n).*pi/180; % Pre-pulse FA in radians

        for IT_FAs_n = 1:length(IT_FAs)
            IT_FA = IT_FAs(IT_FAs_n).*pi/180; % Image Train FA in radians

            for TR_n = 1:length(TRs)
                TR = TRs(TR_n);

                for T1_n = 1:length(T1s)
                    T1 = T1s(T1_n);


                    [IT1(PP_FAs_n,IT_FAs_n,TR_n,T1_n,:),IT2(PP_FAs_n,IT_FAs_n,TR_n,T1_n,:)] = Op2_Seq(PP_FA,IT_FA,TR,T1,T2,IT_TR,IT_TE,Train_Size,phi1,phi2,kg,Diff_co,Velocity,Angle,prep_spoils,train_spoils,noadd,man_spoil,HR_TR,Dummy_Scans,Ejection_Frac,magtrack_flag,Segment_Sizes,nFreqSamples,RF_Phase);


                end
            end
        end
    end

    % ----------------------- T1 Parameter Data Sim ----------------------- %


elseif length(ratios) == 1 && ratios(1) ~=0 && size(T1s,2) >10
    magtrack_flag = ones(1,size(magtrack_flag,2)); % Don't magtrack
    IT_FAs = zeros(length(PP_FAs),length(ratios))';
    IT1 = zeros(length(PP_FAs),1,length(TRs),length(T1s),Train_Size);
    IT2 = zeros(length(PP_FAs),1,length(TRs),length(T1s),Train_Size);
    for PP_FAs_n = 1:length(PP_FAs)
        PP_FA = PP_FAs(PP_FAs_n).*pi/180; % Pre-pulse FA in radians

        IT_FA = (PP_FA./ratios(1)); % Image Train FA in radians

        for TR_n = 1:length(TRs) % TR unusued
            TR = TRs(TR_n);

            for T1_n = 1:length(T1s)
                T1 = T1s(T1_n);

                [IT1(PP_FAs_n,1,TR_n,T1_n,:),IT2(PP_FAs_n,1,TR_n,T1_n,:)] = Op2_Seq(PP_FA,IT_FA,TR,T1,T2,IT_TR,IT_TE,Train_Size,phi1,phi2,kg,Diff_co,Velocity,Angle,prep_spoils,train_spoils,noadd,man_spoil,HR_TR,Dummy_Scans,Ejection_Frac,magtrack_flag,Segment_Sizes,nFreqSamples,RF_Phase);

            end
        end
    end


    % ----------------------- ratio Data Sim ----------------------- %

else
    IT1 = zeros(length(PP_FAs),length(ratios),length(TRs),length(T1s),Train_Size);
    IT2 = zeros(length(PP_FAs),length(ratios),length(TRs),length(T1s),Train_Size);
    IT_FAs = zeros(length(PP_FAs),length(ratios))';
    for PP_FAs_n = 1:length(PP_FAs)
        PP_FA = PP_FAs(PP_FAs_n).*pi/180; % Pre-pulse FA in radians

        for ratio_n = 1:length(ratios)
            IT_FAs(ratio_n,:) = PP_FAs./ratios(ratio_n);  % Image Train Flip Angles in Degrees
            IT_FA = IT_FAs(ratio_n,PP_FAs_n).*pi/180; % Image Train FA in radians

            for TR_n = 1:length(TRs)
                TR = TRs(TR_n);

                for T1_n = 1:length(T1s)
                    T1 = T1s(T1_n);

                    [IT1(PP_FAs_n,ratio_n,TR_n,T1_n,:),IT2(PP_FAs_n,ratio_n,TR_n,T1_n,:),magtrack_flag] = Op2_Seq(PP_FA,IT_FA,TR,T1,T2,IT_TR,IT_TE,Train_Size,phi1,phi2,kg,Diff_co,Velocity,Angle,prep_spoils,train_spoils,noadd,man_spoil,HR_TR,Dummy_Scans,Ejection_Frac,magtrack_flag,Segment_Sizes,nFreqSamples,RF_Phase);

                end
            end
        end
    end
end

% Sequence Function
    function [Train1,Train2,magtrack_flag] = Op2_Seq(PP_FAs,IT_FAs,TR,T1,T2,IT_TR,IT_TE,Train_Size,phi1,phi2,kg,Diff_co,Velocity,Angle,prep_spoils,train_spoils,noadd,man_spoil,HR_TR,Dummy_Scans,Ejection_Frac,magtrack_flag,Segment_Sizes,nFreqSamples,RF_Phase)
        % Function runs inside relevent epg_SCHEME function and simulates a single set of image trains for specified parameters.

        P = [0 0 1]'; % Start magnetisation in equilibrium
        IT_rTE = IT_TR - IT_TE;
        tFreqSample = 0.015893e-3; % Sampling Time (s)
        TotFreqSampleTime = tFreqSample*nFreqSamples;
        phi1_n = 1;
        phi2_n = 1; % Phase increment counter
        prep_spoils_n =1; % pre_pulse spoiler counter
        Train1 = zeros(Train_Size,1,nFreqSamples);
        Train2 = zeros(Train_Size,1,nFreqSamples);
        RF_Time = 1e-3; % Time for RF pulse (just for magtrack)
        Spoil_Time = 1e-3;
        if HR_TR == 1
            seq_TRs = TR.*ones(1,Dummy_Scans+size(Segment_Sizes,2)).*(1.2 + 0.1*randn(1,Dummy_Scans+size(Segment_Sizes,2))); % Variable TRs

        elseif HR_TR == 0
            seq_TRs = TR.*ones(1,Dummy_Scans+size(Segment_Sizes,2));
        end

        %%% START OF SEQUENCE %%%

        [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FAs(1),T1,0,1,1,0,0,magtrack_flag); % Initial magtrack

        if Dummy_Scans ~= 0
            for Dummy_n = 1:Dummy_Scans

                
                    Mode_n = 1; % Set mode for dummy scans

                PP_FA = PP_FAs(Mode_n);
                IT_FA = IT_FAs(Mode_n);

                if Ejection_Frac == 1
                    if T1 == 2.29 % T1 of blood (hence currently simulating blood voxel)
                        P = [0 0 1]'; % basically have manual spoiling due to ejected/refreshed blood
                    end
                end

                % Relaxation between dummy scans
                T_relax = seq_TRs(Dummy_n) - RF_Time - 2.*Segment_Sizes(1).*IT_TR; % Time for Relaxation (s)
                if all(magtrack_flag) % if all flags are 1, true
                    P = epg_grelax(P,T1,T2,T_relax,kg,Diff_co,1,0,Velocity,Angle); % relaxation + spoiling
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,T_relax,magtrack_flag);
                else
                    % Split up relaxation for better plotting
                    for dT_relax = (T_relax/10).*ones(1,10)
                        P = epg_grelax(P,T1,T2,dT_relax,kg,Diff_co,1,0,Velocity,Angle); % relaxation + spoiling
                        [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,dT_relax,magtrack_flag);
                    end
                end

                % Pretend image trains
                for IT_n = 1:Segment_Sizes(1)
                    P = epg_rf(P,IT_FA,RF_Phase(Mode_n)+phi2(phi2_n)); % Imaging RF
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,RF_Time,magtrack_flag);

                    P = epg_grelax(P,T1,T2,IT_TR,kg,Diff_co,0,0,Velocity,Angle); % relaxation no spoiling
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,IT_TR,magtrack_flag);

                    % Dummy Scans not stored
                    for i = 1:train_spoils(IT_n)
                        P = epg_grad(P); % spoiling
                    end
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                    if man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                        [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                    end
                    phi2_n = phi2_n +1; % increment rf spoiling counter
                end

                % Pre-pulse for dummy scans
                P = epg_rf(P,PP_FA,RF_Phase(Mode_n)+phi1(phi1_n));
                [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,RF_Time,magtrack_flag);
                phi1_n = phi1_n +1; % increment rf spoiling counter
                for i = 1:prep_spoils(prep_spoils_n)
                    P = epg_grad(P); % spoiling
                end
                [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                prep_spoils_n = prep_spoils_n +1; % Increment prep spoils counter
                if man_spoil ==1
                    P = [0 0 P(3,1)]'; % manual spoiling
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                end

                % Pretend image trains
                for IT_n = 1:Segment_Sizes(1)
                    P = epg_rf(P,IT_FA,RF_Phase(Mode_n)+phi2(phi2_n)); % Imaging RF
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,RF_Time,magtrack_flag);

                    P = epg_grelax(P,T1,T2,IT_TR,kg,Diff_co,0,0,Velocity,Angle); % relaxation no spoiling
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,IT_TR,magtrack_flag);
                    % Dummy Scans not stored
                    for i = 1:train_spoils(IT_n)
                        P = epg_grad(P); % spoiling
                    end
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                    if man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                        [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                    end
                    phi2_n = phi2_n +1; % increment rf spoiling counter
                end
            end
        end % End of Dummy Scans


        for Mode_n = 1:1 % repeat for different modes of mTx array
            PP_FA = PP_FAs(Mode_n);
            IT_FA = IT_FAs(Mode_n);

            if Ejection_Frac == 1
                if T1 == 2.29 % T1 of blood (hence currently simulating blood voxel)
                    P = [0 0 1]'; % basically have manual spoiling due to ejected/refreshed blood
                end
            end

            for Segment_n = 1:size(Segment_Sizes,2)

                % Relaxation
                T_relax = seq_TRs(Dummy_Scans+Segment_n) - RF_Time - 2.*Segment_Sizes(Segment_n).*IT_TR; % Time for Relaxation (s)
                if all(magtrack_flag)
                    P = epg_grelax(P,T1,T2,T_relax,kg,Diff_co,1,0,Velocity,Angle); % relaxation + spoiling
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,T_relax,magtrack_flag);
                else
                    % Split up relaxation for better plotting
                    for dT_relax = (T_relax/10).*ones(1,10)
                        P = epg_grelax(P,T1,T2,dT_relax,kg,Diff_co,1,0,Velocity,Angle); % relaxation + spoiling
                        [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,dT_relax,magtrack_flag);
                    end
                end
                if man_spoil ==1
                    P = [0 0 P(3,1)]'; % manual spoiling
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                end


                % First quarter of Image Train 1
                for IT_n = sum(Segment_Sizes(1:Segment_n - 1))+1 : sum(Segment_Sizes(1:Segment_n))
                    P = epg_rf(P,IT_FA,RF_Phase(Mode_n)+phi2(phi2_n)); % Imaging RF
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,RF_Time,magtrack_flag);

                    P = epg_grelax(P,T1,T2,IT_TE-(TotFreqSampleTime./2),kg,Diff_co,0,0,Velocity,Angle); % relaxation no spoiling
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,IT_TE-(TotFreqSampleTime./2),magtrack_flag);


                    for nSample = 1:nFreqSamples
                        Train1(IT_n,Mode_n,nSample) = P(1,1).*exp(-1i*phi2(phi2_n)); % Store Phase-Demodulated 'signal' F+0
                        P = epg_grelax(P,T1,T2,tFreqSample,kg,Diff_co,0,0,Velocity,Angle); % relaxation no spoiling
                        [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,tFreqSample,magtrack_flag);
                    end

                    P = epg_grelax(P,T1,T2,IT_rTE-(TotFreqSampleTime./2),kg,Diff_co,0,0,Velocity,Angle); % relaxation no spoiling
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,IT_rTE-(TotFreqSampleTime./2),magtrack_flag);

                    for i = 1:train_spoils(Segment_Sizes(Segment_n))
                        P = epg_grad(P); % spoiling
                    end
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                    if man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                        [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                    end
                    phi2_n = phi2_n +1; % increment rf spoiling counter
                end


                % Pre-pulse
                P = epg_rf(P,PP_FA,RF_Phase(Mode_n)+phi1(phi1_n));
                [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,RF_Time,magtrack_flag);
                phi1_n = phi1_n +1; % increment rf spoiling counter
                for i = 1:prep_spoils(prep_spoils_n)
                    P = epg_grad(P); % spoiling
                end
                [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                prep_spoils_n = prep_spoils_n + 1;
                if man_spoil ==1
                    P = [0 0 P(3,1)]'; % manual spoiling
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                end


                % First quarter Image Train 2
                for IT_n = sum(Segment_Sizes(1:Segment_n - 1))+1 : sum(Segment_Sizes(1:Segment_n))
                    P = epg_rf(P,IT_FA,RF_Phase(Mode_n)+phi2(phi2_n));
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,RF_Time,magtrack_flag);

                    P = epg_grelax(P,T1,T2,IT_TE-(TotFreqSampleTime./2),kg,Diff_co,0,0,Velocity,Angle); % relaxation no spoiling
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,IT_TE-(TotFreqSampleTime./2),magtrack_flag);


                    for nSample = 1:nFreqSamples
                        Train2(IT_n,Mode_n,nSample) = P(1,1).*exp(-1i*phi2(phi2_n)); % Store Phase-Demodulated 'signal' F+0
                        P = epg_grelax(P,T1,T2,tFreqSample,kg,Diff_co,0,0,Velocity,Angle); % relaxation no spoiling
                        [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,tFreqSample,magtrack_flag);

                    end

                    P = epg_grelax(P,T1,T2,IT_rTE-(TotFreqSampleTime./2),kg,Diff_co,0,0,Velocity,Angle); % relaxation no spoiling
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,IT_rTE-(TotFreqSampleTime./2),magtrack_flag);

                    for i = 1:train_spoils(Segment_Sizes(Segment_n))
                        P = epg_grad(P); % spoiling
                    end
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                    if man_spoil ==1
                        P = [0 0 P(3,1)]'; % manual spoiling
                        [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
                    end
                    phi2_n = phi2_n +1; % increment rf spoiling counter
                end

                if Ejection_Frac == 1
                    if T1 == 2.29 % T1 of blood (hence currently simulating blood voxel)
                        P = [0 0 1]'; % basically have manual spoiling due to ejected/refreshed blood
                    end
                end



            end

            % Squeeze in the case of no frequency samples or transmit modes to remove dimensions of 1
            Train1 = squeeze(Train1);
            Train2 = squeeze(Train2);
        end
        %%% END OF SEQUENCE %%%
        [~,~,~,~,magtrack_flag] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,'End',magtrack_flag);
    end

ITs = {IT1,IT2};
end
