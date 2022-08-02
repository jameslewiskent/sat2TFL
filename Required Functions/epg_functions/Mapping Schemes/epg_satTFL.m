function [ITs,PP_FAs,IT_FAs,Segment_Sizes] = epg_satTFL(TRs,T1s,T2,ratios,Velocity,Angle,Diff_co,Train_Size,Segment_Factor)
% EPG function to simulate image trains for reference B1 mapping scheme.
% Controls looping structure over different flip angles, T1s etc.
%
%   INFO:
%       Replicates Chung, S. et al. MRM. 2010. 64(2):439-446. original
%       pre-pulse method of B1 mapping using a 5T1 TR to allow for full
%       relaxation of the longitudinal magnesation. Can turn on optional
%       saturation recovery module and alter the Delay Time (TD_SR).
%
%	INPUT:
%       TRs - Array of Repetition times to simulate (s) (NOT used for Ref sequence)
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
%	OUTPUT:
%       ITs - 1x2 Cell containing Image Train 1 and 2
%       PP_FAs - Array of Pre-pulse Flip Angles used to simulate (degrees)
%       IT_FAs - Array of Image Train Flip Angles used to simulate (degrees)
%       Segment_Sizes - Number of segments acquired and phase encodes in each segment(centrically reordered)
%
% J. Kent. 2021. Using B.Hargreaves EPG Functions.
nFreqSamples = 1;
RF_Spoiling = 1; % RF_Spoiling on (1) or off (0)
PP_FAs = linspace(1,180,180); % Pre-pulse Flip Angles in Degrees
IT_TE = 1.78e-3; % Image Train Echo Time (s)
IT_TR = 4.2e-3; % Image Train Repitition Time (s)
TR_5T1 = 10; % TR time (s)
phi1 = 0; % Pre-pulse RF phase
RF_Spoiling_Increment_Checks = 20; % Accounts for checks run on scanner which increment RF spoiling increment
phi2 = cumsum((0:(ceil(Train_Size/Segment_Factor)*2 +2*Train_Size+Segment_Factor+RF_Spoiling_Increment_Checks))*50).*(pi./180); % Imaging phases (rad) (for RF spoiling- need enough for optional dummy scans too)
phi2 = phi2(RF_Spoiling_Increment_Checks:end);

prep_spoils = 5.* ones(1,Segment_Factor+1); %2*2.^(0:Segment_Factor); % Number of unit gradients to move through after pre-pulse (spoiling- need enough for optional dummy scans too)
train_spoils = 5.* ones(1,ceil(Train_Size/Segment_Factor));%1:1:1*ceil(Train_Size/Segment_Factor);
noadd = 0; % = 1 to NOT add any higher-order states to spoilers. This speeds up simulations, compromises accuracy!
man_spoil = 1; % Sets transverse magnetisation to 0 (if =1) assumes 'perfect' spoiling - useful for debugging
kg = 50e3; % k-space traversal due to gradient (rad/m) for diffusion/flow
SR_Module = 0; % Additional 'PERFECT' Saturation recovery module, destroys all longitudinal magnetisation prior to
TD_SR = 950e-3; % Saturation recovery time delay (s)
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


if SR_Module == 1
    disp('Reference Scheme is using saturation recovery module, the simulation assumes perfect saturation i.e. no longitudinal magnetisation, which could be unrealistic')
end

if RF_Spoiling == 0
    disp('RF spoiling is turned off')
    phi2 = zeros(size(phi2));
end

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
            
            for TR_n = 1:length(TRs) % TR unusued
                
                for T1_n = 1:length(T1s)
                    T1 = T1s(T1_n);
                    
                     [IT1(PP_FAs_n,IT_FAs_n,TR_n,T1_n,:),IT2(PP_FAs_n,IT_FAs_n,TR_n,T1_n,:)] = Ref_Seq(PP_FA,IT_FA,T1,T2,TR_5T1,IT_TR,IT_TE,Train_Size,phi1,phi2,kg,Diff_co,Velocity,Angle,prep_spoils,train_spoils,noadd,man_spoil,SR_Module,TD_SR,magtrack_flag,Segment_Sizes,nFreqSamples,RF_Phase);
            
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
                
                for T1_n = 1:length(T1s)
                    T1 = T1s(T1_n);
                    
                    [IT1(PP_FAs_n,1,TR_n,T1_n,:),IT2(PP_FAs_n,1,TR_n,T1_n,:)] = Ref_Seq(PP_FA,IT_FA,T1,T2,TR_5T1,IT_TR,IT_TE,Train_Size,phi1,phi2,kg,Diff_co,Velocity,Angle,prep_spoils,train_spoils,noadd,man_spoil,SR_Module,TD_SR,magtrack_flag,Segment_Sizes,nFreqSamples,RF_Phase);
                   
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
            
            for TR_n = 1:length(TRs) % TR unusued
                
                for T1_n = 1:length(T1s)
                    T1 = T1s(T1_n);
                    
                    [IT1(PP_FAs_n,ratio_n,TR_n,T1_n,:),IT2(PP_FAs_n,ratio_n,TR_n,T1_n,:),magtrack_flag] = Ref_Seq(PP_FA,IT_FA,T1,T2,TR_5T1,IT_TR,IT_TE,Train_Size,phi1,phi2,kg,Diff_co,Velocity,Angle,prep_spoils,train_spoils,noadd,man_spoil,SR_Module,TD_SR,magtrack_flag,Segment_Sizes,nFreqSamples,RF_Phase);
                    
                end
            end
        end
    end
    
end

% Sequence Function
    function [Train1,Train2,magtrack_flag] = Ref_Seq(PP_FAs,IT_FAs,T1,T2,TR_5T1,IT_TR,IT_TE,Train_Size,phi1,phi2,kg,Diff_co,Velocity,Angle,prep_spoils,train_spoils,noadd,man_spoil,SR_Module,TD_SR,magtrack_flag,Segment_Sizes,nFreqSamples,RF_Phase)
        % Function runs inside relevent epg_SCHEME function and simulates a single set of image trains for specified parameters.
        
        IT_rTE = IT_TR - IT_TE; 
        tFreqSample = 0.015893e-3; % Sampling Time (s)
        TotFreqSampleTime = tFreqSample*nFreqSamples;
        P = [0 0 1]'; % Start magnetisation in equilibrium
        phi_n = 1; % Phase increment counter
        Train1 = zeros(Train_Size,1);
        Train2 = zeros(Train_Size,1);
        TR = TR_5T1; % Time between image trains is fixed to TR_5T1 value, which may not be 5 * T1
        RF_Time = 1e-3; % Time for RF pulse (just for magtrack)
        Spoil_Time = 1e-3;
        [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA(1),T1,0,1,1,0,0,magtrack_flag); % Initial magtrack
        
        
                for Mode_n = 1:1 % repeat for different modes of mTx array
        PP_FA = PP_FAs(Mode_n);
        IT_FA = IT_FAs(Mode_n); 
        
        
        for Segment_n = 1:size(Segment_Sizes,2)
        

        if SR_Module == 1
            P = [0 0 0]'; % Assume 'PERFECT' saturation, all longitudinal signal is nulled
            
            [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,IT_TR,magtrack_flag); % Time assumed as IT_TR..?
            
            if magtrack_flag == 1
                % TD_SR relaxation
            P = epg_grelax(P,T1,T2,TD_SR,kg,Diff_co,0,0,Velocity,Angle); % relaxation no spoiling
            [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,TD_SR,magtrack_flag);
            elseif magtrack_flag == 0
                % break up relaxation for better plotting
                for dTD_SR = (TD_SR./10).*ones(1,10) % break up relaxation for better magtrack plots
                    P = epg_grelax(P,T1,T2,dTD_SR,kg,Diff_co,0,0,Velocity,Angle); % relaxation + spoiling (Period is 5 x T1 seconds)
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,dTD_SR,magtrack_flag);
                end
            end
            TR = 0; % TR set to zero, because it has relaxation during saturation recovery time delay instead
        else
            % Period of relaxation (this does nothing is SR_module is active
            if all(magtrack_flag) % if all flags are 1, true
                P = epg_grelax(P,T1,T2,TR,kg,Diff_co,0,0,Velocity,Angle); % relaxation + spoiling (Period is 5 x T1 seconds)
                [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,TR,magtrack_flag);
            else
                for dTR = (TR./10).*ones(1,10) % break up relaxation for better magtrack plots
                    P = epg_grelax(P,T1,T2,dTR,kg,Diff_co,0,0,Velocity,Angle); % relaxation + spoiling (Period is 5 x T1 seconds)
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,dTR,magtrack_flag);
                end
            end
        end
        
        
        
        % Image Train 1
        for IT_n = sum(Segment_Sizes(1:Segment_n - 1))+1 : sum(Segment_Sizes(1:Segment_n))
            P = epg_rf(P,IT_FA,RF_Phase(Mode_n)+phi2(phi_n)); % Imaging RF
            
            [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,RF_Time,magtrack_flag);
            
            P = epg_grelax(P,T1,T2,IT_TE-(TotFreqSampleTime./2),kg,Diff_co,0,0,Velocity,Angle); % relaxation no spoiling (FpFmZ,T1,T2,RelaxTime,kg,Diff_co,Gon,noadd)
            [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,IT_TE-(TotFreqSampleTime./2),magtrack_flag);
            
            for nSample = 1:nFreqSamples
                Train1(IT_n,Mode_n,nSample) = P(1,1).*exp(-1i*phi2(phi_n)); % Store Phase-Demodulated 'signal' F+0
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
            end
            phi_n = phi_n +1; % increment rf spoiling counter
        end
        
        end
        
        for Segment_n = 1:size(Segment_Sizes,2)
        
        if SR_Module == 1
            P(3,1) = 0 ; % Assume 'PERFECT' saturation, all longitudinal signal is nulled
            
            [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,IT_TR,magtrack_flag);
            
            % TD_SR relaxation
            if magtrack_flag == 1
            P = epg_grelax(P,T1,T2,TD_SR,kg,Diff_co,0,0,Velocity,Angle); % relaxation no spoiling
            [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,TD_SR,magtrack_flag);
            elseif magtrack_flag == 0
                % break up relaxation for better plotting
                for dTD_SR = (TD_SR./10).*ones(1,10) % break up relaxation for better magtrack plots
                    P = epg_grelax(P,T1,T2,dTD_SR,kg,Diff_co,0,0,Velocity,Angle); % relaxation + spoiling (Period is 5 x T1 seconds)
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,dTD_SR,magtrack_flag);
                end
            end
            
            TR = 0; % TR set to zero, because it has relaxation during saturation recovery time delay instead
        else
            % Period of relaxation (this does nothing is SR_module is active
            if all(magtrack_flag) % if all flags are 1, true
                P = epg_grelax(P,T1,T2,TR,kg,Diff_co,0,0,Velocity,Angle); % relaxation + spoiling (Period is 5 x T1 seconds)
                [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,TR,magtrack_flag);
            else
                for dTR = (TR./10).*ones(1,10) % break up relaxation for better magtrack plots
                    P = epg_grelax(P,T1,T2,dTR,kg,Diff_co,0,0,Velocity,Angle); % relaxation + spoiling (Period is 5 x T1 seconds)
                    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,dTR,magtrack_flag);
                end
            end
            
        end
    

% Pre-pulse
P = epg_rf(P,PP_FA,RF_Phase(Mode_n)+phi1);

[Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,RF_Time,magtrack_flag);

for i = 1:prep_spoils
    P = epg_grad(P); % spoiling
end
[Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,Spoil_Time,magtrack_flag);
if man_spoil ==1
    P = [0 0 P(3,1)]'; % manual spoiling
end



% Image Train 2
for IT_n = sum(Segment_Sizes(1:Segment_n - 1))+1 : sum(Segment_Sizes(1:Segment_n))
    P = epg_rf(P,IT_FA,RF_Phase(Mode_n)+phi2(phi_n)); % Imaging RF
    
    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,RF_Time,magtrack_flag);
    
    P = epg_grelax(P,T1,T2,IT_TE-(TotFreqSampleTime./2),kg,Diff_co,0,0,Velocity,Angle); % relaxation
    [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,IT_TE-(TotFreqSampleTime./2),magtrack_flag);
    
    for nSample = 1:nFreqSamples
        Train2(IT_n,Mode_n,nSample) = P(1,1).*exp(-1i*phi2(phi_n)); % Store Phase-Demodulated 'signal' F+0
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
    end
    phi_n = phi_n +1; % increment rf spoiling counter
end
               

            % Just inserted extra relxation for graphic figure of T1 curve
            % Period of relaxation (this does nothing is SR_module is active
%             if all(magtrack_flag) % if all flags are 1, true
%                 P = epg_grelax(P,T1,T2,TR,kg,Diff_co,0,0,Velocity,Angle); % relaxation + spoiling (Period is 5 x T1 seconds)
%                 [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,TR,magtrack_flag);
%             else
%                 for dTR = (TR./10).*ones(1,10) % break up relaxation for better magtrack plots
%                     P = epg_grelax(P,T1,T2,dTR,kg,Diff_co,0,0,Velocity,Angle); % relaxation + spoiling (Period is 5 x T1 seconds)
%                     [Mxy,Mz,M_n,Time] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,dTR,magtrack_flag);
%                 end
%             end
            
        end
[~,~,~,~,magtrack_flag] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,'End',magtrack_flag);
       % Squeeze in the case of no frequency samples to remove dims of 1
       Train1 = squeeze(Train1);
       Train2 = squeeze(Train2);
                end
end

ITs = {IT1,IT2};
end
