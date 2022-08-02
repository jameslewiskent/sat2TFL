% Code to recreate simulation figures for 
% 'Rapid 3D Absolute B1+ Mapping using a Sandwiched Train Pre-Saturated
% TurboFLASH Sequence at 7T for the Brain and Heart' 
%
% James L. Kent (a),Iulius Dragonu (b), Ladislav Valkovič (c,d) and Aaron T. Hess (a)
%
% (a) Wellcome Centre for Integrative Neuroimaging, FMRIB, Nuffield Department
% of Clinical Neurosciences, University of Oxford, Oxford, United Kingdom
% (b) Siemens Healthcare Limited, Frimley, United Kingdom 
% (c) Oxford Centre for Clinical Magnetic Resonance Research (OCMR),
% University of Oxford, Oxford, United Kingdom 
% (d) Department of Imaging Methods, Institute of Measurement Science, 
% Slovak Academy of Sciences, Bratislava, Slovakia
% 
% Correspondence to Aaron T. Hess, Wellcome Centre for Integrative
% Neuroimaging, FMRIB, Nuffield Department of Clinical Neurosciences, John
% Radcliffe Hospital, Hedley Way, Oxford OX3 9DU, UK.
% Email: aaron.hess@ndcn.ox.ac.uk

% Simulations were based on Brian Hargreaves EPG code (http://web.stanford.edu/~bah/software/epg). 

% Repeats reduced to speed up simulations

cd(fileparts(matlab.desktop.editor.getActiveFilename))
addpath(genpath('Required Functions'))

%% Generate and Plot Figure 2

HS8_pulse = GetHSn(8,15);
sinc_pulse = 3*sinc(linspace(-15,15,512)); % SINC ()
rect_pulse = ones(1,50); % RECT (0.5 ms)
dT_s = 10e-6; % 10 us RF sampling time

Energy_sinc = sum(sinc_pulse.*conj(sinc_pulse).*dT_s)/50; % Energy
Energy_HS8 = sum(HS8_pulse.*conj(HS8_pulse).*dT_s./2)/50; % HS8 has half sample time due to twice number of samples (1024)
Energy_rect = sum(rect_pulse.*conj(rect_pulse).*dT_s)/50;

% Scale HS8 to have same energy at Sinc pulse
HS8_scal_fac = sqrt(Energy_sinc./Energy_HS8);
HS8_pulse = HS8_pulse.*HS8_scal_fac;

Energy_sinc = sum(sinc_pulse.*conj(sinc_pulse).*dT_s)/50; % Energy
Energy_HS8 = sum(HS8_pulse.*conj(HS8_pulse).*dT_s./2)/50;

figure()
plot((1:size(HS8_pulse,2)).*dT_s/2,abs(HS8_pulse),'b'); hold on; plot((1:size(sinc_pulse,2)).*dT_s,abs(sinc_pulse),'g'); hold on;  plot((1:size(rect_pulse,2)).*dT_s,rect_pulse,'r'); hold on; 

HS8relPeakB1 = abs(HS8_pulse(size(HS8_pulse,2)/2))./sinc_pulse(size(sinc_pulse,2)/2);

% Bloch Simulate HS8 and Sinc pulses
dB1 = 1; b1_range_Hz = 0:dB1:1500;
freq_range_Hz = -5000:10:5000; [~,iCFreq ] = min(abs(freq_range_Hz)); % Centre frequency
diffFreq = 500; [~,idiffFreq ] = min(abs(freq_range_Hz-diffFreq)); % Centre frequency
[mxy_HS8, mz_HS8] = SimulatePulse(HS8_pulse,b1_range_Hz,freq_range_Hz,dT_s/2); % HS8 (total 5ms)
[mxy_sinc, mz_sinc] = SimulatePulse(sinc_pulse,b1_range_Hz,freq_range_Hz,dT_s); % SINC (5 ms)

figure('Units','normalized','Color','w','Position',[0.2,0.2,0.571354166666667,0.521296296296296]) 
tiledlayout(2,2,'TileSpacing','compact','Padding','none')
nexttile(1); imagesc(freq_range_Hz,b1_range_Hz,mz_sinc,[-1 1])
xlabel('B_0 (Hz)'); ylabel('B_1 (Hz)'); title('a) Sinc');
xticks((-5e3:1e3:5e3))
nexttile(3); imagesc(freq_range_Hz,b1_range_Hz,mz_HS8,[-1 1])
xlabel('B_0 (Hz)'); ylabel('B_1 (Hz)'); title('b) HS8');
xticks((-5e3:1e3:5e3))
nexttile(2,[2,1]); plot(b1_range_Hz,acos(mz_HS8(:,iCFreq)).*180/pi,'b'); hold on; plot(b1_range_Hz,acos(mz_HS8(:,idiffFreq)).*180/pi,'b--');
plot(b1_range_Hz,acos(mz_sinc(:,iCFreq)).*180/pi,'r'); plot(b1_range_Hz,acos(mz_sinc(:,idiffFreq)).*180/pi,'r--');
xlabel('B_1 (Hz)'); ylabel(['Flip angle (',char(176),')']); legend('HS8 (f_0)',['HS8 (f_0',char(177),num2str(diffFreq),' Hz)'],'SINC (f_0)',['SINC (f_0',char(177),num2str(diffFreq),' Hz)'],'Location','southeast','NumColumns',2); ylim([0 180])

%% Generate and Plot Figures 3 and S1

ratios = 0; % Simulate parameter map (0) or specific ratio
T1s = 2; % T1 values to simulate (s)
T2 = 25e-3; % T2 (s)
TRs = 1; % Repetition times (s) (Not used for satTFL Schemes which has hard coded TR of 10 s or 1 s for Short TR)
repeats = 200; % Number of noise repeats to average
Noise_pcs = [0,0.004]; % Simulated Noise Levels 0.004
Scheme = {'satTFL','ShortTR','Sandwich'}; % Simulate Chosen Pulse Sequence(s) 'satTFL', 'ShortTR', 'Sandwich'
Segment_Factor = [1,4,4];
Reordering = 'CentricOut'; % Reordering of phase encodes, 'CentricOut' or 'Linear'.
Phase_Resolution =  0.75;
Partial_Fourier =  0.75; % 4/8, 5/8, 6/8, 7/8 or 1
Matrix_Size = 64; % Zero-filled up to this size if PF <1

if ~exist('Figure3Data', 'dir')
    mkdir('Figure3Data')
else
    addpath('Figure3Data')
end
Kappa_error_results = cell(2,size(Scheme,2));
for scheme_n = 1:size(Scheme,2)
[mean_Kappa_error,std_Kappa_error,PP_FAs,IT_FAs] = GenDataAndAnalyse(Scheme(scheme_n), ratios, T1s, T2, TRs, Noise_pcs, repeats,Reordering,Phase_Resolution, Partial_Fourier, Matrix_Size,Segment_Factor(scheme_n));
Kappa_error_results{1,scheme_n} = mean_Kappa_error;
Kappa_error_results{2,scheme_n} = std_Kappa_error;
end
save('Figure3Data/Figure3DataTopRow','PP_FAs','IT_FAs','Kappa_error_results');

ratios = 10; % Simulate parameter map (0) or specific ratio
T1s = 0.5:0.05:3; % T1 values to simulate (s)
Kappa_error_results = cell(2,size(Scheme,2));
for scheme_n = 1:size(Scheme,2)
[mean_Kappa_error,std_Kappa_error,PP_FAs,IT_FAs] = GenDataAndAnalyse(Scheme(scheme_n), ratios, T1s, T2, TRs, Noise_pcs, repeats,Reordering,Phase_Resolution, Partial_Fourier, Matrix_Size,Segment_Factor(scheme_n));
Kappa_error_results{1,scheme_n} = mean_Kappa_error;
Kappa_error_results{2,scheme_n} = std_Kappa_error;
end
save('Figure3Data/Figure3DataBottomRow','PP_FAs','IT_FAs','Kappa_error_results');

PlotFigure3(T1s)
PlotFigureS1(T1s)

%% Generate and Plot Figure 4 
% This could take a while to run, the number repeats has been reduced to
% reduce this.
% If the Parallel Computing Toolbox is not available, replace parfor loop
% with a regular for loop.

ratios = 1:0.5:40; % Simulate parameter map (0) or specific ratio
T1s = 2; % T1 values to simulate (s)
T2 = 25e-3; % T2 (s)
TRs = 1; % Repetition times (s) (Not used for satTFL Schemes which has hard coded TR of 10 s or 1 s for Short TR)
repeats = 100; % Number of noise repeats to average
Noise_pcs = 0.004; % Simulated Noise Levels 0.004
Scheme = 'Sandwich'; % Simulate Chosen Pulse Sequence(s) 'satTFL', 'ShortTR', 'Sandwich'
Segment_Factor = 4;
Reordering = 'CentricOut'; % Reordering of phase encodes, 'CentricOut' or 'Linear'.
Phase_Resolution =  1;
Partial_Fourier =  1; % 4/8, 5/8, 6/8, 7/8 or 1
Matrix_Sizes = 8:4:256; % Zero-filled up to this size if PF <1
parfor m = 1:size(Matrix_Sizes,2)
    Matrix_Size  = Matrix_Sizes(m);
if ~exist('Figure4Data', 'dir')
    mkdir('Figure4Data')
else
    addpath('Figure4Data')
end
[mean_Kappa_error,std_Kappa_error,~,~,~] = GenDataAndAnalyse(Scheme, ratios, T1s, T2, TRs, Noise_pcs, repeats,Reordering,Phase_Resolution, Partial_Fourier, Matrix_Size,Segment_Factor);
disp(num2str(Matrix_Size))
savedataparfor(Matrix_Size,mean_Kappa_error,std_Kappa_error);
end
PlotFigure4(ratios,Matrix_Sizes)

