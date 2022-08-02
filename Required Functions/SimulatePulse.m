function [mxy, mz] = SimulatePulse(pulse,b1_range_Hz,freq_range_Hz,dT_s)
%% Simulate an RF pulse across 
addpath('C:\Users\jkent\OneDrive - Nexus365\Bloch Simulation\Hargreaves_bloch');

pulse_length = length(pulse);

gmmaHzPerG = 1e2*42.57747892;  %1e2 is MHz/T to Hz/G

nB1 = length(b1_range_Hz);
nF = length(freq_range_Hz);
% Hz to G = b1_Hz/Gamma
mxy = zeros(nB1,nF);
mz = mxy;

for iB1 = 1:nB1
        b1 = pulse*b1_range_Hz(iB1)/gmmaHzPerG; % excitation RF field in units of Hz (flip angle = timestep*b1*360deg)
        gr = zeros(size(pulse)); % gradients (NOT USED)
        tp = ones(size(pulse))*dT_s; % time (in seconds)
        t1 = 2; % T1 (in seconds)
        t2 = 0.05; % T2 (in seconds)
        df = freq_range_Hz; % frequency offset from zero (in Hz)
        dp = 0; % position offset (would combine with gr but not used at the moment)
        mode = 0; % 2 (not used)
    
        [tmx,tmy,tmz]=bloch(b1,gr,tp,t1,t2,df,dp,mode,0,0,1);
        mz(iB1,:) = tmz;
        mxy(iB1,:) = tmx+1i*tmy;
    
end

figure(1141)
imagesc(freq_range_Hz,b1_range_Hz,mz)
xlabel('Frequency (Hz)')
ylabel('B1 (Hz)')