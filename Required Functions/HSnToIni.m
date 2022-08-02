clear all

HS8_pulse = GetHSn(8,15);
sinc_pulse = 3*sinc(linspace(-15,15,512)); % SINC ()
%sinc_pulse = sinc(linspace(-2,2,70)); % SINC ()
%sinc_pulse = 4*sin(linspace(-8.1.*pi,8.1.*pi,512))./linspace(-8.1.*pi,8.1.*pi,512); % SINC ()
rect_pulse = ones(1,50); % RECT (0.5 ms)
dT_s = 10e-6; % 10 us RF sampling time

Energy_sinc = sum(sinc_pulse.*conj(sinc_pulse).*dT_s); % Parseval Energy
Energy_HS8 = sum(HS8_pulse.*conj(HS8_pulse).*dT_s./2);
Energy_rect = sum(rect_pulse.*conj(rect_pulse).*dT_s);

rect_scal_fac = 1%sqrt(Energy_sinc./Energy_rect);
rect_pulse = rect_pulse.*rect_scal_fac;
HS8_scal_fac = sqrt(Energy_sinc./Energy_HS8);
HS8_pulse = HS8_pulse.*HS8_scal_fac;

Energy_sinc = sum(sinc_pulse.*conj(sinc_pulse).*dT_s); % Parseval Energy
Energy_HS8 = sum(HS8_pulse.*conj(HS8_pulse).*dT_s./2);
Energy_rect = sum(rect_pulse.*conj(rect_pulse).*dT_s);

figure()
plot((1:size(HS8_pulse,2)).*dT_s/2,abs(HS8_pulse),'b'); hold on; plot((1:size(sinc_pulse,2)).*dT_s,abs(sinc_pulse),'g'); hold on;  plot((1:size(rect_pulse,2)).*dT_s,rect_pulse,'r'); hold on; 

abs(HS8_pulse(size(HS8_pulse,2)/2))./sinc_pulse(size(sinc_pulse,2)/2)
%%

dB1 = 1; b1_range_Hz = 0:dB1:1000;
freq_range_Hz = -5000:10:5000; [~,iCFreq ] = min(abs(freq_range_Hz)); % Centre frequency
diffFreq = 500; [~,idiffFreq ] = min(abs(freq_range_Hz-diffFreq)); % Centre frequency

[mxy_HS8, mz_HS8] = SimulatePulse(HS8_pulse,b1_range_Hz,freq_range_Hz,dT_s/2); % HS8 (total 5ms)
[mxy_sinc, mz_sinc] = SimulatePulse(sinc_pulse,b1_range_Hz,freq_range_Hz,dT_s); % SINC (5 ms)
%[mxy_rect, mz_rect] = SimulatePulse(rect_pulse,b1_range_Hz,freq_range_Hz,dT_s); % RECT (0.5 ms)

%%
figure('Units','normalized','Color','w','Position',[0.2,0.2,0.571354166666667,0.521296296296296]) 
tiledlayout(2,2,'TileSpacing','compact','Padding','none')
%nexttile(1); imagesc(freq_range_Hz,b1_range_Hz,mz_rect,[-1 1])
%xlabel('B_0 (Hz)'); ylabel('B_1 (Hz)'); title('a) RECT (0.5 ms)');
nexttile(1); imagesc(freq_range_Hz,b1_range_Hz,mz_sinc,[-1 1])
xlabel('B_0 (Hz)'); ylabel('B_1 (Hz)'); title('a) Sinc');
xticks([-5e3:1e3:5e3])
nexttile(3); imagesc(freq_range_Hz,b1_range_Hz,mz_HS8,[-1 1])
xlabel('B_0 (Hz)'); ylabel('B_1 (Hz)'); title('b) HS8');
xticks([-5e3:1e3:5e3])
nexttile(2,[2,1]); plot(b1_range_Hz,acos(mz_HS8(:,iCFreq)).*180/pi,'b'); hold on; plot(b1_range_Hz,acos(mz_HS8(:,idiffFreq)).*180/pi,'b--');
%plot(b1_range_Hz,acos(mz_rect(:,iCFreq)).*180/pi,'r'); plot(b1_range_Hz,acos(mz_rect(:,idiffFreq)).*180/pi,'r--');
plot(b1_range_Hz,acos(mz_sinc(:,iCFreq)).*180/pi,'r'); plot(b1_range_Hz,acos(mz_sinc(:,idiffFreq)).*180/pi,'r--');
xlabel('B_1 (Hz)'); ylabel(['Flip angle (',char(176),')']); legend('HS8 (f_0)',['HS8 (f_0',char(177),num2str(diffFreq),' Hz)'],'SINC (f_0)',['SINC (f_0',char(177),num2str(diffFreq),' Hz)'],'Location','southeast','NumColumns',2); ylim([0 180])

%%

figure('Units','normalized','Color','w','Position',[0.2,0.2,0.571354166666667,0.521296296296296]) 
tiledlayout(3,2,'TileSpacing','compact','Padding','none')
nexttile(1); imagesc(freq_range_Hz,b1_range_Hz,mz_rect,[-1 1])
xlabel('B_0 (Hz)'); ylabel('B_1 (Hz)'); title('a) RECT (0.5 ms)');
nexttile(2); imagesc(freq_range_Hz,b1_range_Hz,mz_sinc,[-1 1])
xlabel('B_0 (Hz)'); ylabel('B_1 (Hz)'); title('a) Sinc');
nexttile(3); imagesc(freq_range_Hz,b1_range_Hz,mz_HS8,[-1 1])
xlabel('B_0 (Hz)'); ylabel('B_1 (Hz)'); title('b) HS8');
