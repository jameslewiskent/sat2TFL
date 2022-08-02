function [Mxy,Mz,M_n,Time,magtrack_flag] = epg_magtrack(P,PP_FA,T1,Mxy,Mz,M_n,Time,AddTime,magtrack_flag)
% Function spread throughout epg_SCHEME to track evolution of longitudinal
% and transverse magnetisation.
%   Just saves taking up space in main sequence

% Values to track

PP_FA_track = [20,60,90,120,160]*(pi/180);
T1_track = 2;%[0.5,1,1.5,2,2.5,3]; 

PP_FA_track = PP_FA_track(1,1:size(magtrack_flag,2)); %resize if flags are smaller
[~,minind] = min(magtrack_flag);

if ((( PP_FA == PP_FA_track(minind) && T1 == T1_track) || PP_FA == 3.142 || ((T1 == 2.29 || T1 == 1.925) &&  (PP_FA - PP_FA_track(minind)) < 0.05 )) && magtrack_flag(minind) == 0) % Could add in other flip angles later

if M_n ~= 1 && ~ischar(AddTime)
Time(M_n) = Time(M_n-1) + AddTime;
Mxy(M_n) = P(1,1);
Mz(M_n) = P(3,1);



elseif ischar(AddTime) % E.g. if AddTime is 'End'
    
    if minind == 1 && magtrack_flag(minind) ~= 1
    %disp(['Plotting magnetisation from epg_magtrack.m for ',num2str(size(PP_FA_track,2)), ' tracked prep-pulse flip angles.'])
    figure('Name','Magtrack')
    end
    
    % These were used to generate data to plot T1 recovery curves for
    % example figure
%     save(['T1Time',num2str(minind),'PrepFA',num2str(PP_FA*180/pi),'Mxy'],'Mxy');
% save(['T1Time',num2str(minind),'PrepFA',num2str(PP_FA*180/pi),'Mz'],'Mz');
% save(['T1Time',num2str(minind),'PrepFA',num2str(PP_FA*180/pi),'Time'],'Time');
    
    subplot(2,3,minind)
    hold on
    plot(Time, abs(Mz),'r')
    plot(Time, abs(Mxy),'b')
    ylabel('Amplitude')
    xlabel('Time (s)');
    title(['Prep-pulse flip angle: ', num2str(PP_FA*180/pi),', T1: ',num2str(T1)])
    axis tight
    ylim([0, 1])
    hold off
    if minind == 1
    legend('Mz','Mxy');
    end
    
    magtrack_flag(minind) = 1; % Sets flag so only ran once per PP_FA
end
M_n = M_n + 1;
end

end

