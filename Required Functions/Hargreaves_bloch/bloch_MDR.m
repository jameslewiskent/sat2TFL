function [outmx,outmy,outmz] = bloch_MDR(b1,gr,tp,T1,T2,df,dp,mode,mx,my,mz)
% bloch_MDR
% this is the same as Hargreaes code, but hopefully this works!
% It is a simple forward simulation of the Bloch equations, nothing clever
% operations assume all RF is instantaneous
%
%   Detailed explanation goes here
%
% b1 = excitation RF field in units of Hz (flip angle = timestep*b1*360deg)
% gr= gradients (NOT USED)
% tp = time (in seconds), this is a list of times and corresponds to the times in the
% b1, output data are output at each of these times
% t1 = T1 (in seconds)
% t2 = T2 (in seconds)
% df = frequency offset from zero (in Hz)
% dp = position offset (would combine with gr but not used at the moment)
%mode = 2 (not used)
% mx = starting mx
% my = starting my
% mz = starting mz
%


outmx = 0.0*tp;
outmy = 0.0*tp;
outmz = 0.0*tp;

outmx(1) = mx;
outmy(1) = my;
outmz(1) = mz;

for counter = 2 : length(tp)-1
    
    m = [outmx(counter-1) outmy(counter-1) outmz(counter-1)];

    delta_t = (tp(counter+1) - tp(counter));
    %first apply RF pulses to rotate the magnetization
    rot_anglex = real(b1(counter))*delta_t*2*pi;
    rot_x = [ 1 0 0 ; 0 cos(rot_anglex) -sin(rot_anglex) ; 0 sin(rot_anglex) cos(rot_anglex)];
    rot_angley = imag(b1(counter))*(tp(counter+1) - tp(counter))*2*pi;
    rot_y = [ cos(rot_angley) 0 sin(rot_angley) ; 0 1 0 ; -sin(rot_angley) 0 cos(rot_angley)];
    m_rot = m * rot_x * rot_y;

    
    %Now apply relaxation

    mxy = sqrt(sum(m_rot(1).^2 + m_rot(2).^2));
    mxy_T2 = mxy * exp(-delta_t/T2);
    mz = 1 - (1 - m_rot(3))*exp(-delta_t/T1);

    %Now rotate due to precession
    outmx(counter) = mxy_T2* sin (atan2(m_rot(1),m_rot(2)) + delta_t*df*2*pi);
    outmy(counter) = mxy_T2* cos (atan2(m_rot(1),m_rot(2)) + delta_t*df*2*pi); %check  sin vs cos here
    outmz(counter) = mz;

end

