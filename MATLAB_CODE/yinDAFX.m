function [PitchContour, pitch] = yinDAFX(x,fs)
% function pitch = yinDAFX(x,fs,f0min,hop)
% Author: Adrian v.d. Knesebeck
% determines the pitches of the input signal x at a given hop size %
% input:
% x         input signal
% fs        sampling frequency
% f0min     minimum detectable pitch
% hop       hop size
%
% output:
% pitch     pitch frequencies in Hz at the given hop size


global config;
% initialization
f0min = 70;
OverlapSize = round(0.01 *fs);
WindowSize = round(0.03 *fs);
input_len=length(x);
NumberOfFrames = floor((input_len - WindowSize) / OverlapSize) + 4;

yinTolerance = 0.22;
taumax = round(1/f0min*fs);
yinLen = 1024;
k = 0;

% frame processing
for i = 1:OverlapSize:(length(x)-(yinLen+taumax))
    k=k+1;                              %this is just a counter, the total number of pitches detected
    xframe = x(i:i+(yinLen+taumax));
    yinTemp = zeros(1,taumax);
    % calculate the square differences, ITO YUNG MATAGAL EH
        for tau=1:taumax
            for j=1:yinLen
                yinTemp(tau) = yinTemp(tau) + (xframe(j) - xframe(j+tau))^2; 
            end
        end
    % calculate cumulated normalization 
        tmp = 0;
        yinTemp(1) = 1;
        for tau=2:taumax
            tmp = tmp + yinTemp(tau);
            yinTemp(tau) = yinTemp(tau) *(tau/tmp); 
        end
    
    % determine lowest pitch
      tau=1;
      while(tau<taumax)
          if(yinTemp(tau) < yinTolerance)
            % search turning point
            while (yinTemp(tau+1) < yinTemp(tau))
                  tau = tau+1;
            end
            pitch(k) = fs/tau;
            break 
          else
                tau = tau+1;
          end
    % if no pitch detected
      pitch(k) = 0;
      end
end


pitch = medfilt1(pitch, 5);
pitch_len = length(pitch);
pitch = [pitch zeros(1,NumberOfFrames-pitch_len)];
pitch = pitch';

%Limiting Pitch Input for better filtering
for i = 2:length(pitch)
    if pitch(i)>config.maxPitchInput 
        pitch(i)=0;
    else
        pitch(i)=pitch(i);
    end
end


% calculate pitch contour
PitchContour = zeros(input_len, 1);
for i = 1 : input_len
    PitchContour(i) = pitch(floor(i / OverlapSize) + 1);
end
end