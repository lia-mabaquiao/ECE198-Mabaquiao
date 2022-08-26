%--------------------------------------------------------------------------
% pitch estimation using autocorrelation method
%--------------------------------------------------------------------------
function [PitchContour,FramePitch] = PitchEstimation(x, fs)
global config;

% Init parameters
WindowLength = round(fs * 0.030); %number of sample in 30 ms
Overlap = round(fs * 0.010);      %number of sample in 10 ms
WaveLength = length(x);
NumberOfFrames = floor((WaveLength - WindowLength) / Overlap) + 2;
FramePitch = zeros(NumberOfFrames + 2, 1);

% Calculate pitch for each window frame
Range = 1 : WindowLength;
for Count = 2 : NumberOfFrames
    %%FramePitch(Count) = PitchDetection(x(Range), fs);
        audio_snippet = x(Range);
        fs = 44100;
        MinLag = round( fs / 1500);
        MaxLag = round( fs / 70);
        
        %Center Clipping
        MaxAmplitude = max( abs(audio_snippet) );
        ClipLevel = MaxAmplitude * 0.3;
        PositiveSet = find( audio_snippet > ClipLevel);
        NegativeSet = find (audio_snippet < -ClipLevel);
        cc = zeros( size(audio_snippet) );
        cc(PositiveSet) = audio_snippet(PositiveSet) - ClipLevel;
        cc(NegativeSet) = audio_snippet(NegativeSet) + ClipLevel;
       
        AutoCorr = xcorr(cc, MaxLag, 'coeff'); 			% normalized ACF (AutoCorrelation Function)
        AutoCorr = AutoCorr(MaxLag + 1 : 2*MaxLag); %take half of ACF
        
        
        [MaxValue, MaxIndex] = max(AutoCorr(MinLag : MaxLag)); %search max value of ACF in search region
        MaxIndex = MaxIndex + MinLag - 1; %this is just to get the real value of the index
        
        HalfIndex = round(MaxIndex/2);
        HalfValue = AutoCorr(HalfIndex);
        
        [MinValue, MinIndex] = min(AutoCorr(1 : MaxIndex));
        
        MeanValue = mean(AutoCorr);
        
        
        if MaxIndex == MinLag || MaxIndex == MaxLag
            v = false; %checks that it is not equal to the lags
        end
        if AutoCorr(MaxIndex) < AutoCorr(MaxIndex-1) || AutoCorr(MaxIndex) < AutoCorr(MaxIndex+1)
            v = false; %checks if it is really the maximum
        end
        v = true; %true means it really is the maximum

        if MaxValue > 0.35 && MinValue < 0 && v
            FramePitch(Count) = fs / (MaxIndex);
            disp(MaxIndex)
        else 
            FramePitch(Count) = 0;
        end

        Range = Range + Overlap;
end

% Using median filter for Post-processing
FramePitch1 = medfilt1(FramePitch, 5);

FramePitch2 = zeros(size(FramePitch1));
%Limiting Pitch Input for better filtering
for i = 2:length(FramePitch1)
    if FramePitch1(i)>650 %650 is an arbitrary number, set as the max pitch
        FramePitch2(i)=0;
    else
        FramePitch2(i)=FramePitch1(i);
    end
end

% calculate pitch contour
PitchContour = zeros(WaveLength, 1);
for i = 1 : WaveLength
    PitchContour(i) = FramePitch2(floor(i / Overlap) + 1);
end

end