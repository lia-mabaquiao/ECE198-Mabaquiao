clear all;
clc;

fileIn = ['/Users/liamabaquiao/Documents/THESIS/DATASET_VOCALS_10_12_SEC/lib-fem1.wav'];
[WaveIn, fs] = audioread(fileIn); %read input signal from file	

global config;
%perfect fifth = 1.5 (3:2)
%perfect fourth = 1.33 (4:3)
%minor third = 1.2 (6:5)
%major third = 1.25 (5:4)
%octave = 2 (2:1)

config.pitchScale           = 1.5;	    %3UP_MAJ
%config.pitchScale           = 0.8;	    %3DOWN_MAJ
%config.pitchScale           = 1.5;	    %5UP
%config.pitchScale           = 0.6667;	%5DOWN

config.timeScale            = 1;	    %time scale ratio 
config.resamplingScale      = 1;		%resampling ratio to do formant shifting
config.reconstruct          = 1;		%if true do low-band spectrum reconstruction
config.displayPitchMarks    = 1;		%if true display pitch mark results
config.maxPitchInput        = 650;      %distortion comes from regular autocorrelation, we just blindly filter it


WaveIn = WaveIn - mean(WaveIn); 	
[PitchContour,pitch] = PitchEstimation(WaveIn,fs);

% Voiced/Unvoiced Detection
[u,v] = UVSplit(PitchContour);
p=PitchContour;
x=WaveIn;

% Pitch Marking for voiced segments
pm = [];
ca = [];
first = 1;
waveOut = []; 
for i = 1 : length(v(:,1))
    range = (v(i, 1) : v(i, 2));
    unvoiced_range =(u(i,1): u(i,2));
    in = x(range);
    [marks, cans] = voiced_segment_ind(in, p(range), fs);

    pm = [pm  (marks + range(1))];
    ca = [ca;  (cans + range(1))];
    
    %ra starts and ends with pitch marks, and it contains the unvoiced
    %samples in between 
    ra = first:marks(1)+range(1)-1;
    
    %printing for error checking
    X = sprintf('first = %d, and end of ra = %d',first,ra(end));
    disp(X)
    first = marks(end)+range(1)+1;
    Y = sprintf('updated first = %d',first);
    disp(Y)
    
    waveOut = [waveOut x(unvoiced_range)];  %Processes the Unvoiced Segments
    waveOut = [waveOut PSOLA(in, fs, marks)];
end

waveOut=waveOut';
out = [waveOut; zeros(length(WaveIn)-length(waveOut),1)];

audiowrite('FOR_ATE.wav',out+WaveIn,fs);
%audiowrite('3DOWN_MAJ_wonderwall_fem2.wav',out,fs);
%audiowrite('5UP_wonderwall_fem2.wav',out,fs);
%audiowrite('5DOWN_wonderwall_fem2.wav',out,fs);