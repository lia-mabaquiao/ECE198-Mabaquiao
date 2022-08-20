classdef Harmonizer198 < audioPlugin
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        PitchShift = 1.5;

    end

    properties (Access = private)
        resamplingScale = 1;
        displayPitchMarks = 1;
        maxPitchInput = 650;
        fs = 44100;
        frameSize = 441;
    end 

    
    properties (Constant)
        PluginInterface = audioPluginInterface(...          %<---
            audioPluginParameter('PitchShift',...           %<---
            'DisplayName','Pitch Shift',...                 %<---
            'Label','factor',...                            %<---
            'Mapping',{'lin',0, 2}))                        %<---
    end
    
%% 
    methods(Access = public)
        function out = PitchEstimation(plugin,input)
            WindowLength = round(plugin.fs * 0.030); %number of sample in 30 ms
            Overlap = round(plugin.fs * 0.010); %number of sample in 10 ms
            WaveLength = length(input);
            NumberOfFrames = floor((WaveLength - WindowLength) / Overlap) + 2;
            if NumberOfFrames < 0
                out = 70*ones(length(input),1); %sets the input as the minimum pitch
                return
            end
            FramePitch = zeros(NumberOfFrames + 2, 1);
            
            % Calculate pitch for each window frame
            Range = 1 : WindowLength;
            for Count = 2 : NumberOfFrames
                    audio_snippet = input(Range);
                    MinLag = round( plugin.fs / 1500);
                    MaxLag = round( plugin.fs / 70);
                    
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
                    MaxIndex = MaxIndex + MinLag - 1;
                    
                    HalfIndex = round(MaxIndex/2);
                    HalfValue = AutoCorr(HalfIndex);
                    
                    [MinValue, ~] = min(AutoCorr(1 : MaxIndex));
                    
                    MeanValue = mean(AutoCorr);
                    
                    
                    if MaxIndex == MinLag || MaxIndex == MaxLag
                        v = false;
                    end
                    if AutoCorr(MaxIndex) < AutoCorr(MaxIndex-1) || AutoCorr(MaxIndex) < AutoCorr(MaxIndex+1)
                        v = false;
                    end
                    v = true;
            
                    if MaxValue > 0.35 && MinValue < 0 && v
                        FramePitch(Count) = plugin.fs / (MaxIndex);
                    else 
                        FramePitch(Count) = 0;
                    end
            
                    Range = Range + Overlap;
            end
            
            
            %Limiting Pitch Input for better filtering
            for i = 2:length(FramePitch)
                if FramePitch(i)>650
                    FramePitch(i)=0;
                else
                    FramePitch(i)=FramePitch(i);
                end
            end
            
            % calculate pitch contour
            out = zeros(WaveLength, 1);
            for i = 1 : WaveLength
                out(i) = FramePitch(floor(i / Overlap) + 1);
            end
        end

%% 
        function [Marks, Candidates] = VoicedSegmentMarking(plugin, x, p)
            counter = 0;
            for t = 2:length(p)
                if p(t) == p(t-1)
                    counter = counter + 1;
                end
            end
            
            if counter+1 == length(p)
                Marks = 1:length(p);
                Candidates = ones(length(p),8);
                return;
            end

            %fs = getSampleRate(plugin);
            MaxCandidateNumber = 3;
            pitch = p;
            [MaxAmp, i] = max(x); % find global maximum of amplitudes in voice segment
            len = length(x);
            
            % first candidate is maximum amplitude sample in voiced segment
            first = ones(1,MaxCandidateNumber + 5);
            first(1, 1:MaxCandidateNumber + 5) = i;
            first(1, 2:MaxCandidateNumber + 1) = 0;
            coder.varsize('RightCandidates');
            coder.varsize('LeftCandidates');
            coder.varsize('Candidates');
            % find marks in right handside
            RightCandidates = plugin.IncreaseMarking(x(i:len), p(i:len), MaxCandidateNumber);
            %find mark in left handside via fliped signal
            LeftCandidates = plugin.IncreaseMarking(flipud(x(1:i)), flipud(p(1:i)), MaxCandidateNumber);
            % combine candidates
            Candidates = [flipud((i + 1) - LeftCandidates); first; (i - 1) + RightCandidates ];
            Candidates( find(Candidates == i + 1 | Candidates == i - 1) ) = 0; %restore zero guards 
            
            % =========================================================================
            % Dynamic programming
            % =========================================================================
            % init
            d = plugin.Pitch2Duration(pitch,plugin.fs);
            
            cost = zeros(len, 1);
            trace = zeros(len, 1);
            len = length(Candidates(:,1));
            curr_val = 0;
            prev_val = 0;
            coder.varsize('curr_val')
            coder.varsize('prev_val')

            imin = Candidates(1, MaxCandidateNumber + 2);
            imax = Candidates(1, MaxCandidateNumber + 3);
            curr_val = plugin.ListCandidates( Candidates(1, 1 : MaxCandidateNumber));
            for w = 1:length(curr_val)
                curr = curr_val(w);
                cost(curr) = log( plugin.StateProb(x(curr), x(imin), x(imax)) );
                trace(curr) = 0;
            end
            % loop
            for k = 2 : len                     %we start at 2 kasi yung 1, yun na yung pitch mark. 
                imin = Candidates(k, MaxCandidateNumber + 2);
                imax = Candidates(k, MaxCandidateNumber + 3);
                curr_val = plugin.ListCandidates(Candidates(k, 1 : MaxCandidateNumber));
                for i=1:length(curr_val)
                    curr = curr_val(i);
                    if trace(curr) ~= 0
                        disp('[error] overlap search region');
                        break;
                    end
                    MaxProb = -999999;
                    prev_val = plugin.ListCandidates(Candidates(k - 1, 1 : MaxCandidateNumber));
                    for q = 1:length(prev_val)
                        prev = prev_val(q);
                        Prob = log(plugin.TransitionProb(prev, curr, d) ) + cost(prev);
                        if Prob > MaxProb
                            MaxProb = Prob;
                            trace(curr) = prev;
                        end
                    end % prev
                    cost(curr) = MaxProb + log( plugin.StateProb(x(curr), x(imin), x(imax)) );
                end % curr
            end % k
            
            % result
            Marks = zeros(1, len);
            last = plugin.ListCandidates(Candidates(len, 1 : MaxCandidateNumber));
            [~, index] = max( cost(last) );
            curr = last(index);
            prev = trace(curr);
            while (prev ~= 0)
                Marks(len) = curr;
                len = len - 1;
                curr = prev;
                prev = trace(curr);
            end
            Marks(len) = curr;
            
            if len ~= 1 
                disp('[error] do not find all pitch marks');
            end
            
            return
         end
%% 
function out = PSOLA(plugin, input, anaPm)
        len = length(anaPm);
        if len <= 3
            out = input';
            return;
        end
        
        pitchPeriod = zeros(length(input),1);
        for i=2:len
            pitchPeriod(anaPm(i-1):anaPm(i)-1) = anaPm(i)-anaPm(i-1);
        end
        pitchPeriod(1:anaPm(1)-1) = pitchPeriod(anaPm(1));
        pitchPeriod(anaPm(len):length(input)) = pitchPeriod(anaPm(len) - 1);
        
        % Find synPm of synthetic pitch marks
        synPm = zeros(1,5);
        coder.varsize('synPm');
        synPm(1) = anaPm(1);
        count = 1;
        index = synPm(count);
        
        len = length(input);
        while index < len
            LHS = 0;
            RHS = 1;
        
            while (LHS < RHS) && (index < len)
                index = index + 1;
                LHS = (index - synPm(count))^2;        
                RHS =  sum(pitchPeriod(synPm(count):index)) / plugin.PitchShift;
            end
            
            if LHS > RHS
                count = count + 1;
                synPm(count) = index;
                index = synPm(count) + 1; 
            end
        end
        
        % UnitWaveforms
        input = input';
        unit_wave = zeros(1,1000);
        len_wave=zeros(1,length(anaPm));
        coder.varsize('unit.wave');
        for i = 2 : length(anaPm) - 1
            left = anaPm(i-1);
            right = anaPm(i+1);
            unit_wave(i,:)  = input(left : right) .* hann(right - left + 1)';
            len_wave(i) = right-left+1;
        end
        
        % Output signal
        outPm = round(synPm);
        out = zeros(1, outPm(count));
        
        maxLen = 0;
        minLen = length(out);
        
        first = 1;
        for j=2:length(synPm)-1
            for i=first:length(anaPm) - 2
                if (synPm(j) >= anaPm(i) && synPm(j) < anaPm(i+1)) 
        
                    first = i;
                    k = plugin.selectCorrectPos(i,anaPm);
                    
                    gamma = (synPm(j) - anaPm(k)) / (anaPm(k+1) - anaPm(k));            
        
                    unitWave1 = unit_wave(k,1:len_wave(i));
                    unitWave2 = unit_wave(k+1,1:len_wave(i+1));

        
                    newUnitWave = plugin.addTogether((1 - gamma)*unitWave1, gamma*unitWave2, 1, maxLen);
                    
                    
                    if (plugin.resamplingScale ~= 1)
                        newUnitWave = plugin.resampling(newUnitWave, plugin.resamplingScale);
                    end                    
                    
                    out = plugin.addTogether(out, newUnitWave, outPm(j-1),maxLen);
                    
                    if outPm(j-1) < minLen
                        minLen = outPm(j-1);
                    end
                    break;
                end                     
            end
        end
        out = out(minLen:outPm(count));
    end

        function ca = IncreaseMarking(plugin, x, p, m)
                leftDuration = round( plugin.fs / p(1)); 
                i = 1 + leftDuration;
                len = length(x);
                Count = 0;
                
                % for each search regions
                % find m pitch mark candidates
                ca = ones(1,8);
                coder.varsize('ca')
                LeftThr = 0;
                while (i < len)
                    rightDuration = round( plugin.fs / p(i));          %these three lines are just allowing you to set a range 
                    leftHalf = floor(leftDuration * 0.3);       %in looking for the pitch marks. so we are -0.3 < x < 0.3 
                    rightHalf = floor(rightDuration * 0.3);
                    Range = (max(i - leftHalf, LeftThr) : min(i + leftHalf, len)); % --------[max --- first --- min]-----
                
                    Offset = Range(1) - 1;
                    c = plugin.FindPeakCandidates( x(Range), m , Offset);
                
                    Count = Count + 1;
                    ca(Count, :) = c;
                
                    i = c(1); %position of max amplitude in current search region
                    leftDuration = round( plugin.fs / p(i));
                    i = i + leftDuration;
                end
        end

        function [u,v] = UVSplit(plugin, p)
            counter = 0;
            for x = 2:length(p)
                if p(x) == p(x-1)
                    counter = counter + 1;
                end
            end
            
            if counter+1 == length(p)
                u = 0;
                v = [1,length(p)];
                return;
            end
            

            max = length(p);
            last = 1;
            count = 1;
            first = 0;
            for i = 1:length(p)
                if p(i) == 0
                    first = first + 1;
                else 
                    break
                end
            end
            
            first = first + last;
            u = zeros(1,2);
            v = zeros(1,2);

            terminator = 0;
            while (terminator ~= 1)
                u(count, 1) = last;
                u(count, 2) = first - 1;
                
                %last = find(p(first : max) == 0, 1, 'first') + first - 1;
                finder = p(first:max);
                counter_last = 0;
                for i = 1:length(finder)
                if finder(i) == 0
                    counter_last = counter_last + 1;
                else 
                    last = counter_last + first - 1;
                    break
                end
                end

                v(count, 1) = first;
                v(count, 2) = last - 1;
            
                %first = find(p(last + 1 : max), 1, 'first') + last;
                finder = p(last + 1:max);
                counter_first = 0;
                for i = 1:length(finder)
                if finder(i) == 0
                    counter_first = counter_first + 1;
                else 
                    break
                end
                end

                if counter_first == length(finder)
                    terminator = 1;
                end

                first = counter_first + last;
                count = count + 1;
            end
            
            u(count, 1) = v(count - 1, 2) + 1;
            u(count, 2) = max;
        end
    end


    
    methods(Static)
        
        %To compute for state probabilities used in pitch marking decision
        %making
          function sc = StateProb(h, min, max)
            % State Probability
                alpha = 1;
            %     disp(sprintf('%d %d %d', h, min, max));
                if min == max %the first pitch mark
                    sc = 1;
                    return;
                end
                sc = ((h - min) / (max - min))^alpha;
            return
            end
            
            % =========================================================================
            function tc = TransitionProb(i, k, d)
            % Transition Probability
                beta = 0.7;
                gamma = 0.6;
                dur = (d(i) + d(k)) / 2; %pitch period
                tc = (1 / (1 - beta * abs(dur - abs(k - i) ) ) )^gamma; 

            end
            
            % =========================================================================
            function d = Pitch2Duration(input_pitch,fs)
            %calculate T0
                d = input_pitch;
                for i = 1:length(input_pitch)
                    if input_pitch(i)~=0 
                        d(i) = fs / input_pitch(i);
                    end
                end
            end
            
            % =========================================================================
            function list = ListCandidates(c)
                list = c( find(c));
            end

            function PeakCandidates = FindPeakCandidates(x, MaxCandidateNumber, Offset)
                    len = length(x);
                    x1 = circshift(x, [1 1]);           %bottom element goes up, the rest go down
                    x2 = circshift(x, [-1 -1]);         %top element goes down, the rest go up
                    PeakIndices = find(x >= x1 & x >= x2); %find peaks
                    % PeakIndices = PeakIndices( find(PeakIndices ~= 1 & PeakIndices ~= len)); % remove two bound peak candidates
                    [~, SortedIndices] = sort( x(PeakIndices), 'descend');  %sort peaks in descending amplitude
                    
                    MinDur = round(len / 7);
                    l = length(SortedIndices);
                    i = 1;
                    % iterative select and remove
                    while (i < l)
                        j = i + 1;
                        % remove peaks located in MinDur range of current selected peak
                        while (j <= l)
                            if abs(PeakIndices(SortedIndices(i)) - PeakIndices(SortedIndices(j))) < MinDur
                                SortedIndices(j) = [];
                                l = l - 1;
                            else
                                j = j + 1;
                            end
                        end
                        i = i + 1;
                    end
                    
                    % basic information
                    PeakCandidates = zeros(1, MaxCandidateNumber);
                    range = 1 : min( MaxCandidateNumber, length(SortedIndices) );
                    PeakCandidates( range ) = PeakIndices( SortedIndices(range) ) + Offset;
                    
                    % added information
                    [~, imin] = min(x);
                    [~, imax] = max(x);
                    PeakCandidates(1, MaxCandidateNumber + (1:5)) = [- Offset; imin; imax; 1; len] + Offset;

            end

            function i = selectCorrectPos(i,anaPm)
                if i == 1
                    i = 2;
                elseif i >= length(anaPm) - 1
                    i = length(anaPm) - 2;
                end
            end   

            function y = addTogether(y, x, t, maxLen)
                len = length(x);
                max = t + len - 1;
                range = t : max; 
        
                if max > maxLen
                    maxLen = max;
                end
        
                len = length(y);
                if len < max
                    y = [y zeros(1, max - len)];
                end
        
                y(range) = y(range) + x;
            end

        %%
            function y = resampling(x, scale)
                m = length(x);     
                range = 1: scale : m; 
                y = interp1(x, range);
            end         
       end

% ||This is the main method|| function invocations
    methods
        function out = process(plugin,in)
            sz = size(in);
            %mono = in;
            if sz(1) > 1
                mono = (in(:,1)+in(:,2))/2;
                in = mono - mean(mono); 
            else
               mono = 0;
               in = in - mean(in);
            end
        	
            PitchContour = plugin.PitchEstimation(in);
            
            %if PitchContour ~= 0
            [u,v] = plugin.UVSplit(PitchContour);
            
            p=PitchContour;
            x=in;
            
            % Pitch Marking for voiced segments
            pm = zeros(1,3);
            ca = zeros(1,8);
            coder.varsize('pm')
            coder.varsize('ca')
            first = 1;
            waveOut = []; 
            coder.varsize('in_x');
            for i = 1 : length(v(:,1))
                range = (v(i, 1) : v(i, 2));
                if u == 0
                    unvoiced_range = 0;
                else
                unvoiced_range =(u(i,1): u(i,2));
                end
                in_x = x(range);
                p_in = p(range);
                [marks, cans] = plugin.VoicedSegmentMarking(in_x, p_in);
            
                pm = [pm  (marks + range(1))];
                ca = [ca;  (cans + range(1))];
                
                %ra starts and ends with pitch marks, and it contains the unvoiced
                %samples in between 
                ra = first:marks(1)+range(1)-1;
                
                coder.varsize('waveOut')
                if unvoiced_range == 0
                    out = in;
                    out = [out out];
                    return
                else
                waveOut = [waveOut x(unvoiced_range)'];  %Processes the Unvoiced Segments
                waveOut = [waveOut plugin.PSOLA(in, marks)];
                out1 = [waveOut; zeros(length(in)-length(waveOut),1)];
                out = out1 + waveOut;
           % else
                end
            end
            out = [out out];
        %VALend
        end


        function reset(plugin)
            plugin.fs = getSampleRate(plugin);
        end
    end


end