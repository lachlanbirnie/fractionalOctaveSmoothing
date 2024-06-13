function [smoothSpectra, binFreqs] = spectrumOctaveSmoothing(spectra, nfft, fs, bandsPerOctave, windowMethod)
% SPECTRUMOCATVESMOOTHING - 1/b fractional octave smoothing for spectra.
% An algorithm to smooth the FFT of an impulse response (or other spectra)
% by averaging over octaves / fractional octaves.
%
% The method follows:
% 
% For each NFFT bin (denoted k) find the smoothed spectrum value.
%
%   Set the frequency of the k-th bin as the octave centre frequency.
%
%   Find the octave's lower and upper edge frequencies using the
%   formula from MATLAB poctave function.
%   (https://au.mathworks.com/help/signal/ref/poctave.html)
%
%   Account for the miss-match between the actual ocatve edge frequency
%   and the frequency value of the nearest NFFT bin. By assigning a 
%   weight to the edge bins for the amount of percentages overlap between
%   the octave and bin edges. Following the description in the MATLAB 
%   poctave function.
%   (https://au.mathworks.com/help/signal/ref/poctave.html)
%
%   Create a window over the octave for obtaining a weighted average.
%   Three options:
%
%       'basic-rectangle' - Rectangle window, apply equal weight to each
%       bin in the octave. This has a bias to the smoothed value where both
%       low and high frequencies are given equal weight even though there
%       are more high frequencies in an octave than low.
%       Similar approach to: AARAE toolbox 'octavesmoothing.m'
%       (https://github.com/densilcabrera/aarae)
%
%       'basic-hann' - Apply a Hann window over the octave. Similar to the
%       rectangle window. The middle of the hann window will not match the
%       center frequency of the octave, so this is probably not a good
%       approach.
%
%       'log-compensated' - Following the approach of Method 3 by J. G.
%       Tylka in: 
%       A Generalized Method for Fractional-Octave Smoothing of
%       Transfer Functions that Perserves Log-Frequency Symmery, AES, 2017.
%       (https://doi.org/10.17743/jaes.2016.0053)
%       This method creates a window over the octave that has symmetric
%       weight to the low and high frequencies in a log-scale.
%
%   Apply the window to the octave and sum to get the smoothed value.
% 
% Syntax:  [smoothSpectra] = spectrumOctaveSmoothing(spec, nfft, fs, bandsPerOctave, windowMethod)
%
% Inputs:
%   spectra - [bins, channels] frequency responses / spectra data.
%   nfft - NFFT size of the frequency response / spectrum.
%   fs - Sampling frequency of the frequency response / spectrum.
%   bandsPerOctave - Number of bands per octave for 1/b-octave smoothing.
%   windowMethod - 'basic-rectrangle' or 'basic-hann' or 'log-compensated'
%
% Outputs:
%   smoothSpectra - [bins, channels] smoothed spectra single sided.
%   binFreqs - Frequency values of the smoothed spectrum bins.
%
% Example: 
% b = 1;  % Bands per octave for smoothing.
% fs = 48000;
% nfft = 4096;
% spec = zeros(2049, 1);
% spec(428,1) = 1;
% recSpec = spectrumOctaveSmoothing(spec,nfft,fs,b,'basic-rectangle');
% [hannSpec, x] = spectrumOctaveSmoothing(spec,nfft,fs,b,'basic-hann');
% logSpec = spectrumOctaveSmoothing(spec,nfft,fs,b,'log-compensated');
% figure('Color',[1 1 1]);
% semilogx(x,spec,'k','DisplayName','Unsmooth spectra');
% hold on;
% semilogx(x,recSpec.*100,'g','Displayname','Rectangle window');
% semilogx(x,hannSpec.*100,'r','Displayname','Hann window');
% semilogx(x,logSpec.*100,'b','DisplayName','Log-Compensated');
% legend('show'); xlabel('Hz'); xlim([0, fs/2]);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: images/generateImages.m
%
% Author: Lachlan Birnie
% Audio & Acoustic Signal Processing Group - Australian National University
% Email: Lachlan.Birnie@anu.edu.au
% Website: https://github.com/lachlanbirnie
% Creation: May 2024
% Last revision: 13-June-2024
%
% References:
%
%   https://au.mathworks.com/matlabcentral/fileexchange/55161-1-n-octave-smoothing
%   https://dsp.stackexchange.com/questions/16635/1-n-octave-complex-smoothing
%   https://doi.org/10.17743/jaes.2016.0053
%   https://au.mathworks.com/help/signal/ref/poctave.html
%


% Default to 1-octave smoothing.
if nargin < 4
    bandsPerOctave = 1;
end

% Default to log-compensated smoothing window.
if nargin < 5
    windowMethod = 'log-compensated';
end

% Check if spectrum is full band otherwise make it half band.
if size(spectra, 1) == nfft
    isFullBandSpectrum = true;
else
    isFullBandSpectrum = false;
    
    % Make spectra half band.
    isNfftEven = ~mod(nfft, 2);
    if isNfftEven
        spectra = spectra(1 : nfft/2 + 1, :);
    else
        spectra = spectra(1 : ceil(nfft/2), :);
    end
end

% Find frequency value of each spectrum bin.
binFreqs = linspace(0, fs, nfft+1).';
binFreqs = binFreqs(1 : end-1);  % Full band frequencies.

if ~isFullBandSpectrum
    % Get half band frequencies.
    if isNfftEven
        binFreqs = binFreqs(1 : nfft/2 + 1);
    else
        binFreqs = binFreqs(1 : ceil(nfft/2));
    end
end

binBandwidth = binFreqs(2);  % Bandwidth of nfft bins.
nBins = numel(binFreqs);  % Total number of nfft bins.

% Smooth each sample in the spectrum.
smoothSpectra = zeros(size(spectra));

% Exit if no smoothing applied.
if isempty(bandsPerOctave)
    smoothSpectra = spectra;
    return;
end

% iBin denotes the index of the octave's centre bin, to be smoothed.
for iBin = (1 : nBins)

    % Find edge frequencies for the fractional octave given the centre.
    octCentreFreq = binFreqs(iBin);
    octLowerFreq = octCentreFreq * 10^(0.3)^(-1 / (2 * bandsPerOctave));
    octUpperFreq = octCentreFreq * 10^(0.3)^(1 / (2 * bandsPerOctave));

    % Truncate the octave band to spectrum bandwidth.
    if octUpperFreq > binFreqs(nBins)
        octUpperFreq = binFreqs(nBins);
    end
    
    % Find bins nearest to the octave lower and upper edges.
    [~, iLowerBin] = min(abs(binFreqs - octLowerFreq));
    [~, iUpperBin] = min(abs(binFreqs - octUpperFreq));

    % Make sure lower bin frequency is below octave's lower frequency.
    if binFreqs(iLowerBin) > octLowerFreq
        iLowerBin = iLowerBin - 1;
    end

    % Make sure upper bin frequency is above octave's upper frequency.
    if binFreqs(iUpperBin) < octUpperFreq
        iUpperBin = min(iUpperBin + 1, nBins);
    end

    % Get total number of bins covering the octave.
    nOctaveBins = iUpperBin - iLowerBin + 1;

    % - Account for octave edge falling in the middle of a bin - 

    % Following MATLAB octave smoothing description.
    % (https://au.mathworks.com/help/signal/ref/poctave.html)
    % Assign percentage weights to the edge bins for how much overlap there
    % is between the bin's bandwidth and the octave.

    % Percent of the lowest bin that contains the octave band.
    lowerEdgeWeight = 1 - (octLowerFreq - binFreqs(iLowerBin)) / binBandwidth;

    % Percent of the upper bin that contains the octave band.
    upperEdgeWeight = 1 - (binFreqs(iUpperBin) - octUpperFreq) / binBandwidth;
    
    % - Make Window Over Octave - 

    % -- Basic Methods - window over the fractional octave --
    
    if strcmp(windowMethod, 'basic-rectangle')
        octWindow = ones(nOctaveBins, 1) ./ nOctaveBins;
    end

    if strcmp(windowMethod, 'basic-hann')
        octWindow = hann(nOctaveBins) ./ ((nOctaveBins-1) / 2);
    end

    % -- AES Paper Method 3 - Logarithmically-compensated window --
    % (https://doi.org/10.17743/jaes.2016.0053)

    if strcmp(windowMethod, 'log-compensated')
        % Need to find the value of the window for every k' centered
        % around the k-th bin by integrating a rectangle window over the 
        % limits of phi_lower and phi_upper.
        kDashes = (iLowerBin : 1 : iUpperBin).';
        
        % --- Old version that shows the process ---
        % wind = zeros(numel(kDashes), 1);
        % for iKDashes = 1 : numel(kDashes)
        %     kd = kDashes(iKDashes);  % Solve the window's value for k'. 
        % 
        %     % Find the limits of integration for k'.
        %     phiLower = log2((kd - 0.5) / k);  % Lower limit Eq.17.
        %     phiUpper = log2((kd + 0.5) / k);  % Upper limit Eq.17.
        % 
        %     % Truncate to only include Phi's that are inside the octave.
        %     % From Eq.15 Wr = 0 if abs(phi) > 1/(2b).
        %     phiLower = max(phiLower, -1/(2*bandsPerOctave));
        %     phiUpper = min(phiUpper, 1/(2*bandsPerOctave));
        % 
        %     % Width is how much of the window overlaps the octave.
        %     % If phu_u - phi_l is negative then there is no overlap.
        %     width = max(phiUpper - phiLower, 0);
        % 
        %     % Height is b for a rectangle window. (Eq.15 in paper).
        %     height = bandsPerOctave;
        % 
        %     % Approx integral of rectangle function Wr as it's area.
        %     W_k_kdash = width * height;  % Eq.16.
        %     wind(iKDashes) = W_k_kdash;
        % end
        
        % --- Faster version to solve all k' at the same time ---
        phiEdges = log2(([kDashes(1)-1; kDashes] + 0.5) / iBin);
        phiLowerEdges = max(phiEdges(1:nOctaveBins), -1/(2*bandsPerOctave));
        phiUpperEdges = min(phiEdges(2:nOctaveBins+1), 1/(2*bandsPerOctave));
        integralWidth = max(phiUpperEdges - phiLowerEdges, 0);
        octWindow = integralWidth .* bandsPerOctave;

    end

    % Apply edge weights to octave window.
    octWindow(1) = octWindow(1) * lowerEdgeWeight;
    octWindow(nOctaveBins) = octWindow(nOctaveBins) * upperEdgeWeight;
    
    % Scale window to keep unit area.
    octWindow = octWindow ./ sum(octWindow);

    % Apply window to the spectrum's octave and sum.
    smoothSpectra(iBin, :) = sum(spectra(iLowerBin:iUpperBin, :) .* octWindow, 1);

end

% Plot spectra and smoothed spectra if no output.
if nargout < 1
    figure('Color', [1 1 1]);
    lineColors = [[0, 0, 0]; [1, 0.6, 0.3]; [0.3, 1, 0.6]; [0.6, 0.3, 1]];
    for iSpec = (1 : size(spectra, 2))
        S1 = semilogx(binFreqs, spectra(:,iSpec), ...
                      '-', ...
                      'DisplayName', 'Original spectra', ...
                      'HandleVisibility', 'off');
        hold on;
        S2 = semilogx(binFreqs, smoothSpectra(:,iSpec), ...
                      '-', ...
                      'LineWidth', 2.5, ...
                      'DisplayName', sprintf('Smoothed Spectra %i', iSpec));

        if iSpec <= size(lineColors, 1)
            set(S1, 'Color', [lineColors(iSpec,:), 0.3]);
            set(S2, 'Color', [lineColors(iSpec,:), 1]);
        else
            col = get(S1, 'Color');
            set(S1, 'Color', [col, 0.3]);
            set(S2, 'Color', [col, 1]);
        end
    end
    if isFullBandSpectrum
        xlim([0, fs]);
    else
        xlim([0, fs/2]);
    end
    xlabel('Freqency (Hz)');
    ylabel('Magnitude');   
    legend('show', 'Location', 'northwest');
    title(sprintf('1/%3.1f Fractional Octave Smoothed Spectra', bandsPerOctave));
    set(gca, 'FontSize', 14, 'XGrid', 'on');
end

end  % spectrumOctaveSmoothing()