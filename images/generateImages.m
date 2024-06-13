%% Plot the smoothed spectra of a unit impulse.
% Lachlan Birnie, 13-June-2024

legibleFactor = 100;  % Scale smoothed spectra for legibility.

bandsPerOctave = 1;

fs = 48000;
impulseFreq = 5000;
nfft = 4096;

impulseBin = ceil(impulseFreq / (fs/nfft)) + 1;

% Make unit impulse spectra.
spectra = zeros(nfft/2 + 1, 1);
spectra(impulseBin, 1) = 1;

% Apply fractional octave smoothing.
[recSpectra, binFreqs] = spectrumOctaveSmoothing(spectra, nfft, fs, bandsPerOctave, 'basic-rectangle');
hannSpectra = spectrumOctaveSmoothing(spectra, nfft, fs, bandsPerOctave, 'basic-hann');
logSpectra = spectrumOctaveSmoothing(spectra, nfft, fs, bandsPerOctave, 'log-compensated');

% Plot results.
hFig = figure('Color', [1, 1, 1]);

hPlot = semilogx(binFreqs, spectra, 'k', 'LineWidth', 2, 'DisplayName', 'Original Spectra');

hold on;

semilogx(binFreqs, recSpectra .* legibleFactor, 'g', 'LineWidth', 2, 'Displayname', 'Rectangle Window');
semilogx(binFreqs, hannSpectra .* legibleFactor, 'c', 'LineWidth', 2, 'Displayname', 'Hann Window');
semilogx(binFreqs, logSpectra .* legibleFactor, 'm', 'LineWidth', 2, 'DisplayName', 'Log-Compensated Window');

legend('show');

xticks([1/4, 1/2, 1, 2, 4] .* impulseFreq);
xlim([min(xticks), max(xticks)]);
xlabel('Frequency (Hz)');

ylim([0, 1]);
ylabel('Magnitude');

title(sprintf('Original and %.1f-octave smoothed spectra for a unit impulse', 1./bandsPerOctave));

hAxis = hPlot.Parent;
set(hAxis, 'FontSize', 14, 'XGrid', 'on');

hFig.Position(3) = hFig.Position(3) * 2;

% saveas(hFig, 'smoothedImpulseSpectra.png');

%% Plot original and smoothed transfer function.
% Lachlan Birnie, 13-June-2024

nfft = 2048;  % Assumed even number.
fs = 48000;

bandsPerOctave = 3;

% Create an example impulse response.
sig = zeros(nfft, 1);
sig([0,1,2] + 10, 1) = [0.4, 1, 0.4];
sig([0,1,2] + 110, 1) = [-0.2, -0.4, -0.2];
sig([0,1,2] + 115, 1) = [0.1, 0.3, 0.1];
sig([0,1,2] + 130, 1) = [-0.05, -0.1, -0.05];

sig = awgn(sig, 20, "measured");  % Add some noise.

% Get the original transfer function.
spectra = 20.*log10(abs(fft(sig)));
spectra = spectra(1:nfft/2 + 1);

% Get the log-compensated fractional octave smoothed spectrum.
[logSpectra, binFreqs] = spectrumOctaveSmoothing(spectra, nfft, fs, bandsPerOctave, 'log-compensated');

% Plot the original and smoothed spectrum.
hFig = figure('Color', [1,1,1]);

hPlot = semilogx(binFreqs, spectra, 'k', 'LineWidth', 2, 'DisplayName', 'Original Spectrum');

hold on;

semilogx(binFreqs, logSpectra, 'm', 'LineWidth', 2, 'DisplayName', 'Log-Compensated Smoothed Spectrum');

legend('Location', 'southwest');

xlim([0, fs/2]);
xlabel('Frequency (Hz)');

ylabel('Magnitude');

title(sprintf('Original and 1/%i-octave smoothed transfer function', bandsPerOctave));

hAxis = hPlot.Parent;
set(hAxis, 'FontSize', 14, 'XGrid', 'on');

hFig.Position(3) = hFig.Position(3) * 2;

% saveas(hFig, 'smoothedTransferFunctionSpectra.png');