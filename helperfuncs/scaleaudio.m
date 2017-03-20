function [outsources outresidue] = scaleaudio(realsources, estsources, estresidue, opt)
% Scale audio from estimated sources in such a way that it
% corresponds to the scale of the real sources.
% This is necessary for good computation of the metrics.
% Also removes any DC offset

% author: Joachim Ganseman


% initialize output
outsources = cell(1, opt.nrsources);
outresidue = [];

% it's possible there's a DC offset in the separated signals, because of
% having tried to unmix the DC component as well. First, make the DC value
% same as the original source for the sources, and equal to 0 for the
% residue. DC value is mean of all samples
for i=1:opt.nrsources
    meanEst = mean(estsources{i});
    meanOrig = mean(realsources{i});
    estsources{i} = estsources{i} - meanEst + meanOrig;
end

if ~isempty(estresidue)
    estresidue = estresidue - mean(estresidue);
end

% now we can scale the audio appropriately. Apply the same scale factor to
% all separated sources, to have a fair judgment of the previously run
% separation routine. Sum all sources, and compare the mean value of
% absolute samples (that's probably a better indicator of overall gain 
% than the maximum sample value)

% first, calculate the sum of all sources and estimates
sourcesum = sum( cell2mat(realsources), 2);
estimates = [cell2mat(estsources) estresidue];

% ratio of mean absolute values will define scaling
ratio = mean(abs(sourcesum))/mean(abs(sum(estimates,2)));

% sanity check: assert this is strictly positive
assert((ratio > 0), 'ERROR: scaling ratio must be positive');

% apply this scaling
for i=1:opt.nrsources
    outsources{i} = estsources{i} .* ratio;
end

if ~isempty(estresidue)
    outresidue = estresidue .* ratio;
end

