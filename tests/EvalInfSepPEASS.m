function [SDR SIR SAR ISR qGlobal qTarget qInterf qArtif OPS TPS IPS APS] = EvalInfSepPEASS(extractSrc, origSrc, options)
%EvalInfSepPEASS Calculate PEASS scores for any separated output
% [SDR SIR SAR ISR qGlobal qTarget qInterf qArtif OPS TPS IPS APS] = EvalInfSepPEASS(extractSrc, origSrc, options)
% 	calculates the PEASS metrics for the extracted sources, given the original sources. 
%	All sources and extracted output must be single-channel and given in column format.
%	Options to the PEASS system can be defined in the options struct.
%
%	Extra: you have the possibility to compute these metrics for every source in parallel. 
%		Do this by replacing every "for" with "parfor" and having a matlabpool open in your session.
%
% input parameters: 
% extractSrc: cell array with in each cell an extracted source, in column format.
% origSrc: cell array with in each cell an original source, in column format (as read by 'wavread')
% options: a struct containing all parameters that need to be kept track of, as defined in runComparison.m (or by yourself)
%	The following fields from this structure are passed to PEASS:
%	* options.segmentationFactor (optional): for low-memory systems or long signals (10 seconds or more) 
%		make this a positive integer. Calculates the metrics in an overlap-add fashion. If nonexisting, put to 0.
%	* options.resultdir (obligatory): the results directory is used as a place for temporary files too. 
% 		In this project, PEASS code has been slightly adapted such that all temporary files are deleted.
%
% output parameters:
% SDR, SIR, SAR, ISR, qGlobal, qTarget, qInterf, qArtif, OPS, TPS, IPS, APS: column arrays which contain the 
%	PEASS metrics for every source. 
% 
% author (and copyright): Joachim Ganseman, 2012 
% some minor bugfixes in preparation for PhD thesis made in 2017
    
    % load PEASS options. Refer to extractDistortionComponents.m for more info.
    % SegFactor: choose this number to define length of segments calculated
    % through overlap-add. Default to 0 for now, but for longer signals, perhaps 
    % put something like "floor(options.datalength / 200000);" in the options struct
    if isfield(options, 'segmentationFactor')
        PEASSopt.segmentationFactor = options.segmentationFactor;
    else
    	PEASSopt.segmentationFactor = 0;
    end
        
    PEASSopt.tmpDir = options.resultdir;
    PEASSopt.destDir = options.resultdir;     %use output dir as temp dir.
    
    % write all sources and extracted files to disk
    [extUnique extFiles srcFiles] = deal( cell(length(extractSrc), 1));
    % to parallellize, replace "for" by "parfor"
    for i = 1:length(extractSrc)
        % write to unique files. Generate filename:
        extUnique{i} = strrep(tempname, tempdir, options.resultdir);
        extFiles{i} = strcat(extUnique{i}, 'ext', num2str(i), '.wav');
        srcFiles{i} = strcat(extUnique{i}, 'src', num2str(i), '.wav');
        wavwrite(extractSrc{i}, options.sr, extFiles{i});
        wavwrite(origSrc{i}, options.sr, srcFiles{i});
    end
    
    % initialize output in column vectors    
    [ SDR SIR SAR, ISR, qGlobal qTarget qInterf qArtif OPS TPS IPS APS ] ...
        = deal( zeros(length(extractSrc), 1));
    
    % calculate for every source the PEASS scores. These are: 
    % 1. traditional SDR, SIR, SAR, ISR
    % 2. PEMO-Q scores: qGlobal,qTarget,qInterf,qArtif
    % 3. Total scores: OPS, TPS, IPS, APS
    % to parallellize, replace "for" by "parfor"
    for i = 1:length(extractSrc)
        
        % reorder the src filenames, so that the source under investigation
        % is first. Consider all sources in sequence: circular shift backwards
        curSrcFiles{i} = circshift(srcFiles, -i+1);
        
        result{i} = PEASS_ObjectiveMeasure(curSrcFiles{i},extFiles{i}, PEASSopt);
        
        SDR(i) = result{i}.SDR;
        SIR(i) = result{i}.SIR;
        SAR(i) = result{i}.SAR;
        ISR(i) = result{i}.ISR;
        qGlobal(i) = result{i}.qGlobal;
        qTarget(i) = result{i}.qTarget;
        qInterf(i) = result{i}.qInterf;
        qArtif(i) = result{i}.qArtif;
        OPS(i) = result{i}.OPS;
        TPS(i) = result{i}.TPS;
        IPS(i) = result{i}.IPS;
        APS(i) = result{i}.APS;               
    end
    
    % delete temporary files.
    % to parallellize, replace "for" by "parfor"
    for i = 1:length(extractSrc)
        delete(srcFiles{i});
        delete(extFiles{i});
        delete(strrep(extFiles{i}, '.wav', '_true.wav'));
        delete(strrep(extFiles{i}, '.wav', '_eArtif.wav'));
        delete(strrep(extFiles{i}, '.wav', '_eInterf.wav'));
        delete(strrep(extFiles{i}, '.wav', '_eTarget.wav'));
    end
