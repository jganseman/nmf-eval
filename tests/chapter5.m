% Chapter 5: creating the images for Chapter 5 of my PhD Thesis.
% (c) Joachim Ganseman, october 2017

% NOTE: run this file while inside the 'tests' folder in Matlab for correct path

% NOTE: AFTER RUNNING, EXTERNALIZE THE DATA TABLE READING in 2 steps:
% 1. add the following code before subimporting the created tex file:
%     \pgfplotstableread{./images/chapter2/audiosignalexcerpt-1.tsv}\loadedtable
% 2. in the tex file itself, change the following:
%     table[row sep=crcr,format=inline]{\loadedtable}; 


clear;
addpath(genpath('../'));     % add subdirectories to path. 

PRINTTOFILE = 0;

%% Part 0: plot the different beta divergences
% use equations in definition, substitute y=1 

deuc = @(x) (x.^2 -2.*x +1);
dkl = @(x) (x.*log(x) -x +1);
dis = @(x) (x -log(x) -1);

figure;
fplot(deuc, [-0.5, 5], 'black');
text(4.1,9,texlabel('d_{EUC}(x|y=1) = x^2 - 2x + 1'));
hold on;
fplot(dkl, [-0.5, 10], 'red');
text(7,7.3,texlabel('d_{KL}(x|y=1) = x log(x) - x + 1'));
hold on;
fplot(dis, [-0.5, 10], 'blue');
text(6,2.8,texlabel('d_{IS}(x|y=1) = x - log(x) - 1'));
hold off;
grid on;

%% Part 1: find a good plot of the plotBetaNMF test (testing optimal beta/power for NMF). 

%load('nmf-beta-test-chord1-base48.mat');
%disp('display bss_eval metrics separating C3-E3-G3 piano chord.')

load('nmf-beta-test-chord1-base60.mat');
disp('display bss_eval metrics separating C4-E4-G4 piano chord.')

%load('nmf-beta-test-chord1-base72.mat');
%disp('display bss_eval metrics separating C5-E5-G5 piano chord.')

% set up as: arpeggiated mix of 3 sources from MAPS dataset
% for each power 0.5:0.1:2.5, for each beta 0.0:0.1:2.5, 3 SDR values [1;2;3]
% try to plot as 3 separate 3D plots

%%
maketable = cell2mat(cat(1,mySDR{:})); % put into tabular format
source1sdr = maketable(1:3:end, :)';
source2sdr = maketable(2:3:end, :)';
source3sdr = maketable(3:3:end, :)';

%only keep values larger than 0
source1sdr = max(0, source1sdr);
source2sdr = max(0, source2sdr);
source3sdr = max(0, source3sdr);

%display max value
%max(source1sdr(:)) 
% keep within the maximum and 5 dB below that?
% source1sdr = max(max(source1sdr(:))-5.0, source1sdr);

%some variables for plot details
nrcontours = 20;

%% plot these as a 3d mesh/contour plot

figure1 = figure;
axes1 = axes('Parent',figure1);
title('BSS_EVAL 3.0 scores, extraction of notes from C4-E4-G4 piano chord');
ylabel('Beta');
xlabel('exponent');
hold(axes1,'on');
box(axes1,'on');
axis(axes1,'tight');
set(axes1,'BoxStyle','full','Layer','top','YMinorTick','on','YTick',...
    [1 6 11 16 21 26],'YTickLabel',{'0.0','0.5','1.0','1.5','2.0','2.5'},...
    'XMinorTick','on','XTick',[1 6 11 16 21],'XTickLabel',...
    {'0.5','1.0','1.5','2.0','2.5'});

myaxes = set(axes1);

%
s1ax = subplot(3, 3, 1);
contourf(source1sdr, nrcontours);
colorbar;
title('C SDR');
grid on;
ylabel('Beta');
% draw line of maxima along x-axis
[~, xmaxes] = max(source1sdr);
hold on; p = plot( 1:1:size(xmaxes,2) , xmaxes, 'red', 'LineWidth', 1.0, 'Marker', '*'); hold off;
% draw line of maxima along y-axis
[~, ymaxes] = max(source1sdr');
hold on; plot( ymaxes , 1:1:size(ymaxes,2), 'green', 'LineWidth', 1.0, 'Marker', 'd'); hold off;

%
s2ax = subplot(3, 3, 2);
contourf(source2sdr, nrcontours);
colorbar;
title('E SDR');
grid on;
% draw line of maxima along x-axis
[~, xmaxes] = max(source2sdr);
hold on; p = plot( 1:1:size(xmaxes,2) , xmaxes, 'red', 'LineWidth', 1.0, 'Marker', '*'); hold off;
% draw line of maxima along y-axis
[~, ymaxes] = max(source2sdr');
hold on; plot( ymaxes , 1:1:size(ymaxes,2), 'green', 'LineWidth', 1.0, 'Marker', 'd'); hold off;

%
s3ax = subplot(3, 3, 3);
contourf(source3sdr, nrcontours);
colorbar;
title('G SDR');
grid on;
% draw line of maxima along x-axis
[~, xmaxes] = max(source3sdr);
hold on; p = plot( 1:1:size(xmaxes,2) , xmaxes, 'red', 'LineWidth', 1.0, 'Marker', '*'); hold off;
% draw line of maxima along y-axis
[~, ymaxes] = max(source3sdr');
hold on; plot( ymaxes , 1:1:size(ymaxes,2), 'green', 'LineWidth', 1.0, 'Marker', 'd'); hold off;

%% plot SIR
maketable = cell2mat(cat(1,mySIR{:})); % put into tabular format
source1sir = maketable(1:3:end, :)';
source2sir = maketable(2:3:end, :)';
source3sir = maketable(3:3:end, :)';

%only keep values larger than 0
source1sir = max(0, source1sir);
source2sir = max(0, source2sir);
source3sir = max(0, source3sir);

%
% figure2 = figure;
% axes1 = axes('Parent',figure1);
% ylabel('Beta');
% xlabel('exponent');
% hold(axes1,'on');
% box(axes1,'on');
% axis(axes1,'tight');
% set(axes1,'BoxStyle','full','Layer','top','YMinorTick','on','YTick',...
%     [1 6 11 16 21 26],'YTickLabel',{'0.0','0.5','1.0','1.5','2.0','2.5'},...
%     'XMinorTick','on','XTick',[1 6 11 16 21],'XTickLabel',...
%     {'0.5','1.0','1.5','2.0','2.5'});
%
s4ax = subplot(3, 3, 4);
contourf(source1sir, nrcontours);
colorbar;
title('C SIR');
grid on;
ylabel('Beta');
% draw line of maxima along x-axis
[~, xmaxes] = max(source1sir);
hold on; p = plot( 1:1:size(xmaxes,2) , xmaxes, 'red', 'LineWidth', 1.0, 'Marker', '*'); hold off;
% draw line of maxima along y-axis
[~, ymaxes] = max(source1sir');
hold on; plot( ymaxes , 1:1:size(ymaxes,2), 'green', 'LineWidth', 1.0, 'Marker', 'd'); hold off;

%
s5ax = subplot(3, 3, 5);
contourf(source2sir, nrcontours);
colorbar;
title('E SIR');
grid on;
% draw line of maxima along x-axis
[~, xmaxes] = max(source2sir);
hold on; p = plot( 1:1:size(xmaxes,2) , xmaxes, 'red', 'LineWidth', 1.0, 'Marker', '*'); hold off;
% draw line of maxima along y-axis
[~, ymaxes] = max(source2sir');
hold on; plot( ymaxes , 1:1:size(ymaxes,2), 'green', 'LineWidth', 1.0, 'Marker', 'd'); hold off;

%
s6ax = subplot(3, 3, 6);
contourf(source3sir, nrcontours);
colorbar;
title('G SIR');
grid on;
% draw line of maxima along x-axis
[~, xmaxes] = max(source3sir);
hold on; p = plot( 1:1:size(xmaxes,2) , xmaxes, 'red', 'LineWidth', 1.0, 'Marker', '*'); hold off;
% draw line of maxima along y-axis
[~, ymaxes] = max(source3sir');
hold on; plot( ymaxes , 1:1:size(ymaxes,2), 'green', 'LineWidth', 1.0, 'Marker', 'd'); hold off;

%% plot SAR
maketable = cell2mat(cat(1,mySAR{:})); % put into tabular format
source1sar = maketable(1:3:end, :)';
source2sar = maketable(2:3:end, :)';
source3sar = maketable(3:3:end, :)';

%only keep values larger than 0
source1sar = max(0, source1sar);
source2sar = max(0, source2sar);
source3sar = max(0, source3sar);

%
% figure3 = figure;
% axes1 = axes('Parent',figure1);
% ylabel('Beta');
% xlabel('exponent');
% hold(axes1,'on');
% box(axes1,'on');
% axis(axes1,'tight');
% set(axes1,'BoxStyle','full','Layer','top','YMinorTick','on','YTick',...
%     [1 6 11 16 21 26],'YTickLabel',{'0.0','0.5','1.0','1.5','2.0','2.5'},...
%     'XMinorTick','on','XTick',[1 6 11 16 21],'XTickLabel',...
%     {'0.5','1.0','1.5','2.0','2.5'});
%
s7ax = subplot(3, 3, 7);
contourf(source1sar, nrcontours);
colorbar;
title('C SAR');
grid on;
ylabel('Beta');
xlabel('exponent');
% draw line of maxima along x-axis
[~, xmaxes] = max(source1sar);
hold on; p = plot( 1:1:size(xmaxes,2) , xmaxes, 'red', 'LineWidth', 1.0, 'Marker', '*'); hold off;
% draw line of maxima along y-axis
[~, ymaxes] = max(source1sar');
hold on; plot( ymaxes , 1:1:size(ymaxes,2), 'green', 'LineWidth', 1.0, 'Marker', 'd'); hold off;

%
s8ax = subplot(3, 3, 8);
contourf(source2sar, nrcontours);
colorbar;
title('E SAR');
grid on;
xlabel('exponent');
% draw line of maxima along x-axis
[~, xmaxes] = max(source2sar);
hold on; p = plot( 1:1:size(xmaxes,2) , xmaxes, 'red', 'LineWidth', 1.0, 'Marker', '*'); hold off;
% draw line of maxima along y-axis
[~, ymaxes] = max(source2sar');
hold on; plot( ymaxes , 1:1:size(ymaxes,2), 'green', 'LineWidth', 1.0, 'Marker', 'd'); hold off;

%
s9ax = subplot(3, 3, 9);
contourf(source3sar, nrcontours);
colorbar;
title('G SAR');
grid on;
xlabel('exponent');
% draw line of maxima along x-axis
[~, xmaxes] = max(source3sar);
hold on; p = plot( 1:1:size(xmaxes,2) , xmaxes, 'red', 'LineWidth', 1.0, 'Marker', '*'); hold off;
% draw line of maxima along y-axis
[~, ymaxes] = max(source3sar');
hold on; plot( ymaxes , 1:1:size(ymaxes,2), 'green', 'LineWidth', 1.0, 'Marker', 'd'); hold off;

%% set axes
set([s1ax, s2ax, s3ax, s4ax, s5ax, s6ax, s7ax, s8ax, s9ax],...
    'BoxStyle','full','Layer','top','YMinorTick','on','YTick',...
    [1 6 11 16 21 26],'YTickLabel',{'0.0','0.5','1.0','1.5','2.0','2.5'},...
    'XMinorTick','on','XTick',[1 6 11 16 21],'XTickLabel',...
    {'0.5','1.0','1.5','2.0','2.5'});

%% if wanted, print to file

%save this plot if you want
if PRINTTOFILE
    filename = '../../../thesis/images/chapter5/betavsexponent.tex';
    matlab2tikz('height', '\figureheight', 'width', '\figurewidth', 'filename', filename, 'relativeDataPath', '.')
end
