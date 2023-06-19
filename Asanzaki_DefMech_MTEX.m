%% calculate grains 
% calculate grains including not indexed with 10deg threshold criteria
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd, 'threshold',10*degree );

% remove grains < 3*step size (1 um)
grains = grains(grains.grainSize>=3);
ebsd = ebsd(grains); 

% repeat reconstruction
[grains,ebsd.grainId,ebsd.mis2mean] = calcGrains(ebsd, 'threshold',10*degree );

% smooth grain boundaries
grains = smooth(grains,3);

figure
plot(grains);

%% calcite grain size analysis 
% subset grains dataset to just calcite
calcite = grains('Calcite');

% remove calcite grains with size (no. of pixels) < 3x step size (1um)
calcite = calcite(calcite.grainSize > 3);

% create a variable of grain diameters to perform calculations on
grain_size = calcite.diameter;

% plot all calcite grains and colour by diameter
figure
plot(calcite,grain_size)
mtexColorbar('title','Grain size (um)')

%% calculate mean diameter of calcite grains
mean_size = mean(calcite.diameter);

% plot a histogram of calcite grain diameter that is normalised (frequency%
% rather than count
figure
hist = histogram(calcite.diameter, 'Normalization','probability');

grid on
hold on

% add a dotted line showing where mean diameter is
xline(mean_size, ':r', round(mean_size, 1), 'linewidth', 2, ...
    'LabelHorizontalAlignment', 'center', 'LabelOrientation', 'horizontal');
xlabel('Grain size (um)')
ylabel('Freq (%)')
hold off

%% Grain average orientation
% calculate orientation using all calcite points 
ori = ebsd(grains('Calcite')).orientations;

plot(grains('Calcite'), grains('Calcite').meanOrientation)

%% Misorientation to grain mean orientation (internal misorientation)
mis2meanG = calcGROD(ebsd, grains);

figure
plot(ebsd('Calcite'), ebsd('Calcite').mis2mean.angle./degree);
mtexColorbar
hold on
plot(grains.boundary, 'edgecolor', 'k', 'linewidth', .5)
hold off

%% Grain orientation spread (GOS)
GOS = ebsd.grainMean(mis2meanG.angle);

figure
plot(grains, GOS ./ degree)
mtexColorbar('title','GOS in degree')

% histogram of GOS 
figure
histogram(ebsd.grainMean(mis2meanG.angle)./degree, 'Normalization','probability')
xlabel('Grain orientation spread (deg)')
ylabel('Freq (%)')
grid on

%% calculate ODF
% crystal symmetry for calcite
cs = ebsd('Calcite').CS;

% Miller indices of calcite planes commonly activated at T < 400 C
h = [Miller(0,0,1,'hkl',cs),Miller(1,0,8,'hkl',cs),Miller(1,0,4,'hkl',cs),...
    Miller(0,1,2,'hkl',cs)];

figure
plotPDF(ebsd('Calcite').orientations,h,'antipodal')

%% orientation distribution
% density distribution of all calcite orientations (CPO)
odf = calcDensity(ebsd('Calcite').orientations);

figure
plotSection(odf,'contourf')
mtexColorbar

%% plot odf - MTEX mislabels Z as Y!!

setMTEXpref('xAxisDirection','east');

plotPDF(odf,h,'antipodal', 'lower', 'silent')
setColorRange('equal') % use common colour bar across all stereonets
mtexColorbar('title', 'm.u.d.')
annotation('textbox',[0 .7 0 .3],'String',...
    sprintf('n = %i',length(ebsd('Calcite').orientations)),'FitBoxToText','on');

%% Grain average misorientation (GAM/KAM)
gam = ebsd.grainMean(ebsd.KAM);

figure
plot(grains,gam./degree)
mtexColorbar('title','GAM in degree')
setColorRange([0,1.5])

% calculate m-index (measure of misorientation strength)
MI = calcMIndex(odf);

% histogram of GAM 
figure
histogram(ebsd.grainMean(ebsd.KAM)./degree, 'Normalization','probability')
xlabel('Grain average misorientation (deg)')
ylabel('Freq (%)')
annotation('textbox',[.73 .6 0 .3],'String',sprintf('m = %0.4f',MI),'FitBoxToText','on');
grid on

%% Calculate CPO using one average orientation from each grain rather than all measured orientations 

figure
plotPDF(grains('Calcite').meanOrientation, h,'contourf', 'lower', 'smooth', 'points', 'all')
setColorRange('equal')
mtexColorbar
