%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% General script for data alignment and warping, putting everything in  %
%%% one Mtex variable on same grid (EBSD data, SEM microstructure images, %
%%% strain data, (hr-ebsd data)                                           %
%%%                                                                       %
%%% IMPORTANT: all dimensions in micrometers                              %    
%%% Authors: Tijmen Vermeij, Gert-Jan Slokker                              %
%%% Eindhoven University of Technology; 2023                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%
%%%%%% THIS SCRIPT NEEDS INPUT AND DEFINITION OF OPTIONS AT MULTIPLE
%%%%%% LOCATIONS, SO BEST TO RUN SECTION BY SECTION IN THE BEGINNING.
%%%%%%%%

%%%%%%%%%%%
%%% Data used from Allan Harte et al., Acta Mater 2020; https://zenodo.org/doi/10.5281/zenodo.3691902
%%%%%%%%%%%


%% INPUT and MTEX Initiation

clearvars; close all; clc;
%--------------------------------------------------------------------------

% add path of src scripts
addpath('src');

%--------------------------------------------------------------------------
% set mtex pref
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');


%%% crystal symmetry for importing EBSD data (not always necessary, depends
% on how ebsd data is imported)

% crystal symmetry. In this case: HCP Zn
CS = {  'notIndexed',...
        crystalSymmetry('m-3m', [3.6 3.6 3.6], 'mineral', 'Ni-superalloy', 'color', [0.53 0.81 0.98])};

%%% Run and Export settings

%%% fnames of selection points EBSD,BSE and DIC [if defined --> automatically use for this run]
SelectionFnameEBSD = 'Selection_Points_Backup_EBSD_02-05-2024_17-28.mat';

exportfig = 0;                          %if 0: no export of overview figure      & 1: export of overview figure
exportdat = 1;                          %if 0: no export of data                 & 1: export of data

%--------------------------------------------------------------------------


%%% Plotting & Colorbars settings
Fontsize = 9;                                          %fontsize titles
Clim = [0.02 0.10];                                    %plotting range Equivalent strain

colorbar_Equi = importdata('cm_haline.mat');          %colorbar Eequivalent Field
colorbar_CI = 'parula';                               %colorbar CI index
%--------------------------------------------------------------------------

%%% Experimental settings
% L0 = 8500;                          %gauge length of dogbone [um]
% incs =[100 200 300 400 500];         %elongation increments [um]
% straininc = (incs./L0)*100;          %"global"strain increments [%]
%--------------------------------------------------------------------------





%% INPUT FOR DATA IMPORT
% % names of types of raw data (filenames, psize, datatypes in STRUCTS)
% Data_names = {'EBSD';'BSE';'DIC_BSE';'DIC_REF';'DIC_BSE1'};
% raw_data.Data_names = Data_names;



%{
%%% for each set of data, define the data name, the filename, datatype, and possible the
%%% pixel size and imregion (the region of the tif which should be taken as
%%% the imagedata. For  ".mat" data, define the name of the variables
%%% of the position fields (in um!) and the image data (data)
%}

dataFolder = 'dataNi';

% EBSD data [input: .crc file] %%% NOTE THAT WE APPLY AN EXTRA ROTATION TO
% PROPERLY IMPORT INTO MTEX. PLEASE ALWAYS CHECK YOURSELF.
DataName = 'EBSD0';                                                             % name of the input data
rawData.(DataName).fname = fullfile(dataFolder, 'r1c3.crc');         % Edax/EMSphinx/oxford output
rawData.(DataName).dtype = 'EBSD';                                              % data type
rawData.(DataName).rot = rotation.byEuler(180*degree,0,0);         % EXTRA ROTATION to assure correct import in Mtex. Please check yourself
rawData.(DataName).CS = CS;                                                     % crystal symmetry information
clear DataName

%%% store list of data names in rawData %%%
rawData.dataNames = fieldnames(rawData);


%--------------------------------------------------------------------------
%% LOAD ALL RAW DATA
rawData = loadRawData(rawData);


%%%%% no DIC, BSE, SE images available. So, DIC displacement data is put on
%%%%% a grid and used for alignment. We do this using a "dummy" EBSD
%%%%% variable
fname_dic = fullfile(dataFolder,'r1c3.txt');
% import DIC data
dic=importdata(fname_dic);
% pixel size
psize = 0.0146;
% make arrays of all data in "prop" struct, which will be input to mtex
prop.x = dic.data(:,1)*psize;
prop.y = dic.data(:,2)*psize;
prop.u = dic.data(:,3)*psize;
prop.v = dic.data(:,4)*psize;

% dummy phase and orientation fields
phase = ones(size(prop.x));
ori = repmat(orientation.byEuler([0 0 0]),size(prop.x));

% create dummy EBSD variable
DIC = EBSD(ori,phase,CS,prop);
DICgrid = DIC.gridify;

% calculate some strains to later use for alignment
data.X = DICgrid.x; data.Y = DICgrid.y; data.U = DICgrid.prop.u; data.V = DICgrid.prop.v;
[~,~,strain] = calcStrainsFromRawDisp(data,1);
%%%%%

%%% add data manually to rawData struct
rawData.DIC.fname = '';
rawData.DIC.dtype = 'mat_data';
rawData.DIC.data = strain.Eequi;
rawData.DIC.data(rawData.DIC.data > 0.2) = 0.2; %to avoid having to change colorbar later
rawData.DIC.X = DICgrid.x;
rawData.DIC.Y = DICgrid.y;
rawData.DIC.psize = DICgrid.dx;
[rawData.DIC.m, rawData.DIC.n] = size(DICgrid);

rawData.DICU = rawData.DIC;
rawData.DICU.data = DICgrid.prop.u;

rawData.DICV = rawData.DIC;
rawData.DICV.data = DICgrid.prop.v;





%%%%%% ALIGNMENTS -- OPTIONS AND INPUT ARE REQUIRED %%%%%%%


%--------------------------------------------------------------------------
%% 1). Alignment of EBSD data to BSE Data
%{
% FOR EACH DATASET WHICH NEEDS TO BE ALIGNED, THE REFERENCE (WITH RESPECT TO 
% WHAT EVERYTHING SHOULD BE ALIGNED TO; FOR NOW SHOULD BE THE SAME IN ALL 
% ALIGNMENTS)
% SOME OPTIONS ARE REQUIRED
%}

% Choose reference data name (should be a field in rawData)
REF = 'DIC';                            % undeformed = DIC strain data

% Choose "deformed" (Warped) data name, which needs to be aligned. (should be a field in rawData)
DEF = 'EBSD0';                           % deformed = EBSD data

%%%%
%%%%%%%%
%%%%

% Setting the point selection and alignment options
%{
% Alignment is done by using coordinates of selected points (in reference and
"deformed/warped" to fit a distortion/warping field. The main input
variables are the number of degrees of freedom (DOFs), in "options.Ndof",
and the number of selection points, in "options.Nsel". 
- Ndof determines how many terms of a 2D polynomial field are used. 
--- Ndof = 6 --> 1st order, Ndof = 12 --> 2nd order, Ndof = 20 --> 3rd
    order (values in between also possible)
- If Ndof == 0, only pure rigid body motion will be considered, (rigid body
translation and rotation) using a special fitting function. This is
recommended to perform alignment between fast SEM scans.
- Nsel should be at least half of Ndof. Higher Nsel means there will be
least square fitting which means lower errors generally.
%}

% Ndof, number of DOFs
options.Ndof = 6;    % might not be enough for EBSD, be careful.

% Nsel, number of selection points
options.Nsel = 4;                            %Nsel = at least Ndof*0.5 (more is better)

% Choose name for alignment, used for filename to save selection points and
% to save alignment data
options.case = 'EBSD';

% Define some information for overlay of EBSD and CI data

%%% in this case, the IPF orientations are added as a third
%%% plot to allow easier point selection
options.selectionIPF = 1;

% choose which EBSD "property" to use for indexing (IQ, CI, BC,...)
% since we use Oxford Instruments data here, we use BC (band contrast)
options.EBSD_selection = 'bc';

% colormap for alignment
options.colormap = 'gray';

% Using pre-defined, stored selection points [only if the fnameEBSD are defined]
if exist('SelectionFnameEBSD','var')
    options.sel_points_fname = SelectionFnameEBSD;
end

% Align position fields of data
alignment_data_EBSD = pointSelectionForAlignment(rawData.(REF),rawData.(DEF),options);

if ~exist('SelectionFnameEBSD','var')
    % Saving alignment_data_EBSD
    save(['alignment_data_EBSD_',datestr(now,'mm-dd-yyyy_HH-MM')],'alignment_data_EBSD')
end

% Remove the sel_points_fname field to be sure
clear options


%--------------------------------------------------------------------------
%% %% Check DIC,EBSD alignment & verification plotting
%{
%first create basic mtex struct, only from warped ebsd, use the specified pixel size 
%in the new dataset, in um.
%}
% Defining a new pixelsize [um]
psize_new = 0.2;

% Defining a new grid, in the same area as the BSE scan (which was used as
% the reference for alignment)
[Xnew,Ynew] = meshgrid(0:psize_new:max(rawData.DIC.X(:)),0:psize_new:max(rawData.DIC.Y(:))); 

% warp EBSD data to new grid
newEBSD = warpCreateEBSDdata(rawData.EBSD0.data_grid,Xnew,Ynew,alignment_data_EBSD);
newEBSD  = newEBSD.gridify;

% include the DIC data also in the new mtex variable (no alignment data
% needed since the reference is equal to the DIC configuration)
newEBSD = warpIncludeDataInMtex(newEBSD,rawData.DIC,'Eequi'); 
newEBSD = warpIncludeDataInMtex(newEBSD,rawData.DICU,'U'); 
newEBSD = warpIncludeDataInMtex(newEBSD,rawData.DICV,'V'); 

%--------------------------------------------------------------------------
%%% Verification plotting EBSD-BSE alignment (with update option)

% Calculate some quick grains to plot boundaries for checking alignment
[grains] = calcGrains(newEBSD,'angle',15*degree); 

figure()
h = newMtexFigure('nrows',1,'ncols',1,'visible','on','Position',get(0,'ScreenSize'));
plot(newEBSD,newEBSD.prop.Eequi,'micronbar','off')                       %DIC_BSE reference
hold on
plot(grains.boundary,'lineColor','red','linewidth',1.5,'Marker','.')                  %warped grain boundaries







%--------------------------------------------------------------------------
%% Cutting image (mapped data)
figure();
h = newMtexFigure();
plot(newEBSD,newEBSD.prop.Eequi)          %Warped Equi. strain data (incr. 2)
hold on
plot(grains.boundary,'lineColor','red','linewidth',1.5)         %Warped grain boundaries
mtexTitle('DIC','Fontsize',Fontsize)
set(gca,'ColorScale','log')                     
setColorRange([10^-2 0.25*max(max(newEBSD.prop.Eequi))]);                       %Setting colorscale (removing noise <0.03 strain) 
colormap(colorbar_Equi)                                %Alternative options: viridis, jet, gray, default
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',3,'XColor','k','YColor','k')

% Selecting relevant area

title('Click upper left and lower right point to cut out a new region')
xycut = ginput(2);
newEBSD_Cut = newEBSD(inpolygon(newEBSD,[xycut(1,1),xycut(1,2),xycut(2,1)-xycut(1,1),xycut(2,2)-xycut(1,2)]));


nextAxis
plot(newEBSD_Cut,newEBSD_Cut.orientations) % Cutted grain boundary data
mtexTitle('EBSD','Fontsize',Fontsize)
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',3,'XColor','k','YColor','k')


%% Overview plotting

% Re-calculating and cleaning the grains for new cutted area
[grains, newEBSD_Cut.grainId] = calcGrains(newEBSD_Cut,'angle',15*degree); % 5-15 (not use below 5)

%
close all
figure()
h = newMtexFigure('nrows',1,'ncols',2,'visible','on','Position',get(0,'ScreenSize'));

newEBSD_Cut= newEBSD_Cut.gridify;

%define the plottingrange of the current figure
plottingrange = [min(newEBSD_Cut('indexed').x) max(newEBSD_Cut('indexed').x) ...
                 min(newEBSD_Cut('indexed').y) max(newEBSD_Cut('indexed').y)];


plot(newEBSD_Cut,newEBSD_Cut.prop.(['Eequi']),'micronbar','off')  
hold on

colormap(h.currentAxes,colorbar_Equi)                                   
set(gca,'ColorScale','log')
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',3,'XColor','k','YColor','k')
h.currentAxes.CLim = Clim;
mtexTitle(['DIC-Incr1'],'Fontsize',Fontsize)


hold on
plot(grains.boundary,'lineColor','red','linewidth',1.5)


cb_equi = colorbar('location','eastoutside');
cb_equi.Label.String = '$E_{eff}$ [-]';
cb_equi.TickLabelInterpreter = 'latex';      
cb_equi.Label.Interpreter = 'latex';
cb_equi.FontSize = 10;   




nextAxis

% plot other things


plot(newEBSD_Cut,newEBSD_Cut.orientations,'micronbar','off')        %Warped MTEX data
hold on
plot(grains.boundary,'lineColor','black','linewidth',1.5)         %Warped grain boundaries
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',3,'XColor','k','YColor','k')
mtexTitle('EBSD','Fontsize',Fontsize)



%plotting an grain overview with current added grainid's 
figure()
h = newMtexFigure();
plot(grains,grains.meanOrientation,'micronbar','off')
hold on
plot(grains.boundary,'lineColor','k','linewidth',1.5)           
text(grains,int2str(grains.id),'color','k')

set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',3,'XColor','k','YColor','k')




 



