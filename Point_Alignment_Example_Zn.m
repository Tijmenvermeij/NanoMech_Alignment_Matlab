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
            crystalSymmetry('6/mmm', [2.6649 2.6649 4.9468], 'X||a', 'Y||b*', 'Z||c*', 'mineral', 'Zinc', 'color', [1 0.6 0])};

%%% Run and Export settings

%%% fnames of selection points EBSD,BSE and DIC [if defined --> automatically use for this run]
% SelectionFnameEBSD = 'Selection_Points_Backup_EBSD_02-01-2024_15-16.mat';
% SelectionFnameBSE = 'selection_points_backup_BSE_04-12-2023_11-32';
% SelectionFnameDIC = 'selection_points_backup_DIC_04-12-2023_09-53';

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
L0 = 8500;                          %gauge length of dogbone [um]
incs =[100 200 300 400 500];         %elongation increments [um]
straininc = (incs./L0)*100;          %"global"strain increments [%]
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

dataFolder = 'dataZn';

% EBSD data [input: .ang file]
DataName = 'EBSD0';                                                             % name of the input data
rawData.(DataName).fname = fullfile(dataFolder, 'N1_20kV_depth_WD20_FOV50_stepsize0_15.ang');         % Edax/EMSphinx output
rawData.(DataName).dtype = 'EBSD';                                              % data type
rawData.(DataName).CS = CS;                                                     % crystal symmetry information
clear DataName

% BSE microstructure [input: .tif file (FOV,Pixels)]
DataName = 'BSE';                                                   % name of the input data
rawData.(DataName).fname = fullfile(dataFolder,'N1_x10.5794_ymin0.5027_FOV50.tif');      % SE, BSE image filename
rawData.(DataName).psize = 50/3072;                                 % pixelsize in um (could be extracted from e.g. hdr file in theory)
rawData.(DataName).dtype = 'SEM';                                   % data type
rawData.(DataName).imregion = [3073 1 3071 3071];                   % startx starty numx-1 numy-1 (in this case the second half of the image, without databar)
clear DataName

% DIC pattern Reference high Kv [input: .tif file (FOV,Pixels)] (BSE)
DataName = 'DIC_BSE';                                                               % name of the input data
rawData.(DataName).fname = fullfile(dataFolder,'N1_x10.5622_ymin0.6514_FOV50.tif');             % SE,BSE image filename
rawData.(DataName).psize = 50/3072;                                        % pixelsize
rawData.(DataName).dtype = 'SEM';                                          % data type
rawData.(DataName).imregion = [3073 1 3071 3071];                          % startx starty numx-1 numy-1
clear DataName

% DIC pattern Reference high Kv [input: .tif file (FOV,Pixels)] (SE) (SAME
% SCAN AS BSE, BUT DIFFERENT QUADRANT
DataName = 'DIC_SE';                                                               % name of the input data
rawData.(DataName).fname = fullfile(dataFolder,'N1_x10.5622_ymin0.6514_FOV50.tif');             % SE,BSE image filename
rawData.(DataName).psize = 50/3072;                                        % pixelsize
rawData.(DataName).dtype = 'SEM';                                          % data type
rawData.(DataName).imregion = [1 1 3071 3071];                          % startx starty numx-1 numy-1
clear DataName

% DIC Reference low Kv, (used to align strains etc) [input: .tif file (FOV,Pixels)]
DataName = 'DIC';
rawData.(DataName).fname = fullfile(dataFolder, 'Square_IBSE_N1_328_5_Elong_0um_Vfield35um_PIX4096_5kV_10Bi_corrected.tif');   %(IB)-SE (low kV) image
rawData.(DataName).psize = 35/4096;                                        %pixelsize
rawData.(DataName).dtype = 'SEM';                                          %data type
rawData.(DataName).imregion = [1 1 4095 4095];                             %startx starty numx-1 numy-1
clear DataName

%%% store list of data names in rawData %%%
rawData.dataNames = fieldnames(rawData);


%--------------------------------------------------------------------------
%% LOAD ALL RAW DATA
rawData = loadRawData(rawData);


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
REF = 'BSE';                            % undeformed = BSE data

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
options.Ndof = 12;

% Nsel, number of selection points
options.Nsel = 10;             %Nsel = at least Ndof*0.5 (more is better)

% Choose name for alignment, used for filename to save selection points and
% to save alignment data
options.case = 'EBSD';

% Define some information for overlay of EBSD and CI data

%%% in this case, the IPF orientations are added as a third
%%% plot to allow easier point selection
options.selectionIPF = 1;

% choose which EBSD "property" to use for indexing (IQ, CI, BC,...)
% since we used EmSphInX for indexing, we want to use CI instead of IQ for
% alignment
options.EBSD_selection = 'ci';

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
%% %% Check BSE,EBSD alignment & verification plotting
%{
%first create basic mtex struct, only from warped ebsd, use the specified pixel size 
%in the new dataset, in um.
%}
% Defining a new pixelsize [um]
psize_new = 0.05;

% Defining a new grid, in the same area as the BSE scan (which was used as
% the reference for alignment)
[Xnew,Ynew] = meshgrid(0:psize_new:max(rawData.DIC_BSE.X(:)),0:psize_new:max(rawData.DIC_BSE.Y(:))); 

% warp EBSD data to new grid
newEBSD = warpCreateEBSDdata(rawData.EBSD0.data_grid,Xnew,Ynew,alignment_data_EBSD);
newEBSD  = newEBSD.gridify;

%include the DIC data also in the new mtex variable
newEBSD = warpIncludeDataInMtex(newEBSD,rawData.BSE,'BSE');  

%--------------------------------------------------------------------------
%%% Verification plotting EBSD-BSE alignment (with update option)

% Calculate some quick grains to plot boundaries for checking alignment
[grains] = calcGrains(newEBSD,'angle',15*degree); 

figure()
h = newMtexFigure('nrows',1,'ncols',1,'visible','on','Position',get(0,'ScreenSize'));
plot(newEBSD,newEBSD.prop.BSE,'micronbar','off')                       %DIC_BSE reference
hold on
plot(grains.boundary,'lineColor','red','linewidth',1.5,'Marker','.')                  %warped grain boundaries





%--------------------------------------------------------------------------
%% 2). Alignment of BSE to BSE-DIC data

% Choose reference data
REF = 'DIC_BSE';                          %reference = BSE-DIC (BSE) data

% Choose "deformed" (Warped) data, which needs to be aligned
DEF = 'BSE';                          %deformed = BSE                          

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

% Ndof, here 0 to use pure RBM
options.Ndof = 0;           
options.Nsel = 4;                %Nsel; 2 should be sufficient but we use 4 for safety and accuracy

% Choose name for alignment, used to save selection points
options.case = 'BSE';
options.colormap = 'gray';

% Using pre-defined, stored selection points [only if the fnameBSE are defined]
if exist('SelectionFnameBSE','var')
    options.sel_points_fname = SelectionFnameBSE;
end

% Align position fields of data
alignment_data_BSE = pointSelectionForAlignment(rawData.(REF),rawData.(DEF),options);

if ~exist('SelectionFnameBSE','var')
    % Saving alignment_data_BSE
    save(['alignment_data_BSE_',datestr(now,'mm-dd-yyyy_HH-MM')],'alignment_data_BSE')
end

% Remove the sel_points_fname field to be sure
clear options




%--------------------------------------------------------------------------
%% 3). Align DIC (low kV) to BSE-DIC (SE)
% Choose reference data
REF = 'DIC_SE';        %reference = BSE-DIC data (SE)

% Choose "deformed" (Warped) data, which needs to be aligned
DEF = 'DIC';        %deformed = DIC (low kv)

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

options.Ndof = 0;           
options.Nsel = 4;               

% Choose name for alignment, used for filename to save selection points
options.case = 'DIC';
options.colormap = 'gray';

% Using pre-defined, stored selection points [only if the fnameEBSD are defined]
if exist('SelectionFnameDIC','var')
   options.sel_points_fname = SelectionFnameDIC;
end

% Align position fields of data
alignment_data_DIC = pointSelectionForAlignment(rawData.(REF),rawData.(DEF),options);

if ~exist('SelectionFnameDIC','var')
    % Saving alignment_data_DIC
    save(['alignment_data_DIC_',datestr(now,'mm-dd-yyyy_HH-MM')],'alignment_data_DIC')
end


clear options

%--------------------------------------------------------------------------
%% Adding BSE and DIC data to the REF data (created earlier)

% warp BSE mapped EBSD dataset to the DIC_BSE dataset            
newEBSD2 = warpCreateEBSDdata(newEBSD,Xnew,Ynew,alignment_data_BSE);   
newEBSD2 = newEBSD2.gridify;
newEBSD2 = warpIncludeDataInMtex(newEBSD2,rawData.DIC_BSE,'DIC_BSE');
% warp the DIC_SEM dataset to BSE,DIC_BSE mapped EBSD dataset
newEBSD2 = warpIncludeDataInMtex(newEBSD2,rawData.DIC,'DIC_SEM',alignment_data_DIC);
% overwrite EBSD dataset
newEBSD = newEBSD2;  

%% Calculate some quick grains to plot boundaries for checking alignment
[grains] = calcGrains(newEBSD,'angle',15*degree);  
%

figure()
h = newMtexFigure('nrows',1,'ncols',3);
plot(newEBSD,newEBSD.prop.BSE)
hold on
plot(grains.boundary,'lineColor','red','linewidth',1.5)  
nextAxis
plot(newEBSD,newEBSD.prop.DIC_BSE)  
hold on
plot(grains.boundary,'lineColor','red','linewidth',1.5)  
nextAxis
plot(newEBSD,newEBSD.prop.DIC_SEM)  
hold on
plot(grains.boundary,'lineColor','red','linewidth',1.5)  
%


%--------------------------------------------------------------------------
%% Add DIC displacement and strain fields to aligned data

% Defining folder with DIC data
dicFolder = fullfile(pwd, dataFolder, 'DIC_DATA_N1_sub33_step3_Vfield35_pix4096');

% Setting pixelsize [um] (= pixelsize of DIC in rawData) 
psize = rawData.DIC.psize;

% Defining the used stepsize [in pixels] within MatchID
stepsize = 3;            

% Settting images (strain incr.) used for correlation (purely to easily get
% the DIC data)

% Number of incs
n = length(incs);
for i = 1:n
    incnames{i} = (['Square_IBSE_N1_328_5_Elong',num2str(incs(i)),'um_Vfield35um_PIX4096_5kV_10Bi_corrected.tif']);
end

% Strain filter, in pixels, for smoothening of strain fields (see Exp.
% Mech. paper for details)
strain_filter = 1;    %since stepsize is 3, filter of 1 is good enough


% Get position field of dic data, in pixels, (should be the same for all incs)
X = importdata(fullfile(dicFolder, 'x_pic', [incnames{1}, '_x_pic.csv']));
Y = importdata(fullfile(dicFolder, 'y_pic', [incnames{1}, '_y_pic.csv']));

[Xp,Yp] = meshgrid(min(X(:)):stepsize:max(X(:)),min(Y(:)):stepsize:max(Y(:)));

% Defining the position fields; stored in struct for easy use
dicData.X = Xp;
dicData.Y = Yp;

% Calculating Equi. strains and align Disp. and strains fields 
for i = 1:n

    fprintf(['Progression: ',num2str(i),' of ',num2str(n),'\n'])
    
    % import displacement field
    dicData.U = importdata(fullfile(dicFolder, 'u', [incnames{i}, '_u.csv']));
    dicData.V = importdata(fullfile(dicFolder, 'v', [incnames{i}, '_v.csv']));
     
    %----------------------------------------------------------------------
    % compute strain field
    [Xraw,Yraw,strainsRaw] = calcStrainsFromRawDisp(dicData,strain_filter);

    % X and Y are in pixels, so convert to microns (to be the same as the ref image X and Y array)
    straindata_raw.X = Xraw*psize;
    straindata_raw.Y = Yraw*psize;
    straindata_raw.data = strainsRaw.Eequi;

    % Displacements
    Udata_raw = straindata_raw; Vdata_raw = straindata_raw;
    Udata_raw.data=dicData.U*psize; Vdata_raw.data=dicData.V*psize;

    % warp and include data in mtex variable
    newEBSD = warpIncludeDataInMtex(newEBSD,straindata_raw,strcat('inc_',num2str(i),'_Eequi'),alignment_data_DIC);
    newEBSD = warpIncludeDataInMtex(newEBSD,Udata_raw,strcat('inc_',num2str(i),'_U'),alignment_data_DIC);
    newEBSD = warpIncludeDataInMtex(newEBSD,Vdata_raw,strcat('inc_',num2str(i),'_V'),alignment_data_DIC);

%     %----------------------------------------------------------------------
%     % Interpolating MISSING Displacement and Strain Field %
%     [dicData.U, dicData.V] = interpDispField(dicData.U,dicData.V,strain_filter,'cubic');        %=Std of applied filter for displacement == strain_filter
% 
%     [X,Y,strains] = calc_strains_from_raw_disp(dicData,strain_filter);
% 
%     % X and Y are in pixels, so convert to microns (to be the same as the ref image X and Y array)
%     straindata.X=X*psize;
%     straindata.Y=Y*psize;
% 
%     % Only save equivalent strains for now
%     straindata.data = strains.Eequi;
% 
%     % Displacements
%     Udata = straindata; Vdata = straindata;
%     Udata.data=dicData.U*psize; Vdata.data=dicData.V*psize;
% 
%     % Rotation
%     rotdata = straindata;
%     rotdata.data = strains.rot;
%     
%     newEBSD = warpIncludeDataInMtex(newEBSD,straindata,strcat('inc_',num2str(i),'_Eequi'),alignment_data_DIC);
%     newEBSD = warpIncludeDataInMtex(newEBSD,Udata,strcat('inc_',num2str(i),'_U'),alignment_data_DIC);
%     newEBSD = warpIncludeDataInMtex(newEBSD,Vdata,strcat('inc_',num2str(i),'_V'),alignment_data_DIC);
%     newEBSD = warpIncludeDataInMtex(newEBSD,rotdata,strcat('inc_',num2str(i),'_Rot'),alignment_data_DIC);
end


%--------------------------------------------------------------------------
%% Cutting image (mapped data)
figure();
h = newMtexFigure();
plot(newEBSD,newEBSD.prop.inc_5_Eequi)          %Warped Equi. strain data (incr. 2)
hold on
plot(grains.boundary,'lineColor','red','linewidth',1.5)         %Warped grain boundaries
mtexTitle('DIC','Fontsize',Fontsize)
set(gca,'ColorScale','log')                     
setColorRange([10^-2 0.25*max(max(newEBSD.prop.inc_1_Eequi))]);                       %Setting colorscale (removing noise <0.03 strain) 
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
[grains, newEBSD_Cut.grainId] = calcGrains(newEBSD_Cut,'angle',5*degree); % 5-15 (not use below 5)

%
close all
figure()
h = newMtexFigure('nrows',2,'ncols',6,'visible','on','Position',get(0,'ScreenSize'));

newEBSD_Cut= newEBSD_Cut.gridify;

%define the plottingrange of the current figure
plottingrange = [min(newEBSD_Cut('indexed').x) max(newEBSD_Cut('indexed').x) ...
                 min(newEBSD_Cut('indexed').y) max(newEBSD_Cut('indexed').y)];

for i = 1:length(incs)

    plot(newEBSD_Cut,newEBSD_Cut.prop.(['inc_' num2str(i) '_Eequi']),'micronbar','off')  
    hold on

    colormap(h.currentAxes,colorbar_Equi)                                   
    set(gca,'ColorScale','log')
    set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',3,'XColor','k','YColor','k')
    h.currentAxes.CLim = Clim;
    mtexTitle(['DIC-Incr ',num2str(i),'. Elong. ',num2str(incs(i)),'$\mathbf{\mu}$m [$\mathbf{\epsilon_{g}}$: ' num2str(straininc(i),'%.2f') '%]'],'Fontsize',Fontsize)


    hold on
    plot(grains.boundary,'lineColor','red','linewidth',1.5)

%     newEBSD_grid_cut_raw = newEBSD_grid_cut_raw.gridify;
% 
%     boundaries = zeros(size(newEBSD_grid_cut_raw)); % make a matrix to define the boundaries
%     % everywhere the strainfield is NAN, the boundary matrix becomes 1
%     boundaries(isnan(newEBSD_grid_cut_raw.prop.(['inc_' num2str(i) '_Eequi_raw']))) = 1;
% 
%     % plot the boundaries of the NAN-zones in pink
%     vb = visboundaries(boundaries,'color','magenta','LineWidth',0.5);
% 
%     % the boundaries are in "point coordinates" and need to be placed on top of
%     % the EBSD map by using the coordinates from the ebsd
%     vb.Children(1).XData = vb.Children(1).XData*psize_new+min(newEBSD_grid_cut_raw.x(:))-psize_new; % 0.02 depends on the pixelsize
%     vb.Children(1).YData = vb.Children(1).YData*psize_new+min(newEBSD_grid_cut_raw.y(:))-psize_new;
%     vb.Children(2).XData = vb.Children(2).XData*psize_new+min(newEBSD_grid_cut_raw.x(:))-psize_new;
%     vb.Children(2).YData = vb.Children(2).YData*psize_new+min(newEBSD_grid_cut_raw.y(:))-psize_new;
% 
%     %plot(grains.boundary,'lineColor','red','linewidth',1.0)                     %Warped grain boundaries 
% 
    if i == length(incs)
        cb_equi = colorbar('location','eastoutside');
        cb_equi.Label.String = '$E_{eff}$ [-]';
        cb_equi.TickLabelInterpreter = 'latex';      
        cb_equi.Label.Interpreter = 'latex';
        cb_equi.FontSize = 10;   
    end
 
    nextAxis

end

nextAxis

% plot other things


plot(newEBSD_Cut,newEBSD_Cut.orientations,'micronbar','off')        %Warped MTEX data
hold on
plot(grains.boundary,'lineColor','black','linewidth',1.5)         %Warped grain boundaries
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',3,'XColor','k','YColor','k')
mtexTitle('EBSD','Fontsize',Fontsize)

nextAxis
plot(newEBSD_Cut,newEBSD_Cut.prop.DIC_BSE,'micronbar','off')        %Reference DIC_BSE image
colormap(h.currentAxes,'gray') 
hold on
plot(grains.boundary,'lineColor','red','linewidth',1.5)     %Warped grain boundaries
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',3,'XColor','k','YColor','k')
mtexTitle('REF','Fontsize',Fontsize)
                  
nextAxis
plot(newEBSD_Cut,newEBSD_Cut.prop.BSE,'micronbar','off')            %Warped BSE image
colormap(h.currentAxes,'gray') 
hold on
plot(grains.boundary,'lineColor','red','linewidth',1.5)     %Warped grain boundaries
set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',3,'XColor','k','YColor','k')
mtexTitle('BSE','Fontsize',Fontsize)

nextAxis

% Exporting overview figure
if exportfig == 1
    exportgraphics(figure(1),fullfile(pwd,['overview_image_',datestr(now,'mm-dd-yyyy_HH-MM'),'.png']),'Resolution',1024)
    savefig(figure(1),fullfile(pwd,['overview_image_',datestr(now,'mm-dd-yyyy_HH-MM'),'.fig']))
end

%plotting an grain overview with current added grainid's 
figure()
h = newMtexFigure();
plot(grains,grains.meanOrientation,'micronbar','off')
hold on
plot(grains.boundary,'lineColor','k','linewidth',1.5)           
text(grains,int2str(grains.id),'color','k')

set(h.currentAxes, 'box', 'on', 'Visible', 'on', 'xtick', [], 'ytick', [],'linewidth',3,'XColor','k','YColor','k')

if exportfig == 1
    exportgraphics(figure(3),fullfile(pwd,['Grainid_Overview.png']),'Resolution',1024)
end




 



