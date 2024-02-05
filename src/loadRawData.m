function rawData = loadRawData(rawData,varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Function for loading of raw data for alignment                        %
%%% Athors: Tijmen Vermeij, Gert-Jan Slokker, Jorn Verstijnen             %
%%% Eindhoven University of Technology; 2023                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% RAW DATA LOADING

% loop over different data sets, load all data and put in one struct
% pick loading indices (of rawData.dataNames)
if nargin > 1
    ind = varargin(1);
else
    ind = 1:length(rawData.dataNames);
end

% loop over data
for i=ind    
    DataName = rawData.dataNames{i};
    % switch between types of data (EBSD or SEM for now)
    switch rawData.(DataName).dtype
        case 'EBSD'
            % load EBSD data, gridify, extract X and Y matrices in um, sizes and
            % X and Y matrices in pixels
            rawData.(DataName).data = EBSD.load(rawData.(DataName).fname, rawData.(DataName).CS, 'convertEuler2SpatialReferenceFrame','setting 2');
            % if given as input, apply rotation to EBSD data.
            if isfield(rawData.(DataName),'rot')
                rawData.(DataName).data = rawData.(DataName).data.rotate(rawData.(DataName).rot);
            end
            
            % get pixelsize
            rawData.(DataName).psize = sum(abs(rawData.(DataName).data.unitCell(1,:)));
            % gridify data
            rawData.(DataName).data_grid = rawData.(DataName).data.gridify;        
            % get x and y grid
            rawData.(DataName).X = rawData.(DataName).data_grid.prop.x;
            rawData.(DataName).Y = rawData.(DataName).data_grid.prop.y;
            % get other stuff
            [rawData.(DataName).m, rawData.(DataName).n] = size(rawData.(DataName).X);
            [rawData.(DataName).Xp,rawData.(DataName).Yp] = meshgrid(1:rawData.(DataName).n,1:rawData.(DataName).m);
                        
        case 'SEM'
            % load SEM data, create X and Y matrices in um, sizes and
            % X and Y matrices in pixels
            rawData.(DataName).data = double(imcrop(imread(rawData.(DataName).fname),rawData.(DataName).imregion));
            [rawData.(DataName).m, rawData.(DataName).n] = size(rawData.(DataName).data);
            [rawData.(DataName).Xp,rawData.(DataName).Yp] = meshgrid(1:rawData.(DataName).n,1:rawData.(DataName).m);
            rawData.(DataName).X = (rawData.(DataName).Xp-1) * rawData.(DataName).psize; % start at zero
            rawData.(DataName).Y = (rawData.(DataName).Yp-1) * rawData.(DataName).psize; % start at zero
        case 'mat_data'
            % load matlab data directly from mat file, using specified
            % filenames
            matdata = load(rawData.(rawData.dataNames{i}).fname,rawData.(rawData.dataNames{i}).xdata_name,rawData.(rawData.dataNames{i}).ydata_name,rawData.(rawData.dataNames{i}).data_name);
            rawData.(DataName).data = matdata.(rawData.(rawData.dataNames{i}).data_name);
            rawData.(DataName).X = matdata.(rawData.(rawData.dataNames{i}).xdata_name);
            rawData.(DataName).Y = matdata.(rawData.(rawData.dataNames{i}).ydata_name);
            
            rawData.(DataName).psize = diff(rawData.(DataName).X(1,1:2));
            [rawData.(DataName).m, rawData.(DataName).n] = size(rawData.(DataName).data);
            
            clear matdata
        otherwise
            error('dtype not supported')
    end
end
