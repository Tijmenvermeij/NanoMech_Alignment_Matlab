function ebsd = warpIncludeDataInMtex(ebsd,rawDataField,dataName,alignmentData)
%%% this function adds data to an Mtex struct, can be sem image, strain,
%%% etc. Alignment Data is used to warp it correctly to the EBSD grid, if
%%% necessary
%%% alignmentData can be used as input if warping is required, if only 3
%%% input arguments, no warping is performed, but only interpolation

%%% author: Tijmen Vermeij, TUe, 2018-2023


% interpmethod = 'nearest';
% interpmethod = 'linear';
interpmethod = 'cubic';


if nargin == 4
    %%%% extend the displacement (distortion) field, if necessary (when the to
    %%%% be used grid is larger than the available grid)
    if min(ebsd.x(:)) < min(alignmentData.X(:)) || min(ebsd.y(:)) < min(alignmentData.Y(:)) || ...
            max(ebsd.x(:)) > max(alignmentData.X(:)) || max(ebsd.y(:)) > max(alignmentData.Y(:))

        %%% recreate displacement fields, based on new positions

        Xf = ebsd.x;
        Yf = ebsd.y;

        % extract displacement fields on original f grid 
        Ux_f = zeros(size(Xf));
        Uy_f = zeros(size(Yf));
        
        if alignmentData.options.Ndof > 0
            options.basis = 'p';
            options.normalized = 0;

            % loop over the shapefunctions
            for kdof = 1:length(alignmentData.u)
                % call each shapefunction
                [Phi dir] = basis(Xf(1,:),Yf(:,1)',alignmentData.phi,kdof,options);
                % rebuild the displacement fields
                if dir == 1
                    Ux_f  = Ux_f + alignmentData.u(kdof) * Phi ;
                elseif dir == 2
                    Uy_f  = Uy_f + alignmentData.u(kdof) * Phi ;
            %             elseif dir == 3
            %                 Uz  = Uz + u(kdof) * Phi ;
            %             else
            %                 error('unexpected value for the direction of phi')
                end
            end
        else
            % generate displacement fields using translation and rotation
            % parameters
            Ux_f = Xf * alignmentData.regParams.R(1,1) + Yf * alignmentData.regParams.R(1,2) + alignmentData.regParams.t(1) - Xf;
            Uy_f = Xf * alignmentData.regParams.R(2,1) + Yf * alignmentData.regParams.R(2,2) + alignmentData.regParams.t(2) - Yf;
        end
        
        alignmentData.X = Xf;
        alignmentData.Y = Yf;

        alignmentData.Ux = Ux_f;
        alignmentData.Uy = Uy_f;
    end
    
    
    % interpolate warping displacement fields to mtex grid
    Ux = interp2(alignmentData.X,alignmentData.Y,alignmentData.Ux,ebsd.prop.x,ebsd.prop.y,interpmethod);
    Uy = interp2(alignmentData.X,alignmentData.Y,alignmentData.Uy,ebsd.prop.x,ebsd.prop.y,interpmethod);
    
    % warp data
    data = interp2(rawDataField.X,rawDataField.Y,rawDataField.data,ebsd.prop.x+Ux,ebsd.prop.y+Uy,interpmethod,0);
else
    % interpolate data to required grid
    data = interp2(rawDataField.X,rawDataField.Y,rawDataField.data,ebsd.prop.x,ebsd.prop.y,interpmethod,0);
end

% input in Mtex variable
ebsd.prop.(dataName) = data;
    


end