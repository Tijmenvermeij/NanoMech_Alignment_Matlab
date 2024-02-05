function newEBSD = warpCreateEBSDdata(oldEBSD,Xnew,Ynew,alignmentData)
%%% this function creates a new Mtex EBSD variable/struct by warping an old one using alignment data created from point selection.
%%% Xnew and Ynew are the aimed grids in the target configuration (i.e.
%%% reference), align_data contains the warping information

%%% author: Tijmen Vermeij, TUe, 2018-2023

if any(size(oldEBSD.prop.x)==1)
    error('ebsd input data must be gridded')
end


if nargin == 4
    % get relevant displacement fields by interpolating from thez original to the new grid of positions in
    % the ref config
    Ux = interp2(alignmentData.X,alignmentData.Y,alignmentData.Ux,Xnew,Ynew);
    Uy = interp2(alignmentData.X,alignmentData.Y,alignmentData.Uy,Xnew,Ynew);
else
    warning('no alignmentData given, so no warping is performed, but only interpolation.')
    Ux = zeros(size(Xnew));
    Uy = zeros(size(Xnew));
end

% use MTEX interp function to interpolate the old EBSD data to the new
newEBSD = interp(oldEBSD,Xnew+Ux,Ynew+Uy,'method','nearest');
newEBSD.prop.x = Xnew(:);
newEBSD.prop.y = Ynew(:);

% create new unitCell variable, to be sure.
newEBSD.unitCell = calcUnitCell([Xnew(:) - mean(Xnew(:)),Ynew(:) - mean(Ynew(:))]);

end