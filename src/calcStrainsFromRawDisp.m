function [X,Y,straindata] = calcStrainsFromRawDisp(data,strain_filter)
%%% INPUT: dic is the struct with the disp and pixel position data, inc is
%%% the increment to calc strains, strain_filter is the smoothing used for
%%% before strain calculation

%%% written by Tijmen Vermeij, TUe, 2019-2024


%%% should have nans already in place

%%% get inc data from dic struct
% take the data from the required increment
X = data.X;
Y = data.Y;

U = data.U;
V = data.V;

% keep track of nan locations
nanU = isnan(U);
nanV = isnan(V);
%%% define filtering option

% filter displacement data to reproduce typical strain filter smoothing, see
% function for more details.

% 'strain window size', specifically filter size, defined in 'data points' (i.e. in terms of the
% stepsize used for DIC)

% present example

% define filter options:

% cutofffraction  , e.g. "90% center-weighted Gaussian filter" used in
% some DIC packages (use 0.9 in this case), basically, the gaussian filter is cutoff
% after decaying 90percent.
options.cutofffraction = 1; 


% define std of filter, only used for option 1
options.filt_std = strain_filter;

%%% execute filtering
data = filterDisplacementsTijmen(U,V,options);

%%% calc strains from (filtered) disp data
data = calcEfromU(X,Y,data.U,data.V);

Exx = data.Exx;
Exy = data.Exy;
Eyy = data.Eyy;

Fxx = data.Fxx;
Fxy = data.Fxy;
Fyx = data.Fyx;
Fyy = data.Fyy;

% put NaN's back in strain field when required
Exx(nanU)=NaN; Eyy(nanU)=NaN; Exy(nanU)=NaN; 

Fxx(nanU)=NaN;
Fxy(nanU)=NaN;
Fyx(nanU)=NaN;
Fyy(nanU)=NaN;


straincutoff = 5;
% remove all unthrusworthy data, e.g. strains above 500%
Exx(Exx<-straincutoff|Exx>straincutoff) = NaN; 
Eyy(Eyy<-straincutoff|Eyy>straincutoff) = NaN; 
Exy(Exy<-straincutoff|Exy>straincutoff) = NaN; 


% equivalent strain, assuming no strains in z-direction
Eequi = sqrt(2)/3 * (((Exx-Eyy).^2+(Exx).^2+(Eyy).^2) + 6*Exy.^2).^(0.5);


% store filtered deformation and strain data
straindata.Exx = Exx;
straindata.Exy = Exy;
straindata.Eyy = Eyy;

straindata.Fxx = Fxx;
straindata.Fxy = Fxy;
straindata.Fyx = Fyx;
straindata.Fyy = Fyy;


straindata.Eequi = Eequi;

end



function [ data ] = filterDisplacementsTijmen(U,V,options)
% filter displacement data to reproduce the usual strain filter smoothing

% created by Tijmen Vermeij & Johan Hoefnagels

%%% create filter



if ~isfield(options,'cutofffraction')
    options.cutofffraction = 1;
end
if ~isfield(options,'filt_std')
    options.filt_std = 3;
    warning('no filter standard deviation given, using 5 pixels');
end

filt_std = options.filt_std;
if options.cutofffraction < 1
    options.cutofffraction = 1;
end
window = ceil(15*filt_std) + 1 - rem(ceil(15*filt_std),2);
imageFilter=fspecial('gaussian',window ,filt_std);



%%% perform filtering

% find nans
nanU = isnan(U);
nanV = isnan(V);

% blur displacement fields (use nanconv script which can handle nan values well)
if options.filt_std > 0
    Uf = nanconv(U,imageFilter, 'nonanout','edge');
    Vf = nanconv(V,imageFilter, 'nonanout','edge');
else
    Uf = U;
    Vf = V;
end

% place NaNs at correct positions
Uf(nanU) = NaN;
Vf(nanV) = NaN;


data.imageFilter = imageFilter;
data.U = Uf;
data.V = Vf;
end

function data = calcEfromU(X,Y,U,V)
%function to calculate green lagrange strain field from a (filtered) displacement field
%   input: X, Y, U and V in matrix form.
%   output: data struct with def grad tensor components and strain components   

%  created by Tijmen Vermeij, as per GDIC GUI by Jan Neggers

% determine spacing of position fields
dx = mean(mean(diff(X,1,2)));
dy = mean(mean(diff(Y,1,1))); % multiply this value with -1 to get equal shear strains as VIC-2D gives, not sure if this is correct


% calculate numerical gradients
[duxdx, duxdy] = gradient(U,dx,dy);
[duydx, duydy] = gradient(V,dx,dy);

% calculate deformation gradient tensor components
Fxx = duxdx + 1;
Fyy = duydy + 1;
Fxy = duxdy;
Fyx = duydx;


% calculate green Lagrange Strain tensor: E = 1/2 ( C - I ), with C = F^T * F
Exx = 0.5 * ( Fxx.*Fxx + Fyx.*Fyx - 1);
Eyy = 0.5 * ( Fxy.*Fxy + Fyy.*Fyy - 1);
Exy = 0.5 * ( Fxx.*Fxy + Fyx.*Fyy );


data.Fxx = Fxx;
data.Fxy = Fxy;
data.Fyx = Fyx;
data.Fyy = Fyy;

data.Exx = Exx;
data.Eyy = Eyy;
data.Exy = Exy;


end
