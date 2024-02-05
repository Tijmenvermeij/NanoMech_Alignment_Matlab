function [Phi dir] = basis(x,y,phi,kdof,options)
% this function calls the correct basis-function type build-function, and
% returns the field for this basis-function and direction

% normalize the domain
m = length(x);
n = length(y);
if options.normalized
    x = linspace(-1,1,m);
    y = linspace(-1,1,n);
end

% store the direction in which the basis operates
dir = phi(kdof,3);

% store the order
a = phi(kdof,1);
b = phi(kdof,2);

% call the correct basis-function
if any(strcmpi(options.basis,{'polynomial','p'}))
    Phi = basis_polynomial(x,y,a,b);
elseif any(strcmpi(options.basis,{'chebyshev','c'}))
    Phi = basis_chebyshev(n,m,a,b);
elseif any(strcmpi(options.basis,{'bspline','b'}))
    % for splines also get the number of (unique) knots and the order
    na = phi(kdof,4);
    nb = phi(kdof,5);
    np = phi(kdof,6);
    Phi = basis_bspline(n,m,a,b,na,nb,np);
else
    error('unknown basis: [%s]',options.basis);
end


% =======================================================================

function Phi = basis_polynomial(x,y,a,b)
% This function creates a polynomial basisfunction
% which is defined at the pixel level

% this basis creates a two dimensional Nth order field of the type
% Phi = x^a * y^b;

% create the shape in each direction
x = x.^a;
y = y.^b;

% expand each direction to a field
[X Y] = meshgrid(x,y);

% multiply the two
Phi = X .* Y;

% Note, the order of the power and meshgrid commands is irrelevant

% =======================================================================

function Phi = basis_chebyshev(n,m,a,b)
% This function creates a chebyshev basisfunction of the first kind
% which is defined at the pixel level. 
% http://mathworld.wolfram.com/ChebyshevPolynomialoftheFirstKind.html

% chebyshev functions only make sense on a normalized base
x = linspace(-1,1,m);
y = linspace(-1,1,n);

% For the x direction
% =====================

% initial the first two chebyshev polynomials
Tn = ones(1,m);
Tk = x;

% Recursively build the polynomials from the first two
if a == 0
    Ta = Tn;
elseif a == 1
    Ta = Tk;
else
    for k = 2:a
        
        % calculate the new polynomial
        Ta = 2*x.*Tk - Tn;
        
        % update the previous two
        Tn = Tk;
        Tk = Ta;
    end
end

% For the y direction
% =====================

% initial the first two chebyshev polynomials
Tn = ones(1,n);
Tk = y;

% Recursively build the polynomials from the first two
if b == 0
    Tb = Tn;
elseif b == 1
    Tb = Tk;
else
    for k = 2:b
        
        % calculate the new polynomial
        Tb = 2*y.*Tk - Tn;
        
        % update the previous two
        Tn = Tk;
        Tk = Tb;
    end
end

% expand each direction to a field
[X Y] = meshgrid(Ta,Tb);

% combine the direction to one basis-function
Phi = X.*Y ;


% =======================================================================


function     Phi = basis_bspline(n,m,a,b,na,nb,np)
% nk = number of (unique) knots
% np = polynomial degree

% forloop for x and y directions
for xy = 1:2
    % space the knots uniformly
    if xy == 1
        nk = na;
    else
        nk = nb;
    end
    knots = linspace(-1,1,nk);
    
    % Make the B-Spline Open: repeat end knots p times
    knots = [knots(1)*ones(1,np) knots knots(end)*ones(1,np)];
    % new number of knots
    Ni = length(knots);
    
    % Initiate the parametric space
    if xy == 1
        Nzeta = m;
    else
        Nzeta = n;
    end
    zeta = linspace(knots(1),knots(end),Nzeta);
    
    % Zero order
    % =============
    
    % degree
    p = 0;
    
    % number of basis-functions
    Nn = Ni-(p+1);
    
    % initiate matrix for basis functions
    N0 = zeros(Nn,Nzeta);
    
    for i = 1:Nn
        if i == Nn-np
            % only for the right most (single) knot
            I = find(zeta >= knots(i) & zeta <= knots(i+1));
        else
            % for all other knots
            I = find(zeta >= knots(i) & zeta < knots(i+1));
        end
        
        % set the function to zero
        N0(i,I) = 1;
    end
    % copy the zero basis for later use
    N = N0;
    
    % Subsequent orders
    % =============
    for p = 1:np
        % calculate the number of shape functions for this degree
        Nn = Ni-(p+1);
        
        % store the previous N as N1
        N1 = N;
        % initiate the current N
        N  = zeros(Nn,Nzeta);
        
        % forloop over the knots
        for i = 1:Nn
            if knots(i+p) == knots(i)
                % if first term == zero (double knot on right side)
                N(i,:) = N1(i+1,:).* (knots(i+p+1) - zeta ) ./ ( knots(i+p+1) - knots(i+1) );
            elseif knots(i+p+1) == knots(i+1)
                % if second term is zero (double knot on left side)
                N(i,:) = N1(i,:)  .* (   zeta - knots(i)   ) ./ (  knots(i+p)   - knots(i)  ) ;
            else
                % for all other knots
                N(i,:) = N1(i,:)  .* (   zeta - knots(i)   ) ./ (  knots(i+p)   - knots(i)  ) + ...
                    N1(i+1,:).* (knots(i+p+1) - zeta ) ./ ( knots(i+p+1) - knots(i+1) );
            end
        end
        
    end
    
    % Select the combination of shapes
    if xy == 1
        Nx = N(a,:);
    else
        Ny = N(b,:);
    end
    
end

% expand each direction to a field
[X Y] = meshgrid(Nx,Ny);

% combine the direction to one basis-function
Phi = X.*Y ;
