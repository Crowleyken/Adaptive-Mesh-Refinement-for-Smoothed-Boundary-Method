function Fx = central_diff( F, x )
% Central Difference Gradient for unevenly spaced univariate data
% in the interior and second-order forward/backward differences 
% at the left/right ends
%
% usage:    gradient = central_diff( F, x )
%
% inputs:   F - Values of a function evaluated at x 
%               to be differentiated with respect to x
%               (matrix: number of rows = length of x, unless x is scalar)
%           x - Monotonically increasing coordinate values 
%               where F is evaluated (vector, length=number rows of F),
%               or dx spacing for evenly spaced coordinates (scalar)
%
% output:  gradient - numerically evaluated gradient by:
%                     forward difference at the left end;
%                     backward difference at the right end;
%                     central difference in the interior
%                     (matrix, same size as F)
%
%  Written by:   Robert A. Canfield
%  email:        bob.canfield@vt.edu
%  Version:      2.0
%
%  Created:      10/19/00
%  Modified:     10/01/15
%
%  Description: The central_diff function calculates a numeric gradient
%  using second-order accurate difference formula for evenly or unevenly
%  spaced coordinate data. It operates in a similar fashion to the MATLAB
%  function, gradient, except that it permits only one independent
%  variable, x, and correctly handles unevenly spaced values of the
%  x-coordinate data. Accuracy is increased at the ends relative to the
%  MATLAB gradient function, which uses only first-order forward or
%  backward differences at the ends, by instead using second-order forward
%  difference at the left end and second-order backward difference at the
%  right end.

%  MATLAB's gradient function is incorrect for unevenly spaced coordinates.
%  This central_diff function uses the correct formula.
%  Tested under MATLAB versions 5.2, 5.3.1, and 8.3.
%  (logical operators & and | replaced with && and || for 8.3)
%
% Alternatively, you may patch MATLAB's gradient function 
% to make unevenly spaced interior points second-order accurate
% (leaving left and right ends first-order accurate)
% by replacing the following lines...
%
%>  % Take centered differences on interior points
%>  if n > 2
%>     h = h(3:n) - h(1:n-2);
%>     g(2:n-1,:) = (f(3:n,:)-f(1:n-2,:))./h(:,ones(p,1));
%>  end
%
% with...
%   % Take centered differences on interior points
%   if n > 2
%      if all(abs(diff(h,2)) < eps) % only use for uniform h (RAC)
%         h = h(3:n) - h(1:n-2);
%         g(2:n-1,:) = (f(3:n,:)-f(1:n-2,:))./h(:,ones(p,1));
%      else   % new formula for un-evenly spaced coordinates (RAC)
%         h = diff(h); h_i=h(1:end-1,ones(p,1)); h_ip1=h(2:end,ones(p,1));
%         g(2:n-1,:) =  (-(h_ip1./h_i).*f(1:n-2,:) + ...
%                         (h_i./h_ip1).*f(3:n,:)   )./ (h_i + h_ip1) + ...
%                         ( 1./h_i - 1./h_ip1 ).*f(2:n-1,:);
%      end
%   end

%--Modifications
%  10/23/00 
%  10/01/01 - Copyright (c) 2001 Robert A. Canfield (BSD License)
%  10/01/15 - Second-order accurate at left and right ends

%% Ensure compatible vectors and x monotonically increasing or decreasing
if nargin<1
   disp('usage:  gradient = central_diff( F, x )')
   return
elseif nargin<2
   m = 1;
   x = 1;
else
   m = length(x);
end

Tflag = 0;
if ismatrix(F) && size(F,1)==1 % Treat row vector as a column vector
   F = F.';
   Tflag = 1;
end;
[n,p] = size(F);
if m==0
   error('x cannot be null')
elseif m~=1 && m~=n
   error('First dimension of F and x must be the same.')
elseif m>1 && ~( all(diff(x)>0) || all(diff(x)<0) )
   error('Vector x must be monotonically increasing or decreasing.')
elseif n<=1
   Fx = F / x;
   return
end

Fx = zeros(size(F));
x  = x(:);

%% Forward difference at left end, and Backward difference at right end
if m>1
   H = x(2) - x(1);
else
   H = x;
end
if n==2 % First-order difference for a single interval with end values
   Fx(1,:) = ( F(2,:) - F(1,:) ) / H;
   Fx(2,:) = Fx(1,:);
else    % Second-order differences
   % Left end forward difference
   if m==1 || abs(diff(x(2:3))-H)<=eps % evenly spaced
      Fx(1,:) = ( [-3, 4, -1]/(2*H)*F(1:3,:) ).';
   else                              % unevenly spaced
      h   = diff(x(1:3));
      hph = sum(h); % h_1 + h+2
      Fx(1,:) = hph/h(1)/h(2)*F(2,:) ...
        - ((2*h(1)+h(2))/h(1)*F(1,:) ...
        +           h(1)/h(2)*F(3,:))/hph;
   end
   % Right end backward difference
   if m==1 || abs(diff(diff(x(m-2:m))))<eps % evenly spaced
      Fx(end,:) = ( [1, -4, 3]/(2*H)*F(m-2:m,:) ).';
   else                                    % unevenly spaced
      h   = diff(x(m-2:m));
      hph = sum(h);
      Fx(end,:) = ( h(2)/h(1)*F(end-2,:) ...
         + (h(1)+2*h(2))/h(2)*F(end,:) )/hph ...
              - hph/h(1)/h(2)*F(end-1,:);
   end
end

%% Central Difference in interior (second-order)
if n > 2
   if m==1 || all(abs(diff(x)-H)<=eps)
      % Evenly spaced formula used in MATLAB's gradient routine
      Fx(2:n-1) = ( F(3:n,:) - F(1:n-2,:) ) / (2*H);
   else
      % Unevenly spaced central difference formula
      h = diff(x); h_i=h(1:m-2,ones(p,1)); h_ip1=h(2:m-1,ones(p,1));
      Fx(2:n-1,:) =  (-(h_ip1./h_i).*F(1:n-2,:) + ...
                       (h_i./h_ip1).*F(3:n,:)   )./ (h_i + h_ip1) + ...
                       ( 1./h_i - 1./h_ip1 ).*F(2:n-1,:);
   end
end
if Tflag, Fx=Fx.'; end
end