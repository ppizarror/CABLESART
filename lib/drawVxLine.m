function obj = drawVxLine(x, style, lw, color)
% drawVxLine: This function draw a vertical line on x value.
%
% Author: Pablo Pizarro @ppizarror.com, 2017.
%
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

if ~exist('lw', 'var')
    lw = 1;
end
y = get(gca, 'ylim');
if isempty(y) || length(y) == 1
    error('No hay limite en y');
end
if ~exist('color', 'var')
    obj = plot([x, x], y, style, 'Linewidth', lw);
else
    obj = plot([x, x], y, style, 'Linewidth', lw, 'color', color);
end

end % drawVxLine function