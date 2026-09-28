function s = pStars(p, levels)
% pStars  The asterisk for a p value.
%
%   s = pStars(p)           uses alpha = 0.05
%   s = pStars(p, levels)   levels(1) is alpha; any further entries are ignored
%
%   '*'  p < alpha
%   ''   otherwise, and for NaN or a missing p
%
% One asterisk, one threshold, used by every figure.
%
% levels may be a vector (e.g. [alpha 0.01 0.001]); only the first entry is used.

    if nargin < 2 || isempty(levels), levels = 0.05; end

    alpha = levels(1);

    if isempty(p) || ~isfinite(p) || p >= alpha
        s = '';
    else
        s = '*';
    end
end
