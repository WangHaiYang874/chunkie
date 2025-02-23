function obj = lap2d(type, eta)
%KERNEL.LAP2D   Construct the Laplace kernel.
%   KERNEL.LAP2D('s') or KERNEL.LAP2D('single') constructs the single-layer
%   Laplace kernel.
%
%   KERNEL.LAP2D('d') or KERNEL.LAP2D('double') constructs the double-layer
%   Laplace kernel.
%
%   KERNEL.LAP2D('sp') or KERNEL.LAP2D('sprime') constructs the derivative
%   of the single-layer Laplace kernel.
%
%   KERNEL.LAP2D('c', ETA) or KERNEL.LAP2D('combined', ETA) constructs the
%   combined-layer Laplace kernel with parameter ETA, i.e.,
%   KERNEL.LAP2D('d') + eta*KERNEL.LAP2D('s').
%
% See also CHNK.LAP2D.KERN.

if ( nargin < 1 )
    error('Missing Laplace kernel type.');
end

obj = kernel();
obj.name = 'laplace';
obj.opdims = [1 1];

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 's');
        obj.fmm  = @(eps,s,t,sigma,pgt) chnk.lap2d.fmm(eps, s, t, 's', sigma, pgt);
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'd');
        obj.fmm  = @(eps,s,t,sigma,pgt) chnk.lap2d.fmm(eps, s, t, 'd', sigma, pgt);
        obj.sing = 'smooth';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'sprime');
        % TODO: FMM for sprime.
        obj.sing = 'smooth';

    case {'st', 'stau'}
        obj.type = 'st';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'stau');
        % TODO: FMM for stau.
        obj.sing = 'pv';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'dprime');
        % TODO: FMM for sprime.
        obj.sing = 'hs';

    case {'sg', 'sgrad'}
        obj.type = 'sg';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'sgrad');
        obj.fmm = @(eps,s,t,sigma,pgt) chnk.lap2d.fmm(eps, s, t, 'sgrad', sigma, pgt);
        obj.sing = 'pv';
        obj.opdims = [2,1];

    case {'dg', 'dgrad'}
        obj.type = 'dg';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'dgrad');
        obj.fmm = @(eps,s,t,sigma,pgt) chnk.lap2d.fmm(eps, s, t, 'dgrad', sigma, pgt);
        obj.sing = 'hs';
        obj.opdims = [2,1];

    case {'c', 'combined'}
        if ( nargin < 2 )
            warning('Missing combined layer parameter eta. Defaulting to 1.');
            eta = 1;
        end
        obj.type = 'c';
        obj.params.eta = eta;
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'c', eta);
        obj.fmm  = @(eps,s,t,sigma,pgt) chnk.lap2d.fmm(eps, s, t, 'c', sigma, pgt, eta);
        obj.sing = 'log';

    otherwise
        error('Unknown Laplace kernel type ''%s''.', type);

end

end
