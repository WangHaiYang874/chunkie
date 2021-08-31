function in = chunkerinterior(chnkr,pts,opts)
%CHUNKERINTERIOR returns an array indicating whether each point in
% pts is inside the domain. Assumes the domain is closed.
%
% Syntax: in = chunkerinterior(chnkr,pts,opts)
%
% Input:
%   chnkr - chunker object describing geometry
%   pts - (chnkr.dim,:) array of points to test
%
% Optional input:
%   opts - options structure with entries:
%       opts.flam = boolean, use FLAM routines (true)
%
% Output:
%   in - logical array, if in(i) is true, then pts(:,i) is inside the
%       domain
%
% Examples:
%   chnkr = chunkerfunc(@(t) starfish(t));
%   pts = 2*randn(2,100);
%   in = chunkerinterior(chnkr,pts);
%

% author: Travis Askham (askhamwhat@gmail.com)

assert(chnkr.dim == 2,'interior only well-defined for 2D');

if nargin < 3
    opts = [];
end

useflam = true;
if isfield(opts,'flam')
    useflam = opts.flam;
end

kernd = @(s,t) chnk.lap2d.kern(s,t,'d');
dens1 = ones(chnkr.k,chnkr.nch);
wts = weights(chnkr);

opdims = [1 1];

if useflam
    xflam1 = chnkr.r(:,:);
    matfun = @(i,j) chnk.flam.kernbyindexr(i,j,pts,chnkr,wts,kernd,opdims);
    [pr,ptau,pw,pin] = chnk.flam.proxy_square_pts();

    pxyfun = @(rc,rx,cx,slf,nbr,l,ctr) chnk.flam.proxyfunr(rc,rx,slf,nbr,l, ...
        ctr,chnkr,wts,kernd,opdims,pr,ptau,pw,pin);
    F = ifmm(matfun,pts,xflam1,200,1e-14,pxyfun);
    vals1 = ifmm_mv(F,dens1(:),matfun);
else
    optskerneval = []; optskerneval.usesmooth = 1;
    vals1 = chunkerkerneval(chnkr,kernd,dens1,pts,optskerneval);
end

in = abs(vals1+1) < abs(vals1);

% find nearest neighbors at certain level of refinement (here chosen
% uniformly, for simplicity)

nt = size(pts,2);
xx = zeros(2,chnkr.npt + nt);
xx(:,1:chnkr.npt) = chnkr.r(:,:);
xx(:,chnkr.npt+1:chnkr.npt+nt) = pts;

pt2chnk = repmat(1:chnkr.nch,chnkr.k,1);

chunklens = zeros(chnkr.nch,1);
ws = weights(chnkr);
chunklens(:) = sum(ws,1);
lmax = max(chunklens)*3/chnkr.k;

T = hypoct_uni(xx,lmax);
targ_dists = Inf(nt,1);
itarg_dists = (1:nt).';

rnorm = normals(chnkr);
normdist_flag = false(size(targ_dists));

for i = 1:length(T.nodes)
    xi = T.nodes(i).xi;
    if nnz(xi <= chnkr.npt) > 0
        isrc = xi(xi <= chnkr.npt);
        inbor = [T.nodes(T.nodes(i).nbor).xi T.nodes(i).xi];
        if nnz(inbor > chnkr.npt) > 0
            
            % get all pairwise distances
            itarg = inbor(inbor > chnkr.npt)-chnkr.npt;
            srcs = reshape(chnkr.r(:,isrc),chnkr.dim,1,length(isrc));
            targs = reshape(pts(:,itarg),chnkr.dim,length(itarg),1);
            dists = reshape(sqrt(sum((bsxfun(@minus,targs,srcs)).^2,1)), ...
                length(itarg),length(isrc));
            [mindists,inds] = min(dists,[],2);
            isrc2 = isrc(inds);
            
            % if minimum distance for any targ is smaller than
            % previously seen, we update
            ifnew = mindists < targ_dists(itarg);
            targ_dists(itarg(ifnew)) = mindists(ifnew);
            itarg_dists(itarg(ifnew)) = isrc2(ifnew);
            
            % if vector from point to boundary point aligns 
            % with normal, probably on inside.
            % if angle is near right angle or distance is small,
            % flag it for more precise testing
            rdiff = chnkr.r(:,isrc2(ifnew)) - pts(:,itarg(ifnew));
            rnorms = rnorm(:,isrc2(ifnew));
            lens = chunklens(pt2chnk(isrc2(ifnew)));
            
            rdrn = sum(rdiff.*rnorms,1); rdrn = rdrn(:);
            lens = lens(:);
            in(itarg(ifnew)) = rdrn > 0;
            normdist_flag(itarg(ifnew)) = abs(rdrn) < 0.1*lens;            
        end
    end
end

% more precise testing

iflagged = find(normdist_flag);
[~,~,u] = lege.exps(chnkr.k);
tover = lege.exps(2*chnkr.k);
ainterpover = lege.matrin(chnkr.k,tover);
for i = 1:length(iflagged)
    ii = iflagged(i);
    ich = pt2chnk(itarg_dists(ii));
    iadj = chnkr.adj(:,ich);
    adj1 = iadj(iadj > 0);
    iverts = -iadj(iadj < 0);
    adj2 = [chnkr.vert{iverts}];
    ich2 = [ich; adj1(:); adj2(:)];
    [rn,dn] = nearest(chnkr,pts(:,ii),ich2,[],u,tover,ainterpover);
    in(ii) = (rn(:)-pts(:,ii)).'*[dn(2);-dn(1)] > 0;
end




end