function npxy = nproxy_square(kerns,width,opts)
%CHNK.FLAM.NPROXY_SQUARE find the number of points to place on a
% side of a proxy surface for given kernel, width, and tolerance
%
% Syntax: npxy = chnk.flam.nproxy_square(kern,width,eps)
%
% Input:
%   kern - kernel function of the form kern(srcinfo,targinfo)
%   width - box/cube width (of center box, proxy surface at
%                1.5*width)
%   opts - options structure
%       opts.rank_or_tol - tolerance or maximum rank (default 1e-13)
%       opts.eps - alias for rank_or_tol for backward compat (default 1e-13)
%       opts.nsrc - number of sources to use in test (default 200)  
%
% Output:
%   npxy - number of points on perimeter of box to be sent to proxy routine
%

  nsrc = 400;
  rank_or_tol = 1e-13;

  if nargin < 3
    opts = [];
  end

  %if isfield(opts,'eps')
  %  rank_or_tol = opts.eps;
  %end
  %if isfield(opts,'rank_or_tol')
  %  rank_or_tol = opts.rank_or_tol;
  %end
  %if isfield(opts,'nsrc')
  %  nsrc = opts.nsrc;
  %end

		      % set up sources with randomly oriented d and d2

  srcinfo = [];
  srcinfo.r = [-0.5;-0.5]*width + rand(2,nsrc)*width;
  
  srcinfo.d = randn(2,nsrc);
  srcinfo.n = randn(2,nsrc); srcinfo.n = ... 
      (srcinfo.n)./(sum((srcinfo.n).^2,1));
  srcinfo.d2 = randn(2,nsrc);

  npxy = 16;
  numel(srcinfo.r)
				% double until you get enough
  for i = 1:11
    [pr,ptau] = chnk.flam.proxy_square_pts(npxy);
    targinfo.r = pr*width;
    targinfo.d = ptau;
    targinfo.d2 = zeros(2,npxy);
    targinfo.n = [-ptau(2,:);ptau(1,:)] ./ sqrt(sum(ptau.^2,1));
    
    if length(kerns) == 1
      mat = kerns(srcinfo,targinfo);
      mat = permute(mat,[2,1,3]);
      nr = size(mat,1);
      nc = numel(mat)/nr;
      mat = reshape(mat,[nr,nc]);
    else
      mat = cell(length(kerns),1);
      
      for j = 1:length(kerns)
        mat_j = kerns{j}(srcinfo,targinfo);
        mat_j = permute(mat_j,[2,1,3]);
        nr = size(mat_j,1);
        nc = numel(mat_j)/nr;
        mat_j = reshape(mat_j,[nr,nc]);
        mat{j} = mat_j;
      end
      mat = cat(1,mat{:});
    end

    [sk,~] = id(mat,rank_or_tol);
    
    if length(sk) < min(nsrc,npxy)
      npxy = (floor((length(sk)-1)/4+0.1)+1)*4;
      break;
    end
    npxy = 2*npxy;
  end
  
end
