
function submat= kern(zk,srcinfo,targinfo,type,varargin)
%CHNK.HELM2D.KERN standard Helmholtz layer potential kernels in 2D
% 
% Syntax: submat = chnk.heml2d.kern(zk,srcinfo,targingo,type,varargin)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%  
% Kernels based on G(x,y) = i/4 H_0^{(1)}(zk |x-y|)
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% S'(x,y) = \nabla_{n_x} G(x,y)
%
% Input:
%   zk - complex number, Helmholtz wave number
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   targinfo - description of targets in ptinfo struct format,
%                if info not relevant (d/d2) it doesn't need to
%                be provided. sprime requires tangent info in
%                targinfo.d
%   type - string, determines kernel type
%                type == 'd', double layer kernel D
%                type == 's', single layer kernel S
%                type == 'sprime', normal derivative of single
%                      layer S'
%                type == 'c', combined layer kernel D + i eta S
%   varargin{1} - eta in the combined layer formula, otherwise
%                does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% see also CHNK.HELM2D.GREEN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

if strcmpi(type,'d')
  srcnorm = srcinfo.n;
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'sprime')
  targnorm = targinfo.n;
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat((targnorm(1,:)).',1,ns);
  ny = repmat((targnorm(2,:)).',1,ns);

  submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
end

if strcmpi(type,'sdtau')
  targtan = targinfo.d;
  [~,grad] = chnk.helm2d.green(zk,src,targ);
  dx = repmat((targtan(1,:)).',1,ns);
  dy = repmat((targtan(2,:)).',1,ns);
  dn = sqrt(dx.*dx+dy.*dy);
  submat = (grad(:,:,1).*dx./dn + grad(:,:,2).*dy)./dn;
end

if strcmpi(type,'s')
  submat = chnk.helm2d.green(zk,src,targ);
end

if strcmpi(type,'dprime')
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  [~,~,hess] = chnk.helm2d.green(zk,src,targ);
  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
end

if strcmpi(type,'c')
  srcnorm = srcinfo.n;
  coef = varargin{1};
  [submats,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
  submat = coef(1)*submatd + coef(2)*submats;
end

if strcmpi(type,'all')
  %targnorm = chnk.normal2d(targinfo);
  %srcnorm = chnk.normal2d(srcinfo);
  % added by Shidong Jiang to avoid O(N^2) calculation of normals
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  cc = varargin{1};
  
  submat = zeros(2*nt,2*ns);
  % S
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);
  %[submat(1:2:2*nt,2:2:2*ns),grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  % D'
  %submat(2:2:2*nt,1:2:2*ns) = -(hess(:,:,1).*nxsrc.*nxtarg ... 
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
%  submat(2:2:2*nt,2:2:2*ns) = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
%  submat(1:2:2*nt,1:2:2*ns)  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  
  
  submat(1:2:2*nt,1:2:2*ns) = submatd*cc(1,1);
  submat(1:2:2*nt,2:2:2*ns) = submats*cc(1,2);
  submat(2:2:2*nt,1:2:2*ns) = submatdp*cc(2,1);
  submat(2:2:2*nt,2:2:2*ns) = submatsp*cc(2,2);
end

if strcmpi(type,'eval trans')

  srcnorm = srcinfo.n;
  [submats,grad] = chnk.helm2d.green(zk,src,targ);
  nx = repmat(srcnorm(1,:),nt,1);
  ny = repmat(srcnorm(2,:),nt,1);
  submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);

  submat = zeros(1*nt,2*ns);
  submat(1:1:1*nt,1:2:2*ns) = submatd;
  submat(1:1:1*nt,2:2:2*ns) = submats;
end

if strcmpi(type,'c and cprime')
  targnorm = targinfo.n;
  srcnorm = srcinfo.n;
  coef = varargin{1};

  submat = zeros(2*nt,1*ns);
  % S
  [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);
  %[submat(1:2:2*nt,2:2:2*ns),grad,hess] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  nxtarg = repmat((targnorm(1,:)).',1,ns);
  nytarg = repmat((targnorm(2,:)).',1,ns);
  % D'
  %submat(2:2:2*nt,1:2:2*ns) = -(hess(:,:,1).*nxsrc.*nxtarg ... 
  submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
      + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);
  % S'
%  submat(2:2:2*nt,2:2:2*ns) = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
  % D
%  submat(1:2:2*nt,1:2:2*ns)  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
  
  submatc      = coef(1,1) * submatd  + coef(1,2) * submats;
  submatcprime = coef(2,1) * submatdp + coef(2,2) * submatsp;
  
  submat(1:2:2*nt,1:1:1*ns) = submatc;
  submat(2:2:2*nt,1:1:1*ns) = submatcprime;
end

if strcmpi(type,'eval')
  coef = varargin{1};
  
  srcnorm = srcinfo.n;
  
  submat = zeros(nt,2*ns);
  % S
  [submats,grad] = chnk.helm2d.green(zk,src,targ);

  nxsrc = repmat(srcnorm(1,:),nt,1);
  nysrc = repmat(srcnorm(2,:),nt,1);
  % D
  submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
    
  submat(:,1:2:2*ns) = coef(1)*submatd;
  submat(:,2:2:2*ns) = coef(2)*submats;
end


if strcmpi(type,'evalg')
  coef = varargin{1};
  
  [~,ns] = size(src);
  [~,nt] = size(targ);
  
  if strcmpi(type,'d')
    srcnorm = srcinfo.n;
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
  end
  
  if strcmpi(type,'sprime')
    targnorm = targinfo.n;
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);
  
    submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);
  end
  
  if strcmpi(type,'sdtau')
    targtan = targinfo.d;
    [~,grad] = chnk.helm2d.green(zk,src,targ);
    dx = repmat((targtan(1,:)).',1,ns);
    dy = repmat((targtan(2,:)).',1,ns);
    dn = sqrt(dx.*dx+dy.*dy);
    submat = (grad(:,:,1).*dx./dn + grad(:,:,2).*dy)./dn;
  end
  
  if strcmpi(type,'s')
    submat = chnk.helm2d.green(zk,src,targ);
  end
  
  if strcmpi(type,'dprime')
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    [~,~,hess] = chnk.helm2d.green(zk,src,targ);
    nxsrc = repmat(srcnorm(1,:),nt,1);
    nysrc = repmat(srcnorm(2,:),nt,1);
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
        + hess(:,:,3).*nysrc.*nytarg);
  end
  
  if strcmpi(type,'c')
    srcnorm = srcinfo.n;
    eta = varargin{1};
    [s,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    submat = -(grad(:,:,1).*nx + grad(:,:,2).*ny) + 1i*eta*s;
  end
  
  if strcmpi(type,'all')
    %targnorm = chnk.normal2d(targinfo);
    %srcnorm = chnk.normal2d(srcinfo);
    % added by Shidong Jiang to avoid O(N^2) calculation of normals
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    cc = varargin{1};
    
    submat = zeros(2*nt,2*ns);
    % S
    [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);
    %[submat(1:2:2*nt,2:2:2*ns),grad,hess] = chnk.helm2d.green(zk,src,targ);
  
    nxsrc = repmat(srcnorm(1,:),nt,1);
    nysrc = repmat(srcnorm(2,:),nt,1);
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    % D'
    %submat(2:2:2*nt,1:2:2*ns) = -(hess(:,:,1).*nxsrc.*nxtarg ... 
    submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
        + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
        + hess(:,:,3).*nysrc.*nytarg);
    % S'
  %  submat(2:2:2*nt,2:2:2*ns) = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
    submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
    % D
  %  submat(1:2:2*nt,1:2:2*ns)  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
    submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
    
    
    submat(1:2:2*nt,1:2:2*ns) = submatd*cc(1,1);
    submat(1:2:2*nt,2:2:2*ns) = submats*cc(1,2);
    submat(2:2:2*nt,1:2:2*ns) = submatdp*cc(2,1);
    submat(2:2:2*nt,2:2:2*ns) = submatsp*cc(2,2);
  end

  if strcmpi(type,'cprime')
    %targnorm = chnk.normal2d(targinfo);
    %srcnorm = chnk.normal2d(srcinfo);
    % added by Shidong Jiang to avoid O(N^2) calculation of normals
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    eta = varargin{1};
    
    submat = zeros(2*nt,2*ns);
    % S
    [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);
    %[submat(1:2:2*nt,2:2:2*ns),grad,hess] = chnk.helm2d.green(zk,src,targ);
  
    nxsrc = repmat(srcnorm(1,:),nt,1);
    nysrc = repmat(srcnorm(2,:),nt,1);
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    % D'
    %submat(2:2:2*nt,1:2:2*ns) = -(hess(:,:,1).*nxsrc.*nxtarg ... 
    submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
        + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
        + hess(:,:,3).*nysrc.*nytarg);
    % S'
  %  submat(2:2:2*nt,2:2:2*ns) = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
    submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
    % D
  %  submat(1:2:2*nt,1:2:2*ns)  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
    submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
    
    submat = submatdp + 1i*eta*submatsp;
    
  end

  if strcmpi(type,'eval trans')
  
    srcnorm = srcinfo.n;
    [submats,grad] = chnk.helm2d.green(zk,src,targ);
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    submatd = -(grad(:,:,1).*nx + grad(:,:,2).*ny);
  
    submat = zeros(1*nt,2*ns);
    submat(1:1:1*nt,1:2:2*ns) = submatd;
    submat(1:1:1*nt,2:2:2*ns) = submats;
  end
  
  if strcmpi(type,'c and cprime')
    % opdim = (2,1)
    %targnorm = chnk.normal2d(targinfo);
    %srcnorm = chnk.normal2d(srcinfo);
    % added by Shidong Jiang to avoid O(N^2) calculation of normals
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    eta = varargin{1};
    cc = varargin{2};

    submat = zeros(2*nt,1*ns);

    % S
    [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);
  
    nxsrc = repmat(srcnorm(1,:),nt,1);
    nysrc = repmat(srcnorm(2,:),nt,1);
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    % D'
    submatdp = -(hess(:,:,1).*nxsrc.*nxtarg ...
        + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
        + hess(:,:,3).*nysrc.*nytarg);
    % S'
    submatsp = (grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);
    % D
    submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);
    
    submatc = submatd + 1i*eta*submats;
    submatcp = submatdp + 1i*eta*submatsp;

    submat(1:2:2*nt,1:1:1*ns) = cc(1) * submatc;
    submat(2:2:2*nt,1:1:1*ns) = cc(2) * submatcp;
  end
  
  if strcmpi(type,'eval')
    coef = varargin{1};
    srcnorm = srcinfo.n;
    
    submat = zeros(nt,2*ns);
    % S
    [submats,grad] = chnk.helm2d.green(zk,src,targ);
  
    nxsrc = repmat(srcnorm(1,:),nt,1);
    nysrc = repmat(srcnorm(2,:),nt,1);
    % D
    submatd  = -(grad(:,:,1).*nxsrc + grad(:,:,2).*nysrc);

    submat(:,1:2:2*ns) = coef(1)*submatd;
    submat(:,2:2:2*ns) = coef(2)*submats;
  end
  
  
  if strcmpi(type,'evalg')
    coef = varargin{1};
    
    srcnorm = srcinfo.n;
    
    submat = zeros(nt,ns,6);
    % S
    [submats,grad,hess] = chnk.helm2d.green(zk,src,targ);
  
  submat(1:2:2*nt,1:2:2*ns) = alpha1.*(c2.*submatd2-c1.*submatd1);
  submat(1:2:2*nt,2:2:2*ns) = alpha1.*(submats2-submats1);
  submat(2:2:2*nt,1:2:2*ns) = -alpha2.*(submatdp2-submatdp1);
  submat(2:2:2*nt,2:2:2*ns) = -alpha2.*(1./c2.*submatsp2-1./c1.*submatsp1);
end
