% modified based on DistMesh by Per-Olof Persson
% Function: create a mesh of a cell and save it as cell.xml 
% h0 controls how refined the mesh is, the smaller h0 the more refined the mesh is
close all
h0 = 25;
pfix =[y_fix',x_fix'];
dptol=.001; ttol=.1; Fscale=1.2; deltat=.2; geps=.001*h0; deps=sqrt(eps)*h0;
densityctrlfreq=30;

% 1. Create initial distribution in bounding box (equilateral triangles)
[x,y]=meshgrid(bbox(1,1):h0:bbox(2,1),bbox(1,2):h0*sqrt(3)/2:bbox(2,2));
x(2:2:end,:)=x(2:2:end,:)+h0/2;                      % Shift even rows
p=[x(:),y(:)];                                       % List of node coordinates

% 2. Remove points outside the region, apply the rejection method
distance = fd(p,boundary_points);
p=p(distance<geps,[1 2]);                 % Keep only d<0 points
h_size = fh(p);
r0=1./h_size.^2;                    % Probability to keep point
p=p(rand(size(p,1),1)<r0./max(r0),:);                % Rejection method
if ~isempty(pfix), p=setdiff(p,pfix,'rows'); end     % Remove duplicated nodes
pfix=unique(pfix,'rows'); nfix=size(pfix,1);
p=[pfix; p];                                         % Prepend fix points
N=size(p,1);                                         % Number of points N

count=0;
pold=inf;                                            % For first iteration
clf,view(2),axis equal,axis off
while 1
  count=count+1;
  % 3. Retriangulation by the Delaunay algorithm
  if max(sqrt(sum((p-pold).^2,2))/h0)>ttol           % Any large movement?
    pold=p;                                          % Save current positions
    t=delaunayn(p);                                  % List of triangles
    pmid=(p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;    % Compute centroids
    distance_p = fd(pmid,boundary_points);
    t=t(distance_p<-geps,[1 2 3]);         % Keep interior triangles
    % 4. Describe each bar by a unique pair of nodes
    bars=[t(:,[1,2]);t(:,[1,3]);t(:,[2,3])];         % Interior bars duplicated
    bars=unique(sort(bars,2),'rows');                % Bars as node pairs
    % 5. Graphical output of the current mesh
    cla,patch('vertices',p,'faces',t,'edgecol','k','facecol',[.8,.9,1]);
    drawnow
  end

  % 6. Move mesh points based on bar lengths L and forces F
  barvec=p(bars(:,1),:)-p(bars(:,2),:);              % List of bar vectors
  L=sqrt(sum(barvec.^2,2));                          % L = Bar lengths
  hbars=fh((p(bars(:,1),:)+p(bars(:,2),:))/2);
%   hbars=feval(fh,(p(bars(:,1),:)+p(bars(:,2),:))/2,varargin{:});
  L0=hbars*Fscale*sqrt(sum(L.^2)/sum(hbars.^2));     % L0 = Desired lengths
  
  % Density control - remove points that are too close
  if mod(count,densityctrlfreq)==0 & any(L0>2*L)
      p(setdiff(reshape(bars(L0>2*L,:),[],1),1:nfix),:)=[];
      N=size(p,1); pold=inf;
      continue;
  end
  
  F=max(L0-L,0);                                     % Bar forces (scalars) 
  Fvec=F./L*[1,1].*barvec;                           % Bar forces (x,y components)
  Ftot=full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
  Ftot(1:size(pfix,1),:)=0;                          % Force = 0 at fixed points
  p=p+deltat*Ftot;                                   % Update node positions

  % 7. Bring outside points back to the boundary
  distance_p = fd(p,boundary_points);
  ix=distance_p>0;                 % Find points outside (d>0)
  dgradx=(fd([p(ix,1)+deps,p(ix,2)],boundary_points)-distance_p(ix))/deps;
%   dgradx=(feval(fd,[p(ix,1)+deps,p(ix,2)],boundary_points,varargin{:})-d(ix))/deps; % Numerical
  dgrady=(fd([p(ix,1),p(ix,2)+deps],boundary_points)-distance_p(ix))/deps;
%   dgrady=(feval(fd,[p(ix,1),p(ix,2)+deps],boundary_points,varargin{:})-d(ix))/deps; %    gradient
  dgrad2=dgradx.^2+dgrady.^2;
  p(ix,:)=p(ix,:)-[distance_p(ix).*dgradx./dgrad2,distance_p(ix).*dgrady./dgrad2];    % Project

  % 8. Termination criterion: All interior nodes move less than dptol (scaled)
  if max(sqrt(sum(deltat*Ftot(distance_p<-geps, [1 2]).^2,2))/h0)<dptol, break; end
 end

% Clean up and plot final mesh
[p,t]=fixmesh(p,t);
simpplot(p,t)
hold on
scatter(y_fix,x_fix,'r')

% save mesh "cell.xml"
vrts = p.*1.613.*10.^(-1);
tets = t;
% vrts = round(p);
% tets = round(t);
% xmlmesh(vrts,tets,'cell.xml')