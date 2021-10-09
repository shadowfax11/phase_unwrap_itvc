function [map] = integrate_derivatives(dx_map,dy_map,strategy,init_val,init_pos)
%INTEGRATE_DERIVATIVES Given the x- and y-derivative 2D maps, and an
%integration strategy, the output of this function is the integrated map
%from the derivative maps
%   Strategy = 1 implies a strategy of finding a minimum spanning tree
%   based on edge weights equal to the magnitude of the derivative
%   corresponding to that edge
%   Strategy = 2 implies a strategy of calculating the curl (2x2 loop
%   summations of the derivatives) and then removing edges corresponding to
%   a curl value > 1e-3
%   init_val refers to the starting integration constant that is added to
%   the map. If not specified, init_val = 0
%   init_pos refers to the position at which init_val is the integration
%   constant. If not specified, init_pos = [1,1]

% function argument preprocessing and checks
if nargin<5
    init_pos = [1,1];
    if nargin < 4
        init_val = 0;
    end
end
H = size(dx_map,1);
W = size(dy_map,2);
% [H,W] = size(dx_map);
% assert(~any([H,W]~=size(dy_map)));

% create graph instance
v = reshape(1:H*W,H,W);
s1 = v(1:H-1,:); t1 = v(2:H,:);
s2 = v(:,1:W-1); t2 = v(:,2:W);
if strategy==2
    curl = calculate_curl(dx_map,dy_map);
    I = find(abs(curl)>1e-3);
    I(I==sub2ind([H,W],init_pos(1),init_pos(2))) = [];
    [r,c] = ind2sub([H-1,W-1],I);
    Ix = sub2ind([H,W-1],r,c);
    Iy = sub2ind([H-1,W],r,c);
    dx_map(Ix) = inf;
    dy_map(Iy) = inf;
    dx_map(Ix+1) = inf;
    dy_map(Iy+H-1) = inf;
    % refine this maybe ??
end
wt1 = abs(dy_map(1:H-1,:)); wt2 = abs(dx_map(:,1:W-1));
wt = [wt1(:);wt2(:)];
s = [s1(:);s2(:)]; s(isinf(wt)) = [];
t = [t1(:);t2(:)]; t(isinf(wt)) = [];
wt(isinf(wt)) = [];
s(isnan(wt)) = [];
t(isnan(wt)) = [];
wt(isnan(wt)) = [];
G = graph(s, t, wt);

% find minimum spanning tree beginning at init_pos
init_node = sub2ind([H,W],init_pos(1),init_pos(2));
[T,pred] = minspantree(G,'Root',findnode(G,init_node));

% integrate derivatives based on the minimum spanning tree
map = Inf*ones(H,W);
node_ind = pred~=0 & ~isnan(pred);
rootedTree = digraph(pred(node_ind),find(node_ind));
e_new = bfsearch(rootedTree,init_node,'edgetonew');
v_curr = e_new(1,1);
[y,x] = ind2sub([H,W],v_curr);
assert(y==init_pos(1));
assert(x==init_pos(2));
map(y,x) = init_val;
i=1;
while i<=size(e_new,1)
    v_curr = e_new(i,2);
    v_prev = e_new(i,1);
    [yp, xp] = ind2sub([H,W],v_prev);
    [y,x] = ind2sub([H,W],v_curr);
    if (x-xp)==1
        map(y,x) = map(yp,xp) + dx_map(yp,xp);
    end
    if (xp-x)==1
        map(y,x) = map(yp,xp) - dx_map(y,x);
    end
    if (y-yp)==1
        map(y,x) = map(yp,xp) + dy_map(yp,xp);
    end
    if (yp-y)==1
        map(y,x) = map(yp,xp) - dy_map(y,x);
    end
    i = i+1;
end
map(isinf(map)) = NaN;
end

