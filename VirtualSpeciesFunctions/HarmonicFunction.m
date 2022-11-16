function Harmonic = HarmonicFunction(samples, limdef, plotting)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

if nargin <2
    limdef = 3;
end

if nargin<3
    plotting = false;
end

tic

harm = Dalfa(1:limdef);%creación de armónicos hasta 2 grado
ddefM = cell(1,limdef);
ddefM{1} = eye(3);
for i = 2:limdef
    len = size(harm{i},2);
    ddefM{i} = rand(3,len);
end
gN1 = mapping(harm,ddefM);%Mapeo de deformación
vars1 = symvar(gN1);%Variables del mapeo
assume(vars1(1),'Real')
assume(vars1(2),'Real')
syms t
assume(t,'Real')
if plotting
    figure(1)
    fsurf(symfun(gN1(1),vars1),symfun(gN1(2),vars1),symfun(gN1(3),vars1),[0 1 0 2*pi])
    savefig('Harmonic.fig')
    %hold on
    r = t*samples(1,1)+(1-t)*samples(2,1);
    theta = t*samples(1,2)+(1-t)*samples(2,2);
    gN = subs(gN1,{vars1(1),vars1(2)},{r,theta});
    vars2 = symvar(gN);
    %fplot3(symfun(gN(1),vars2),symfun(gN(2),vars2),symfun(gN(3),vars2),[0 1],'--or','LineWidth',2,'MarkerFaceColor','auto')
end

funr = simplify(diff(gN1,vars1(1)));
funtheta = simplify(diff(gN1,vars1(2)));
g11 = simplify(dot(funr,funr));
g12 = simplify(dot(funr,funtheta));
g22 = simplify(dot(funtheta,funtheta));

r = @(t) t*samples(1,1)+(1-t)*samples(:,1);
theta = @(t) t*samples(1,2)+(1-t)*samples(:,2);
dr = samples(1,1)-samples(:,1);
dtheta = samples(1,2)-samples(:,2);
fg11 = symfun(g11,[vars1(1),vars1(2)]);
fg11 = matlabFunction(fg11);

fg12 = symfun(g12,[vars1(1),vars1(2)]);
fg12 = matlabFunction(fg12);
fg22 = symfun(g22,[vars1(1),vars1(2)]);
fg22 = matlabFunction(fg22);

prefun = @(t) sqrt(fg11(r(t),theta(t)).*(dr).^2 + 2*fg12(r(t),theta(t)).*dr.*dtheta + fg22(r(t),theta(t)).*(dtheta).^2);
distances = integral(prefun,0,1,'ArrayValued',true);
distances = (distances-max(distances))/(min(distances)-max(distances));

if plotting
    figure(2)
    scatter(samples(:,1),samples(:,2),[],distances,'filled')
    savefig('ColoredSpace.fig')
    colormap jet
end
%OUTPUT STORAGE
Harmonic.distances = distances;

end


function inds = multiind(d,n) %multiindices
    P = partitions(d,1:d,n);
    N = size(P,1);
    C = cell(1,N);
    for k = 1:N
        tmp = repelem(1:d,P(k,:));
        tmp(end+1:n) = 0;
        C{k} = unique(perms(tmp),'rows');
    end
    inds = vertcat(C{:});
end
function res = Dalfa(range)%armónicos para los valores en range
%esta es la función principal para extraer las funciones de deformación
    syms x_1 x_2 x_3
    x = [x_1,x_2,x_3];
    expre = 1/sqrt(x_1^2+x_2^2+x_3^2);
    len2 = length(range);
    res = cell(1,len2);
    
    for j=1:len2
        inds=multiind(range(j),3);
        len = size(inds,1);
        harmon=sym(zeros(1,len));
        for i=1:len
            harmon(i)=harmdiff(expre,inds(i,:),x);
        end
        res{j}=harmon;
    end
end

function expre = harmdiff(expre,inds,x)%armónico asociado a un multiíndice
    expre = diff(expre,x(1),inds(1));
    expre = diff(expre,x(2),inds(2));
    expre = diff(expre,x(3),inds(3));
    expre = -numden(expre);
end

function [hemN,hemS] = hemisferios(expr)%funciones para cada hemisferio
    syms u v r theta
    vars=symvar(expr);
    exprN = subs(expr,vars,{2*u/(u^2+v^2+1),2*v/(u^2+v^2+1),(1-u^2-v^2)/(u^2+v^2+1)});
    exprS = subs(expr,vars,{2*u/(u^2+v^2+1),2*v/(u^2+v^2+1),-(1-u^2-v^2)/(u^2+v^2+1)});
    hemN = simplify(subs(exprN,{u,v},{r*cos(theta),r*sin(theta)}));
    hemS = simplify(subs(exprS,{u,v},{r*cos(theta),r*sin(theta)}));
end

function [hemN,hemS] = mapping(harm,matrix)%creación de las 3 funciones esfera deformada
    mapp=0;
    for i=1:length(matrix)
        mapp=mapp+harm{i}*matrix{i}';
    end
    [hemN,hemS] = hemisferios(mapp);
end

function plist = partitions(total_sum,candidate_set,max_count,fixed_count)

% default for candidate_set
if (nargin<2) || isempty(candidate_set)
  candidate_set = 1:total_sum;
end

% how many candidates are there
n = length(candidate_set);

% error checks
if any(candidate_set<0)
  error('All members of candidate_set must be >= 0')
end
% candidates must be sorted in increasng order
if any(diff(candidate_set)<0)
  error('Efficiency requires that candidate_set be sorted')
end
% check for a max_count. do we supply a default?
if (nargin<3) || isempty(max_count)
  % how high do we need look?
  max_count = floor(total_sum./candidate_set);
elseif length(max_count)==1
  % if a scalar was provided, then turn it into a vector
  max_count = repmat(max_count,1,n);
end

% check for a fixed_count
if (nargin<4) || isempty(fixed_count)
  fixed_count = [];
elseif (fixed_count<0) || (fixed_count~=round(fixed_count))
  error('fixed_count must be a positive integer if supplied')
end

% check for degenerate cases
if isempty(fixed_count)
  if total_sum == 0
    plist = zeros(1,n);
    return
  elseif (n == 0)
    plist = [];
    return
  elseif (n == 1)
    % only one element in the set. can we form
    % total_sum from it as an integer multiple?
    p = total_sum/candidate_set;
    if (p==fix(p)) && (p<=max_count)
      plist = p;
    else
      plist = [];
    end
    return
 end
else
  % there was a fixed_count supplied
  if (total_sum == 0) && (fixed_count == 0)
    plist = zeros(1,n);
    return
  elseif (n == 0) || (fixed_count <= 0)
    plist = [];
    return
  elseif (n==1)
    % there must be a non-zero fixed_count, since
    % we did not trip the last test. since there
    % is only one candidate in the set, will it work?
    if ((fixed_count*candidate_set) == total_sum) && (fixed_count <= max_count)
      plist = fixed_count;
    else
      plist = [];
    end
    return
  end
end

% finally, we can do some work. start with the
% largest element and work backwards
m = max_count(end);
% do we need to back off on m?
c = candidate_set(end);
m = min([m,floor(total_sum/c),fixed_count]);

plist = zeros(0,n);
for i = 0:m
  temp = partitions(total_sum - i*c, ...
      candidate_set(1:(end-1)), ...
      max_count(1:(end-1)),fixed_count-i);
  plist = [plist;[temp,repmat(i,size(temp,1),1)]];  %#ok
end
end
