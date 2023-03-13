function Harmonic = HarmonicFunction(samples, deformations_limit, plotting)

% This function generates a harmonic function, maps it onto a surface, and
% calculates distances from given samples to the surface. 

% INPUTS:
%   samples: a 2x2 matrix of [r,theta] values representing two points in 
%            polar coordinates
%   deformations_limit: an integer representing the degree of harmonic functions to 
%           create, default is 3
%   plotting: a boolean to indicate whether or not to create plots, default is false

% OUTPUTS:
%   Harmonic.distances: an array of calculated distances

% If limdef is not specified, set it to 3
if nargin < 2
    deformations_limit = 3;
end

% If plotting is not specified, set it to false
if nargin<3
    plotting = false;
end

% Create harmonic functions up to degree 2
harm = Dalfa(1:deformations_limit);

% Create cell array to store deformation matrices for each degree
deformation_matrices = cell(1,deformations_limit);
deformation_matrices{1} = eye(3);
for i = 2:deformations_limit
    len = size(harm{i},2);
    deformation_matrices{i} = rand(3,len);
end

% Map deformation
gN1 = mapping(harm,deformation_matrices);

% Get variables used in the mapping
vars1 = symvar(gN1);
assume(vars1(1),'Real')
assume(vars1(2),'Real')
syms t
assume(t,'Real')

if plotting
    % Create surface plot of mapped deformation
    figure(1)
    fsurf(symfun(gN1(1),vars1),symfun(gN1(2),vars1),symfun(gN1(3),vars1),[0 1 0 2*pi])
    savefig('Harmonic.fig')
    %hold on
    % Create a line plot of the samples in the same plot
    r = t*samples(1,1)+(1-t)*samples(2,1);
    theta = t*samples(1,2)+(1-t)*samples(2,2);
    gN = subs(gN1,{vars1(1),vars1(2)},{r,theta});
    vars2 = symvar(gN);
    %fplot3(symfun(gN(1),vars2),symfun(gN(2),vars2),symfun(gN(3),vars2),[0 1],'--or','LineWidth',2,'MarkerFaceColor','auto')
end

% Calculate derivatives and metric tensor components
funr = simplify(diff(gN1,vars1(1)));
funtheta = simplify(diff(gN1,vars1(2)));
g11 = simplify(dot(funr,funr));
g12 = simplify(dot(funr,funtheta));
g22 = simplify(dot(funtheta,funtheta));

% Define functions for r and theta
r = @(t) t*samples(1,1)+(1-t)*samples(:,1);
theta = @(t) t*samples(1,2)+(1-t)*samples(:,2);

% Calculate differences in r and theta
dr = samples(1,1)-samples(:,1);
dtheta = samples(1,2)-samples(:,2);

% Convert symbolic expressions to anonymous functions for integration
fg11 = symfun(g11,[vars1(1),vars1(2)]);
fg11 = matlabFunction(fg11);
fg12 = symfun(g12,[vars1(1),vars1(2)]);
fg12 = matlabFunction(fg12);
fg22 = symfun(g22,[vars1(1),vars1(2)]);
fg22 = matlabFunction(fg22);

% Define a function handle for the integrand
prefun = @(t) sqrt(fg11(r(t),theta(t)).*(dr).^2 + 2*fg12(r(t),theta(t)).*dr.*dtheta + fg22(r(t),theta(t)).*(dtheta).^2);

% Integrate the function over [0,1] to get the distances between samples
distances = integral(prefun,0,1,'ArrayValued',true);

% Normalize the distances between [0,1] and store the result
distances = (distances-max(distances))/(min(distances)-max(distances));

% Plot the results if plotting is enabled
if plotting
    % Create a scatter plot of the samples with different colors based on the calculated distances
    figure(2)
    scatter(samples(:,1),samples(:,2),[],distances,'filled')
    colormap jet
    savefig('ColoredSpace.fig')
end

% Store the distances into a structure
Harmonic.Distances = distances;

end

function inds = multiind(d,n)
% This function generates all possible indices of a multidimensional array of size d with n elements in each dimension.
    % Generate all partitions of the set {1,2,...,d} into n nonempty subsets.
    P = partitions(d,1:d,n);
    
    % Determine the number of partitions generated.
    N = size(P,1);
    
    % Initialize a cell array to store unique permutations of each partition.
    C = cell(1,N);
    
    % Generate all unique permutations of each partition and store them in C.
    for k = 1:N
        % Repeat each element of the partition P(k,:) P(k,i) times.
        tmp = repelem(1:d,P(k,:));
        % Add zeros to the end of the array to make it of size n.
        tmp(end+1:n) = 0;
        % Generate all unique permutations of the array and store them in C{k}.
        C{k} = unique(perms(tmp),'rows');
    end
    
    % Concatenate all permutations in C into a single matrix and return it as inds.
    inds = vertcat(C{:});

end

function res = Dalfa(range)
% This function computes harmonic functions for a given range of multi-indices.
    % Define symbolic variables and expression for harmonic functions
    syms x_1 x_2 x_3
    x = [x_1, x_2, x_3];
    expre = 1 / sqrt(x_1^2 + x_2^2 + x_3^2);
    
    len2 = length(range);
    res = cell(1, len2); % initialize a cell array to store results
    
    % Loop over the values in the input argument "range"
    for j = 1:len2
        inds = multiind(range(j), 3); % call the "multiind" function to generate multi-indices
        len = size(inds, 1);
        harmon = sym(zeros(1, len)); % initialize a symbolic array to store harmonic functions
        
        % Loop over the multi-indices and compute the corresponding harmonic functions
        for i = 1:len
            harmon(i) = harmdiff(expre, inds(i, :), x); % call the "harmdiff" function to compute harmonic functions
        end
        
        res{j} = harmon; % store the computed harmonic functions in the cell array "res"
    end
end

function expre = harmdiff(expre, inds, x)
    % Compute the harmonic function associated with a multi-index
    % by taking derivatives of the expression "expre" with respect to
    % the variables in "x" and the corresponding indices in "inds".
    
    % Take the first partial derivative of "expre" with respect to the
    % first variable in "x" and the corresponding index in "inds".
    expre = diff(expre, x(1), inds(1));
    
    % Take the second partial derivative of "expre" with respect to the
    % second variable in "x" and the corresponding index in "inds".
    expre = diff(expre, x(2), inds(2));
    
    % Take the third partial derivative of "expre" with respect to the
    % third variable in "x" and the corresponding index in "inds".
    expre = diff(expre, x(3), inds(3));
    
    % Negate the numerator and denominator of "expre".
    expre = -numden(expre);
end

function [hemN, hemS] = hemisferios(expr)
    % Compute the harmonic functions for each hemisphere of the unit sphere.
    
    % Define symbolic variables for the spherical coordinates.
    syms u v r theta
    
    % Get the symbolic variables used in the input expression "expr".
    vars = symvar(expr);
    
    % Compute the expression for the northern hemisphere.
    exprN = subs(expr, vars, {2*u/(u^2+v^2+1), 2*v/(u^2+v^2+1), (1-u^2-v^2)/(u^2+v^2+1)});
    
    % Compute the expression for the southern hemisphere.
    exprS = subs(expr, vars, {2*u/(u^2+v^2+1), 2*v/(u^2+v^2+1), -(1-u^2-v^2)/(u^2+v^2+1)});
    
    % Convert the expressions to spherical coordinates and simplify.
    hemN = simplify(subs(exprN, {u,v}, {r*cos(theta), r*sin(theta)}));
    hemS = simplify(subs(exprS, {u,v}, {r*cos(theta), r*sin(theta)}));
end

function [hemN, hemS] = mapping(harm, matrix)
    % Create three deformed spherical harmonic functions by mapping the unit sphere using a given matrix.
    
    mapp = 0;
    % Loop over the elements in the input matrix.
    for i = 1:length(matrix)
        % Compute the dot product of the i-th element of the input harmonic functions with the i-th row of the input matrix.
        mapp = mapp + harm{i} * matrix{i}';
    end
    
    % Compute the harmonic functions for each hemisphere of the deformed sphere.
    [hemN, hemS] = hemisferios(mapp);
end

function plist = partitions(total_sum, candidate_set, max_count, fixed_count)
%PARTITIONS Generate all possible partitions of a total sum using a set of candidate integers.

%   Inputs:
%   total_sum       - The integer to partition
%   candidate_set   - A vector of non-negative integers that are candidates for partitioning
%   max_count       - A vector of non-negative integers that restrict the number of times each 
%                     candidate integer can appear in a partition
%   fixed_count     - A non-negative integer representing the fixed number of a particular 
%                     candidate integer to include in the partition
%
%   Outputs:
%   plist           - A matrix containing all possible partitions of total_sum using candidate_set,
%                     max_count, and fixed_count.
%
%   Example:
%       plist = partitions(4) returns the matrix [4,0,0; 3,1,0; 2,2,0; 2,1,1; 1,1,1,1].
%       This represents all possible ways to partition the integer 4 using integers 1 to 4.
%
%   Notes:
%   This function uses a recursive algorithm to generate all possible partitions of the input integer.
%   The function generates partitions in descending order of size and uses candidate_set to restrict 
%   the possible partitions.

  % If candidate_set is not provided, set it to 1:total_sum
  if (nargin<2) || isempty(candidate_set)
    candidate_set = 1:total_sum;
  end
  
  % Determine the number of candidates
  n = length(candidate_set);
  
  % Error checks
  if any(candidate_set<0)
    error('All members of candidate_set must be >= 0')
  end
  
  % Sort candidates in increasing order for efficiency
  if any(diff(candidate_set)<0)
    error('Efficiency requires that candidate_set be sorted')
  end
  
  % If max_count is not provided, set it to the highest possible count for each candidate
  if (nargin<3) || isempty(max_count)
    max_count = floor(total_sum./candidate_set);
    % If max_count is a scalar, convert it to a vector of length n
  elseif length(max_count)==1
    max_count = repmat(max_count,1,n);
  end
  
  % If fixed_count is not provided, set it to an empty array
  if (nargin<4) || isempty(fixed_count)
    fixed_count = [];
    % If fixed_count is not a positive integer, throw an error
  elseif (fixed_count<0) || (fixed_count~=round(fixed_count))
    error('fixed_count must be a positive integer if supplied')
  end
  
  % Check for degenerate cases
  if isempty(fixed_count)
    % If total_sum is 0, return a row of zeros
    if total_sum == 0
      plist = zeros(1,n);
      return
      % If there are no candidates, return an empty array
    elseif (n == 0)
      plist = [];
      return
      % If there is only one candidate, check if total_sum is a multiple of that candidate
    elseif (n == 1)
      p = total_sum/candidate_set;
      if (p==fix(p)) && (p<=max_count)
        plist = p;
      else
        plist = [];
      end
      return
    end
  else
    % If a fixed_count is provided
    if (total_sum == 0) && (fixed_count == 0)
      % If both total_sum and fixed_count are 0, return a row of zeros
      plist = zeros(1,n);
      return
      % If there are no candidates or fixed_count is non-positive, return an empty array
    elseif (n == 0) || (fixed_count <= 0)
      plist = [];
      return
      % If there is only one candidate, check if fixed_count and total_sum are multiples of that candidate
    elseif (n==1)
      if ((fixed_count*candidate_set) == total_sum) && (fixed_count <= max_count)
        plist = fixed_count;
      else
        plist = [];
      end
      return
    end
  end
  
  % Start with the largest candidate and work backwards
  m = max_count(end);
  c = candidate_set(end);
  m = min([m,floor(total_sum/c),fixed_count]);
  
  % Initialize plist
  plist = zeros(0,n);
  
  % Recursively generate partitions for each possible count of the largest candidate
  for i = 0:m
    temp = partitions(total_sum - i*c, ...
      candidate_set(1:(end-1)), ...
      max_count(1:(end-1)),fixed_count-i);
    plist = [plist;[temp,repmat(i,size(temp,1),1)]];
  end

end
