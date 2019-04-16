% script for homework 1 problem 
% solve -(au')' + u = f on [ 0 , 1 ] 
% BCs: u(x=0) = alpha and u(x=1) = beta
% Use CG finite element method 
% Author: Jieun Lee
% Date: 

% flags 
errorOn=1;  % calculate errors (see end of script; set to 0 to turn off) 
testCase=3; % which test case to run 

% parameters 
nElements= 10;   % elements on first grid (each refinement doubles again)
nQuadrature= 2;       % quadrature order for error calculation
numGrids=6;                  % number of grids to run on
L2Errors=zeros(numGrids,1);  % store errors in L2-norm on each grid
H1Errors=zeros(numGrids,1);  % store errors in H1-norm on each grid

if (testCase==1)
 alpha= 1;
 beta= 2;
elseif (testCase==2)
 alpha= 1;
 beta= 2;
elseif (testCase==3)
 alpha= 1;
 beta= exp(1);
else
 error('hw1 :: invalid test case... exiting.');
end

for gridIndex=1:numGrids

% generate grid information; all nodes and boundary flags 
numNodes=1+nElements;
grid=linspace(0,1,numNodes);  % mesh nodes 
dof=zeros(1,numNodes);
% set the Dirichlet boundary node flags
% In 2D you could test if the coordinates are on the bndy
dof(1)=-1;  
dof(end)=-2;  
bndyData=[alpha;beta];

nFree=numNodes-2;               % number of degrees of freedom
dof(2:(numNodes-1))= (1:nFree); % Simply maps to index of a node within array 
                                % of free nodes

% generate elements 
element=zeros(nElements,2);
for n=1:nElements
    element(n,:)=[n (n+1)]; % each row "n" stores node indices for element "n"
end

% build mass matrix 
M=buildMass(grid, dof, element, nElements, nFree);

% build stiffness matrix 
K=buildStiffness(grid, dof, element, nElements, nFree, testCase);

% NOTE: mass and stiffness need not be stored separately, but that's OK ;-)

% build rhs vector (forcing and BCs) 
rhs=buildRHS(grid, dof, element, nElements, nFree, bndyData, testCase);

% solve the linear system (call UMFPACK)
uI=(M+K)\rhs;

% build global solution vector (add left boundary value)
u=zeros(numNodes,1);
u=[alpha;uI;beta];
% get errors 
if errorOn
if testCase<4 % true solution known

   L2Err=0;
   H1Err=0;

   for i=1:nElements

       [x, w]=gl_weight(grid(i), grid(i+1), nQuadrature);
       % L2 error
       trueSoln=trueU(x, testCase);
       appxSoln=feEval(x, grid, dof, element, i, u);
       L2Err=L2Err+sum((trueSoln - appxSoln).^2.*w); 
       % H1 semi-norm
       trueSoln=trueUDerivative(x, testCase);
       appxSoln=feEvalDerivative(x, grid, [dof(i), dof(i+1)], element, i, u);
       H1Err=H1Err+sum((trueSoln - appxSoln).^2.*w);

   end % elements
   
   H1Err=H1Err+L2Err; % convert to full H1-norm
   L2Errors(gridIndex)=sqrt(L2Err);
   H1Errors(gridIndex)=sqrt(H1Err);

else 
    error('invalid testCase');
end % test case
end % flag for errors 


%plot testCase (c) for n = 80, 160, 320
if testCase == 3
    if gridIndex == 4
        plot(grid, u - transpose(trueU(grid, testCase)), 'red')
        hold on
    end
    
    if gridIndex == 5
        plot(grid, u - transpose(trueU(grid, testCase)), 'blue')
        hold on
    end
    
    if gridIndex == 6
         plot(grid, u - transpose(trueU(grid, testCase)), 'green')
    end            
   
end

nElements=2*nElements;

end % loop over grids

if errorOn
if testCase<4 % true solution known
disp('L2 errors:')
L2Errors
disp(' ')
if (numGrids>1) 
disp('L2 rates:')
log(L2Errors(1:numGrids-1)./L2Errors(2:numGrids))/log(2.0)
disp(' ')
end
disp('H1 errors:')
H1Errors
if (numGrids>1)
disp(' ')
disp('H1 rates:')
log(H1Errors(1:numGrids-1)./H1Errors(2:numGrids))/log(2.0)
end
end
end

