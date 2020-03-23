
%% forward - backward sweep of radial distribution network
function [ node_V,line_I ] = FBS_LF( Branches,RootNode,Vbase ,baseSpotLoads,newSpotLoads )

%   Inputs : Branch,RootNode,Voltage base, SpotLoadNode, DGNode]
%   Outputs: node_V line_I

assert(isRadial(Branches),'The grid is not radial.')

%% adding new DG loads
if (nargin == 6) % combine new load with base load
    for i = 1:size(newSpotLoads,1)
        NotInList = true;
        for j = 1:size(baseSpotLoads,1)
            if (baseSpotLoads(j,1) == newSpotLoads (i,1)) % check if there is a load in the node the new load is added
                baseSpotLoads(j,2) = baseSpotLoads(j,2) + newSpotLoads(i,2);
                NotInList = false;
                break;
            end
        end
        if (NotInList) % add new load to end of spotload if there is not any load in this new load node
            baseSpotLoads(size(baseSpotLoads,1)+1,:) = [newSpotLoads(i,1) newSpotLoads(i,2)];
        end
    end
end

%% Initialization of variables
Branches = CRN(Branches, RootNode); % Root Node check
Nodes = unique(Branches(:,1:2)); % extract nodes
Path = PathFinder(Branches,Nodes); % path matrix
Nx = size(Nodes,1); 
node_Vprev = [Nodes, zeros(Nx,1)]; % could be any thing, it does not matter
node_V = [Nodes, Vbase * ones(Nx,1)]; % init voltage

%% FBS LoadFlow

while any(any (abs(node_V(:,2) - node_Vprev(:,2)) > 1e-6))
    node_Vprev = node_V; % set the last computed voltage to be init value
    Load_I = LoadCurrent(node_Vprev,baseSpotLoads); % compute load currents based on old voltage set
    % Since this is a radial network, end node of each line is unique a
    % line can be represented by end node of it.
    line_I = BranchCurrent(Path,Load_I); % compute line currents based on old voltage set
    for i = 1:Nx % compute new voltages
        V1 = node_V(i,2); %set the strat node voltage
        V1n = node_V(i,1); %set the strat node voltage number
        for j = 1: Nx-1
            if (V1n == Branches(j,1))
                V2n = Branches(j,2); %set the end node voltage number
                Z = Branches(j,3); %set the end node voltage
                for x = 1:Nx
                    if (V2n == abs(line_I(x,1)))
                        I = line_I(x,2); % set the line current
                        V2 = V1 - Z * I; % compute end node new voltage
                        for k = 1: Nx
                            if (abs(node_V(k,1)) == V2n)
                                node_V(k,2) = V2; %put computed voltage to new voltage mat
                                break;
                            end
                        end
                     end
                end
            end
        end
    end
end
end
%% sub functions

% checking if the network is radial
function y = isRadial(Branches)
%isRadial checks whether a grid defined by Branch is radial
%   The idea is checking if there exist two lines ending in one node.
%   if there are more than one line ending in one node the the unique
%   function will get rid of it! and number of unique ending are less than
%   all of them.
y=length(unique(Branches(:,2))) == length(Branches(:,2));
end

% check rootnode
function [ Branch ] = CRN( Branch,RootNode)
%CRN changes the direction of branches to make the input node be the root of
%the tree
Bx = size (Branch,1);
for i = 1:Bx
    if (Branch(i,2) == RootNode)
        NewRoot = Branch(i,1);
        Branch(i,2) = NewRoot;
        Branch(i,1) = RootNode;
        Branch = CRNi(Branch,NewRoot,RootNode);
    elseif (Branch(i,1) == RootNode)
        NewRoot = Branch(i,2);
        Branch = CRNi(Branch,NewRoot,RootNode);
    end
end
end
function [ Branch ] = CRNi( Branch,RootNode,PreviuseRoot )
%CRNi Internal part of CRN function
Bx = size (Branch,1);
for i = 1:Bx
    if (Branch(i,2) == RootNode) && (Branch(i,1) ~= PreviuseRoot)
        NewRoot = Branch(i,1);
        Branch(i,2) = NewRoot;
        Branch(i,1) = RootNode;
        Branch = CRNi(Branch,NewRoot,RootNode);
    elseif (Branch(i,1) == RootNode)
        NewRoot = Branch(i,2);
        Branch = CRNi(Branch,NewRoot,RootNode);
    end
end
end

% path matrix formation
function [ Path ] = PathFinder( Branch,Nodes )
%PathFinder will find all sub nodes of a node
Bx = size(Branch,1);
Path = zeros(Bx+2);
Path(:,1) = [0; Nodes];
Path(1,:) = [0, Nodes.'];
for i = 2:Bx+2
    Path(i,i) = 1;
    A = Path(i,1);
    j = 1;
    while (j <= Bx)
        if (A == Branch(j,2))
            B = Branch(j,1);
            for k = 2:Bx+2
                if(Path(1,k) == B)
                    Path(k,i) = 1;
                    A = B;
                    j = 0;
                    break;
                end
            end
        end
        j = j + 1;
    end
end
end

% load current computation
function [ LoadCurrent ] = LoadCurrent( NodeVoltage,SpotLoadNode )
%LoadCurrent computes loads current
Current = 0;
Lx = size(SpotLoadNode,1);
LoadCurrent = zeros(Lx,2);
for i = 1:Lx
    for j = 1:size(NodeVoltage,1)
        if (abs(NodeVoltage(j,1)) == abs(SpotLoadNode(i,1)))
            Current = conj(SpotLoadNode(i,2)./NodeVoltage(j,2));
            %break;
        end
    end
    LoadCurrent(i,:) = [SpotLoadNode(i,1) Current];
end
end

%line current computation
function [ BranchCurrent ] = BranchCurrent( Path,LoadCurrent )
%BranchCurrent computes current in each line
Lx = size (LoadCurrent,1);
Px = size(Path,1);
Current = 0;
BranchCurrent = zeros(Px-1,2); % create the output
for i = 2:Px
    for j = 2:Px
        if Path(i,j)
            for k = 1:Lx
                if Path(1,j) == LoadCurrent(k,1)
                    Current = Current + LoadCurrent(k,2);
                    break;
                end
            end
        end
    end
    BranchCurrent(i-1,:) = [Path(i,1) Current];
    Current = 0;
end
end