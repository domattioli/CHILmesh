classdef CHILmesh
    %CHILMESH Triangulation, Quadrangulation, or Mixed-Elements in 2-D.
    %   CHILMESH supports topological and geometric queries for 2D
    %   triangular, quadrangular, and mixed finite element meshes. CHILMESH
    %   is a subclass to the parent CHILgraph superclass, where CHILMESH is
    %   an object of a nodal graph assembled into properties and methods
    %   for analyzing the element nodal connectivity list. Various CHILMESH
    %   methods exist that allow for queries involving the geometry and
    %   topology of vertices, edges, and elements ("vertex", "node", and
    %   "point" are interchangeable throughout this class' documentation.
    %   Additional methods allow for plotting the mesh.
    %   Note: maybe CHILmesh should be rebuilt with inheritance as a
    %   subclass of the built-in MATLAB classes triangulation/graph.
    %   
    %   CHILMESH creates an object for the mesh from select input file
    %   formats or converts the mesh from a 2D triangular mesh of the
    %   MATLAB triangulation class. CHILMESH employs many methods and
    %   structure similar to MATLAB triangulation.
    %
    %   CM	= CHILMESH(FILENAME) creates a 2D mesh from the inputted mesh
    %   FILENAME. Some acceptable formats include GMESH and FORT14
    %   (extentions .msh and .14, respectively). The resulting mesh object
    %   structure CM contains fields for the points and connectivity list
    %   that are extracted from the inputted file.
    %
    %   CM = CHILMESH(C,P) creates a 2D mesh from the connectivity
    %   list/matrix C and the points list/matrix P, where C is an Mx3 or
    %   Mx4 matrix with each row listing the indices that define the
    %   element with respect to P, M is the number of triangles or
    %   quadrilaterals, and P is an NxD matrix, where N is the number of
    %   vertices in the mesh and D is the dimension. Note that while
    %   CHILMESH is for 2D meshes, D may be equal to 3, however, no
    %   alterations will every occur to the z-coordinates during any of the
    %   CHILMESH methods.
    %
    %   CM = CHILMESH(C,P,BC) creates a 2D mesh from the connectivity list
    %   C and the points list P, where CONN is an Mx3 or Mx4 matrix, P is a
    %   Nx2 or Nx3 matrix, and BC is a structure defining the boundary
    %   conditions of the mesh.
    %   
    %   Note: The following examples are not yet built-in (future-work).
    %   Example 1: Load a 2D mesh from a .14 file and query, plot the
    %       centroids of each element as red asterisks.
    %       load testCHILmesh2D                 % Loads the mesh.
    %       CM	= CHILMESH(FILENAME)            % Creates mesh object.
    %       Ct	= CM.centroid;                  % Shows diagonals as matrix.
    %       CM.plot;                            % Plots mesh.
    %       hold on; plot(Ct(:,1),Ct(:,2),'r*');
    %
    %   Example 2: Create a delanuay triangulation of random points,
    %   convert it to a CHILMESH, then query and plot the boundary edges.
    %       P	= rand(100,2);
    %       DT	= delaunayTriangulation(P);
    %       CM  = CHILmesh(DT,P)                % Convert to CHILMESH.
    %       fe  = CM.boundaryEdges'             % Boundary edges.
    %       CM.plot;                            % Plots mesh.
    %       plot(P(fe,1),P(fe,2),'r','LineWidth',2); hold off;
    %
    %   Example 3: Load a 2D mesh from a .msh file and query, plot the
    %   layers of the mesh via contour map.
    %
    %   CHILMESH methods:
    %       BOUNDARYEDGES       - Edges that define the boundary of the mesh.
    %       BUILDADJACENCIES	- Builds the 6 adjacency lists of the mesh.
    %       CENTROID            - Geometric centroid of triangle or quadrilateral in the mesh.
    %       CHILMESH            - Creates a mesh object class for an inputted mesh.
    %       DIAGONALS        	- Vertices defining the 2 diagonals of a quadrilateral.
    %       EDGELENGTH         	- Length of edge in euclidian distance.
    %       EDGEMIDPOINT     	- Midpoint coordinates of edges in the mesh.
    %       ELEM2ELEM           - Neighbors to a triangle or quadrilateral.
    %       ELEMQUALITY         - Geometric quality of triangle or quadrilateral in the mesh.
    %       ELEMTYPE            - Identifies triangles, quadrilaterals, and 1D elements in the mesh.
    %       INTERIORANGLES    	- Interior angles of triangle or quadrilateral in the mesh.
    %       ISPOLYCCW           - Sets nodal conectivity to a counter-clockwise orientation.
    %       MEDIANS         	- Non-adjacent edges of each element's vertices in the mesh.
    %       MESHLAYERS          - Discretize the elements of the mesh into layers.
    %       SIGNEDAREA         	- Signed area of polygonal element.
    %
    %	CHILMESH properties:
    %       GRIDNAME            - The name of the mesh.
    %     	POINTS              - The coordinates of the mesh's points.
    %     	CONNECTIVITYLIST	- The polygonal vertex-connectivity list.
    %       BOUNDARYCONDITION	- The boundary condition(s) of the mesh (future work).
    %     	MAPPING           	- The current geographical mapping system (future work).
    %       PLOT                - Ploting object for the mesh (future work).
    %
    %
    %   CHILMESH was developed by Ethan J. Kubatko, Dominik Mattioli,
    %   Dustin J. West, Dylan Wood, and Omar El Khoury in the Computational
    %   Hydrodynamics & Informatics Laboratory (CHIL) at The Ohio State
    %   University in 2017, with funding in part by Aquaveo and U.S. Army
    %   Corps of Engineers.
    %   Some aspects of this object class are inspired by:
    %       1. The MathWorks, Inc. built-in function for triangulations.
    %       2. ***Source for adjacencies***
    %       3. ***Other sources***
    %
    %   See also DELAUNAYTRIANGULATION, TRIANGULATION, GRAPH.
    %======================================================================
    
    %{
properties
    %GRIDNAME The name of the mesh
    %   GRIDNAME is the name of the mesh retrieved from the input file, if
    %   any. If the input is a triangulation, GRIDNAME must be entered
    %   manually.
    %
    %   See also CHILMESH.
    GridName;
    
    %POINTS The coordinates of the points in the mesh.
    %	POINTS is a matrix of size nVerts-by-ndim, where nVerts is the
    %   number of points and ndim is the number of dimensions (2 or 3). If
    %   column vectors of X,Y or X,Y,Z coordinates are used, the data is
    %   consolidated into a single matrix.
    %
    %   See also CHILMESH, TRIANGULATION/POINTS.
    Points;
    
    %CONNECTIVITYLIST The polgyonal-element vertex-connectivity list.
    %	CONNECTIVITYLIST is a matrix is of size nElems-by-nv, where nElems
    %   is the number of polygonal elements and nv is the number of
    %   vertices (3 or 4). If the mesh contains both triangular and
    %   quadrilateral elements, nv = 4 and all triangular elements contain
    %   a redundant point/vertex/node - usually in the first column. This
    %   is for plotting purposes.
    %
    %   See also CHILMESH, TRIANGULATION/CONNECTIVITYLIST.
    ConnectivityList;
    
    %BOUNDARYCONDITION Boundary condition(s) of the mesh.
    %	BOUNDARYCONDITION is a structure with fields containing information
    %   in regards to the boundary condition(s) of the mesh.
    %
    %   See also CHILMESH, CHILMESH/READGRIDFILE.
    BC;
end
    %}
    
    properties (SetAccess = protected)              % Protected properties.
        GridName
        Type
    end
    
    properties (Hidden = false)                 	% Visible properties.
        Points
        ConnectivityList
        BoundaryCondition
    end
    
    properties (SetAccess = protected)
        Mapping
        Plot
    end
    
    properties (Hidden = true)                      % Hidden properties.
        Adjacencies
        cpp
        dim
        FileName
        Layers
        nEdges
        nElems
        nLayers
        nVerts
    end
    
    methods
        %% ==========================Instantiation=========================
        % Construct CHILmesh Class Object
        function CM	= CHILmesh(varargin)            % Class Constructor.
            %CHILMESH Constructs object for 2D polygonal (tri or quad) meshes.
            %   CM = CHILMESH() ******** explain possible inputs and the
            %   outputted data structure.
            %
            %   See also CHILMESH, TRIANGULATION.
            %==============================================================
            
            % Check I/O argument(s).
            CM.GridName = [];                       % Initialize object.
            CM.Points   = [];
            CM.ConnectivityList = [];
            CM.BoundaryCondition = [];
            switch nargin
                case {0,1}                          % Grid or triangulation given.
                    % Create an input if none exist.
                    if nargin == 0
                        % Ask for a file or triangulation to read.
                        askInput  = questdlg('Do you have an input for CHILmesh?',...
                            'Zero Input Dialog','Yes','No','Yes');
                        noInput	= true;
                        attempt = 1;
                        while noInput           	% Ask for an input.
                            if strcmp('Yes',askInput) 	% "Yes" input.
                                inputType	= questdlg('What type of input?','Input Type',...
                                    'Grid File','Connectivity & Point Lists','Grid File');
                                
                                % Retrieve an input to assign.
                                switch inputType
                                    case 'Grid File'% Ask for file.
                                        noFileChosen	= true;
                                        fileIter	= 1;
                                        while noFileChosen
                                            [FN,PN]	= uigetfile({'*.14;*.gmsh'},...
                                                'Select A Grid File.');
                                            if FN == 0  % No file selected.
                                                warndlg('Please select a file.','No File Selected');
                                                
                                            elseif fileIter	== 3    % Max # of chances.
                                                warndlg('Empty CHILmesh object created','No Input.');
                                                return
                                                
                                            else    	% Confirm that file exists.
                                                if exist([PN,FN],'file') ~= 2
                                                    warndlg('File doesn''t exist; select another.','No File Found');
                                                else    % File chosen, exit loop.
                                                    noFileChosen	= false;
                                                end
                                            end
                                            fileIter	= fileIter + 1;
                                        end
                                        % Assign file and path.
                                        varargin{1}	= [PN,FN];
                                        
                                    case 'Connectivity & Point List'
                                        % Ask for lists.
                                        varargin{1}.ConnectivityList	= inputdlg('Input Connectivity List','Input Lists');
                                        varargin{1}.Points	= inputdlg('Input Points List','Input Lists');
                                        
                                        % Confirm validity of inputs.
                                        
                                    otherwise
                                        attempt	= attempt + 1;
                                        continue
                                end
                                noInput	= false;        % End input search.
                                
                            else                        % "No" input.
                                warndlg('No input given; CHILmesh requires an input.','Warning: Zero Inputs');
                                if attempt == 3         % Max of 3 attempts.
                                    warndlg('Empty CHILmesh object created','No Input.');
                                    return
                                end
                                attempt	= attempt + 1;  % Chronicle attempt.
                            end
                        end
                    end
                    
                    % Read input (either triangulation or grid file).
                    switch class(varargin{1})       % Now 1 input exists.
                        case {'triangulation',...	% Given grid as Triangulation.
                                'delaunayTriangulation','struct'}
                            CM.GridName	= ['MATLAB ', class(varargin{1})];
                            CM.Points(:,1)	= double(varargin{1}.Points(:,1));
                            CM.Points(:,2)	= double(varargin{1}.Points(:,2));
                            CM.ConnectivityList	= uint32(varargin{1}.ConnectivityList);
                            CM.BoundaryCondition.id	= [];
                            
                        case 'char'                 % Given grid as a file.
                            CM	= readGridFile(CM,varargin{1});
                    end
                    
                case {2,3,4}% Given points, connectivity lists & boundary conditions.
                    % Assign connectivity list and points list inputs.
                    CM.Points	= varargin{2};
                    CM.ConnectivityList	= uint32(varargin{1});
                    if nargin >= 3
                        iGridName   = false(nargin,1);
                        for idx = 1:nargin          % Find input for GridName.
                            iGridName(idx) = ischar(varargin{idx});
                        end
                        if all(iGridName == false)
                            CM.GridName	= 'Custom Input';
                            CM.BoundaryCondition    = varargin{3};
                        else
                            CM.GridName     = varargin{iGridName};
                            if nargin == 4
                                CM.BoundaryCondition     = varargin{setdiff(3:4,find(iGridName))};
                            else
                                CM.BoundaryCondition	= [];
                            end
                        end
                    else
                        CM.GridName	= 'Custom Input';
                        CM.BoundaryCondition	= [];
                    end
                    
                otherwise                           % Unsupported input.
                    errordlg('Too many inputs.','Unsupported Input');
            end
            
            % Assign other necessary properties.
            CM.nVerts   = size(CM.Points,1);
            if size(CM.Points,2) == 2               % Add z-coord.
                CM.Points   = [CM.Points,zeros(CM.nVerts,1)];
            end
            CM.nElems   = size(CM.ConnectivityList,1);
            [~,CM]  = CM.isPolyCCW;                 % Ensure CCW orientation.
            CM  = CM.buildAdjacencies;              % Build mesh adjacency lists.
            CM  = CM.meshLayers;                    % Identify mesh layers.
        end
        
        % Read Grid file
        function CM	= readGridFile(CM,FILENAME)
            %==============================================================
            %READGRIDFILE Identify file format and read in mesh.
            %   CM = ReadGridFile(CM,FILENAME) extracts the file extension
            %   from an inputted mesh grid file (FILENAME), reads the file,
            %   and then constructs a CHILmesh object CM. The method of
            %   reading depends on the extension type of the input file.
            %
            %   See also CHILMESH, READFORT14, READGMSH.
            %==============================================================
            % Identify grid file format.
            [~,gn,ext]	= fileparts(FILENAME);
            
            % Read in mesh file
            switch ext
                case '.14'                          % Fort14 (ADCIRC) form.
                    [CM,bc]	= readFort14(CM,FILENAME);
                    
                case '.msh'                         % Gmsh file format.
                    [CM,bc]	= readGMSH(CM,FILENAME);% Get connect., points.
                    
                otherwise                           % Unknown file format.
                    warndlg('Mesh file format not supported; Mesh will be empty.','Unsupported input');
            end
            CM.GridName	= gn;                       % Assign grid name.
            CM.BoundaryCondition	= bc;        	% Assign boundary conditions.
        end
        
        % Import Fort14 Grid File
        function [CM,BoundaryCondition]	= readFort14(CM,FILENAME)
            %==============================================================
            %READFORT14 Read in FE mesh from the ADCIRC text file format.
            %   CM = READFORT14(CM,FILENAME) *** explain the requirements
            %   for an ADCIRC text file format for a mesh.
            %
            %   See also CHILMESH, READGRIDFILE, READGMSH.
            %==============================================================
            
            % Read in the number of elements & number of grid points.
            FID	= fopen(FILENAME, 'r');             % Open grid file.
            CM.GridName	= fgetl(FID);             	% Read grid name.
            INFO= textscan(FID,'%f %f %*[^\n]',1);  % Read # elements, pts.
            P	= textscan(FID,...                  % Read coordinates of
                '%*f %f %f %f %*[^\n]', INFO{2});	% each point in mesh.
            C	= textscan(FID,...                  % Read connectivity.
                '%*f %f %f %f %f %f %*[^\n]', INFO{1},'EmptyValue',0);
            
            % Read in the Number of Open Ocean Segments.
            nOOS	= cell2mat(textscan(FID, '%f %*[^\n]', 1));
            
            % Skip over the total number of elevation specified boundary nodes
            textscan(FID, '%*f %*[^\n]', 1,'EmptyValue',0);
            
            % Accomodate for any open ocean segments.
            if isempty(nOOS) || nOOS == 0         	% No open ocean segments.
                nOOS	= 0;
                
            else                                    % Open ocean segments exist.
                % Intitialize boundary condition structure.
                BoundaryCondition(nOOS,1)	= struct('id',[],'des',[],'nodes',[],'att',[]);
                
                for k	= 1:nOOS
                    % Get boundary condition ID and description.
                    BoundaryCondition(k).id    = -1;
                    BoundaryCondition(k).des   = 'Open Ocean';
                    
                    % Read number of nodes in open ocean segment k.
                    NVDLL	= textscan(FID, '%f %*[^\n]', 1);
                    
                    % Node numbers in boundary segment k.
                    Nodes	= textscan(FID, '%f %*[^\n]', NVDLL{1});
                    
                    % Check nodes for any errors.
                    if any(isnan(Nodes{:}))
                        errordlg('Missing nodes in open ocean boundary segments','Missing Nodes');
                        return
                    end
                    
                    % Store node string of Open Ocean Boundary.
                    BoundaryCondition(k).nodes	= Nodes{:};
                    BoundaryCondition(k).att	= 0;% No attributes.
                end
            end
            
            % Read in the Number of Normal Flow Boundary Segments.
            NNFBS	= cell2mat(textscan(FID, '%f %*[^\n]', 1,'EmptyValue',0));
            
            % Scan over the number of normal flow specified boundary nodes.
            textscan(FID, '%*f %*[^\n]', 1);
            
            % Accomodate for any normal flow specified boundary segments.
            if isempty(NNFBS) || NNFBS == 0     	% No normal flow boundary segments.
                NNFBS	= 0;
                
            else                                    % Store normal flow boundary coordinates.
                % Append to boundary condition structure.
                BoundaryCondition(nOOS+NNFBS) = struct('id',[],'des',[],'nodes',[],'att',[]);
                
                for k = (nOOS+1):(nOOS+NNFBS)
                    % Read in the number of nodes and the boundary type
                    INFO	= textscan(FID, '%f %f %*[^\n]', 1);
                    [NVELL, BoundType]	= deal(INFO{:});
                    
                    % Get Node numbers in boundary segment k, assign boundary
                    % condition(s) ID, description, node string, and attributes.
                    BoundaryCondition(k).id	= BoundType;
                    if any(BoundType == [0 2 10 12 20 22])	% No flow boundary conditions (mostly)
                        Nodes	= textscan(FID, '%f %*[^\n]', NVELL);
                        BoundaryCondition(k).des   = 'External Boundary';
                        BoundaryCondition(k).nodes = Nodes{1};
                        BoundaryCondition(k).att   = 0; % No attributes.
                        
                    elseif any(BoundType == 6)      % Radiation bondary condition.
                        Nodes	= textscan(FID, '%f %*[^\n]', NVELL);
                        BoundaryCondition(k).des   = 'Radiation Boundary';
                        BoundaryCondition(k).nodes = Nodes{1};
                        BoundaryCondition(k).att   = 0; % No attributes
                        
                    elseif any(BoundType == [1 11 21])	% Internal boundary condition.
                        % Node numbers in boundary segment k
                        Nodes	= textscan(FID, '%f %*[^\n]', NVELL);
                        BoundaryCondition(k).des   = 'Internal Boundary';
                        BoundaryCondition(k).nodes = Nodes{1};
                        BoundaryCondition(k).att   = 0;% No attributes
                        
                    elseif any(BoundType == [3 13 23])	% Barrier boundary condition.
                        % Node numbers in boundary segment k
                        Nodes	= textscan(FID, '%f %f %f %*[^\n]', NVELL);
                        BoundaryCondition(k).des   = 'External Barrier';
                        BoundaryCondition(k).nodes = Nodes{1};
                        BoundaryCondition(k).att   = horzcat(Nodes{2:end});
                        
                    elseif any(BoundType == 18) 	% Channel boundary condition.
                        Nodes	= textscan(FID, '%f %f %f %f %f %f %*[^\n]', NVELL);
                        BoundaryCondition(k).des   = 'Line';
                        BoundaryCondition(k).nodes = Nodes{1};
                        BoundaryCondition(k).att   = horzcat(Nodes{2:end});
                        
                    elseif any(BoundType == [4 24])	% Internal barrier boundary condition.
                        Nodes	= textscan(FID, '%f %f %f %f %f %*[^\n]', NVELL);
                        BoundaryCondition(k).des   = 'Internal Barrier Type 1';
                        BoundaryCondition(k).nodes = [Nodes{1}, Nodes{2}];
                        BoundaryCondition(k).att   = horzcat(Nodes{3:end});
                        
                    elseif any(BoundType == [5 25])	 % Internal barrier boundary condition
                        Nodes	= textscan(FID, '%f %f %f %f %f %f %f %f %*[^\n]', NVELL);
                        BoundaryCondition(k).des   = 'Internal Barrier Type 2';
                        BoundaryCondition(k).nodes = [Nodes{1}, Nodes{2}];
                        BoundaryCondition(k).att   = horzcat(Nodes{3:end});
                        
                    elseif any(BoundType == [102 112 122])
                        errordlg('need to account for these boundary types still');
                    end
                end
            end
            
            % If there are no boundary conditions:
            if nOOS == 0 && NNFBS == 0              % Set BC as empty.
                BoundaryCondition.id	= [];	BoundaryCondition.des	= [];
                BoundaryCondition.nodes	= [];	BoundaryCondition.att	= [];
            end
            fclose(FID);                         	% Close file.
            
            % Get mesh dimension.
            idx1D   = C{4} == 0;                    % Identify 1D elements.
            do1Dexist	= any(idx1D);
            if C{1}(1) == 2 || do1Dexist
                CM.dim	= 1;
            elseif C{1}(1) > 2 && ~any([BoundaryCondition.id] == 18)
                CM.dim	= 2;
            elseif C{1}(1) > 2 && any([BoundaryCondition.id] == 18) && do1Dexist
                CM.dim	= [1 2];
            end
            
            % Initialize mesh structure array and store components of mesh.
            switch C{1}(1)
                case 2                              % 1D elements.
                    CM.Points        	= horzcat(P{:});
                    CM.ConnectivityList	= uint32(horzcat(C{:,2:3}));
                    CM.cpp            	= [nan nan];
                    
                case 3                              % 2D triangles.
                    CM.Points          	= horzcat(P{:});
                    CM.ConnectivityList	= uint32(horzcat(C{:,2:4}));
                    CM.cpp             	= [nan nan];
                    
                case 4                              % 2D quads/mixed-elements.
                    CM.Points         	= horzcat(P{:});
                    CM.ConnectivityList	= uint32(horzcat(C{:,2:end}));
                    CM.cpp            	= [nan nan];
                    
                otherwise                           % No 3D elements.
                    errordlg('Element type not supported','Invalid Element Type')
            end
            
            % Isolate 1D elements, remove them from the connectivity list.
            if do1Dexist                            % 1D elements exist.
                IDX1DELEM	= nOOS+NNFBS+1;
                BoundaryCondition(IDX1DELEM) = struct('id',19,...
                    'des','1D Elements',...
                    'nodes',CM.ConnectivityList(idx1D,1:2),'att',[]);
                CM.ConnectivityList(idx1D,:)    = [];
            end
        end
        
        % Import Gmsh Grid File
        function [CM,BoundaryCondition]	= readGMSH(CM,FILENAME)
            %==============================================================
            %READGMSH Read in FE mesh from the Gmsh ASCII text file format.
            %   CM = READGMSH(CM,FILENAME) *** explain the requirement(s)
            %   for an Gmsh text file format for a mesh. converts to .14***
            %
            %   READGMSH is written to reflect the MSH ASCII file format
            %   described on the Gmsh database website at
            %   http://gmsh.info/doc/texinfo/gmsh.html#File-formats. The
            %   last revision to this code occurred on: 9/13/2017.
            %
            %   See also INITMESH, READFORT14, READGRIDFILE.
            %==============================================================
            
            % FILENAME = 'LakeGeorgeBasinGeometry_frontal1.msh'
            % Open doc, skip over version-number, file-type, & data-size.
            FID	= fopen(FILENAME,'r+');              % Open grid file.
            numPhysName	= textscan(FID,'%q %*[^\n]',5);
            
            % Get physical dimension {1}, number {2}, and name{3}.
            physicalNames	= textscan(FID,'%d %d %q',...
                str2double(numPhysName{1}{5}));   	% (2) may not be monotonically increasing.
            %             physicalDimension	= physicalNames{1};
            %             physicalNumber	= physicalNames{2};
            %             physicalName	= physicalNames{3};
            
            % Get nodes of mesh.
            numNodes	= textscan(FID,'%q %*[^\n]',3);
            Nodes	= textscan(FID,'%d %f %f %f',...
                str2double(numNodes{1}{3}));        % X-Y-Z coordinates.
            %             nodeNumber  = Nodes{1};
            
            % Get elements of mesh.
            numElems    = textscan(FID,'%q %*[^\n]',3);
            Elements	= textscan(FID,'%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d',...
                str2double(numElems{1}{3}));        % Connectivity info.
            %             elementNumber	= Elements{1};
            %             elementType 	= Elements{2};
            %             numberOfTags	= Elements{3};
            fclose(FID);                            % Close doc.
            
            % Identify max number of columns with text for Elements.
            dontStop	= true;
            idx = size(Elements,2);                 % Last column.
            while dontStop
                if any(Elements{idx} > 0)           % Find entries in column.
                    lastCol	= idx;                  % Last column w text found.
                    dontStop	= false;            % Stop loop.
                else
                    idx	= idx - 1;
                end
            end
            
            % Create mesh class by assigning points, connectivity lists.
            CM.Points   = [Nodes{2:4}];             % Points list.
            idx2NodeLine	= Elements{2} == 1;     % Index 1D line segments.
            idx3NodeTri 	= Elements{2} == 2;     % Index 2D triangles.
            idx4NodeQuad	= Elements{2} == 3;     % Index 2D quadrilaterals.
            if any(idx4NodeQuad)                    % Determine max type number.
                maxTypeNum  = 4;
            elseif any(idx3NodeTri)
                maxTypeNum  = 3;
            elseif any(idx2NodeLine)
                maxTypeNum  = 4;
            end
            CM.ConnectivityList	= ...               % 1D line segments;
                [Elements{lastCol-3}(idx2NodeLine),...
                Elements{lastCol-2}(idx2NodeLine),...
                nan(sum(idx2NodeLine),maxTypeNum-2);...
                Elements{lastCol-3}(idx3NodeTri),...% 2D triangles;
                Elements{lastCol-2}(idx3NodeTri),...
                Elements{lastCol-1}(idx3NodeTri),...
                nan(sum(idx3NodeTri),maxTypeNum-3)];
            if maxTypeNum == 4                      % 2D quadrilaterals.
                CM.ConnectivityList = [CM.ConnectivityList;
                    Elements{lastCol-3}(idx4NodeQuad),...
                    Elements{lastCol-2}(idx4NodeQuad),...
                    Elements{lastCol-1}(idx4NodeQuad),...
                    Elements{lastCol}(idx4NodeQuad)];
            end
            
            % Get boundary condition information.
            BoundaryCondition  = physicalNames;
        end
        
        %% ========================Getters & Setters=======================
        %%% Note: I didn't know what these were when I first wrote it; the
        %%% class gets along (inefficiently) without them, should be
        %%% rewritten to include them.
%         % Getters.
%         function gotPoints = get.Points(CM)
%             % Don't know how to write this correctly.
%             gotPoints   = size(CM.Points,1);
%         end
%         
%         % Setters.
%         function CM	= set.Points(CM,newPoints)
%             % Don't know how to write this correctly.
%             CM.Points   = [CM.Points; newPoints];
%         end
%         

        %% ==========================Define Mesh===========================
        % Build Adjacency Lists
        function CM	= buildAdjacencies(CM)
            %BUILDADJACENCIES Obtain adjacency lists for a 2D mesh.
            %   CM = BUILDADJACENCIES(CM) computes all 6 adjacency lists
            %   for the mesh (stored in hidden field M.Adjacencies). Half
            %   of the adjacency lists are "downward", indicating that they
            %   provide the lower dimensional entities of the input entity.
            %   The remaining adjacency lists are "upwards" because they
            %   provide the input entitiy's higher dimensional entities.
            %   The adjacency lists are stored as fields in the struct CM.
            %
            %   The "downward" adjacencies:
            %   1. Elem2Vert - Identifies VERTICES of each ELEMENT.
            %                - A triangulation's element-connectivity list.
            %   2. Edge2Vert - Identifies VERTICES of each EDGE.
            %                - A triangulation's edge-connectivity list.
            %   3. Elem2Edge - Identifies EDGES of each ELEMENT. A value of
            %   zero indicates a singularity, i.e. redundant vertices that
            %   allows quad connectivity lists to accomodate tris.
            %
            %	The "upward" adjacencies:
            %   4. Vert2Edge - Identifies EDGES attached to each VERTEX.
            %   5. Vert2Elem - Identifies ELEMENTS attached to each VERTEX.
            %                - A triangulation's vertex-attachments list.
            %   6. Edge2Elem - Identifies ELEMENTS attached to each EDGE.
            %                - A triangulation's edge-attachments list.
            %
            %   It is helpful to read the name of an individual field in
            %   reverse order to discern the meaning of the contents, e.g.
            %   CM.Elem2Vert contains the vertices of each element.
            %
            %   CM = BUILDADJACENCIES(CM,TR) builds the 6 adjacency lists
            %   for the connectivity list TR.
            %
            %   Example:
            %   % Points list.
            %   CM.Points = [1 1 0; 2 1 0; 3 1 0; 4 1 0; 1 2 0; 2 2 0; 3 2
            %   0; 4 2 0; 1 4 0; 2 3 0; 3 4 0; 4 3 0];
            %   % Connectivity list.
            %   CM.ConnectivityList = [1 2 6 5; 2 3 7 6; 3 4 8 7; 5 6 10 9;
            %   6 6 7 10; 7 8 12 10; 10 12 11 9];
            %   % This is a mixed-element mesh, where element/face number 4
            %   is a triangle and the remaining elements/faces are quads.
            %   The triangle is given a redundant vertex (the edge is 1D)
            %   in order assimilate it's connectivity in with the quads.
            %   CM = BuildAdjacencies(M);
            %   % Returns:
            %   CM.Adjacencies.Elem2Vert = [1 2 6 5; 2 3 7 6; 3 4 8 7;
            %   5 6 10 0; 6 6 7 10; 7 8 12 10; 10 12 11 9]
            %   CM.Adjacencies.Edge2Vert = [5 1; 1 2; 6 2;2 3; 7 3; 3 4;
            %   9 5; 5 6; 10 6; 6 7; 10 7; 4 8; 7 8; 11 9; 9 10; 12 11;
            %   8 12; 10 12]
            %	CM.Adjacencies.Elem2Edge = [2 3 8; 4 5 10; 6 12 13; 8 9 15;
            %	0 10 11; 13 17 18; 18 16 14]
            %   CM.Adjacencies.Elem2Edge = [5 1; 1 2; 6 2; 2 3; 7 3; 3 4;
            %   9 5; 5 6; 10 6; 6 7; 10 7; 4 8; 7 8; 11 9; 9 10; 12 11;
            %   8 12; 10 12]
            %	CM.Adjacencies.Vert2Edge = see CHILMESH\VERT2EDGE.
            %	CM.Adjacencies.Vert2Elem = see CHILMESH\VERT2ELEM.
            %	CM.Adjacencies.Edge2Elem = [1 0; 1 0; 2 1; 2 0; 3 2; 3 0;
            %	4 0; 4 1; 5 4; 5 2; 6 5; 3 0; 6 3; 7 0; 7 4; 7 0; 6 0; 7 6]
            %
            %   See also CHILMESH, TRIANGULATION, CHILMESH/ELEM2VERT,
            %   CHILMESH/EDGE2VERT, CHILMESH/ELEM2EDGE, CHILMESH/VERT2EDGE,
            %   CHILMESH/VERT2ELEM, CHILMESH/EDGE2ELEM, CHILMESH/ELEM2ELEM.
            %==============================================================
            
            % Check I/O arguments.
            [iT,iQ]	= CM.elemType;                  % Identify tris, quads.
            if isempty(iQ)                          % Triangular mesh.
                CM.Type	= 'Triangular';
            elseif isempty(iT);
                CM.Type = 'Quadrangular';
            else                                    % Tris may exist in the
                CM.Type = 'Mixed-Element';       	% quad/mixed-elem mesh.
                CM.ConnectivityList(iT,end)	= CM.ConnectivityList(iT,1);
            end
            if ~isfloat(CM.ConnectivityList)        % Ensure data type.
                CM.ConnectivityList	= double(CM.ConnectivityList);
            end
            
            % Identify edges of the mesh.
            edges  = cell(1,size(CM.ConnectivityList,2));
            col = [1:length(edges),1];
            for idx = 1:length(col)-1
                edges{idx}	= [CM.ConnectivityList(:,col(idx)),...
                    CM.ConnectivityList(:,col(idx+1))];
            end
            
            % Invert element-to-vertex table then vertex-to-element table.
            numElems	= (1:CM.nElems)';        	% Row-vector for speed.
            vert2elem0	= ...        	% Vert2Elem table.
                sparse(edges{1}(:,1),edges{1}(:,2),numElems,CM.nVerts,CM.nVerts) + ...
                sparse(edges{2}(:,1),edges{2}(:,2),numElems,CM.nVerts,CM.nVerts) + ...
                sparse(edges{3}(:,1),edges{3}(:,2),numElems,CM.nVerts,CM.nVerts);
            if ~strcmp(CM.Type,'Triangular')         % Account for 4th vertex.
                vert2elem0  = vert2elem0 +...
                    sparse(edges{4}(:,1),edges{4}(:,2),numElems,CM.nVerts,CM.nVerts);
            end
            vert2Elem	= fliplr(sort(vert2elem0,2));%#ok<FLPST>
            [~,vert2Elemj]	= find(vert2Elem);
            
            % Construct edge-to-vertex table.
            [edge2Vert1,edge2Vert2]	= find((vert2elem0 - vert2elem0') > 0);
            CM.nEdges	= size(edge2Vert1,1);
            
            % Construct element-to-edge table.
            elem2Edge   = zeros(size(CM.ConnectivityList),'uint32');
            for idx	= 1:size(elem2Edge,2)
                [~,iEdges,iElems]	= intersect(... % Find common edges.
                    [edge2Vert1,edge2Vert2],edges{idx},'rows');
                elem2Edge(iElems,idx)	= iEdges;
                [~,iEdges,iElems]	= intersect(...
                    [edge2Vert2,edge2Vert1],edges{idx},'rows');
                elem2Edge(iElems,idx)	= iEdges;
            end
            
            % Invert edge-to-vertex table to get vertex-to-edge table.
            numEdges	= (1:CM.nEdges)';
            vert2Edge0	= ...
                sparse(edge2Vert1,edge2Vert2,numEdges,CM.nVerts,CM.nVerts) +....
                sparse(edge2Vert2,edge2Vert1,numEdges,CM.nVerts,CM.nVerts);
            vert2Edge	= fliplr(sort(vert2Edge0,2)); %#ok<FLPST>
            [~,vert2Edgej]	= find(vert2Edge);
            
            % Intersect edge2vert and vert2elem tables to get edge2elem table.
            edge2Elem	= [...
                vert2elem0(sub2ind(size(vert2elem0),edge2Vert1,edge2Vert2)),...
                vert2elem0(sub2ind(size(vert2elem0),edge2Vert2,edge2Vert1))];
            
            % Store downward adjacencies.
            CM.Adjacencies.Elem2Vert	= uint32(CM.ConnectivityList);
            CM.Adjacencies.Edge2Vert	= uint32([edge2Vert1,edge2Vert2]);
            CM.Adjacencies.Elem2Edge	= elem2Edge;
            
            % Store upward adjacencies.
            CM.Adjacencies.Vert2Edge	= vert2Edge(:,1:max(vert2Edgej));
            CM.Adjacencies.Vert2Elem	= vert2Elem(:,1:max(vert2Elemj));
            CM.Adjacencies.Edge2Elem	= uint32(full(edge2Elem));
            
            % Convert connectivity list to unsigned int.
            CM.ConnectivityList	= CM.Adjacencies.Elem2Vert;
            if size(CM.ConnectivityList,2) == 4     % Remove redundant
                CM.ConnectivityList(iT,4)   = 0;    % vertices in tri(s).
            end
        end
        
        % Identify Layers of Mesh
        function CM	= meshLayers(CM)
            %==============================================================
            %MESHLAYERS Discretizes mesh into layers.
            %   CM = MESHLAYERS(CM) for mesh class CM returns the layers of
            %   the mesh, where CM.L is a structure containing fields for
            %   the face (element) and vertex details of each mesh layer.
            %
            %   The fields of structure L are:
            %       L.OE: Outer elements in each layer of the mesh.
            %       L.IE: Inner elements in each layer of the mesh.
            %       L.OV: Outer vertices in each layer of the mesh.
            %       L.IV: Inner vertices in each layer of the mesh.
            %
            %   MESHLAYERS discretizes the layers of M by first beginning
            %   with the boundaries edge segments and iteratively advancing
            %   away from these boundaries toward the innermost regions of
            %   the mesh. MESHLAYERS defines the vertices of the mesh's
            %   boundaries as the initial "Outer Vertices" of the first
            %   layer of the mesh. All elements who's connectivity
            %   (Elem2Vert) is comprised by at least one of these vertices
            %   are collectively define the first layer of the mesh. Given
            %   the vertices that define these elements, the "Inner
            %   Vertices" are identified from a set difference with the
            %   Outer Vertices. The Inner Vertices of some layer L then
            %   become the Outer Vertices of layer L-1 and the procedure
            %   repeats itself until no new mesh layers can be defined.
            %
            %   Elements with at least one edge defined by two outer
            %   vertices are called "Outer Elements". The set difference
            %   between all elements of a mesh layer and its Outer Elements
            %   gives the "Inner Elements" of the layer, i.e., the
            %   designations of Outer Elements and Inner Elements of the
            %   same layer are mutually exclusive.
            %
            %   See also CHILMESH.
            %==============================================================
            
            % Get edges and edge neighbors (Edge2Verts and Edges2Elems).
            Edge2VertIDs	= CM.edge2Vert;
            Edge2ElemIDs	= CM.edge2Elem;
            
            %%%  Future edit: Check if open ocean boundary(ies) exist.
%             open_bound_exists	= false;          	% Default: false.
%             if length(CM.BE) == 4              	% Open boundary exists.
%                 open_bound_exists	= true;
%                 iE_open	= sum(ismember(CM.E,CM.BE{4}),2) == 0;
%             end
            
            % Identify layers of CM via inner/outer elements, vertices.
            VertIDs	= (1:size(CM.Points,1))';   	% Index vertices, elems.
            ElemIDs = (1:size(CM.ConnectivityList,1))';
            CM.Layers	= struct('OE',[],'IE',[],'OV',[],'IV',[],'bEdgeIDs',[]);
            iL  = 1;                                % Current Layer.
            while any(Edge2ElemIDs(:) > 0)      	% Until all elems flagged.
                % Identify outer (boundary) edges of layer iL.
                if iL == 1	% Boundary of CM defined by 1st layer's OV.
                    % Check if open ocean boundary(ies) exist.
                    iLbEdgeIDs	= CM.boundaryEdges;	% Index boundary edges.
                    CM.Layers.OV	= [CM.Layers.OV,{VertIDs(ismember(...
                        VertIDs,Edge2VertIDs(iLbEdgeIDs,:)))}];
                    
                else	% OV of iL defined by edges with only 1 neighbor.
                    % Check if open ocean boundary(ies) exist.
                    iLbEdgeIDs	= sum(Edge2ElemIDs > 0,2) == 1;
                    CM.Layers.OV	= [CM.Layers.OV,...
                        unique(Edge2VertIDs(iLbEdgeIDs,:))];
                end
                CM.Layers.bEdgeIDs  = [CM.Layers.bEdgeIDs,{iLbEdgeIDs}];
% CM.plotPoint(CM.Layers.OV{iL},'color','g');pause
                % Outer Elements of iL are the neighbors of iLbEdgeIDs.
                CM.Layers.OE	= [CM.Layers.OE,{ElemIDs(ismember(...
                    ElemIDs,Edge2ElemIDs(iLbEdgeIDs,:)))}];
% CM.plotElem(CM.Layers.OE{iL},'color','m');pause
                % Flag edges in Edge2ElemIDs used by OE.
                Edge2ElemIDs(ismember(Edge2ElemIDs,CM.Layers.OE{end}))	= 0;
                
                % Get all edges associated with iL's OV.
                iLbEdgeIDs	= sum(ismember(Edge2VertIDs,CM.Layers.OV{end}),2) > 0;
                
                % Inner Elements of iL are the neighbors of iLbEdgeIDs.
                CM.Layers.IE	= [CM.Layers.IE,{ElemIDs(ismember(...
                    ElemIDs,Edge2ElemIDs(iLbEdgeIDs,:)))}];
% CM.plotElem(CM.Layers.IE{iL},'color','c');pause
                % Flag used edges of Edge2VertIDs and Edge2ElemIDs.
                Edge2VertIDs(ismember(Edge2VertIDs,CM.Layers.OV{end}))	= 0;
                Edge2ElemIDs(ismember(Edge2ElemIDs,CM.Layers.IE{end}))	= 0;
                
                % Inner Vertices of iL are the vertices of edges defined by iLbEdgeIDs.
                CM.Layers.IV	= [CM.Layers.IV,...
                    {setdiff(CM.ConnectivityList([CM.Layers.OE{iL};...
                    CM.Layers.IE{iL}],:),CM.Layers.OV{iL})}];
% CM.plotPoint(CM.Layers.IV{iL},'color','b');pause

                % Update current layer.
                iL	= iL + 1;
            end
            CM.nLayers	= iL - 1;                   % Number of layers.
        end
        
        %% ==========================Mesh Queries==========================
        % Boundary Edges
        function [EdgeIDs,iBC]	= boundaryEdges(CM,BCIDs)
            %BOUNDARYEDGES Triangulation facets referenced by only 1 face.
            %   EdgeIDs = BOUNDARYEDGES(CM) for mesh class CM returns the
            %   edges that define the free boundary(ies) of the mesh, where
            %   EdgeIDs is an Ex1 vector of indices or logicals with
            %   respect to all of the edges in the mesh, i.e. 1:CM.nEdges.
            %
            %   [EdgeIDs,iBC] = BOUNDARYEDGES(CM) also returns the indices
            %   to specific boundary conditions for each boundary edge in
            %   the mesh, where iBC is an Ex1 vector.
            %
            %   [EdgeIDs,iBC] = BOUNDARYEDGES(CM,BCIDs) returns the
            %   boundary edges corresponding to the specified boundary
            %   conditions BCIDs, where BCIDs is....
            %
            %   Example:
            %
            %   See also TRIANGULATION/FREEBOUNDARY, CHILmesh.
            %==============================================================
            
            % Check I/O arguments.
            %             if nargin == 1                          % If no BC were
            %                 BCIDs	= zeros(CM.nEdges,1,'uint32');% inputted, select
            %                 BCIDs(:)	= 1:CM.nElems;          % all BC.
            %             end
            
            % Identify each edge with only one neighbor.
            EdgeIDs	= find(sum(CM.Adjacencies.Edge2Elem == 0,2) > 0);
            
            % Identify boundary condition of each boundary edge.
            if nargout == 2
                iBC = 1;
            end
        end
        
        % Boundary Vertices
        function [VertIDs,iBC]  = boundaryVerts(CM,BCIDs)
            %BOUNDARYEDGES Triangulation facets referenced by only 1 face.
            %   EdgeIDs = BOUNDARYEDGES(CM) for mesh class CM returns the
            %   edges that define the free boundary(ies) of the mesh, where
            %   EdgeIDs is an Ex1 vector of indices or logicals with
            %   respect to all of the edges in the mesh, i.e. 1:CM.nEdges.
            %
            %   [EdgeIDs,iBC] = BOUNDARYEDGES(CM) also returns the indices
            %   to specific boundary conditions for each boundary edge in
            %   the mesh, where iBC is an Ex1 vector.
            %
            %   [EdgeIDs,iBC] = BOUNDARYEDGES(CM,BCIDs) returns the
            %   boundary edges corresponding to the specified boundary
            %   conditions BCIDs, where BCIDs is....
            %
            %   Example:
            %
            %   See also TRIANGULATION/FREEBOUNDARY, CHILmesh.
            %==============================================================
            
            % Check I/O arguments.
            %             if nargin == 1                          % If no BC were
            %                 BCIDs	= zeros(CM.nEdges,1,'uint32');% inputted, select
            %                 BCIDs(:)	= 1:CM.nElems;          % all BC.
            %             end
            
            % Identify each vert with only 2 incident edges.
            VertIDs	= unique(CM.edge2Vert(CM.boundaryEdges));
            
            % Identify boundary condition of each boundary edge.
            if nargout == 2
                iBC = 1;
            end
        end
        
        % Element Centroid
        function [X,Y,Z]	= centroid(CM,ElemIDs)
            %CENTROID XYZ-Coordinates of element's centroid in the mesh.
            %   [X,Y,Z] = CENTROID(CM) returns the euclidean coordinates
            %   of every polygonal element's centroid in the mesh. The
            %   centroid is computed as the arithmetic mean of the
            %   coordinates that define the polygon's nodal connectivity.
            %
            %   [X,Y,Z] = CENTROID(CM,ElemIDs) returns the coordinates of
            %   the centroids identified by ElemIDs, which is a Nx1 or 1xN
            %   list of indices or logicals with respect to
            %   CM.ConnectivityList and N the number of analyzed elements
            %
            %   CENTROID(MESH) assumes that any mixed-element mesh defined
            %   by CM is defined by a 4-node connectivity, where triangular
            %   elements have psuedo-quadrilateral connectivity, i.e., the
            %   first and second nodes are redundant.
            %
            %   See also CHILMESH, TRIANGULATION\INCENTER,
            %   TRIANGULATION\CIRCUMCENTER.
            %==============================================================
            
            % Check I/O argument(s).
            if nargin == 1                          % If no elements were
                ElemIDs = zeros(CM.nElems,1,'uint32');% inputted, select
                ElemIDs(:)	= 1:CM.nElems;          % all elements.
            end
            
            % Discretize mesh into triangles and quadrilaterals.
            [iT,iQ]	= CM.elemType(ElemIDs);         % Tris, quads in mesh.
            iTb	= ismember(ElemIDs,iT);             % Boolean wrt ElemIDs.
            iQb	= ismember(ElemIDs,iQ);             % Boolean wrt ElemIDs.
            
            % Arithmetic mean position of all elements' vertices' coords.
            [x,y,z]	= CM.elemCoordinates(ElemIDs);  % Coordinates of elems.
            X   = zeros(length(ElemIDs),1);         % Compute X-coordinates.
            X(iTb)   = mean(x(iTb,1:3),2);         	% ^ for triangles.
            X(iQb)   = mean(x(iQb,:),2);            % ^ for quadrilaterals.
            if nargout >= 2                         % Compute Y-coordinates.
                Y   = X;
                Y(iTb)   = mean(y(iTb,1:3),2);    	% ^ for triangles.
                Y(iQb)   = mean(y(iQb,:),2);       	% ^ for quadrilaterals.
                if nargout == 3
                    Z   = Y;
                    Z(iTb)	= mean(z(iTb,1:3),2); 	% ^ for triangles.
                    Z(iQb)	= mean(z(iQb,:),2);   	% ^ for quadrilaterals.
                end
            end
        end
        
        % Diagonals of Quadrilateral                                        * Doesn't accomodate for zero-indices yet.
        function [VertIDs,ElemIDs]	= diagonals(CM,varargin)
            %DIAGONALS Diagonal of quadrilateral element in the mesh.
            %   VertIDs = DIAGONALS(CM) for mesh class CM returns a Nx2
            %   cell array VertIDs containing the diagonals of every
            %   element in the mesh, where N is the number of elements in
            %   the mesh. There are no diagonals for triangles, instead the
            %   corresponding rows of VertIDs contain cells of [0 0].
            %
            %   VertIDs = DIAGONALS(CM,ElemIDs) returns the diagonals of
            %   the elements identified by ElemIDs, which is a Nx1 or 1xN
            %   list of indices or logicals with respect to the
            %   CM.ConnectivityList and N is the number of analyzed elements.
            %
            %   Example:
            %   CM.ConnectivityList = [1 1 2 3; 4 5 6 7] % 1 tri, 1 quad.
            %   VertIDs = DIAGONALS(CM,2)
            %   % returns
            %   VertIDs = {[4 6], [5 7]};
            %
            %   See also CHILMESH, CHILMESH/MEDIAN.
            %==============================================================
            
            % Check I/O arguments.
            p   = inputParser;
            defaultElemIDs	= zeros(CM.nElems,1,'uint32');
            defaultElemIDs(:)	= 1:CM.nElems;   	% All elems.
            allowableStore	= {'Cell','cell','cel','ce','c',...
                'matrix','matri','matr','mat','ma','m'};
            addRequired(p,'CM',@isCHILmesh);
            addParameter(p,'Store','matrix',@(x) any(validatestring(x,allowableStore)));
            addParameter(p,'ElemIDs',defaultElemIDs,@isnumeric);
            parse(p,CM,varargin{:});            	% Parse inputs.
            
            % Diagonals of quadrilateral(s) only.
            iT	= CM.elemType(p.Results.ElemIDs);
            if length(iT) == CM.nElems              % Mesh is only tris.
                return                              % End function.
            else
                % Adjust indices to tris.
                iT  = (iT.*2)';
                iT  = reshape([iT-1; iT],length(iT)*2,1);
            end
            
            % Identify diagonals.
            idx0    = p.Results.ElemIDs == 0;     	% Search for zeros.
            VertIDs	= zeros(length(p.Results.ElemIDs)*2,2,'uint32');
            VertIDs(1:2:end,:)	= CM.ConnectivityList(p.Results.ElemIDs(~idx0),[1 3]);
            VertIDs(2:2:end,:)	= CM.ConnectivityList(p.Results.ElemIDs(~idx0),[2 4]);
            VertIDs(iT,:)	= 0;                 	% Tris dont have diags.
            
            % Adjust output, if necessary.
            if any(strcmpi(p.Results.Store,{'Cell','cell','cel','ce','c'}))
                % Convert to a cell array output.
                VertIDs	= mat2cell([VertIDs(1:2:end,:) VertIDs(2:2:end,:)],...
                    ones(length(VertIDs)/2,1),[2 2]);
            end
            
            % Fulfill second output of indices to elems matching VertIDs.
            if nargout == 2
                ElemIDs  = p.Results.ElemIDs';
                ElemIDs  = reshape([ElemIDs; ElemIDs],length(ElemIDs)*2,1);
            end
        end
        
        % Retrieve Edge-to-Element Adjacency
        function ElemIDs	= edge2Elem(CM,EdgeIDs)
            %EDGE2ELEM Elements attached to edge(s).
            %   ElemIDs = EDGE2ELEM(CM,EdgeIDS) gives all elements (0 to 4)
            %   of the mesh attached to the edge(s) specified by EdgeIDs.
            %   EdgeIDs can be an 1xE or Ex1 vector of indices or logicals,
            %   where E is the number of edges. ElemIDs is an Ex2 matrix
            %   of indices to the M.Connectivity list, where E is the
            %   number of edge indices (size(EdgeIDs,1)). A value of "0"
            %   indicates no attached element on that side of the edge.
            %
            %   See also CHILMESH/BUILDADJACENCIES, CHILMESH/EDGE2VERT,
            %   CHILMESH.CONNECTIVITYLIST, TRIANGULATION/EDGEATTACHMENTS.
            %==============================================================
            
            if nargin == 1                          % If no edges were
                EdgeIDs = zeros(CM.nEdges,1,'uint32');% inputted, select
                EdgeIDs(:)	= 1:CM.nEdges;          % all edges.
            end
            idx0	= EdgeIDs == 0;                 % Search for zeros.
            ElemIDs = zeros(numel(EdgeIDs),2,'uint32');
            ElemIDs(~idx0,:)	= CM.Adjacencies.Edge2Elem(EdgeIDs(~idx0),:);
        end
        
        % Retrieve Edge-to-Vertex Adjacency
        function VertIDs	= edge2Vert(CM,EdgeIDs)
            %EDGE2VERT Points defining an edge(s).
            %   VertIDs = EDGE2VERT(CM,EdgeIDs) gives the 2 vertices of the
            %   mesh that define the edge(s) specified by EdgeIDs. EdgeIDs
            %   can be an 1xE or Ex1 vector of indices or logicals, where E
            %   is the number of edge indices (size(EdgeIDs,1)). VertIDs is
            %   a Ex2 matrix of indices to the CM.Points list.
            %
            %   See also CHILMESH/BUILDADJACENCIES, CHILMESH/EDGE2ELEM,
            %   CHILMESH/POINTS, TRIANGULATION/EDGES.
            %==============================================================
            
            if nargin == 1                          % If no edges were
                EdgeIDs = zeros(CM.nEdges,1,'uint32');% inputted, select
                EdgeIDs(:)	= 1:CM.nEdges;          % all edges.
            end
            idx0	= EdgeIDs == 0;                 % Search for zeros.
            VertIDs = zeros(numel(EdgeIDs),2,'uint32');
            VertIDs(~idx0,:)	= CM.Adjacencies.Edge2Vert(EdgeIDs(~idx0),:);
        end
        
        % Coordinates of Edge
        function [X,Y,Z]	= edgeCoordinates(CM,EdgeIDs)
            %EDGECOORDINATES Coordinates of edge's vertices in the mesh.
            %   [X,Y,Z] = EDGECOORDINATES(CM,EdgeIDs) returns the euclidean
            %   coordinates of edges in the mesh.
            %
            %   See also CHILMESH, CHILMESH/EDGEMIDPOINT,
            %   CHILMESH/VERTCOORDINATES, CHILMESH/ELEMCOORDINATES.
            %==============================================================
            
            if nargin == 1                          % If no edges were
                EdgeIDs = zeros(CM.nEdges,1,'uint32');% inputted, select
                EdgeIDs(:)	= 1:CM.nEdges;          % all edges.
            end
            idx0    = EdgeIDs == 0;                 % Search for zeros.
            VertIDs = zeros(numel(EdgeIDs),2);      % Get verts of edges.
            VertIDs(~idx0,:)	= CM.Adjacencies.Edge2Vert(EdgeIDs(~idx0),:);
            [X,Y,Z]	= CM.vertCoordinates(VertIDs);  % Get coordinates.
            
            % Reshape according to edges.
            X	= reshape(X,numel(EdgeIDs),2);
            if nargout >= 2
                Y	= reshape(Y,numel(EdgeIDs),2);
                if nargout == 3
                    Z	= reshape(Z,numel(EdgeIDs),2);
                end
            end
        end
        
        % Length of Edge
        function EL	= edgeLength(CM,EdgeIDs)
            %EDGELENGTH Length of edge in the mesh.
            %   EL = EDGELENGTH(CM) for mesh class CM returns a vector of
            %   the length of each edge in the mesh.
            %
            %   EL = EDGELENGTH(CM,EdgeIDs) returns the length the indexed
            %   edges EdgeIDs in the mesh.
            %
            %   Example:
            %
            %   See also CHILMESH, CHILMESH/EDGECOORDINATES.
            %==============================================================
            
            % Assign edges for analysis.
            if nargin == 1                          % If no edges were
                EdgeIDs = zeros(CM.nEdges,1,'uint32');% inputted then
                EdgeIDs(:)	= 1:CM.nEdges;          % select all of them.
            end
            
            % Compute edge lengths.
            [X,Y]	= CM.edgeCoordinates(EdgeIDs(:));
            EL  = hypot(X(:,1) - X(:,2),Y(:,1) - Y(:,2));
        end
        
        % Midpoint of Edge
        function [X,Y,Z]    = edgeMidpoint(CM,EdgeIDs)
            %EDGEMIDPOINT Coordinates of an edge's midpoint in mesh.
            %   [X,Y,Z] = EDGEMIDPOINT(CM,EdgeIDs) returns the euclidean
            %   coordinates of the midpoint of edges in the mesh, where
            %   EdgeIDs is an Ex1 or 1xE vector of indices or logicals with
            %   respect to CM.Adjacencies.Edge2Vert.
            %
            %   See also CHILMESH, CHILMESH/EDGELENGTH,
            %   CHILMESH/EDGECOORDINATES.
            %==============================================================
            
            % Check I/O argument(s).
            if nargin == 1                          % If no edges were
                EdgeIDs = zeros(CM.nEdges,1,'uint32');% inputted, select
                EdgeIDs(:)	= 1:CM.nEdges;          % all edges.
            end
            [X,Y,Z]   = CM.edgeCoordinates(EdgeIDs);% Edge's coordinates.
            
            % Average the coordinates.
            X   = mean(X,2);    Y   = mean(Y,2);    Z   = mean(Z,2);
        end
        
        % Order of Edges Around Vertex
        function EdgeIDs    = edgeOrder(CM,VertIDs)
            %EDGEORDER Coordinates of an edge's midpoint in mesh.
            %   [X,Y,Z] = EDGEORDER(M,EdgeIDs) returns the euclidean
            %   coordinates of the midpoint of edges in the mesh, where
            %   EdgeIDs is an Ex1 or 1xE vector of indices or logicals with
            %   respect to CM.Adjacencies.Edge2Vert.
            %
            %   See also CHILMESH, CHILMESH/VERT2EDGE.
            %==============================================================
            
            % Check I/O argument(s).
            if nargin == 1                          % If no edges were
                VertIDs = zeros(CM.nVerts,1,'uint32');% inputted, select
                VertIDs(:)	= 1:CM.nVerts;          % all edges.
                
            elseif size(VertIDs,2) ~= 1             % VertIDs must be row.
                VertIDs = VertIDs(:);
            end
            VertIDsXY   = CM.Points(VertIDs,1:2);   % Datum coordinates.
            
            % Get respective edges of each VertID.
            EdgeIDs	= CM.vert2Edge('vertids',VertIDs,'store','m');
            
            % Get verts of each EdgeID.
            polyVertIDs	= reshape(CM.edge2Vert(EdgeIDs')',size(EdgeIDs,2)*2,size(EdgeIDs,1))';
            polyVertIDs(double(polyVertIDs) - double(VertIDs) == 0)	= 0;
            
            % Define the polygon of vertices surrounding each VertID.
            S	= size(polyVertIDs);                % Identify non-unique
            [B,C]	= sort(polyVertIDs,2);          % values in each row.
            D	= [false(S(1),1),diff(B,1,2)==0];
            R	= (1:S(1))'*ones(1,S(2));
            ipolyVertIDs	= sub2ind(S,R(D),C(D));
            polyVertIDs(ipolyVertIDs)	= 0;    	% Flag non-uniques.
            idx0= polyVertIDs == 0;                 % Identify all zeros.
            
            % Sort angles (edges) of each vert wrt polar coordinates.
            [X,Y]	= CM.vertCoordinates(polyVertIDs);
            X   = reshape(X,S(1),S(2));             % Match shape of polyVertIDs.
            Y   = reshape(Y,S(1),S(2));
            X(idx0)	= NaN;                          % Flag invalid inputs.
            Y(idx0) = NaN;
            X   = X - VertIDsXY(:,1);               % Shift to datum(s)
            Y   = Y - VertIDsXY(:,2);
            theta   = reshape(cart2pol(X,Y),S(1),S(2));
            [~,itheta]  = sort(theta,2,'ascend');	% Ensures CCW orientation.
            
            % Adjust theta wrt polyEdgeIDs (convert to linear indexing).
            itheta(:,S(2)/2+1:end)	= [];           % Remove redundant indices.
            for idx = 1:S(2)/2
                itheta(:,idx)	= sub2ind(S,(1:length(VertIDs))',itheta(:,idx));
            end
            
            % Reorder EdgeIDs (output).
            polyEdgeIDs	= polyVertIDs;
            polyEdgeIDs(:,1:2:end)  = EdgeIDs;
            polyEdgeIDs(:,2:2:end)  = EdgeIDs;
            polyEdgeIDs(idx0)   = 0;
            EdgeIDs	= polyEdgeIDs(itheta);
        end
        
        % Retrieve Element-to-Edge Adjacency
        function EdgeIDs	= elem2Edge(CM,ElemIDs)
            %ELEM2EDGE Edges defining an element(s).
            %   EdgeIDs = ELEM2EDGE(CM,ElemIDs) gives all edges (3 or 4) of
            %   the mesh that are attached to the element(s) specified by
            %   ElemIDs. ElemIDs can be an 1xC or Cx1 vector of indices or
            %   logicals, where C is the number of element indices
            %   (size(ElemIDs,1)). EdgeIDs is a CxD matrix of indices to
            %   the Edge2Vert list of CM, where D is the max number of
            %   sides for any polygon in the mesh (3 -> tris, 4 -> quads).
            %
            %   See also CHILMESH/BUILDADJACENCIES, CHILMESH/ELEM2VERT,
            %   CHILMESH/EDGE2VERT.
            %==============================================================
            
            if nargin == 1                          % If no elements were
                ElemIDs = zeros(CM.nElems,1,'uint32');% inputted, select
                ElemIDs(:)	= 1:CM.nElems;          % all elements.
            end
            
            % Identify number of edges per element.
            nE  = 3;                                
            if size(CM.ConnectivityList,2) == 4     % Mixed/quad mesh.
                nE  = 4;
            end
            
            idx0	= ElemIDs == 0;                 % Search for zeros.
            EdgeIDs = zeros(numel(ElemIDs),nE,'uint32');
            EdgeIDs(~idx0,:)	= CM.Adjacencies.Elem2Edge(ElemIDs(~idx0),:);
        end
        
        % Neighbors of Element                                              *Doesn't accomodate for zero-indices yet.
        function ElemIDs	= elem2Elem(CM,ElemIDs)
            %ELEM2ELEM Neighbor(s) of an element in the mesh.
            %   ElemIDs = ELEM2ELEM(CM) for mesh class CM returns the
            %   neighbors of every element in the mesh as an Nx3 matrix
            %   ElemIDs, where N is the number of elements in the mesh and
            %   E is the maximum number of edges an element in the mesh can
            %   have (N = 3 when the mesh is all triangles, N = 4 when the
            %   mesh is all quadrilaterals or mixed-element).
            %
            %   ElemIDs = ELEM2ELEM(CM,ElemIDs) returns the neighbors of
            %   elements specifed by input ElemIDs. For sake of convention,
            %   the Ex1 or 1xE input of indices or logicals is overwritten
            %   to become the ExN output ElemIDs.
            %
            %   Example:
            %   load _____ % Mixed-Element test case.
            %   CM = buildAdjacencies(CM);
            %   ElemIDs = elem2Elem(CM,1:2:5);
            %   % returns
            %   ElemIDs = [4 0 0 2; 6 2 0 0; 4 0 2 6;
            %
            %   See also CHILMESH, CHILMESH/ELEM2EDGE, CHILMESH/EDGE2ELEM,
            %   TRIANGULATION/NEIGHBORS.
            %==============================================================
            
            % Check I/O argument(s).
            if nargin == 1                          % If no elements were
                ElemIDs = zeros(CM.nElems,1,'int32');% inputted, select
                ElemIDs(:)	= 1:CM.nElems;          % all elements.
            elseif ~isa(ElemIDs,'int32')
                ElemIDs	= int32(ElemIDs);        	% Convert to int32
            end
            if size(ElemIDs,2) ~=1                  % Make sure ElemIDs is
                ElemIDs	= ElemIDs';                 % a column vector.
            end
            
            % Get edges of each element, then get elements of those edges.
            maxNumEdges	= size(CM.ConnectivityList(ElemIDs,:),2);
            edgeOfElem	= CM.elem2Edge(ElemIDs);	% Each element's edges.
            neighbors   = zeros(length(ElemIDs),maxNumEdges*2,'int32');
            icol	= 1:2;                          % Column indices.
            for idx	= 1:maxNumEdges                 % edge2Elem of each edge.
                neighbors(:,icol) = CM.edge2Elem(edgeOfElem(:,idx));
                icol	= icol + 2;                 % Next 2 columns.
            end
            
            % Flag, remove self-indices (elem #1 isn't its own neighbor).
            neighbors(neighbors - repmat(ElemIDs,1,maxNumEdges*2) == 0)	= 0;
            
            % Rotate nonzeros to the first (3 or 4) columns.
            stop = false;
            col = maxNumEdges*2;
            while ~stop
                idx	= neighbors(:,col) == 0;
                if all(idx)
                    col = col-1;
                    if col == maxNumEdges
                        stop    = true;
                    end
                else
                    neighbors(~idx,:)	= neighbors(~idx,...
                        [col,setdiff(1:maxNumEdges*2,col)]);
                end
            end
            ElemIDs	= uint32(neighbors(:,1:maxNumEdges));
        end
        
        % Retrieve Element-to-Vertex Adjacency
        function VertIDs	= elem2Vert(CM,ElemIDs)
            %ELEM2VERT Vertices that define the element(s).
            %   VertIDS = ELEM2VERT(CM,ElemIDS) gives the (3 or 4) vertices
            %   of the mesh that define the element(s) specified by ElemIDs.
            %   ElemIDs can be an 1xN or Nx1 vector of indices or logicals,
            %   where N is the number of elements. VertIDs is an NxP matrix
            %   of indices to the M.Points list, where P is the number of
            %   points defining the elements.
            %
            %   If VertIDs contains triangular and quadrangular elements,
            %   P = 4 and the triangle elements will have a redundant node
            %   in the first column.
            %
            %   See also CHILMESH/BUILDADJACENCIES, CHILMESH/EDGE2ELEM,
            %   CHILMESH/POINTS, TRIANGULATION/CONNECTIVITYLIST.
            %==============================================================
            
            if nargin == 1                          % If no elements were
                ElemIDs = zeros(CM.nElems,1,'uint32');% inputted, select
                ElemIDs(:)	= 1:CM.nElems;          % all elements.
            end
            idx0	= ElemIDs == 0;                 % Search for zeros.
            VertIDs = zeros(numel(ElemIDs),size(CM.ConnectivityList,2),'uint32');
            VertIDs(~idx0,:)	= CM.Adjacencies.Elem2Vert(ElemIDs(~idx0),:);
        end
        
        % Coordinates of Element
        function [X,Y,Z]	= elemCoordinates(CM,ElemIDs)
            %ELEMCOORDINATES Coordinates of element's connectivity in mesh.
            %   [X,Y,Z] = ELEMCOORDINATES(CM,ElemIDS) returns the euclidean
            %   coordinates of elements in the mesh, where (X,Y,Z) are the
            %   respective coordinates of each element's vertices and
            %   size(X) = size(Y) = size(Z) =
            %   size(CM.ConnectivityList(ElemIDs,:)).
            %
            %   See also CHILMESH, CHILMESH/CENTROID,
            %   CHILMESH/VERTCOORDINATES, CHILMESH/EDGECOORDINATES.
            %==============================================================
            
            if nargin == 1                          % If no elements were
                ElemIDs = zeros(CM.nElems,1,'uint32');% inputted, select
                ElemIDs(:)	= 1:CM.nElems;          % all elements.
                
            elseif size(ElemIDs,2) ~= 1
                ElemIDs	= ElemIDs(:);
            end
            
            % Get coordinates of each VertIDs of ElemIDs.
            idx0    = ElemIDs == 0;                 % Search for zeros.
            VertIDs = zeros(numel(ElemIDs),size(CM.ConnectivityList,2));
            VertIDs(~idx0,:)	= CM.Adjacencies.Elem2Vert(ElemIDs(~idx0),:);
            [X,Y,Z]	= CM.vertCoordinates(VertIDs);  % Get coordinates.
            
            % Reshape according to connectivity.
            X	= reshape(X,numel(ElemIDs),size(CM.ConnectivityList,2));
            if nargout >= 2
                Y	= reshape(Y,numel(ElemIDs),size(CM.ConnectivityList,2));
                if nargout == 3
                    Z	= reshape(Z,numel(ElemIDs),size(CM.ConnectivityList,2));
                end
            end
        end
        
        % Quality of Element
        function [Quality,Theta]	= elemQuality(CM,varargin)
            %ELEMQUALITY Quality measurement of element in the mesh.
            %   Quality = ELEMQUALITY(CM) for mesh class CM returns a Nx1
            %   vector Quality containing the quality of every element in
            %   the mesh, where Quality is a quantification of the largest
            %   angular skewness within each element and N is the number of
            %   elements in the mesh.
            %
            %   [Quality,Theta] = ELEMQUALITY(CM) also returns the interior
            %   angles of each element.
            %
            %   Quality = ELEMQUALITY(CM,MEASUREMENT) measures the quality
            %   of every element with respect to the quantification
            %   MEASUREMENT. The default quality measurement for
            %   ELEMQUALITY is angular skewness. Other allowable inputs for
            %   MEASUREMENT include:
            %       1. 'skewness'   - measures largest deviation of
            %       interior angles from an ideal angle (60 tris, 90 quads)
            %       2.
            %
            %   Quality = ELEMQUALITY(CM,ElemIDs) returns the element
            %   quality of the elements identified by ElemIDs, which is a
            %   Nx1 or 1xN list of indices or logicals with respect to
            %   CM.ConnectivityList and N is the number of elements to be
            %   analyzed.
            %
            %   See also CHILMESH, CHILMESH/INTERIORANGLES,
            %==============================================================
            
            % Check I/O argument(s).
            p   = inputParser;
            defaultElemIDs  = zeros(CM.nElems,1,'uint32');
            defaultElemIDs(:)	= 1:CM.nElems;   	% All elements.
            allowableQualities	= {'skew','skewness','angular skewness'};
            addRequired(p,'CM',@isCHILmesh)
            addOptional(p,'ElemIDs',defaultElemIDs,@(x) isnumeric(x))
            addOptional(p,'Quality','skew',@(x) any(validatestring(x,allowableQualities)))
            parse(p,CM,varargin{:})             	% Parse inputs.
            
            % Discretize elements in mesh by element by type.
            [iT,iQ] = CM.elemType(p.Results.ElemIDs);
            
            % Compute quality.
            switch p.Results.Quality
                case {'skew','skewness','angular skewness'}
                    % Compute the maximum and minimum interior angles.
                    Theta	= CM.interiorAngles(p.Results.ElemIDs);
                    iTris	= ismember(p.Results.ElemIDs,iT);
                    Tmax    = max(Theta(iTris,1:3),[],2);
                    Tmin	= min(Theta(iTris,1:3),[],2);
                    
                    % Compute the equiangular skew of all triangles.
                    Quality = zeros(length(p.Results.ElemIDs),1);
                    Quality(iTris)	= 1-max((Tmax-60)./(180-60),(60-Tmin)./60);
                    
                    % Repeat for all quadrilaterals.
                    iQuads	= ismember(p.Results.ElemIDs,iQ);
                    Qmax    = max(Theta(iQuads,:),[],2);
                    Qmin	= min(Theta(iQuads,:),[],2);
                    Quality(iQuads)	= 1-max((Qmax-90)./(180-90),(90-Qmin)./90);
                    
                    % Account for poor interior angle calculation (e.g concave quads).
                    iZero   = sum(Theta,2) <= 179.99 & iTris;
                    iZero   = iZero | (sum(Theta,2) <= 359.99 & iQuads);
                    Quality(iZero)  = 0;
                    
                otherwise                           % No other quality methods yet.
            end
        end
        
        % Determine Element Type
        function [iT,iQ]	= elemType(CM,ElemIDs)
            %ELEMTYPE Identifies triangular and quadrilateral elements.
            %	iT = ELEMTYPE(CM) for mesh class CM returns indices to
            %	all triangular elements iT in the mesh.
            %
            %   iT = ELEMTYPE(CM,ElemIDs) returns indices iT to the
            %   triangular elements in the mesh that are also in the Nx1 or
            %   1xN vector ElemIDs, where ElemIDs are indices with respect
            %   to CM.ConnectivityList and N is the number of elements.
            %
            %   [iT,iQ] = ELEMTYPE(CM) also returns indices to all
            %   quadrilateral elements iQ in the mesh.
            %
            %	Example:
            %   ElemIDs = [1 3 4 10]
            %   CM.ConnectivityList(iC,:) = [1 2 2 3; 4 5 6 4;
            %   1 2 8 9; 7 7 8 9]; % Three triangles, one quadrilateral.
            %   [iT,iQ] = ELEMTYPE(MESH,ElemIDs);
            %   % returns
            %   iT = 1 2 4; % With respect to ElemIDs.
            %   iQ = 3;
            %   ElemIDs(iT), ElemIDs(iQ) % gives:
            %   [1 3 10], [4]
            %
            %   See also CHILMESH.
            %==============================================================
            
            % Assign elements for analysis.
            if nargin == 1
                ElemIDs = zeros(CM.nElems,1,'uint32');
                ElemIDs(:)	= 1:CM.nElems;
            
            elseif size(ElemIDs,2) ~= 1             % Ensure column array.
                ElemIDs = ElemIDs(:);
            end
            
            % Initialize outputs.
            if size(CM.ConnectivityList,2) == 3     % Mesh is all tris.
                iT	= ElemIDs;
                iQ	= [];                           % No quads in the mesh.
                return
            end
            
            % Identify all triangular faces of ElemIDs.
            idx0    = ElemIDs == 0;
            idxT    = false(numel(ElemIDs),1);
            idxT(~idx0)	= CM.ConnectivityList(ElemIDs(~idx0),1) == CM.ConnectivityList(ElemIDs(~idx0),2) |...
                CM.ConnectivityList(ElemIDs(~idx0),2) == CM.ConnectivityList(ElemIDs(~idx0),3) |...
                CM.ConnectivityList(ElemIDs(~idx0),3) == CM.ConnectivityList(ElemIDs(~idx0),4) |...
                CM.ConnectivityList(ElemIDs(~idx0),4) == CM.ConnectivityList(ElemIDs(~idx0),1) |...
                CM.ConnectivityList(ElemIDs(~idx0),4) == 0;
            iT  = ElemIDs(idxT & ~idx0);
            
            % Index to all quadrilaterals in the mesh.
            if nargout == 2
                iQ  = ElemIDs(~idxT & ~idx0);
            end
        end
        
        % Interior Angles of Element
        function Theta	= interiorAngles(CM,ElemIDs)
            %INTERIORANGLES Interior angles of polygonal element in the mesh.
            %   Theta = INTERIORANGLES(CM) for mesh class CM returns a Nx3
            %   or Nx4 matrix Theta containing the interior angles of every
            %   element in the mesh, where N is the number of elements and
            %   the number of columns is 3 for completely triangular and
            %   4 for completely quadrangular/mixed-element meshes.
            %
            %   Theta = INTERIORANGLES(CM,ElemIDs) returns the interior
            %   angles of the elements identified by ElemIDs, which is a
            %   Nx1 or 1xN list of indices or logicals with respect to
            %   CM.ConnectivityList and N is the number of elements to be
            %   analzed.
            %
            %   See also CHILMESH, CHILMESH/ELEMQUALITY.
            %==============================================================
            
            % Check I/O argument(s).
            if nargin == 1                          % If no elements were
                ElemIDs = zeros(CM.nElems,1,'uint32');% inputted, select
                ElemIDs(:)	= 1:CM.nElems;          % all elements.
                
            elseif size(ElemIDs,2) ~= 1
                ElemIDs = ElemIDs(:);
            end
            
            % Get elem types for ElemIDs.
            [iT,iQ]	= CM.elemType(ElemIDs);
            
            % For triangles: Compute edge lengths then interior angles.
            if ~isempty(iT)                         % Tris exist in mesh.
                % Compute edge lengths.
                EL	= reshape(CM.edgeLength(CM.elem2Edge(iT)'),[],length(iT));
                
                % Compute interior angles from edge lengths.
                Theta   = zeros(numel(ElemIDs),size(CM.ConnectivityList,2));
                iElemIDs= ismember(ElemIDs,iT);
                Theta(iElemIDs,1:3)	= [...
                    acosd((EL(1,:).^2 + EL(3,:).^2 - EL(2,:).^2.)...
                    ./(2.*EL(1,:).*EL(3,:)));...    % Law of cosines.
                    acosd((EL(1,:).^2 + EL(2,:).^2 - EL(3,:).^2.)...
                    ./(2.*EL(1,:).*EL(2,:)));...
                    acosd((EL(2,:).^2 + EL(3,:).^2 - EL(1,:).^2.)...
                    ./(2.*EL(2,:).*EL(3,:)))]';     % Transpose back.
            end
            
            % For quadrilaterals: Create 2 tris then compute compute angles.
            if ~isempty(iQ)
                % Get vectors of each quad.
                %%% This sometimes gives incorrect solutions - see concave  quads.
                xQ  = zeros(length(iQ),6);          % X-coords of quads.
                xQ(:,2:5)	= reshape(CM.Points(...
                    CM.ConnectivityList(iQ,:)',1),4,length(iQ))';
                xQ(:,[1 6])	= [xQ(:,5), xQ(:,2)];   % Clever appending stuff.
                yQ  = xQ;                           % Y-coords of quads.
                yQ(:,2:5)	= reshape(CM.Points(...
                    CM.ConnectivityList(iQ,:)',2),4,length(iQ))';
                yQ(:,[1 6])	= [yQ(:,5), yQ(:,2)];
                
                % Compute angles of each vertex. The way I am computing the
                % angle is by combining two sides to make a triangle and then
                % computing the vertex angle using cosine relationship,
                % therefore elemental orientation (CW or CCW) is irrelevant.
                % Compute lengths of edges 1, 2, and 3.
                iElemIDs(:)	= ismember(ElemIDs,iQ);
                for idx	= 2:5
                    L1	= sqrt((xQ(:,idx-1) - xQ(:,idx)).^2 +...
                        (yQ(:,idx-1) - yQ(:,idx)).^2);
                    L2	= sqrt((xQ(:,idx) - xQ(:,idx+1)).^2 +...
                        (yQ(:,idx) - yQ(:,idx+1)).^2);
                    L3	= sqrt((xQ(:,idx-1) - xQ(:,idx+1)).^2 + ...
                        (yQ(:,idx-1) - yQ(:,idx+1)).^2);
                    
                    % Compute the interior angles of the quads.
                    Theta(iElemIDs,idx-1)	= real(acosd((L1.^2 +...
                        L2.^2 - L3.^2)./(2*L1.*L2)));
                end
            end
        end
        
        % Check Class
        function isCHILmesh(CM)
            %ISCHILMESH True for CHILmesh object class.
            %   ISCHILMESH(CM) returns true if CM is a CHILmesh object mesh
            %   class and false otherwise.
            %
            %   See also CHILMESH, ISA.
            %==============================================================
            
            isa(CM,'CHILmesh');
        end
        
        % Enforce Counter-Clockwise Vertex Orientation
        function [iCW,CM,A]	= isPolyCCW(CM,varargin)
            %ISPOLYCCW Set connectivity to a counter-clockwise orientation.
            %   iCW = ISPOLYCCW(CM) for mesh class CM returns a vector
            %   of indices iCW indicating which elements in the mesh have
            %   clockwise connectivity orientation(s).
            %
            %   iCW = ISPOLYCCW(CM,INDEX) specified which index type the
            %   output iCW is, where the default value of INDEX = 'index'
            %   and the alternative value is INDEX = 'logical'.
            %
            %   iCW = ISPOLYCCW(CM,ElemIDs) identifies elements of index
            %   vector ElemIDs have clockwise oriented connectivities,
            %   where ElemIDs is a 1xN or Nx1 vector of indices and N is
            %   the number of elements that will be analyzed. ISPOLYCCW
            %   assumes that ElemIDs are with respect to
            %   CM.ConnectivityList and values of iCW will be as well.
            %
            %   [iCW,CM] = ISPOLYCCW(CM,ElemIDs) also returns the updated
            %   connectivity list of CM such that all elements in CM
            %   indexed by ElemIDs will be oriented counter-clockwise.
            %
            %   [iCW,CM,A] = ISPOLYCCW(CM,ElemIDs) also returns the signed
            %   aread of each element indexed by ElemIDs. ISPOLYCCW
            %   determines the orientation of polygons by computing the
            %   signed area, where positive values of A correspond to a
            %   counter-clockwise orientation and vice versa.
            %
            %   Example:
            %   CM.ConnectivityList = [1 2 3 4; 5 6 7 8; 9 10 11 12]
            %   % The signed areas of these three elements are 1, -2, and
            %   3, respectively (arbitrary for explanatory purposes).
            %   [iCW,CM,A] = ISPOLYCCW(CM)
            %   % returns
            %   iCW = 2
            %   CM.ConnectivityList = [1 2 3 4; 5 8 7 6; 9 10 11 12];
            %   A = [1 -2 3];
            %   CCWC = [1 2 3 4; 5 8 7 6; 9 10 11 12]
            %   I = 2
            %   A = [1; -5; 6]
            %
            %   See also CHILMESH, CHILMESH/AREA.
            %==============================================================
            
            % Check I/O argument(s).
            p   = inputParser;
            defaultElemIDs  = zeros(CM.nElems,1,'uint32');
            defaultElemIDs(:)	= 1:CM.nElems;   	% All elements.
            allowableIndex	= {'index','logical','Index','Logical','Conn','Connectivity'};
            addRequired(p,'CM',@isCHILmesh);
            addOptional(p,'Index','Index',@(x) any(validatestring(x,allowableIndex)));
            addOptional(p,'ElemIDs',defaultElemIDs,@(x) isnumeric(x));
            parse(p,CM,varargin{:});            	% Parse inputs.
            
            % Compute signed areas.
            A	= CM.signedArea('signed',p.Results.ElemIDs);
            
            % Identify any clockwise (CW) connectivity orientation(s).
            if strcmpi(p.Results.Index,{'logical','Logical'})
                iCW = A < 0;                        % Output logical.
            else                                    % Output indices.
                iCW = find(A < 0);
            end
            
            % Reorient any CW connectivity orientation to CCW.
            if nargout > 1                          % Rotate about node 1.
                CM.ConnectivityList(p.Results.ElemIDs(iCW),:)	=...
                    [CM.ConnectivityList(p.Results.ElemIDs(iCW),1),...
                    fliplr(CM.ConnectivityList(p.Results.ElemIDs(iCW),2:end))];
            end
        end
        
        % Medians of Edge
        function VertIDs	= medians(CM,EdgeIDs)
            %MEDIANS Median(s) of an edge in the mesh.
            %   VertIDs = MEDIANS(CM) for mesh class CM returns the VertIDs
            %   of vertices opposite of each edge in the mesh, where
            %   VertIDs is an Ex2 cell array containing the vertices
            %   opposite of the edge on each of its side. A median for an
            %   edge of a triangle has only one opposite vertex within that
            %   triangle, whereas a median for an edge of quadrilateral has
            %   two opposite vertices within that quadrilateral.
            %
            %   VertIDs = MEDIANS(CM,EdgeIDs) returns the vertices opposite
            %   of each edge specified by EdgeIDs, where EdgeIDs is an Ex1
            %   or 1xE vector of indicies or logicals.
            %
            %   Example:
            %   load ___ % Load example mixed-element mesh.
            %   EdgeIDs = [1;3;5];              % Define edges for analysis
            %   VertIDs = MEDIANS(CM,[1 3 5]);  % Medians of EdgeIDs.
            %   %returns
            %   VertIDs = {};
            %   figure; plot(CM);               % Plot mesh.
            %   [X,Y] = midpoint(CM,EdgeIDs);   % Coordinates of midpoints.
            %   for idx = 1:3
            %   plotEdge(CM,EdgeIDs(i),'Color','r');
            %   plotPoint(CM,'Vertex',VertIDs{idx,1},'Color','g');% Side 1.
            %   plotPoint(CM,'Vertex',VertIDs{idx,2},'Color','b');% Side 2.
            %   plot([X(idx);CM.Points(VertIDs{idx,1},1)]],[Y(idx);CM.Points(VertIDs{idx,1},2)]],'g--');
            %   plot([X(idx);CM.Points(VertIDs{idx,2},1)]],[Y(idx);CM.Points(VertIDs{idx,2},2)]],'b--');
            %   pause                           % Plot each median.
            %   end
            %
            %   See also CHILMESH, CHILMESH/DIAGONALS.
            %==============================================================
            
            % Check I/O argument(s).
            if nargin == 1                          % If no elements were
                EdgeIDs = zeros(CM.nEdges,1,'uint32');% inputted, select
                EdgeIDs(:)	= 1:CM.nEdges;          % all elements.
            end
            nEdgeIDs	= length(EdgeIDs);       	% Number of edges.
            
            % Vertices defining EdgeIDs.
            vertsOfEdgeIDs = CM.edge2Vert(EdgeIDs);
            
            % Elements attached to EdgeIDs.
            ElemIDs	= CM.edge2Elem(EdgeIDs);
            
            % Connectivity of ElemIDs on each side of EdgeIDs.
            connSide1	= CM.ConnectivityList(ElemIDs(:,1),:);
            connSide2   = zeros(nEdgeIDs,size(CM.ConnectivityList,2));
            nonZeros    = ElemIDs(:,2) > 0;         % Non-boundary edges.
            connSide2(nonZeros,:)	= CM.ConnectivityList(ElemIDs(nonZeros,2),:);
            oppositeVertex = [connSide1 ~= vertsOfEdgeIDs(:,1)...
                & connSide1 ~= vertsOfEdgeIDs(:,2) & connSide1 > 0;...
                connSide2 ~= vertsOfEdgeIDs(:,1)...
                & connSide2 ~= vertsOfEdgeIDs(:,2) & connSide2 > 0];
            conn    = [connSide1;connSide2];
            
            % Rotate nonzeros to first two columns.
            [sBoolean,idx]	= sort(oppositeVertex,2,'descend');
            newConn = zeros(size(conn));
            for jdx = 1:2                           % For first two columns.
                for kdx	= 1:size(conn,2)
                    iCol    = idx(:,jdx) == kdx;
                    newConn(iCol,jdx) = conn(iCol,kdx);
                    iZero   = sBoolean(:,jdx) == 0;
                    if any(iZero)
                        newConn(iZero,jdx)	= 0;
                    end
                end
            end
            
            % Assign to output.
            VertIDs = [mat2cell(newConn(1:nEdgeIDs,1:2),ones(nEdgeIDs,1),2)...
                mat2cell(newConn(nEdgeIDs+1:end,1:2),ones(nEdgeIDs,1),2)];
        end
        
        % Get Layers of Mesh
        function [ElemIDs,VertIDs]   = layers(CM,LayerIDs)
            %LAYERS Get layers of mesh.
            %   LayerIDS = layers(CM) returns stuff.
            %
            %   See also CHILMESH/meshLayers.
            %==============================================================
            
            %Check I/O argument(s).
            if nargin == 1                          % If no layers were
                LayerIDs	= zeros(CM.nLayers,1,'uint32');% inputted, select
                LayerIDs(:)	= 1:CM.nLayers;      	% all layers.
            end
            
            ElemIDs	= struct('OE',CM.Layers.OE(LayerIDs),...
                'IE',CM.Layers.IE(LayerIDs));
            if nargout == 2
                VertIDs	= struct('OV',CM.Layers.OV(LayerIDs),...
                    'IV',CM.Layers.IV(LayerIDs));
            end
        end
        
        % Number of Edges
        function nEdges	= numEdges(CM)
            %NUMEDGES Total number of edges in the mesh.
            %   nEdges = NUMVERTS(CM) returns the number of edges in the
            %   mesh.
            %
            %   See also CHILMESH, CHILMESH/NUMELEMS, CHILMESH/NUMVERTS.
            %==============================================================
            
            nEdges	= CM.nEdges;
        end
        
        % Number of Elements
        function nElems	= numElems(CM)
            %NUMELEMS Total number of vertices in the mesh.
            %   nElems = NUMELEMS(CM) returns the number of elements in the
            %   mesh.
            %
            %   See also CHILMESH, CHILMESH/NUMVERTS, CHILMESH/NUMEDGES,
            %   TRIANGULATION/SIZE.
            %==============================================================
            
            nElems  = CM.nElems;
%             nElems	= size(CM.ConnectivityList,1);
%             set(CM,'nElems',nElems);                % Update CM accordingly.
        end
        
        % Number of Layers
        function nLayers	= numLayers(CM)
            %NUMLAYERS Total number of layers in the mesh.
            %   nLayers = NUMLAYERS(CM) returns the number of layers in the
            %   mesh.
            %
            %   See also CHILMESH, CHILMESH/NUMVERTS, CHILMESH/NUMEDGES,
            %   CHILMESH/MESHLAYERS.
            %==============================================================
            
            nLayers	= CM.nLayers;
        end
        
        % Number of Vertices
        function nVerts	= numVerts(CM)
            %NUMVERTS Total number of vertices in the mesh.
            %   nVerts = NUMVERTS(CM) returns the number of vertices in the
            %   mesh.
            %
            %   See also CHILMESH, CHILMESH/NUMEDGES, CHILMESH/NUMELEMS.
            %==============================================================
            
            nVerts  = CM.nVerts;
%             nVerts	= size(CM.Points,1);
%             set(CM,'nVerts',nVerts);                % Update CM accordingly.
        end
        
        % Signed Area
        function Area	= signedArea(CM,varargin)
            %SIGNEDAREA Signed area of triangular or quad element(s).
            %   Area = SIGNEDAREA(CM) for mesh class CM returns the signed
            %   area of each element in the mesh. Area is a Cx1 vector,
            %   where C is the number of elements in the mesh.
            %
            %   NOTE:
            %       The signed area of a planar non-self-intersecting,
            %       convex polygon is computed from its vertices.
            %       A negative value of Area indicates a clockwise (CW)
            %       orientation of the polygon's vertices while a positive
            %       value indicates counter-clockwise (CCW) orientation.
            %
            %   Area = SIGNEDAREA(CM,ElemIDs) returns the areas of elements
            %   in the mesh corresponding to the Cx1 or 1xC vector ElemIDs.
            %
            %   Area = SIGNEDAREA(CM,SIGN) specifies whether Area is signed
            %   or unsigned, where the default value of SIGN = 'signed' and
            %   the alternative value is SIGN = 'unsigned'.
            %
            %   NOTE: If CM contains both triangular and quadrilateral
            %   elements, the connectivity list for triangular elements are
            %   expected to be in the format of an Cx4 matrix, where the
            %   first column is identical to the second column.
            %
            %   See also CHILMESH, CHILMESH/INTERIORANGLES, POLYAREA.
            %==============================================================
            
            % Check I/O argument(s).
            p   = inputParser;
            defaultElemIDs  = zeros(CM.nElems,1,'uint32');
            defaultElemIDs(:)	= 1:CM.nElems;   	% All elements.
            allowableSign	= {'signed','unsigned','sign','unsign'};
            addRequired(p,'CM',@isCHILmesh);
            addOptional(p,'Sign','signed',@(x) any(validatestring(x,allowableSign)));
            addOptional(p,'ElemIDs',defaultElemIDs,@(x) isnumeric(x));
            parse(p,CM,varargin{:})             	% Parse inputs.
            
            % Assign connectivity list and identify element types.
            [iT,iQ]	= CM.elemType(p.Results.ElemIDs);
            
            % Compute determinant (x1y2 - x2y1 + ..) and signed area (0.5*det).
            Area	= zeros(numel(p.Results.ElemIDs),1);
            iElemIDs= ismember(p.Results.ElemIDs,iT);
            Area(iElemIDs) = .5*sum([...
                CM.Points(CM.ConnectivityList(iT,1),1).*CM.Points(CM.ConnectivityList(iT,2),2) -...
                CM.Points(CM.ConnectivityList(iT,2),1).*CM.Points(CM.ConnectivityList(iT,1),2) +...
                CM.Points(CM.ConnectivityList(iT,2),1).*CM.Points(CM.ConnectivityList(iT,3),2) -...
                CM.Points(CM.ConnectivityList(iT,3),1).*CM.Points(CM.ConnectivityList(iT,2),2) +...
                CM.Points(CM.ConnectivityList(iT,3),1).*CM.Points(CM.ConnectivityList(iT,1),2) -...
                CM.Points(CM.ConnectivityList(iT,1),1).*CM.Points(CM.ConnectivityList(iT,3),2)],2);
            if ~isempty(iQ)
                iElemIDs= ismember(p.Results.ElemIDs,iQ);
                Area(iElemIDs) = .5*sum([...
                    CM.Points(CM.ConnectivityList(iQ,1),1).*CM.Points(CM.ConnectivityList(iQ,2),2) -...
                    CM.Points(CM.ConnectivityList(iQ,2),1).*CM.Points(CM.ConnectivityList(iQ,1),2) +...
                    CM.Points(CM.ConnectivityList(iQ,2),1).*CM.Points(CM.ConnectivityList(iQ,3),2) -...
                    CM.Points(CM.ConnectivityList(iQ,3),1).*CM.Points(CM.ConnectivityList(iQ,2),2) +...
                    CM.Points(CM.ConnectivityList(iQ,3),1).*CM.Points(CM.ConnectivityList(iQ,4),2) -...
                    CM.Points(CM.ConnectivityList(iQ,4),1).*CM.Points(CM.ConnectivityList(iQ,3),2) +...
                    CM.Points(CM.ConnectivityList(iQ,4),1).*CM.Points(CM.ConnectivityList(iQ,1),2) -...
                    CM.Points(CM.ConnectivityList(iQ,1),1).*CM.Points(CM.ConnectivityList(iQ,4),2)],2);
            end
            
            % Output unsigned area, if specifed.
            if strcmpi(p.Results.Sign,{'unsign','unsigned'})
                Area   = abs(Area);                   	% Unsigned area.
            end
        end
        
        % Retrieve Vertex-to-Edge Adjacency
        function EdgeIDs	= vert2Edge(CM,varargin)
            %VERT2EDGE Edges attached to point(s).
            %   EdgeIDs = VERT2EDGE(CM,VertIDs) gives all edges of the mesh
            %   that are attached to the point(s) (or vertices) specified
            %   by VertIDs. VertIDs can be an 1xP or Px1 vector of indices
            %   of logicals, where P is the number of vertice indices
            %   (size(VertIDs,1)).
            %
            %   EdgeIDs = VERT2EDGE(CM,STORE) specifies how EdgeIDs is
            %   stored as an output. The default value is STORE = 'sparse',
            %   which returns EdgeIDs as a sparse matrix that is directly
            %   accessed from the Vert2Elem adjacency list of CM, whereas
            %   STORE = 'cell' returns EdgeIDs as a Px1 cell array with
            %   each cell containing the edges attached to the vertex.
            %
            %   See also CHILMESH/BUILDADJACENCIES, CHILMESH/VERT2ELEM,
            %   CHILMESH/EDGE2VERT, TRIANGULATION/ISCONNECTED.
            %==============================================================
            
            % Check I/O arguments.
            p   = inputParser;
            defaultVertIDs	= zeros(CM.nVerts,1,'uint32');
            defaultVertIDs(:)	= 1:CM.nVerts;   	% All vertices.
            allowableStore	= {'Sparse','sparse','s','Cell','cell','c',...
                'Matrix','Mat','M','matrix','mat','m'};
            addRequired(p,'CM',@isCHILmesh);
            addParameter(p,'Store','sparse',@(x) any(validatestring(x,allowableStore)));
            addParameter(p,'VertIDs',defaultVertIDs,@isnumeric);
            parse(p,CM,varargin{:});            	% Parse inputs.
            
            % Get edges of each vertex from adjacency list.
            idx0    = p.Results.VertIDs == 0;       % Search for zeros.
            EdgeIDs = cell(numel(p.Results.VertIDs),1);
            edgeIDs = CM.Adjacencies.Vert2Edge(p.Results.VertIDs(~idx0),:);
            
            % Adjust output, if necessary.
            if any(~strcmpi(p.Results.Store,{'Sparse','sparse','s'}))
                % Identify edge neighbors of each vertex.
                edgeIDs = sort(edgeIDs.');          % Presort, work by col to improve speed.
                [~,Vert,Edge]	= find(edgeIDs); 	% Find nonzeros.
                zeroAttachmentsFound	= setdiff(1:max(Vert),Vert);
                if any(zeroAttachmentsFound)
                    idx0(zeroAttachmentsFound)  = true;
                end
                
                % Identify unique vertex element-neighbors rows.
                uVA	= unique([Vert,Edge],'rows');
                
                % Group neighbors of each vertex.
                numOfEdges	= diff([0; find(diff(uVA(:,1))); size(uVA,1)]);
                EdgeIDs(~idx0)	= mat2cell(uVA(:,2),numOfEdges);
                
                % Convert to matrix, if necessary.
                if any(strcmpi(p.Results.Store,{'Matrix','Mat','M','matrix','mat','m'}))
                    EdgeIDs	= NDcell2mat(EdgeIDs,'zero');
                end
                
            else
                EdgeIDs = edgeIDs;                  % Set equal to sparse output.
            end
        end
        
        % Retrieve Vertex-to-Element Adjacency
        function ElemIDs	= vert2Elem(CM,varargin)
            %VERT2ELEM Elements attached to point(s).
            %   ElemIDs = VERT2ELEM(CM,VertIDs) gives all elements of the
            %   mesh that are attached to the point(s) (or vertices)
            %   specified by VertIDs. VertIDs can be an 1xP or Px1 vector
            %   of indices of logicals, where P is the number of vertice
            %   indices (size(VertIDs,1)).
            %
            %   ElemIDs = VERT2ELEM(CM,STORE) specifies how ElemIDs is
            %   stored as an output. The default value is STORE = 'sparse',
            %   which returns ElemIDs as a sparse matrix that is directly
            %   accessed from the Vert2Elem adjacency list of CM, whereas
            %   STORE = 'cell' returns ElemIDs as a Px1 cell array with
            %   each cell containing the element attached to the vertex.
            %
            %   See also CHILMESH/BUILDADJACENCIES, CHILMESH/VERT2EDGE,
            %   CHILMESH/CONNECTIVITYLIST, TRIANGULATION/VERTEXATTACHMENTS
            %==============================================================
            
            % Check I/O arguments.
            p   = inputParser;
            defaultVertIDs	= zeros(CM.nVerts,1,'uint32');
            defaultVertIDs(:)	= 1:CM.nVerts;   	% All vertices.
            allowableStore	= {'Sparse','sparse','s','Cell','cell','c',...
                'Matrix','Mat','M','matrix','mat','m'};
            addRequired(p,'CM',@isCHILmesh)
            addParameter(p,'Store','sparse',@(x) any(validatestring(x,allowableStore)))
            addParameter(p,'VertIDs',defaultVertIDs,@isnumeric);
            parse(p,CM,varargin{:})             	% Parse inputs.
            
            % Get elems of each vertex from adjacency list.
            idx0    = p.Results.VertIDs == 0;       % Search for zeros.
            ElemIDs = cell(numel(p.Results.VertIDs),1);
            elemIDs = CM.Adjacencies.Vert2Elem(p.Results.VertIDs(~idx0),:);
            
            % Adjust output, if necessary.
            if any(~strcmpi(p.Results.Store,{'Sparse','sparse','s'}))
                % Identify element neighbors of each vertex.
                elemIDs = sort(elemIDs.');          % Presort, work by col to improve speed.
                [~,Vert,Elem]	= find(elemIDs); 	% Find nonzeros.
                zeroAttachmentsFound	= setdiff(1:max(Vert),Vert);
                if any(zeroAttachmentsFound)
                    idx0(zeroAttachmentsFound)  = true;
                end
                
                % Identify unique vertex element-neighbors rows.
                uVA	= unique([Vert,Elem],'rows');
                
                % Group neighbors of each vertex.
                numOfElems	= diff([0; find(diff(uVA(:,1))); size(uVA,1)]);
                ElemIDs(~idx0)	= mat2cell(uVA(:,2),numOfElems);
                
                % Convert to matrix, if necessary.
                if any(strcmpi(p.Results.Store,{'Matrix','Mat','M','matrix','mat','m'}))
                    ElemIDs	= NDcell2mat(ElemIDs,'zero');
                end
                
            else
                ElemIDs = elemIDs;                  % Set equal to sparse output.
            end
        end
        
        % Coordinates of Vertices
        function [X,Y,Z]	= vertCoordinates(CM,VertIDs)
            %VERTCOORDINATES Coordinates of vertex in mesh.
            %   [X,Y] = VERTCOORDINATES(CM,VertIDs) returns the euclidean
            %   coordinates of vertices (or points or nodes) in the mesh.
            %
            %   See also CHILMESH, CHILMESH/POINTS.
            %==============================================================
            
            if nargin == 1                          % If no vertices were
                VertIDs = zeros(CM.nVerts,1,'uint32');% inputted, select
                VertIDs(:)	= 1:CM.nVerts;          % all vertices.
            end
            idx0    = VertIDs == 0;                 % Search for zeros.
            X   = zeros(numel(VertIDs),1);          % Initialize [X,Y,Z].
            X(~idx0)	= CM.Points(VertIDs(~idx0),1); 	% Retrieve coordinates.
            if nargout >= 2
                Y   = X;
                Y(~idx0)	= CM.Points(VertIDs(~idx0),2);
                if nargout == 3
                    Z   = X;
                    Z(~idx0)	= CM.Points(VertIDs(~idx0),3);
                end
            end
        end
        
        %% ============================Mapping=============================
        % Convert to Latitude-Longitude Coordinates.
        function CM	= convert2LL(CM)
            
        end
        
        % Convert to SI Coordinates.
        function CM	= convert2SI(CM)
            
        end
        
        % Convert to USC Coordinates.
        function CM	= convert2USC(CM)
            
        end
        
        % Reverse Elevation Sign.
        function CM	= reverseZdir(CM)
            %REVERSEZDIR Flip Sign Of Mesh Points' Z-Coordinate.
            %   CM = reverseZdir(CM) for mesh class object CM returns
            %
            %   See also: CHILmesh.
            %==============================================================
            
            CM.Points(:,3)	= CM.Points(:,3).*-1;   % Reverse sign.
        end
        
        %% ===========================Plotting=============================
        % Compute Axis Attributes.
        function Ax	= axisCHILmesh(CM)
            %AXISCHILMESH Compute axis limits for mesh plotting.
            %   Ax = AXISCHILMESH(CM) returns axis properties appropriate
            %   for plotting CHILmesh mesh object CM.
            %
            %   See also CHILMESH, CHILMESH/PLOT.
            %==============================================================
            
            % If figure exists, do not change axis; otherwise, create fig.
            Ax	= gca;  hold on;                 	% Get current axis.
            if ~isempty(get(Ax,'Children'))         % Axis is open.
                xlim	= Ax.XLim;                  % Keep current axis.
                ylim	= Ax.YLim;
                offset  = 0;
            else                                    % No figure open.
                % Compute axis limits.
                iVerts  = CM.ConnectivityList > 0;
                xlim(1) = min(CM.Points(CM.ConnectivityList(iVerts),1));
                xlim(2) = max(CM.Points(CM.ConnectivityList(iVerts),1));
                ylim(1) = min(CM.Points(CM.ConnectivityList(iVerts),2));
                ylim(2)	= max(CM.Points(CM.ConnectivityList(iVerts),2));
                offset	= min(abs(xlim(1)-xlim(2)),abs(xlim(1)-ylim(2)))*.01;
            end
            
            if sum(Ax.DataAspectRatio ~= 3)         % Figure may already be adjusted.
                % Adjust figure for better aesthetic.
                axis([xlim(1)-offset xlim(2)+offset...
                    ylim(1)-offset ylim(2)+offset]);
                daspect([1 1 1]);                   % Correct data aspect.
                set(gcf,'Color',[1 1 1]);           % White background.
                box on;                             % Make axis pretty:
%                 Ax.Color    = [.5 1 .5];            % Nice green axis background.
                Ax.FontSize	= 14;                   % Big axis labels.
                Ax.FontWeight	= 'bold';           % Bold ^.
                Ax.LineWidth	= 2;                % Big axis lines.
                Ax.Units	= 'normalized';
                Ax.Position	= [0.05 0.05 0.9 0.9];  % Make axis larger.
            end
        end
        
        % Plot Mesh Elements.
        function ph	= plot(CM,varargin)
            %PLOT Plot elements of a mesh (tri, quad, or  mixed).
            %   ph = PLOT(CM) plots every 2D element in the mesh as a
            %   patched face with no color ('none') and every 1D element as
            %   a red line with each node colored red. ph is the returned
            %   plot handle.
            %
            %   ph = PLOT(CM,'ElemColor',ELEMCOLOR) colors every 2D element
            %   in the mesh with respect to the specified ELEMCOLOR. Any 1D
            %   elements are colored red. If ELEMCOLOR is red, then any 1D
            %   elements are colored blue. The edges are colored black by
            %   default.
            %
            %   ph = PLOT(CM,'EdgeColor',EDGECOLOR) plots every element in
            %   the mesh and colors every edge in respect to EDGECOLOR. The
            %   elements patched with no color by default.
            %
            %   ph = PLOT(CM,'LineWidth',LINEWIDTH) plots every element in
            %   the mesh and the edge linewidths are with respect to
            %   LINEWIDTH. The elements patched with no color and the edges
            %   are colored black by default. The default linewidth is 1.
            %
            %   ph = PLOT(CM,'LineStyle',LINESTYLE) plots every element in
            %   the mesh and the edge linestyles are with respect to
            %   LINESTYLE. The elements patched with no color and the edges
            %   are colored black by default. The default line style is '-'.
            %
            %   ph = PLOT(CM,ElemIDs) plots elements identified by an Mx1
            %   or 1xM index vector ElemIDs, where M is the number of
            %   elements for plotting.
            %
            %   See also CHILMESH, TRIPLOT, CHILMESH/EDGE, CHILMESH/PLOTELEM,
            %   CHILMESH/PLOTLAYER, CHILMESH/PLOTPOINT.
            %==============================================================
            
            % Check I/O argument(s).
            p   = inputParser;
            defaultElemIDs  = zeros(CM.nElems,1,'uint32');
            defaultElemIDs(:)	= 1:CM.nElems;   	% All elements.
            allowableColors	= {'b','g','r','c','m','y','k','w','none',...
                'blue','green','red','cyan','magenta','yellow','black','white'};
            allowableStyles	= {'-',':','-.','--','none'};
            allowableLabels	= {'Edges','Elems','Nodes'};
            addRequired(p,'CM',@isCHILmesh);
            addOptional(p,'ElemIDs',defaultElemIDs,@isnumeric);
            addParameter(p,'ElemColor','none',@(x) (isnumeric(x) &...
                numel(x) == 3 & min(size(x)) == 1) | any(strcmpi(x,allowableColors)));
            addParameter(p,'EdgeColor','k',@(x) any(validatestring(x,allowableColors)));
            addParameter(p,'Label','-',@(x) any(validatestring(x,allowableLabels)));
            addParameter(p,'LineWidth',1,@(x) x > 0);
            addParameter(p,'LineStyle','-',@(x) any(validatestring(x,allowableStyles)));
            parse(p,CM,varargin{:});             	% Parse inputs.
            
            % Plot 2D elements: if no ElemColor specified, plot as edges.
            CM.axisCHILmesh;                        % Compute offset.
            if any(strcmp(p.Results.ElemColor,{'none','w'})) ||...
                    sum(p.Results.ElemColor) == 0   % Plot edges.
                % Retrieve edges of each ElemIDs; plot edges.
                EdgeIDs	= CM.elem2Edge(p.Results.ElemIDs);
                ph	= CM.plotEdge(EdgeIDs(:),...
                    'LineStyle',p.Results.LineStyle,...
                    'LineWidth',p.Results.LineWidth,...
                    'Color',p.Results.EdgeColor); %%%%% Maybe asign to CM.Plot?
                
            else                                    % Plot faces.
                % Remove duplicate and empty (0 or NaN) elems.
                ElemIDs	= unique(p.Results.ElemIDs);  	% Assign ElemIDs.
                ElemIDs(ElemIDs == 0 | isnan(ElemIDs))	= [];
                
                % Identify 2D triangular elements.
                [iT,iQ]	= CM.elemType(ElemIDs);
                
                % Create patch for mesh data structure and plotting inputs.
                patchobj.Faces	= CM.ConnectivityList(iT,[1,1:3]);
                if ~isempty(iQ)                         % Elems == faces.
                    patchobj.Faces = [patchobj.Faces; CM.ConnectivityList(iQ,:)];
                end
                patchobj.Vertices	= CM.Points(:,[1 2]);
                patchobj.FaceColor	= p.Results.ElemColor;
                patchobj.EdgeColor	= p.Results.EdgeColor;
                patchobj.LineWidth	= p.Results.LineWidth;
                patchobj.LineStyle	= p.Results.LineStyle;
                ph  = patch(patchobj);
%                 ph  = CM.plotElem(p.Results.ElemIDs,...
%                     'color',p.Results.ElemColor,'EdgeColor',p.Results.EdgeColor,...
%                     'LineWidth',p.Results.LineWidth,'LineStyle',p.Results.LineStyle);
            end
            
            % Plot any 1D elements.
            if exist('CM.BoundaryCondition.id','var')
                seg	= find(ismember([CM.BoundaryCondition.id],[18 19]));
                for idx = 1:length(seg)             % Plot all 1D line segments.
                    edges
                    
                end
                if strcmpi('r',p.Results.ElemColor) % If elem color is red,
                    patchobj.EdgeColor	= 'b';	% change 1Delem color to blue.
                else
                    patchobj.EdgeColor	= 'r';	% Default 1Delem color.
                end
                
                % Insert NaNs where discontinuity is needed.
                x	= NaN(sum(idx1D),3);	y	= x;% X and Y coordinates.
                x(:,1:2)	= CM.ConnectivityList(p.Results.ElemIDs(idx1D),1:2);
                x(:,1)  = CM.Points(x(:,1),1);	x(:,2)  = CM.Points(x(:,2),1);
                y(:,1:2)	= CM.ConnectivityList(p.Results.ElemIDs(idx1D),1:2);
                y(:,1)  = CM.Points(y(:,1),2);	y(:,2)  = CM.Points(y(:,2),2);
                
                % Plot lines for 1D elements.
                ph	= [ph; {plot(x',y','-','LineWidth',p.Results.LineWidth,'Color',...
                    patchobj.EdgeColor,'Marker','.','MarkerSize',p.Results.LineWidth*5),'1D Element Lines'}];
            end
            %
            % Plot any channel segments.
            %                 loc	= find([CM.BoundaryConditions.id] == 18);
            %                 for k	= 1:length(loc)
            %                     nodes   = CM.BoundaryConditions(loc(k)).nodes;
            %                     x	= CM.Points(nodes,1);	y	= CM.Points(nodes,2);
            %                     if all(all(CM.BoundaryConditions(loc(k)).att(:,[1 2 3]) == 0))
            %                         plot(x,y,'color',a.ColorOrder(2,:),'linewidth',2);
            %                     else
            %                         plot(x,y,'color',a.ColorOrder(1,:),'linewidth',2);
            %                     end
            %                     plot(x([1,end]),y([1,end]),'ko','markersize',3,...
            %                         'markerfacecolor',a.ColorOrder(3,:));
            %                 end
            %             end
            
            % Correct axis lims.
            CM.axisCHILmesh;    hold off;
        end
        
        % Plot Mesh Edges.
        function ph	= plotEdge(CM,varargin)
            %PLOTEDGE Plot edges of the mesh.
            %   ph = PLOTEDGE(CM) for mesh class CM plots discontinuous
            %   line segments for each unique edge in the mesh, where an
            %   edge is defined as a 1x2 vector of vertices with respect to
            %   the points list of CM. ph is a line object.
            %
            %   ph = PLOTEDGE(CM,EdgeIDs) plots each edge specified by
            %   EdgeIDs, where EdgeIDs is a 1xE or Ex1 vector of indices or
            %   logicals with respect to CM.Points and E is the number of
            %   edge indices for plotting that is equal to size(EdgeIDs,1).
            %
            %   See also CHILMESH, CHILMESH/PLOT, CHILMESH/PLOTELEM,
            %   CHILMESH/PLOTPOINT, CHILMESH/PLOTLAYER.
            %==============================================================
            
            % Check I/O argument(s).
            p   = inputParser;  p.FunctionName  = 'plotEdge';
            defaultEdgeIDs  = zeros(CM.nEdges,1,'uint32');
            defaultEdgeIDs(:)	= 1:CM.nEdges;   	% All edges.
            allowableColors	= {'b','g','r','c','m','y','k','w','none'};
            allowableStyles	= {'-',':','-.','--','none'};
            addRequired(p,'CM',@isCHILmesh)
            addOptional(p,'EdgeIDs',defaultEdgeIDs,@(x) isnumeric(x) && min(size(x)) == 1)
            addParameter(p,'Color','k',@(x) any(validatestring(x,allowableColors)))
            addParameter(p,'LineWidth',1,@(x) x > 0)
            addParameter(p,'LineStyle','-',@(x) any(validatestring(x,allowableStyles)))
            parse(p,CM,varargin{:})             	% Parse inputs.
            
            % Remove duplicate and empty (0 or NaN) edges.
            EdgeIDs	= unique(p.Results.EdgeIDs);	% Assign EdgeIDs.
            EdgeIDs(EdgeIDs == 0 | isnan(EdgeIDs))	= [];
            
            % Index EdgeIDs into points list to retrive vertices of edge.
            VertIDs	= CM.edge2Vert(EdgeIDs)';
            
            % Create x, y coordinate vectors with NaNs.
            CM.axisCHILmesh;                        % Compute offset, limits.
            emptyCoords	= NaN(1,size(VertIDs,2));
            x   = [CM.Points(VertIDs(1,:),1)';...   % X-Coordinates.
                CM.Points(VertIDs(2,:),1)'; emptyCoords];
            y   = [CM.Points(VertIDs(1,:),2)';...   % Y-Coordinates.
                CM.Points(VertIDs(2,:),2)'; emptyCoords];
            ph  = plot(x(:),y(:),...                % Plot edges.
                'LineStyle',p.Results.LineStyle,...
                'LineWidth',p.Results.LineWidth,...
                'Color',p.Results.Color);
            CM.axisCHILmesh;   hold off;            % Daspect, etc.
        end
        
        % Plot Mesh Elems.
        function ph	= plotElem(CM,varargin)
            %PLOTELEM Plot elems of the mesh.
            %   ph = PLOTELEM(CM) for mesh class CM plots patches for
            %   each unique elem in the mesh.
            %
            %   ph = PLOTELEM(CM,ElemIDs) plots each elem specified by
            %   ElemIDs, where ElemIDs is a 1xE or Ex1 vector of indices or
            %   logicals with respect to CM.ConnectivityList and E is the
            %   number of elem indices for plotting that is equal to
            %   size(ElemIDs,1).
            %
            %   See also CHILMESH, CHILMESH/PLOT, CHILMESH/PLOTEDGE.
            %   CHILMESH/PLOTPOINT, CHILMESH/PLOTLAYER.
            %==============================================================
            
            % Check I/O argument(s).
            p   = inputParser;  p.FunctionName  = 'plotElem';
            defaultElemIDs  = zeros(CM.nElems,1,'uint32');
            defaultElemIDs(:)	= 1:CM.nElems;   	% All elems.
            allowableColors	= {'b','g','r','c','m','y','k','w','none'};
            allowableStyles	= {'-',':','-.','--','none'};
            addRequired(p,'CM',@isCHILmesh)
            addOptional(p,'ElemIDs',defaultElemIDs,@(x) isnumeric(x) && min(size(x)) == 1)
            addParameter(p,'Color','b',@(x) (isnumeric(x) &...
                numel(x) == 3 & min(size(x)) == 1) | any(strcmpi(x,allowableColors)));
            addParameter(p,'EdgeColor','k',@(x) any(validatestring(x,allowableColors)));
            addParameter(p,'LineWidth',1,@(x) x > 0)
            addParameter(p,'LineStyle','-',@(x) any(validatestring(x,allowableStyles)))
            parse(p,CM,varargin{:})             	% Parse inputs.
            
            % Remove duplicate and empty (0 or NaN) elems.
            ElemIDs	= unique(p.Results.ElemIDs);  	% Assign ElemIDs.
            ElemIDs(ElemIDs == 0 | isnan(ElemIDs))	= [];
            
            % Identify 2D triangular elements.
            [iT,iQ]	= CM.elemType(ElemIDs);
            
            % Create patch for mesh data structure and plotting inputs.
            patchobj.Faces	= CM.ConnectivityList(iT,[1,1:3]);
            if ~isempty(iQ)                         % Elems == faces.
                patchobj.Faces = [patchobj.Faces; CM.ConnectivityList(iQ,:)];
            end
            patchobj.Vertices	= CM.Points(:,[1 2]);
            patchobj.FaceColor	= p.Results.Color;
            patchobj.EdgeColor	= p.Results.EdgeColor;
            patchobj.LineWidth	= p.Results.LineWidth;
            patchobj.LineStyle	= p.Results.LineStyle;
            ph  = patch(patchobj);
        end
        
        % Plot Mesh Facets Labels.
        function ph	= plotLabel(CM,varargin)
            %PLOTLABEL Label entities of the mesh.
            %   ph = PLOTLABEL(CM) plots
            %
            %   See also CHILMESH, CHILMESH/PLOT, CHILMESH/EDGE,
            %   CHILMESH/PLOTPOINT.
            %==============================================================
            
            % Check I/O argument(s).
            p   = inputParser;
            allowableLabels	= {'Edge','edge','All','all',...
                'Element','element','Elem','elem','Face','face',...
                'Node','node','Vertex','vertex','Point','point'};
            addRequired(p,'CM',@isCHILmesh);
            addOptional(p,'Label','All',@(x) any(validatestring(x,allowableLabels)));
            addOptional(p,'IDs','All',@(x) isnumeric(x));
            parse(p,CM,varargin{:});             	% Parse inputs.
            
            % Plot facet labels.
            gcf;    CM.axisCHILmesh;
            if any(strcmp(p.Results.Label,{'Node','node','Vertex','vertex','Point','point','All','all'}))
                for ID	= 1:length(p.Results.IDs)   % Plot Vert(s).
                    VertID = p.Results.IDs(ID);     % Entity ID.
                    text(CM.Points(VertID,1),CM.Points(VertID,2),...
                        ['Node ',num2str(VertID)],...
                        'Color','red','HorizontalAlignment','center')
                end
            end
            if any(strcmp(p.Results.Label,{'Edge','edge','All','all'}))
                for ID	= 1:length(p.Results.IDs)   % Plot Edge(s).
                    EdgeID = p.Results.IDs(ID);     % Entity ID.
                    VertIDs = CM.Adjacencies.Edge2Vert(EdgeID,:);
                    text(mean(CM.Points(VertIDs,1)),mean(CM.Points(VertIDs,2)),...
                        ['Edge ',num2str(EdgeID)],...
                        'Color','green','HorizontalAlignment','center')
                end
            end
            if any(strcmp(p.Results.Label,{'Element','element','Elem','elem','Face','face','All','all'}))
                for ID	= 1:length(p.Results.IDs)   % Plot Elem(s).
                    ElemID = p.Results.IDs(ID);     % Entity ID.
                    VertIDs = CM.Adjacencies.Elem2Vert(ElemID,:);
                    text(mean(CM.Points(VertIDs,1)),mean(CM.Points(VertIDs,2)),...
                        ['Elem ',num2str(ElemID)],...
                        'Color','blue','HorizontalAlignment','center')
                end
            end
            CM.axisCHILmesh;
            ph  = gca;  ph	= ph.Children;
        end
        
        % Plot Mesh Layers.
        function ph = plotLayer(CM,varargin)
            %PLOTLAYER Plot contour map of the mesh layers.
            %   pL = PLOTLAYER(CM) for mesh class CM plots a contour mesh
            %   representing all of the elements in the mesh defined by CM
            %   as layers. The default colormap used is parula and overlays
            %   any currently open plot.
            %
            %   PLOTLAYER(CM,cmap) adjusts the colormap used according to
            %   the char input cmap.
            %
            %   PLOTLAYER(CM,iL) plots only the specified layers iL of CM.
            %
            %   See also CHILMESH, CHILMESH/PLOT, CHILMESH/PLOTELEM,
            %   MESHLAYERS, CONTOUR.
            %==============================================================
            
            % Check I/O argument(s).
            p   = inputParser;
            defaultLayers	= zeros(CM.nLayers,1,'uint32'); % All Layers.
            defaultLayers(:)	= 1:CM.nLayers;     % Layers for plotting.
            allowableCmaps	= {'parula','jet','hsv','hot','cool','spring',...
                'summer','autumn','winter','gray','bone','copper','pink',...
                'lines','colorcube','prism','flag','white'};
            addRequired(p,'CM',@isCHILmesh);
            addOptional(p,'Layers',defaultLayers,@isnumeric);
            addOptional(p,'Cmap','parula',@(x) ischar(x) | (isnumeric(x) & min(size(x)) == 1 & max(size(x)) == 3));
            parse(p,CM,varargin{:});                % Parse inputs.
            if any(strcmpi(p.Results.Cmap,allowableCmaps))
                handleCmap	= p.Results.Cmap;     	% Assign colormap.
            else
                handleCmap  = 'parula';             % Default to parula.
            end
            
            % No zero-indices in connectivity list.
            if size(CM.ConnectivityList,2) == 4
                iT  = CM.elemType;
                CM.ConnectivityList(iT,4)	= CM.ConnectivityList(iT,1);
            end
            
            % Create patch for mesh data structure and plotting inputs.
            CM.axisCHILmesh;
            patchobj.Vertices  = CM.Points(:,[1 2]);
            patchobj.EdgeColor = 'k';
            patchobj.LineWidth = .1;
            patchobj.LineStyle = '-';
            
            % Plot contour map.
            nL	= length(p.Results.Layers);      	% # of layers to plot.
            ph  = cell(nL,1);
            handleCmap	= str2func(handleCmap);
            Cmap	= colormap(handleCmap(nL));
            for idx	= 1:nL
                patchobj.Faces	= CM.ConnectivityList(...
                    [CM.Layers.OE{p.Results.Layers(idx)};...
                    CM.Layers.IE{p.Results.Layers(idx)}],:);
                patchobj.FaceColor  = Cmap(idx,:);
                ph{idx}	= patch(patchobj);
            end
            CM.axisCHILmesh;
            
            % Edit colorbar - indexed to layers.
            cB	= colorbar;	cB.TickLabels	= num2cell((1:nL)');
            cB.Ticks	= (1/(2*nL)):(1/nL):(1-(1/(2*nL)));
            cB.Label.String	= 'Mesh Layers';
            cB.Label.FontSize	= 15;   cB.Label.FontWeight	= 'bold';
        end
        
        % Plot Mesh Points
        function ph	= plotPoint(CM,varargin)
            %PLOTPOINT Plot points in the mesh.
            %   ph = PLOTPOINT(CM) for mesh class CM plots a marker at each
            %   unique vertex in the mesh, where a vertex is defined as a
            %   1x3 coordinate vector. The default marker is 'ro' with a
            %   shaded face and ph is the returns a line object.
            %
            %   ph = PLOTPOINT(CM,TYPE) specifies the point to be plotted,
            %   where the default value of TYPE is 'Vertex'. Other values
            %   of TYPE include: 'Centroid' for centroid of an element,
            %   'Midpoint', for the midpoint of edges.
            %
            %   ph = PLOTPOINT(CM,POINT) specifies entities in the mesh for
            %   plotting and corresponding to TYPE, i.e. VertIDs, ElemIDs,
            %   and EdgeIDs. The default parameter plots mesh vertices.
            %
            %   See also CHILMESH, CHILMESH/PLOT, CHILMESH/PLOTEDGE.
            %   CHILMESH/PLOTELEM, CHILMESH/PLOTLAYER.
            %==============================================================
            
            % Check I/O argument(s).
            p   = inputParser;  p.FunctionName  = 'plotPoint';
            allowableTypes  = {'Centroid','centroid','Cent','cent','Centroids','centroids','Cents','cents'...
                'Element','element','Elements','elements','Elem','elem','Elems','elems',...
                'Midpoint','midpoint','Midpoints','midpoints','Midpt','midpt',...
                'Edge','edge','Edges','edges',...
                'Vertex','vertex','Vertices','vertices','Vert','vert','Verts','verts',...
                'Point','point','Points','points','PT','pt','PTs','pts','PTS'};
            allowableColors	= {'b','g','r','c','m','y','k','w','none'};
            allowableMarkers	= {'.','o','x','+','*','s','d','v','^','<','>','p','h'};
            addRequired(p,'CM',@isCHILmesh);
            addOptional(p,'IDs',[],@(x) isnumeric(x));
            addOptional(p,'PointType','Vertex',@(x) any(validatestring(x,allowableTypes)));
            addParameter(p,'Color','r',@(x) any(validatestring(x,allowableColors)));
            addParameter(p,'Marker','o',@(x) any(validatestring(x,allowableMarkers)));
            addParameter(p,'MarkerSize',5,@(x) x > 0);
            parse(p,CM,varargin{:})             	% Parse inputs.
            
            % Assign x,y,z coordinate(s) of plotted point(s).
            switch p.Results.PointType
                case {'Element','element','Elements','elements','Elem',...
                        'elem','Elems','elems','Centroid','centroid',...
                        'Cent','cent','Centroids','centroids','Cents','cents'}
                    % Assign ElemIDs.
                    if isempty(p.Results.IDs)    	% Plot for all elems.
                        ElemIDs	= zeros(CM.nElems,1,'uint32');
                        ElemIDs(:)	= 1:CM.nElems;
                    else                            % Extract specific elems.
                        ElemIDs = p.Results.IDs;
                    end
                    
                    % Remove duplicates and zeros; compute centroid.
                    ElemIDs	= unique(ElemIDs,'rows');
                    ElemIDs(ElemIDs == 0 | isnan(ElemIDs))	= [];
                    if strcmpi(p.Results.PointType,{'Element','element',...
                            'Elements','elements','Elem','elem','Elems','elems'})
                        [X,Y]   = CM.elemCoordinates(ElemIDs);
                        
                    else                            % Plot Centroid.
                        [X,Y]	= CM.centroid(ElemIDs);
                    end
                    
                case {'Edge','edge','Edges','edges','Midpoint','midpoint',...
                        'Midpoints','midpoints','Midpt','midpt'}
                    % Assign EdgeIDs.
                    if isempty(p.Results.IDs)    	% Plot for all edges.
                        EdgeIDs	= zeros(CM.nEdges,1,'uint32');
                        EdgeIDs(:)	= 1:CM.nEdges;
                    else                            % Extract specific edges.
                        EdgeIDs = p.Results.IDs;
                    end
                    
                    % Remove duplicates and zeros; compute centroid.
                    EdgeIDs	= unique(EdgeIDs,'rows');
                    EdgeIDs(EdgeIDs == 0 | isnan(EdgeIDs))	= [];
                    if strcmpi(p.Results.PointType,{'Edge','edge','Edges','edges'})
                        [X,Y]   = CM.edgeCoordinates(EdgeIDs);
                        
                    else                            % Plot Midpoint.
                        [X,Y]	= CM.edgeMidpoint(EdgeIDs);
                    end
                    
                case {'Vertex','vertex','Vertices','vertices','Vert','vert','Verts','verts',...
                        'Point','point','Points','points','PT','pt','PTs','pts','PTS'}
                    % Assign VertIDs.
                    if isempty(p.Results.IDs)    	% Plot for all vertices.
                        VertIDs	= zeros(CM.nVerts,1,'uint32');
                        VertIDs(:)	= 1:CM.nVerts;
                    else                            % Extract specific verts.
                        VertIDs = p.Results.IDs;
                    end
                    
                    % Remove duplicates and zeros; compute centroid.
                    VertIDs	= unique(VertIDs,'rows');
                    VertIDs(VertIDs == 0 | isnan(VertIDs))	= [];
                    [X,Y]	= CM.vertCoordinates(VertIDs);
                    
                otherwise                           % Error will occur by now.
            end
            
            % Plot points within mesh.
            CM.axisCHILmesh;                        % Compute offset, limits.
            ph  = plot(X,Y,...                      % Plot points in mesh.
                'Color',p.Results.Color,'LineStyle','none',...
                'Marker',p.Results.Marker,'MarkerEdgeColor',p.Results.Color,...
                'MarkerFaceColor',p.Results.Color,'MarkerSize',p.Results.MarkerSize);
            CM.axisCHILmesh;                        % Daspect, etc.
        end
        
        % Plot Mesh Quality.
        function ph = plotQuality(CM,varargin)
            %PLOTQUALITY Plot contour map of the mesh layers.
            %   pL = PLOTQUALITY(CM) for mesh class CM plots a contour mesh
            %   representing the quality of all elements in the mesh
            %   defined by CM. The default colormap used is cool and
            %   overlays any currently open plot.
            %
            %   PLOTQUALITY(CM,cmap) adjusts the colormap used according to
            %   the char input cmap.
            %
            %   PLOTQUALITY(CM,ElemIDs) plots only the elements specified
            %   by ElemIDs of CM.
            %
            %   See also CHILMESH, CHILMESH/PLOT, CHILMESH/PLOTELEM,
            %   ELEMQUALITY, CONTOUR.
            %==============================================================
            
            % Check I/O argument(s).
            p   = inputParser;
            defaultElemIDs	= zeros(CM.nElems,1,'uint32');
            defaultElemIDs(:)	= 1:CM.nElems;      % Elems for plotting.
            allowableCmaps	= {'parula','jet','hsv','hot','cool','spring',...
                'summer','autumn','winter','gray','bone','copper','pink',...
                'lines','colorcube','prism','flag','white'};
            addRequired(p,'CM',@isCHILmesh);
            addOptional(p,'ElemIDs',defaultElemIDs,@isnumeric);
            addOptional(p,'Cmap','cool',@(x) ischar(x) | (isnumeric(x) & min(size(x)) == 1 & max(size(x)) == 3));
            parse(p,CM,varargin{:});                % Parse inputs.
            if any(strcmpi(p.Results.Cmap,allowableCmaps))
                handleCmap0	= p.Results.Cmap;     	% Assign colormap.
            else
                handleCmap0  = 'cool';               % Default to parula.
            end
            handleCmap	= str2func(handleCmap0);
            
            % No zero-indices in connectivity list.
            if size(CM.ConnectivityList,2) == 4
                iT  = CM.elemType;
                CM.ConnectivityList(iT,4)	= CM.ConnectivityList(iT,1);
            end
            
            % Compute and bin element qualities; plot contour map.
            q   = CM.elemQuality(p.Results.ElemIDs);
            binEdges= 0:.05:1;
            nBins   = length(binEdges)-1;
            ph  = cell(nBins,1);
            Cmap	= colormap(handleCmap(nBins));
            if strcmp(handleCmap0,'cool')
                Cmap	= flipud(Cmap);
            end
            CM.axisCHILmesh;
            for idx = 1:nBins
                iElemIDs = p.Results.ElemIDs(q >= binEdges(idx) & q < binEdges(idx+1));
                if isempty(iElemIDs)
                    continue
                end
                ph{idx}	= CM.plotElem(iElemIDs,'color',Cmap(idx,:));
                drawnow;
            end
            CM.axisCHILmesh;
            
            % Edit colorbar - indexed to layers.
            cB	= colorbar('Ticks',binEdges,...
                'TickLabels',flipud(num2cell(round(binEdges,3)')));
            cB.Label.String	= 'Mesh Quality';
            cB.Label.FontSize	= 15;   cB.Label.FontWeight	= 'bold';
        end
        
    end
    
    %     events
    %
    %     end
    %
    %     enumeration
    %
    %     end
end
