function M = NDcell2mat(C,FILLER)
    %NDCELL2MAT Cell2mat for an array-cell structure of any length.
    %   M = MYCELL2MAT(C) converts a unidimensional cell array CELL with
    %   contents of the same data type into a single matrix M. There is no
    %   need for the dimensions of the cell's contents to match, however C
    %   must be of dimension Nx1 or 1xN and all cell contents are expected
    %   to be vectors.
    %
    %   The dimensionality of M will be NxL or LxN, depending on the
    %   dimensionality of C, where L corresponds to the maximum length of
    %   all cell contents of C. Given that not all contents of C may be the
    %   same length, any "extra" space will be replaced with a NaN in M.
    %
    %   MYCELL2MAT is not supported for cell arrays containing objects.
    %
    %   Example:
    %   C = {[0]; [0 1]; [0 1 2]; [0 1 2 3]; [4 0]};
    %   M = MYCELL2MAT(C)
    %   %returns
    %   M = [0 NaN NaN NaN; 0 1 NaN NaN; 0 1 2 NaN; 0 1 2 3; 4 0 NaN NaN]
    %
    %   See also cell2mat.
    %======================================================================
    
    % Check I/O argument(s).
    defaultFILLER    = 'nan';
    if nargin == 2
        nanStrings  = {'nan','nans','n','ns'};
        zeroStrings = {'zero','zeros','z','zs'};
        acceptableInput2    = [nanStrings;zeroStrings];
        checkFILLER = find(sum(strcmpi(acceptableInput2,FILLER),2) == 1);
        if any(checkFILLER)
            defaultFILLER	= acceptableInput2{checkFILLER,1};
        end
    end
    
    % Determine largest length of cell contents.
    sC  = cell2mat(cellfun(@(x) size(x),C,'UniformOutput',0));
    
    % Check input.
    if sum(size(C) == 1) == 0 || all(sC(:,1) == sC(1,1))
        if sC(1,1) > sC(1,2)                    % Make sure they don't concactenate.
            C = cellfun(@transpose,C,'UniformOutput',0);
        end
        M	= cell2mat(C);                   	% Built-in cell2mat.
        return                                 
    elseif size(C,1) == 1
        C   = transpose(C);                     % Accomodate for 1xN input.
        transpose_M	= true;                     % Output will be LxN.
    else
        transpose_M	= false;                    % Output will be NxN.
    end
    
    % Step through cell contents via their lengths; extract entries to M.
    Lmax= max(sC(:,1));                        	% Cell contents' max length.
    if strcmp('nan',defaultFILLER)              % Fill with NaNs.
        M	= NaN(length(C),Lmax);
    else                                        % Fill with zeros.
        M	= zeros(length(C),Lmax);
    end
    for L	= 1:Lmax
        % Insert a matrix of the cells of length L into indexed rows of M.
        iM	= sC(:,1) == L;                   	% Index cells of length L.
        try
            M(iM,1:L)	= cell2mat(C(iM));
        catch   % Input cells might be row vectors.
            M(iM,1:L)	= reshape(cell2mat(C(iM)),L,sum(iM))';
        end
    end
    
    % Transpose output M to reflect dimensionality of C, if needed.
    if transpose_M
        M   = transpose(M);
    end
end