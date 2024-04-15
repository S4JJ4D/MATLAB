classdef signalMask < handle & matlab.mixin.CustomDisplay &  matlab.mixin.SetGet 
%signalMask Modify and convert signal masks and extract signal regions of interest 
%   M = signalMask(SOURCE) creates a signal mask object, M. Input SOURCE
%   contains a source mask in the form of a table, a categorical vector
%   sequence, or a matrix of binary sequences. The input SOURCE mask
%   defines the locations of regions of interest of a signal together with
%   the label or category values for each region.
%
%   The output object, M, allows you to represent the input signal mask,
%   SOURCE, in any of its three possible forms: table, categorical
%   sequence, or matrix of binary sequences. Based on the properties of M,
%   you can modify regions of interest of the SOURCE mask by extending or
%   shortening their duration, merging same-category regions that are
%   sufficiently close, or removing regions that are not long enough.
%   Further, you can use the signalMask object to extract signal regions of
%   interest from a signal vector, and to plot a signal with color-coded
%   regions of interest.
%   
%    - When SOURCE is a table, it must contain two variables. The first
%      variable is a two-column matrix. Each row of the matrix contains the
%      beginning and end sample indices of a signal region of interest. If
%      the matrix does not contain integer indices, signalMask rounds them
%      to the nearest integer. The second variable contains the region
%      labels specified as a categorical array or a string array.
%
%    - When SOURCE is a categorical sequence, groups of contiguous
%      same-value category elements indicate a signal region of interest
%      labeled with that particular category. Elements in the signal mask
%      that belong to no category (and hence have no label value) should be
%      specified as the missing categorical, displayed as <undefined>.
%
%    - When SOURCE is a matrix of binary sequences with P columns, each
%      column is interpreted as a signal mask with true elements marking
%      regions of interest for each of P different categories. In this
%      scenario, you can specify a list of P category names corresponding
%      to each column using the 'Categories' property.
%
%   M = signalMask(...,'SampleRate',Fs) specifies a numeric, positive
%   scalar sample rate value, Fs. In this case, when input SOURCE is a
%   table, signalMask assumes that the table contains region limits in
%   seconds. When performing indexing operations, signalMask converts time
%   values to sample indices rounding to the nearest integer. When
%   'SampleRate' is omitted, all region limits are treated as sample
%   indices.
%
%   M = signalMask(...,'Categories',CAT) specifies a P-element string
%   array, CAT, with category names. This property can be set only when
%   SOURCE is a matrix of binary sequences with P columns. signalMask
%   interprets the i-th column in SOURCE as a signal mask corresponding to
%   the i-th category in array CAT. signalMask sets categories to ["1" "2"
%   ... "P"] when 'Categories' is omitted and input SOURCE is a matrix of
%   binary sequences. For any other input SOURCE type, categories are
%   inferred directly from the SOURCE values.
%
%   M = signalMask(...,'LeftExtension',N) modifies mask by extending
%   regions to the left by N samples, where N is a nonnegative integer. The
%   number of extended samples is truncated when the beginning of the
%   sequence is reached. When not specified, 'LeftExtension' defaults to 0
%   samples.
%
%   M = signalMask(...,'RightExtension',N) modifies mask by extending
%   regions to the right by N samples, where N is a nonnegative integer.
%   When not specified, 'RightExtension' defaults to 0 samples.
%
%   M = signalMask(...,'LeftShortening',N) modifies mask by shortening
%   regions from the left by N samples, where N is a nonnegative integer.
%   Regions are removed if shortened by a number of samples equal to or
%   larger than their length. When not specified, 'LeftShortening' defaults
%   to 0 samples.
%
%   M = signalMask(...,'RightShortening',N) modifies mask by shortening
%   regions from the right by N samples, where N is a nonnegative integer.
%   Regions are removed if shortened by a number of samples equal to or
%   larger than their length. When not specified, 'RightShortening'
%   defaults to 0 samples.
%
%   M = signalMask(...,'MergeDistance',N) modifies mask by merging regions
%   of the same-category value that are separated by N samples or less,
%   where N is a nonnegative integer. When not specified, 'MergeDistance'
%   defaults to 0 samples.
%
%   M = signalMask(...,'MinLength',N) modifies mask by removing regions
%   that are shorter than N samples, where N is a positive integer. When
%   not specified, 'MinLength' defaults to 1 sample.
%
%   The region limits of the input SOURCE, for each category, are converted
%   to a matrix of region limits. Modifications based on the
%   'LeftExtension', 'RightExtension', 'LeftShortening', 'RightShortening',
%   'MergeDistance', and 'MinLength' properties are then applied to the
%   resulting region limits matrices in this order:
%   
%      1) Extend regions to the left or right based on 'LeftExtension' and 
%         'RightExtension'.
%      2) Shorten regions from the left or right based on 'LeftShortening' 
%         and 'RightShortening'.
%      3) Merge close-enough regions based on 'MergeDistance'. signalMask
%         always merges same-category regions that overlap, are separated
%         by zero samples, or are repeated.
%      4) Remove short regions based on 'MinLength'.
%
%   signalMask properties:
%
%   SourceType               - Type of input source mask. Can be one of
%                              'roiTable', 'categoricalSequence', or
%                              'binarySequences'. (read-only)
%
%   SampleRate                - Sample rate value. Applies only when 
%                               'SampleRate' is specified when calling
%                               signalMask. (read-only)
%
%   Categories                - List of categories in the SOURCE mask. It
%                               is read-only unless 'SourceType' is
%                               'binarySequences', in which case it can be
%                               set to a string vector with P elements, one
%                               category name for each of the P columns in
%                               the matrix of binary sequences source mask.
%
%   SpecifySelectedCategories - True if you want to select a subset of
%                               categories from the 'Categories' list. If
%                               set to false, all categories in 'Categories'
%                               are selected. Default is false.
%
%   SelectedCategories        - Vector with index values pointing to the 
%                               category elements in the 'Categories'
%                               property that you want selected. Categories
%                               not listed in this property are filtered
%                               out from the input SOURCE mask when calling
%                               the different signalMask object functions.
%                               Category indices must be sorted in
%                               ascending order. This property applies
%                               only when 'SpecifySelectedCategories' is
%                               true.
%
%   LeftExtension             - Number of samples to extend regions to the
%                               left. Default is 0.
%
%   RightExtension            - Number of samples to extend regions to the
%                               right. Default is 0.
%
%   LeftShortening            - Number of samples to shorten regions from 
%                               the left. Default is 0.
%
%   RightShortening           - Number of samples to shorten regions from
%                               the right. Default is 0.
%
%   MergeDistance             - Merge distance in samples. signalMask merges
%                               same-category contiguous regions separated
%                               by MergeDistance samples or less. Default
%                               is 0.
%
%   MinLength                 - Minimum region length in samples. signalMask
%                               removes regions shorter than MinLength.
%                               Default is 1.
%
%   signalMask methods:
%   
%   roimask       - Get ROI table mask
%   catmask       - Get categorical sequence mask
%   binmask       - Get matrix of binary sequences mask
%   extractsigroi - Extract regions of interest based on signal mask
%   plotsigroi    - Plot signal regions based on signal mask
%
%   See also binmask2sigroi, sigroi2binmask, extractsigroi, extendsigroi,
%   shortensigroi, mergesigroi, removesigroi.

%   Copyright 2020 MathWorks, Inc.
    
properties (SetAccess = private)
    %SourceType Type of input source mask. Can be one of
    %'roiTable', 'categoricalSequence', 'binarySequences' (read-only)
    SourceType
end
properties
    %SpecifySelectedCategories True if you want to specify selected
    %categories. False if you want all categories selected
    SpecifySelectedCategories
    %LeftExtension Number of samples to extend regions to the left
    LeftExtension
    %RightExtension Number of samples to extend regions to the right
    RightExtension
    %LeftShortening Number of samples to shorten regions from the left
    LeftShortening
    %RightShortening Number of samples to shorten regions from the right
    RightShortening
    %MergeDistance Merge distance in samples. signalMask merges
    %same-category contiguous regions separated by MergeDistance samples or
    %less
    MergeDistance
    %MinLength Minimum region length in samples. signalMask removes regions
    %shorter than MinLength
    MinLength
end

properties (Dependent)
    %SampleRate Sample rate value. Applies only when 'SampleRate' is
    %specified when calling signalMask (read-only)
    SampleRate
    %Categories List of region of interest categories
    Categories
    %SelectedCategories Indices of selected categories
    SelectedCategories
end

properties (Access = private)
    %pSampleRateSpecified True if sample rate is specified
    pSampleRateSpecified
    % Store dependent property values
    pSampleRate
    pCategories
    pSelectedCategories
    % Number of categories for input source
    pNumCategories
    
    %Source mask - original input mask without any modifications
    pSourceMask
    
    % True if input source mask has inherent categories, in which case
    % users cannot set them
    pIsInherentCategories;
            
    % Property list used to display the object
    pPropertyList = [...
        "SourceType",...        
        "SampleRate",...
        "Categories",...
        "SpecifySelectedCategories",...
        "SelectedCategories",...
        "LeftExtension",...
        "RightExtension",...
        "LeftShortening",...
        "RightShortening",...        
        "MergeDistance",...
        "MinLength"];        
end

methods
    function obj = signalMask(src,varargin)
        %signalMask Construct a signalMask object
        validateSourceMask(obj,src);
        iParseNameValues(obj,varargin{:});
    end
    
    %----------------------------------------------------------------------
    function [tbl,numRegions,cats] = roimask(obj)
        %roimask Get ROI table mask
        %   TBL = roimask(M) returns an ROI table mask, TBL, based on the
        %   input source mask and the properties specified in signalMask
        %   object, M. When 'RightExtension' is nonzero, resulting region
        %   limits can go beyond the sequence lengths when 'SourceType' is
        %   'categoricalSequence' or 'binarySequences'. When 'SampleRate' is
        %   specified, region limits in the ROI table, TBL, are in seconds.
        %
        %   [TBL,NUMROI,CATS] = roimask(M) returns NUMROI, a vector
        %   containing the number of regions found for each of the
        %   categories listed in string vector CATS.
        % 
        %   % EXAMPLE 1:
        %      % Convert a mask of binary sequences to an ROI table mask.
        %      % Extend the regions of interest by one sample to the left 
        %      % and right.
        %      binSeqs = logical([ ...
        %         0 0 0 0 1 1 1 0 0 0 0 0 0 0 1 1 1 1;
        %         1 1 0 0 0 0 0 0 0 1 1 1 0 0 0 0 0 0]');
        %      m = signalMask(binSeqs,'Categories',["A" "B"]);
        %      m.RightExtension = 1;
        %      m.LeftExtension = 1;
        %      roiTbl = roimask(m)
        %   
        %   % EXAMPLE 2:
        %      % Convert a categorical mask to an ROI table mask.
        %      % Merge same-category regions separated by only one sample.
        %      catSeq = categorical(["A" "A" "A" missing missing "B" "B" ...
        %           missing missing "B" "B" missing "B" "A" "A" "A"]);
        %      m = signalMask(catSeq);
        %      m.MergeDistance = 1;
        %      roiTbl = roimask(m)
        
        cats = getSelectedCategoryStrings(obj);        
        
        if obj.SourceType == "roiTable"
            [roiTblVarNames,roiMatrix,roiLabels] = parseTable(obj);
            [tbl,numRegions,cats] = tableMask2table(obj,roiTblVarNames,roiMatrix,roiLabels,cats);
            
        elseif obj.SourceType == "categoricalSequence"
            seq = obj.pSourceMask;
            [tbl,numRegions,cats] = catseqMask2table(obj,seq,cats);
            
        else
            % Select the sequences that belong to selected categories only
            catsIdx = getSelectedCategories(obj);
            seqs = obj.pSourceMask(:,catsIdx);
            [tbl,numRegions,cats] = binseqsMask2table(obj,seqs,cats);
        end
    end
    
    %----------------------------------------------------------------------
    function [seq,numRegions,cats] = catmask(obj,varargin)
        %catmask Get categorical sequence mask
        %   SEQ = catmask(M) returns a categorical sequence mask, SEQ,
        %   based on the input source mask and the properties specified in
        %   the signalMask object, M. The length of the output sequence is
        %   equal to the length of the source mask sequences when
        %   'SourceType' property is 'categoricalSequence' or
        %   'binarySequences'. When 'SourceType' is 'roiTable', you must
        %   specify a sequence length.
        %
        %   Samples in the output sequence, SEQ, that do not belong to a
        %   region of interest (and hence have no label value) are set to
        %   the missing categorical value <undefined>.
        %
        %   SEQ = catmask(M,L) specifies the length of output sequence,
        %   SEQ, as an integer scalar, L. Regions outside of length L are
        %   ignored. Output sequence is be padded with missing values when
        %   L is greater than the source sequence (in the case where
        %   'SourceType' is 'categoricalSequence' or 'binarySequences') or
        %   greater than the maximum region index (in the case where
        %   'SourceType' is 'roiTable').
        %
        %   When 'RightExtension' is nonzero and 'SourceType' is
        %   'categoricalSequence' or 'binarySequences', signalMask extends
        %   regions possibly beyond the sequence length, applies all other
        %   modifications based on 'LeftExtension', 'RightShortening',
        %   'LeftShortening', 'MergeDistance', and 'MinLength', and then
        %   truncates resulting sequence to the original sequence length,
        %   or to the specified length L.
        %
        %   SEQ = catmask(...,'OverlapAction',ACTION) specifies how
        %   signalMask deals with regions having different category values
        %   that overlap. 'OverlapAction' can be set to 'prioritizeByList'
        %   or 'error'. When set to 'prioritizeByList', you can also
        %   specify a priority list using the 'PriorityList' parameter. In
        %   this case, the first category in the list has the highest
        %   priority and all its samples are kept in cases of overlap. Then
        %   the next priority category region is kept in the remaining
        %   non-overlapping samples, and so on. When set to 'error', the
        %   function throws an error if there are overlaps between regions
        %   with different categories. If 'PriorityList' is not specified,
        %   catmask uses the list in the 'Categories' property of the
        %   signalMask object. The default is 'error'.
        %
        %   signalMask modifies sequences based on the 'LeftExtension',
        %   'RightExtension', LeftShortening', 'RightShortening'
        %   'MergeDistance', and 'MinLength' properties and then manages
        %   overlap based on the 'OverlapAction' value.
        %
        %   SEQ = catmask(...,'OverlapAction','prioritizeByList','PriorityList',IDXLIST)
        %   specifies a vector of integer indices, IDXLIST, pointing to the
        %   'Categories' elements and ordered by priority with which they
        %   should be treated when regions with different category values
        %   overlap. IDXLIST must contain indices for all the elements in
        %   'Categories'. The first category in the list has the highest
        %   priority and all its samples are kept in cases of overlap. Then
        %   the next priority category region is kept in the remaining
        %   non-overlapping samples, and so on. If 'PriorityList' is not
        %   specified, catmask prioritizes categories following the order
        %   of the elements in the 'Categories' property.
        %
        %   [SEQ,NUMROI,CATS] = catmask(...) returns NUMROI, a vector
        %   containing the number of regions found for each of the
        %   categories listed in string vector CATS.
        % 
        %   % EXAMPLE 1:
        %      % Convert a mask of binary sequences to a categorical mask.
        %      % Remove the regions of interest shorter than 3 samples.
        %      binSeqs = logical([ ...
        %         0 0 0 1 1 1 0 0 0 1 1 0 0 1 1 1 1 0;
        %         1 1 0 0 0 0 0 0 1 1 1 0 0 0 0 0 1 1]');
        %      m = signalMask(binSeqs,'Categories',["A" "B"]);
        %      m.MinLength = 3;
        %      catSeq = catmask(m)'
        %
        %   % EXAMPLE 2:
        %      % Convert an ROI table to a categorical mask of length 35.
        %      roiTbl = table([2 5; 7 10; 15 25; 28 30],["A" "B" "B" "A"]');
        %      m = signalMask(roiTbl);
        %      catSeq = catmask(m,35)'

        persistent inpPCatMask;
        if isempty(inpPCatMask)
            inpPCatMask = inputParser;
            addOptional(inpPCatMask, 'L',[]);
            addParameter(inpPCatMask, 'OverlapAction',"error");
            addParameter(inpPCatMask, 'PriorityList',[]);
        end
        
        parse(inpPCatMask, varargin{:});
        parsedStruct = inpPCatMask.Results;
                
        L = [];
        if ~isempty(parsedStruct.L)
            validateattributes(parsedStruct.L,{'numeric'},{'integer','positive','scalar'},'catmask','L');
            L = parsedStruct.L;
        end        
                
        overlapAction = validatestring(parsedStruct.OverlapAction,["error","prioritizeByList"],'catmask','OverlapAction');
        
        priorityList = [];
        if ~isempty(parsedStruct.PriorityList)            
            if (overlapAction == "error")
                error(message('signal:internal:segmentation:PriorityListNotApplies'));
            end                                                      
            
            validateattributes(parsedStruct.PriorityList,{'numeric'},{'vector','positive','numel',obj.pNumCategories},'catmask','PriorityList');
            if (numel(parsedStruct.PriorityList) ~= numel(unique(parsedStruct.PriorityList))) || any(parsedStruct.PriorityList > obj.pNumCategories)
                error(message('signal:internal:segmentation:PriorityListInvalid'));
            end
           priorityList =  parsedStruct.PriorityList;
        end
        
        % Get selected categories
        cats = getSelectedCategoryStrings(obj);
        
        if isempty(priorityList)
            priorityListStrings = cats;
        else
            priorityListStrings = obj.pCategories(priorityList);
            % Keep only strings in cats in the priority order
            priorityListStrings = priorityListStrings(ismember(priorityListStrings,cats));
        end
        
        getNumRegionsFlag = (nargout > 1);
        
        if obj.SourceType == "roiTable"
            if isempty(L)
                error(message('signal:internal:segmentation:MustSpecifySequenceLength','catmask'))
            end            
            % Parse table content
            [~,roiMatrix,roiLabels] = parseTable(obj);                                   
            [seq,numRegions,cats] = tableMask2catseq(obj,L,roiMatrix,roiLabels,cats,overlapAction,priorityListStrings,getNumRegionsFlag);
            
        elseif obj.SourceType == "categoricalSequence"            
            if isempty(L)
                 L = numel(obj.pSourceMask);
            end
            % Get reversed priority list as this is more useful to create
            % prioritized sequence
            revPriorityList = flipud(priorityListStrings); 
            [seq,numRegions,cats] = catseqMask2catseq(obj,L,obj.pSourceMask,cats,overlapAction,revPriorityList,getNumRegionsFlag);
            
        else            
            if isempty(L)
                L = size(obj.pSourceMask,1);
            end
            % Put sequences in order of reversed priority list and reverse
            % priorityList
            allCats = obj.Categories;            
            revPriorityList = flipud(priorityListStrings);            
            [~,I] = ismember(revPriorityList,allCats);           
            revSeqs = obj.pSourceMask(:,I);               
            [seq,numRegions,cats] = binseqsMask2catseq(obj,L,revSeqs,cats,overlapAction,revPriorityList,getNumRegionsFlag);
        end
    end
    
    %----------------------------------------------------------------------
    function [seqs,numRegions,cats] = binmask(obj,L)
        %binmask Get matrix of binary sequences mask
        %   SEQS = binmask(M) returns a P-column matrix of binary sequences
        %   mask, SEQS, based on the input source mask and the properties
        %   specified in the signalMask object, M. When
        %   'SpecifySelectedCategories' is false, P is the number of
        %   categories in the 'Categories' property of M. In this case the
        %   i-th column in SEQS contains a binary mask sequence for the
        %   i-th category listed in the 'Categories' property. When
        %   'SpecifySelectedCategories' is true, P is the number of
        %   categories specified in the 'SelectedCategories' property of M.
        %   In this case the i-th column in SEQS contains a binary mask
        %   sequence for the i-th category listed in the
        %   'SelectedCategories' property. The length of the output
        %   sequences is equal to the length of the source mask sequences
        %   when 'SourceType' property is 'categoricalSequence' or
        %   'binarySequences'. When 'SourceType' is 'roiTable', you must
        %   specify a sequence length.
        %
        %   SEQS = binmask(M,L) specifies the length of output sequences,
        %   SEQS, as an integer scalar L. Regions outside of length L are
        %   ignored. The output sequences are padded with missing values
        %   when L is greater than the source sequence (in the case where
        %   'SourceType' is 'categoricalSequence' or 'binarySequences') or
        %   greater than the maximum region index (in the case where
        %   'SourceType' is 'roiTable').
        %
        %   When 'RightExtension' is nonzero and 'SourceType' is
        %   'categoricalSequence' or 'binarySequences', signalMask extends
        %   regions possibly beyond the sequence length, applies all other
        %   modifications based on 'LeftExtension', 'RightShortening',
        %   'LeftShortening', 'MergeDistance', and 'MinLength', and then
        %   truncates the resulting sequences to the original sequence
        %   length, or to the specified length L.
        %
        %   [SEQS,NUMROI,CATS] = binmask(...) returns NUMROI, a vector
        %   containing the number of regions found for each of the
        %   categories listed in string vector CATS.
        %
        %   % EXAMPLE:
        %      % Convert an ROI table to binary mask sequences of length 35.
        %      % Shorten regions by one sample from the right.
        %      roiTbl = table([2 5; 7 10; 15 25; 28 30],["A" "B" "B" "A"]');
        %      m = signalMask(roiTbl);
        %      m.RightShortening = 1;
        %      binSeqs = binmask(m,35)' 
        
        if nargin == 2
            validateattributes(L,{'numeric'},{'integer','positive','scalar'},'binmask','L');
            L = double(L);
        else
            L = [];            
        end
                
        cats = getSelectedCategoryStrings(obj);        
        
        if obj.SourceType == "roiTable"
            % Parse table content
            if isempty(L)
                error(message('signal:internal:segmentation:MustSpecifySequenceLength','binmask'))
            end
            
            [~,roiMatrix,roiLabels] = parseTable(obj);            
            [seqs,numRegions,cats] = tableMask2binseqs(obj,L,roiMatrix,roiLabels,cats);
            
        elseif obj.SourceType == "categoricalSequence"
            if nargin < 2
                L = numel(obj.pSourceMask);
            end
            [seqs,numRegions,cats] = catseqMask2binseqs(obj,L,obj.pSourceMask,cats);
        else
            if nargin < 2
                L = size(obj.pSourceMask,1);
            end     
            catsIdx = getSelectedCategories(obj);
            seqs = obj.pSourceMask(:,catsIdx);            
            [seqs,numRegions,cats] = binseqsMask2binseqs(obj,L,seqs,cats);
        end
    end
       
    %----------------------------------------------------------------------
    function [sigROI, limits, numRegions, cats] = extractsigroi(obj,x,varargin)
        %extractsigroi Extract regions of interest based on signal mask
        %   SIGROI = extractsigroi(M,X) extracts regions of input signal
        %   vector X based on the input source mask and the properties
        %   specified in the signalMask object, M. Output SIGROI is a cell
        %   array. Each element of the cell array contains a cell array
        %   with signal regions extracted for each category in the signal
        %   mask object, M.
        %
        %   extractsigroi modifies and sorts the regions for each category
        %   based on the property settings in object M before extracting
        %   the signal samples.
        %
        %   SIGROI = extractsigroi(...,'ConcatenateRegions',FLAG) specifies
        %   whether to concatenate extracted signal regions for each
        %   category. If FLAG is set to false, output is a cell array with
        %   cell arrays of individual signal regions. When FLAG is set to
        %   true, output is a cell array with each element containing a
        %   vector of concatenated extracted signal regions for each
        %   category in signalMask object, M. If omitted,
        %   'ConcatenateRegions' defaults to false.
        %
        %   SIGROI = extractsigroi(...,'SelectedRegions',SELROIS) specifies
        %   selected regions as a vector of integers. If SELROIS equals 1,
        %   then only the first region of each category is extracted and
        %   returned in SIGROI. If SELROIS equals [i,j,k,...] then i-th,
        %   j-th, k-th, and so on, regions of each category are extracted
        %   and returned in SIGROI. Indices larger than the number of
        %   regions available for a category are ignored. If
        %   'SelectedRegions' is omitted, extractsigroi extracts all signal
        %   regions.
        %
        %   [SIGROI,LIMITS] = extractsigroi(...) returns a cell array,
        %   LIMITS, with each element containing a two-column matrix with
        %   the region limit indices corresponding to the extracted signal
        %   regions. If 'SampleRate' is specified in the signalMask object,
        %   LIMITS is in seconds, otherwise, LIMITS contains integers
        %   pointing to signal sample indices.
        %
        %   [SIGROI,LIMITS,NUMROI,CATS] = extractsigroi(...) returns vector
        %   NUMROI, containing the number of regions found for each of the
        %   categories listed in string vector CATS.
        %
        %   % EXAMPLE 1:
        %      % Extract signal regions based on mask defined by an ROI
        %      % table. Ignore regions shorter than three samples.
        %      roiTbl = table([2 5; 7 10; 12 13; 15 25; 28 30],["A" "B" "A" "B" "A"]');
        %      m = signalMask(roiTbl);
        %      m.MinLength = 3;
        %      x = 1:30;
        %      [sigs,lims] = extractsigroi(m,x);
        %      % Concatenate extracted regions
        %      sigsConcat = extractsigroi(m,x,'ConcatenateRegions',true);
        %      sigsConcat{1}'
        %      sigsConcat{2}'
        %
        %   % EXAMPLE 2:
        %      % Extract signal regions based on a binary mask.
        %      binSeq = logical([0 0 0 1 1 1 1 0 0 0 1 1 1 0 0 0 0 1 1 0 0 1 1 0 0]);
        %      m = signalMask(binSeq);
        %      x = randn(25,1);
        %      [sigs,lims] = extractsigroi(m,x);
        %      figure
        %      plotsigroi(m,x)

        if ~isempty(x)
            validateattributes(x,{'numeric'},{'vector'},'extractsigroi','X');
        end
        persistent inpPExtractSig;
        if isempty(inpPExtractSig)
            inpPExtractSig = inputParser;            
            addParameter(inpPExtractSig, 'ConcatenateRegions',false);
            addParameter(inpPExtractSig, 'SelectedRegions',[]);
        end
        
        parse(inpPExtractSig, varargin{:});
        parsedStruct = inpPExtractSig.Results;
        
        validateattributes(parsedStruct.ConcatenateRegions,{'logical','numeric'},{'scalar','real','finite'},'extractsigroi','ConcatenateRegions');        
        concatRegions = logical(parsedStruct.ConcatenateRegions);
        selectedRegions = [];
        if ~isempty(parsedStruct.SelectedRegions)
            validateattributes(parsedStruct.SelectedRegions,{'numeric'},{'vector','integer','positive'},'extractsigroi','SelectedRegions');
            selectedRegions = unique(parsedStruct.SelectedRegions,'sorted');
        end

        % Get table masks with transformed regions based on properties of
        % signalMask object.
        [roiTbl,numRegions,cats] = roimask(obj);  
        [~,roiMatrix,roiLabels] = parseTableWithTime(obj,roiTbl);
                
        sigROI = cell(numel(cats),1);
        limits = cell(numel(cats),1);
        
        for idx = 1:numel(cats)     
            % Get the region limit rows corresponding to the current
            % category of interest
            lblIdx = (roiLabels == cats(idx));
            if ~any(lblIdx)
                limits{idx} = zeros(0,2,'like',roiMatrix);
                if concatRegions
                    sigROI{idx} = zeros(0,1,'like',x);
                else
                    sigROI{idx} = cell(0,1);
                end
                continue;
            end
            currentROIMatrix = roiMatrix(lblIdx,:);
                                    
            if ~isempty(selectedRegions)
                currentSelectedRegions = selectedRegions(selectedRegions <= size(currentROIMatrix,1));
                if isempty(currentSelectedRegions)
                    numRegions(idx) = 0;
                    limits{idx} = zeros(0,2,'like',roiMatrix);
                    if concatRegions
                        sigROI{idx} = zeros(0,1,'like',x);
                    else
                        sigROI{idx} = cell(0,1);
                    end
                    continue;
                else
                    currentROIMatrix = currentROIMatrix(currentSelectedRegions,:);
                    numRegions(idx) = size(currentROIMatrix,1);
                end
            end

            % Extract signal regions
            sigROI{idx} = reshape(signal.internal.segmentation.extractsigroi(x,currentROIMatrix,concatRegions),[],1);               
            limits{idx} = (currentROIMatrix - obj.pSampleRateSpecified)/getSampleRate(obj);            
        end        
    end
    
    %----------------------------------------------------------------------
    function varargout = plotsigroi(obj,x,patchFlag)  
        %plotsigroi Plot signal regions based on signal mask
        %   plotsigroi(M,X) plots signal vector X with color-coded regions
        %   based on the signalMask object, M. If X is complex-valued,
        %   plotsigroi plots its magnitude.
        %
        %   plotsigroi(M,X,PATCHFLAG) plots regions of interest using
        %   rectangular patches when PATCHFLAG is true. If omitted,
        %   PATCHFLAG defaults to false. Use this option if you want to see
        %   region overlaps.
        %
        %   H = plotsigroi(...) returns the axes handle, H, of the
        %   color-coded plot.
        %
        %   % EXAMPLE:
        %      % Plot signal regions defined by an ROI table mask.
        %      roiTbl = table([2 5; 7 10; 12 13; 15 25; 28 30],["A" "B" "A" "B" "A"]');
        %      m = signalMask(roiTbl);
        %      x = randn(35,1);
        %      figure
        %      plotsigroi(m,x)
        %      figure
        %      plotsigroi(m,x,true)

        validateattributes(x,{'numeric'},{'vector'},'plotsigroi','X');
        if isrow(x)
            x = x(:);
        end
        if ~isreal(x)
            x = abs(x);
        end
        
        if nargin < 3
            patchFlag = false;
        else
            validateattributes(patchFlag,{'logical','numeric'},{'scalar','real','finite'},'plotsigroi','PATCHFLAG');
            patchFlag = logical(patchFlag);
        end
                
        if patchFlag
            [tbl,~,cats] = roimask(obj);
            [~,roiMatrix,roiLabels] = parseTable(obj,tbl);
            [roiMatrix, I] = signal.internal.segmentation.truncateROIMatrix(roiMatrix,(length(x) - obj.pSampleRateSpecified)/obj.getSampleRate);
            roiLabels = roiLabels(I);
        else
            [logicalSeqs,~,cats] = binmask(obj,length(x));
        end
        % Time vector starts at 0 seconds if sample rate specified,
        % otherwise it starts at sample index 1.
        tVect = ((1:length(x))' - obj.pSampleRateSpecified) /obj.getSampleRate;                   
        
        hLine = plot(tVect,x,'k','LineWidth',0.5);    
        if obj.pSampleRateSpecified
            xlabel(getString(message('signal:internal:segmentation:Seconds')))
        else
            xlabel(getString(message('signal:internal:segmentation:Samples')))
        end
      
        hAx = ancestor(hLine,'axes');
        if tVect(1) ~= tVect(end)
            hAx.XLim = [tVect(1) tVect(end)];
        end
        if min(x) ~= max(x)
            hAx.YLim = [min(x) max(x)];
        end
        clrs = lines(numel(cats));
        currentYLimits = hAx.YLim;
        
        if patchFlag
            for idx = 1:numel(cats)
                lblIdx = (roiLabels == cats(idx));
                if ~any(lblIdx)
                    continue;
                end
                currentROIMatrix = roiMatrix(lblIdx,:);
                for kk = 1:size(currentROIMatrix,1)
                    aRegion = currentROIMatrix(kk,:);
                    xCoords = [aRegion(1) aRegion(2) aRegion(2) aRegion(1)];                    
                    yCoords = [currentYLimits(1) currentYLimits(1) currentYLimits(2) currentYLimits(2)];
                    patch(hAx,xCoords,yCoords,clrs(idx,:),'FaceAlpha',.3,'EdgeColor',clrs(idx,:),'LineWidth',1);                    
                end
            end                                                
        else
            for idx = 1:numel(cats)
                y = x;
                y(~logicalSeqs(:,idx)) = NaN;
                line(hAx,tVect,y,'Color',clrs(idx,:));
            end
        end
      
        grid on;
        
        % Create color bar with first category at the top of the colorbar.
        % Use class names for tick marks.
        cmap = flipud(clrs);
        tickLbls = flipud(cats);
        colormap(hAx,cmap)
        c = colorbar(hAx);
        c.TickLabels = tickLbls;
        % Center tick labels
        numCategories = numel(cats);
        c.Ticks = 1/(numCategories*2):1/numCategories:1;
        % Remove tick mark
        c.TickLength = 0;
        
        if nargout > 0
            varargout{1} = hAx;
        end
    end
end % end methods

%--------------------------------------------------------------------------
% Setter/getter methods
%--------------------------------------------------------------------------
methods
    function val = get.SampleRate(obj)
        if ~obj.pSampleRateSpecified
            error(message('signal:internal:segmentation:SampleRateNotApplies'));
        end        
        val = obj.pSampleRate;
    end
    
    function set.SampleRate(obj,~)  
        % Sample rate only applies when it was specified at constructor
        % time. If it applies, then it is read-only.
        if ~obj.pSampleRateSpecified
            error(message('signal:internal:segmentation:SampleRateNotApplies'));
        end        
        error(message('signal:internal:segmentation:SampleRateReadOnly'));
    end    
    
    function set.LeftExtension(obj,val)
        validateattributes(val,{'numeric'},{'scalar','integer','nonnegative'},'signalMask','LeftExtension');
        obj.LeftExtension = val;
    end
    
    function set.RightExtension(obj,val)
        validateattributes(val,{'numeric'},{'scalar','integer','nonnegative'},'signalMask','RightExtension');
        obj.RightExtension = val;
    end
    
    function set.LeftShortening(obj,val)
        validateattributes(val,{'numeric'},{'scalar','integer','nonnegative'},'signalMask','LeftShortening');
        obj.LeftShortening = val;
    end
    
    function set.RightShortening(obj,val)
        validateattributes(val,{'numeric'},{'scalar','integer','nonnegative'},'signalMask','RightShortening');
        obj.RightShortening = val;
    end
        
    function set.MergeDistance(obj,val)
        validateattributes(val,{'numeric'},{'scalar','integer','nonnegative'},'signalMask','MergeDistance');
        obj.MergeDistance = val;
    end
    
    function set.MinLength(obj,val)
        validateattributes(val,{'numeric'},{'scalar','integer','positive'},'signalMask','MinLength');
        obj.MinLength = val;
    end
           
    function set.Categories(obj,val)
        val = validateSetCategories(obj,val);
        obj.pCategories = val;
    end
    function val = get.Categories(obj)
        val = obj.pCategories;
    end
    
    function set.SpecifySelectedCategories(obj,val)
        validateattributes(val,{'numeric','logical'},{'scalar','real','finite'},'signalMask','SpecifySelectedCategories');
        obj.SpecifySelectedCategories = logical(val);
    end
    
    function set.SelectedCategories(obj,val)
        if ~obj.SpecifySelectedCategories
            error(message('signal:internal:segmentation:SelectedCategoriesNotApplies'));
        end
        validateSelectedCategories(obj,val);
        obj.pSelectedCategories = val(:);
    end
    
    function val = get.SelectedCategories(obj)
        if ~obj.SpecifySelectedCategories
            error(message('signal:internal:segmentation:SelectedCategoriesNotApplies'));
        end        
        val = obj.pSelectedCategories;
    end
end
methods (Access = private)    
    function val = getSampleRate(obj)
        % Get sample rate value based on wheter it was specified or not
        if obj.pSampleRateSpecified
            val = obj.pSampleRate;
        else
            val = 1;
        end
    end
    
    function val = getSelectedCategoryStrings(obj)
        % Get selected category strings based on specify selected
        % categories setting
        if obj.SpecifySelectedCategories
            val = obj.Categories(obj.pSelectedCategories);
        else
            val = obj.Categories;            
        end
    end
    
    function val = getSelectedCategories(obj)
        % Get selected category strings based on specify selected
        % categories setting
        if obj.SpecifySelectedCategories
            val = obj.pSelectedCategories;
        else
            val = (1:obj.pNumCategories)';
        end
    end
end % end set/get methods

%--------------------------------------------------------------------------
% Conversion private methods
%--------------------------------------------------------------------------
methods (Access = private)
    
    function [tbl,numRegions,cats] = tableMask2table(obj,roiTblVarNames,roiMatrix,roiLabels,cats)
        % Convert source table to new table based on properties of the
        % signalMask object. Input cats is already the pruned category list
        % based on selected categories value.
        
        %num regions per category
        numRegions = zeros(numel(cats),1);
        accumROIMatrix = zeros(0,2);
        accumCats = categorical(strings(0,1),cats);
        
        % Loop through each category and its regions
        for idx = 1:numel(cats)
            
            % Get the region limit rows corresponding to the current
            % category of interest
            lblIdx = (roiLabels == cats(idx));
            if ~any(lblIdx)
                continue;
            end
            currentROIMatrix = roiMatrix(lblIdx,:);
            
            % Transform the regions based on object settings
            currentROIMatrix = transformRegionsInROIMatrix(obj,currentROIMatrix);
            foundNumRegions = size(currentROIMatrix,1);
            
            % Concatenate regions and categories
            accumROIMatrix = [accumROIMatrix; currentROIMatrix]; %#ok<AGROW>
            accumCats = [accumCats; categorical(repmat(cats(idx),foundNumRegions,1))]; %#ok<AGROW>
            numRegions(idx) = foundNumRegions;
        end
        
        % Create a table after modifications, set roi limits to time if
        % sample rate was specified, and sort table by regions        
        tbl = table((accumROIMatrix - obj.pSampleRateSpecified)/getSampleRate(obj),accumCats,'VariableNames',roiTblVarNames);
        tbl = sortROITable(obj,tbl);
    end
    
    %----------------------------------------------------------------------
    function [tbl,numRegions,cats] = binseqsMask2table(obj,seqs,cats)
        % Convert source binary sequences to a table based on properties
        % of the signalMask object.
                
        numRegions = zeros(numel(cats),1);
        accumROIMatrix = zeros(0,2);
        accumCats = categorical(strings(0,1),cats);
                        
        % Loop through each category and its regions
        for idx = 1:numel(cats)
            
            currentSeq = seqs(:,idx);
            
            % Get regions matrix from sequence            
            currentROIMatrix = signal.internal.segmentation.binmask2sigroi(currentSeq);            
            currentROIMatrix = transformRegionsInROIMatrix(obj,currentROIMatrix);            
            foundNumRegions = size(currentROIMatrix,1);
            
            % Concatenate regions and categories
            accumROIMatrix = [accumROIMatrix; currentROIMatrix]; %#ok<AGROW>
            accumCats = [accumCats; categorical(repmat(cats(idx),foundNumRegions,1))]; %#ok<AGROW>
            numRegions(idx) = foundNumRegions;
        end
        
        % Create a table after modifications, set roi limits to time if
        % sample rate was specified, and sort table by regions
        roiTblVarNames = ["ROILimits","Value"];
        tbl = table((accumROIMatrix - obj.pSampleRateSpecified)/getSampleRate(obj),accumCats,'VariableNames',roiTblVarNames);
        tbl = sortROITable(obj,tbl);
    end
    
    %----------------------------------------------------------------------
    function [tbl,numRegions,cats] = catseqMask2table(obj,seq,cats)
        % Convert source categorical sequence to a table based on
        % properties of the signalMask
        L = numel(seq);
        
        % Convert categorical sequence to binary sequences
        logicalSeqs = catSeq2binSeqs(obj,L,seq,cats);
        
        % Convert to table and apply region modifications
        [tbl,numRegions,cats] = binseqsMask2table(obj,logicalSeqs,cats);        
    end
    
    %----------------------------------------------------------------------
    function [seqs,numRegions,cats] = tableMask2binseqs(obj,L,roiMatrix,roiLabels,cats)
        % Convert source table to binary sequences based on properties of
        % the signalMask object
        
        seqs = false(L,numel(cats));
        numRegions = zeros(numel(cats),1);
        
        for idx = 1:numel(cats)
            
            lblIdx = (roiLabels == cats(idx));
            if ~any(lblIdx)
                continue;
            end
            % Get the region limit rows corresponding to the current
            % category of interest            
            currentROIMatrix = roiMatrix(lblIdx,:);
                                                           
            % Transform regions and truncate to sequence length afterwards
            currentROIMatrix = transformRegionsInROIMatrix(obj,currentROIMatrix);
            currentROIMatrix = signal.internal.segmentation.truncateROIMatrix(currentROIMatrix,L);   
            foundNumRegions = size(currentROIMatrix,1);
            numRegions(idx) = foundNumRegions;
            % Convert region matrix to binary sequence
            if ~isempty(currentROIMatrix)
                seqs(:,idx) = signal.internal.segmentation.sigroi2binmask(currentROIMatrix,L);
            end
        end        
    end
    
    %----------------------------------------------------------------------
    function [seqs,numRegions,cats] = catseqMask2binseqs(obj,L,catSeq,cats)
        % Convert categorical sequence to binary sequences based on
        % properties of the signalMask

        % Get binary sequences from the categorical sequence - truncate or
        % append missing to catSeq to match size of L
        binSeqs = catSeq2binSeqs(obj,L,catSeq,cats);
        
        % Convert binary masks to new binary sequences based on object
        % properties
        [seqs,numRegions,cats] = binseqsMask2binseqs(obj,L,binSeqs,cats);        
    end
    
    %----------------------------------------------------------------------
    function [seqs,numRegions,cats] = binseqsMask2binseqs(obj,L,logicalSeqs,cats)
        % Convert binary sequences to binary sequences based on properties
        % of the signalMask
        
        seqs = false(L,numel(cats));
        numRegions = zeros(numel(cats),1);
        
        % Truncate or append false values to logicalSeqs to match size of L
        logicalSeqs = adjustBinSeqLength(obj,L,logicalSeqs);
        
        % Convert logical sequences to roi matrices and apply conversions,
        % then convert back to logical sequences
        for idx = 1:numel(cats)
            currentLogicalSeq = logicalSeqs(:,idx);            
            currentROIMatrix = signal.internal.segmentation.binmask2sigroi(currentLogicalSeq); 
            currentROIMatrix = transformRegionsInROIMatrix(obj,currentROIMatrix);
            currentROIMatrix = signal.internal.segmentation.truncateROIMatrix(currentROIMatrix,L);
            foundNumRegions = size(currentROIMatrix,1);
            numRegions(idx) = foundNumRegions;
            if ~isempty(currentROIMatrix)
                seqs(:,idx) = signal.internal.segmentation.sigroi2binmask(currentROIMatrix,L);
            end
        end
    end
    
    %----------------------------------------------------------------------
    function [seq,numRegions,cats] = tableMask2catseq(obj,L,roiMatrix,roiLabels,cats,overlapAction,priorityList,getNumRegionsFlag)
        % Convert a table to a categorical sequence based on properties of
        % the signalMask. Priority list has already been truncated to
        % contain only the elements in selected categories cats.
                       
        % Return empty numRegions unless it was requested
        numRegions = [];
        
        % Convert table to binary sequences with all modifications. Get
        % sequences in reversed order of priorityList.
        revPriorityList = flipud(priorityList);
        revLogicalSeqs = tableMask2binseqs(obj,L,roiMatrix,roiLabels,revPriorityList);
                                               
        if overlapAction == "error" && checkOverlapInLogicalSeqs(obj,revLogicalSeqs)
            error(message('signal:internal:segmentation:OverlapError'));
        end
        
        % Convert logical sequences to a cat sequence with the order
        % defined in revPriorityList
        seq = logicalSeqs2catSeq(obj,revLogicalSeqs,revPriorityList,cats,L);
                        
        if getNumRegionsFlag
            % Count number of resulting regions for each category in cats
            numRegions = countRegionsOnCatSeq(obj,seq,cats);
        end                
    end

    %----------------------------------------------------------------------
    function [seq,numRegions,cats] = catseqMask2catseq(obj,L,catSeq,cats,overlapAction,revPriorityList,getNumRegionsFlag)
        % Convert a categorical sequence to a categorical sequence based on
        % properties of the signalMask. revPriorityList is assumed to
        % contain highest priority in the last element (i.e. we assume
        % reversed order here with respect to priorityList input).
        
        % Return empty numRegions unless it was requested
        numRegions = [];
        
        % Convert to binary sequences with all modifications
        revLogicalSeqs = catseqMask2binseqs(obj,L,catSeq,revPriorityList);
        
        if overlapAction == "error" && checkOverlapInLogicalSeqs(obj,revLogicalSeqs)
            error(message('signal:internal:segmentation:OverlapError'));
        end       
        
        % Convert logical sequences to a cat sequence with the order
        % defined in revPriorityList
        seq = logicalSeqs2catSeq(obj,revLogicalSeqs,revPriorityList,cats,L);
        
        if getNumRegionsFlag
            % Count number of resulting regions for each category in cats
            numRegions = countRegionsOnCatSeq(obj,seq,cats);
        end
    end
    
    %----------------------------------------------------------------------
    function [seq,numRegions,cats] = binseqsMask2catseq(obj,L,revLogicalSeqs,cats,overlapAction,revPriorityList,getNumRegionsFlag)
        % Convert a binary sequence to a categorical sequence based on
        % properties of the signalMask. This function assumes that input
        % revLogicalSequences are already in the order of the list in
        % revPriorityList. revPriorityList is assumed to contain highest
        % priority in the last element (i.e. we assume reversed order
        % here with respect to priorityList input).
        
        % Return empty numRegions unless it was requested
        numRegions = [];
        
        % Convert to binary sequences with all modifications
        revLogicalSeqs = binseqsMask2binseqs(obj,L,revLogicalSeqs,revPriorityList);
               
        if overlapAction == "error" && checkOverlapInLogicalSeqs(obj,revLogicalSeqs)
            error(message('signal:internal:segmentation:OverlapError'));
        end
        
        % Convert logical sequences to a cat sequence with the order
        % defined in revPriorityList
        seq = logicalSeqs2catSeq(obj,revLogicalSeqs,revPriorityList,cats,L);
        
        if getNumRegionsFlag
            % Count number of resulting regions for each category in cats
            numRegions = countRegionsOnCatSeq(obj,seq,cats);
        end
    end   
end

%--------------------------------------------------------------------------
% Conversion helper functions
%--------------------------------------------------------------------------
methods (Access = private)
    
    function [tblVarNames,roiMatrix,roiLabels] = parseTable(obj,tbl)
        % Parse roi table contents. Return an roiMatrix in samples (i.e.
        % converted to samples if Fs was specified) and labels as a
        % categorical array.
        
        % obj.pSourceMask contains an ROI table or error out
        if nargin < 2
            tbl = obj.pSourceMask;                    
        end
        
        % obj.pSourceMask contains an ROI table
        tblVarNames = tbl.Properties.VariableNames;
        
        % Get roi matrix
        roiMatrix = tbl{:,1};
        
        % Get labels
        roiLabels = tbl{:,2};        
    end
    
    %----------------------------------------------------------------------
    function [tblVarNames,roiMatrix,roiLabels] = parseTableWithTime(obj,tbl)
        % Parse table and convert time values to samples if working with
        % sample rate                  
         [tblVarNames,roiMatrix,roiLabels] = parseTable(obj,tbl);
         roiMatrix = round(roiMatrix*getSampleRate(obj)) + obj.pSampleRateSpecified;
    end
    
    %----------------------------------------------------------------------
    function roiMatrix = transformRegionsInROIMatrix(obj,roiMatrix,maxIdx)
        % Transform regions of an roi matrix based on object properties.
        % Extend to left and right, merge, and remove short regions. Treat
        % regions of each category independently.

        if nargin < 3
            maxIdx = realmax;
        end
        
        % Extend regions to left and right
        roiMatrix = signal.internal.segmentation.extendsigroi(roiMatrix,...
            obj.LeftExtension,obj.RightExtension,maxIdx);       

        roiMatrix = signal.internal.segmentation.shortensigroi(roiMatrix,...
            obj.LeftShortening,obj.RightShortening);       
        
        % Merge regions based on MergeDistance and MergeRule properties
        roiMatrix = signal.internal.segmentation.mergesigroi(roiMatrix,obj.MergeDistance);  
                
        % Remove regions based on MinLength property - removesigroi removes
        % regions of length L or smaller so subtract 1 from MinLength as we
        % do not want to remove MinLength long regions only regions shorter
        % than that.
        roiMatrix = signal.internal.segmentation.removesigroi(roiMatrix,obj.MinLength-1);
    end        
    
    %----------------------------------------------------------------------
    function tbl = sortROITable(~,roiTbl)
        % Sort table rows based on the first column of the roi index matrix
        roiMatrix = roiTbl{:,1};
        [~,sortI] = sortrows(roiMatrix,1);        
        tbl = roiTbl(sortI,:);
    end
    
    %----------------------------------------------------------------------
    function logicalSeqs = catSeq2binSeqs(obj,L,seq,cats)
        % Convert categorical sequence to multiple binary sequences
        
        % Truncate or append missing values to categorical sequence to
        % match size of L  
        seq = adjustCatSeqLength(obj,L,seq);
        
        logicalSeqs = false(numel(seq),numel(cats));
        for idx = 1:numel(cats)
            logicalSeqs(:,idx) = (seq == cats(idx));
        end
    end   
    
    %----------------------------------------------------------------------
    function catSeq = logicalSeqs2catSeq(~,logicalSeqs,cats,valueSetCats,L)
        % Convert logical seqs to cat seq with priority in the reverse
        % order of cats. Order for logicalSeqs and cats is assumed to be
        % the same. NOTICE that last element in cats has the highest
        % priority as the regions for that category will overwrite any
        % other overlapping regions.
        
        catSeq = repmat(categorical(missing,valueSetCats),L,1);        
        
        for idx = 1:numel(cats)
            currentCat = cats(idx);
            currentSeq = logicalSeqs(:,idx);  
            catSeq(currentSeq) = currentCat;            
        end
    end
   
    %----------------------------------------------------------------------
    function logicalSeqs = adjustBinSeqLength(~,L,logicalSeqs)
        % Truncate or append false values to logicalSeqs to match size of L
         N = size(logicalSeqs,1);
         
        if N > L
            logicalSeqs = logicalSeqs(1:L,:);
        end        
        if N < L
            logicalSeqs(end+1:end+L-N,:) = false;
        end
    end
    
    %----------------------------------------------------------------------
    function seq = adjustCatSeqLength(~,L,seq)
        % Truncate or append missing values to categorical sequence to
        % match size of L        
        N = numel(seq);
        if N > L
            seq = seq(1:L,1);
        end
        if N < L
            seq(end+1:end+L-N,1) = categorical(missing);
        end
    end
    
    %----------------------------------------------------------------------
    function overlapFlag = checkOverlapInLogicalSeqs(~,logicalSeqs)
        % Check if regions on different logical sequences overlap - output
        % is true if any region overlaps. This method assumes oinput
        % sequences are of type logical.
        overlapFlag = any( sum(logicalSeqs,2) > 1);         
    end
end

%--------------------------------------------------------------------------
% Count regions
%--------------------------------------------------------------------------
methods (Access = private)
    function numRegions = countRegionsOnBinSeqs(~,seqs)
        % Count number of regions on P binary sequences        
        P = size(seqs,2);
        numRegions = zeros(P,1);
        for idx = 1:P            
            currentLogicalSeq = seqs(:,idx);
            roiMatrix = signal.internal.segmentation.binmask2sigroi(currentLogicalSeq);
            numRegions(idx) = size(roiMatrix,1);
        end
    end
    
    %--------------------------------------------------------------------------
    function numRegions = countRegionsOnCatSeq(obj,seq,cats)
        % Count number of regions on a cat sequence
        L = numel(seq);
        logicalSeqs = catSeq2binSeqs(obj,L,seq,cats);
        numRegions = countRegionsOnBinSeqs(obj,logicalSeqs);
    end
end

%--------------------------------------------------------------------------
% Parsing and validation private methods
%--------------------------------------------------------------------------
methods (Access = private)
    
    function validateSourceMask(obj,src)
        %Validate source mask
        obj.pIsInherentCategories = false;
        
        if iscategorical(src)
            if ~isvector(src)
                error(message('signal:internal:segmentation:InvalidSourceType'));
            end
            src = src(:);
            obj.pIsInherentCategories = true;
            obj.pCategories = string(categories(src));
            obj.pNumCategories = numel(obj.pCategories);
            obj.pSelectedCategories = (1:obj.pNumCategories)';
            obj.SourceType = "categoricalSequence";
            obj.pSourceMask = src;
            
            return;
        end
        
        if istable(src)
            % Always store a table mask using categorical labels even if
            % input table has string labels
            
            labelDataTypeError = false;
            
            if size(src,2) ~= 2
                error(message('signal:internal:segmentation:InvalidROITable'));
            end            
            validateattributes(src{:,1},{'numeric'},{'ncols',2,'real','finite','nonnegative'},'signalMask','First variable of input table');
            % Verify matrix of indices has non-decreasing intervals
            if any(diff(src{:,1},1,2) < 0)                
                error(message('signal:internal:segmentation:InvalidROIMatrix'));
            end
            labels = src{:,2};
            if iscategorical(labels) || isstring(labels)
                if ~isvector(labels)
                    labelDataTypeError = true;
                end
                if isstring(labels)
                    labels = categorical(labels);
                    src{:,2} = labels;
                end
            else
                labelDataTypeError = true;
            end
            
            if  labelDataTypeError
                error(message('signal:internal:segmentation:InvalidROITblCats'))
            end
                        
            obj.pIsInherentCategories = true;
            obj.pCategories = string(categories(labels));
            obj.pNumCategories = numel(obj.pCategories);
            obj.pSelectedCategories = (1:obj.pNumCategories)';
            obj.SourceType = "roiTable";
            obj.pSourceMask = src;
            
            return;
        end
        
        if ismatrix(src) && (islogical(src) || isnumeric(src))
            if isnumeric(src) && ~all(all(src == 0 | src == 1))
                error(message('signal:internal:segmentation:InvalidSourceType'));
            end
            src = logical(src);                
            if isvector(src)
                src = src(:);
            end
            obj.pIsInherentCategories = false;
            obj.pCategories = string(categories(categorical(1:size(src,2))));
            obj.pNumCategories = numel(obj.pCategories);
            obj.pSelectedCategories = (1:obj.pNumCategories)';
            obj.SourceType = "binarySequences";
            obj.pSourceMask = src;
            
            return;
        end
        
        error(message('signal:internal:segmentation:InvalidSourceType'));
    end
    
    %----------------------------------------------------------------------
    function parsedStruct = iParseNameValues(obj,varargin)
        % Parse name value inputs
        persistent inpP;
        if isempty(inpP)
            inpP = inputParser;
            addParameter(inpP, 'SampleRate',[]);
            addParameter(inpP, 'SelectedCategories',[]);  
            addParameter(inpP, 'Categories',[]);            
            addParameter(inpP, 'LeftExtension',0);
            addParameter(inpP, 'RightExtension',0);
            addParameter(inpP, 'LeftShortening',0);
            addParameter(inpP, 'RightShortening',0);            
            addParameter(inpP, 'MinLength',1);            
            addParameter(inpP, 'MergeDistance',0);
        end
        
        parse(inpP, varargin{:});
        parsedStruct = inpP.Results;
        
        if ~isempty(parsedStruct.SelectedCategories)
            error(message('signal:internal:segmentation:SelectedCategoriesAtConstruction'));
        end
        if isempty(parsedStruct.SampleRate)
            obj.pSampleRateSpecified = false;
            obj.pSampleRate = 1;
            % If no sample rate has been specified, round roi limits to
            % closest integer >= 1
            if obj.SourceType == "roiTable"
                if any(any(obj.pSourceMask{:,1} == 0))
                    error(message('signal:internal:segmentation:InvalidROIMatrixNoSampleRate'));
                end
                obj.pSourceMask{:,1} = max(1,round(obj.pSourceMask{:,1}));
            end
        else
            validateattributes(parsedStruct.SampleRate,{'numeric'},{'scalar','real','finite','positive'},'signalMask','SampleRate')
            obj.pSampleRateSpecified = true;
            obj.pSampleRate = parsedStruct.SampleRate;            
            if obj.SourceType == "roiTable"
                % If sample rate has been specified, store table with index
                % roi limits. Time zero is sample number 1.
                obj.pSourceMask{:,1} = round(obj.pSourceMask{:,1}*obj.pSampleRate) + 1;
            end
        end
        
        obj.LeftExtension = parsedStruct.LeftExtension;
        obj.RightExtension = parsedStruct.RightExtension;
        obj.LeftShortening = parsedStruct.LeftShortening;
        obj.RightShortening = parsedStruct.RightShortening;        
        obj.MinLength = parsedStruct.MinLength;        
        obj.MergeDistance = parsedStruct.MergeDistance;        
        
        if ~isempty(parsedStruct.Categories)
            obj.Categories = parsedStruct.Categories;
        end
        
        obj.SpecifySelectedCategories = false;
    end
    
    %----------------------------------------------------------------------
    function val = validateSetCategories(obj,val)
        if obj.pIsInherentCategories
            error(message('signal:internal:segmentation:CategoriesDoNotApply'));
        end
        
        val = validateStringVect(obj,val,'Categories');
        
        if (numel(unique(val))~= numel(val))
            error(message('signal:internal:segmentation:MustHaveUniqueElements','Categories'));
        end
        if numel(val) ~= obj.pNumCategories
            error(message('signal:internal:segmentation:InvalidCategoriesLength'));
        end
    end
    
    %----------------------------------------------------------------------
    function validateSelectedCategories(obj,val)
        % Validate that selected categories has valid sorted indices that
        % point to the Categories list
        validateattributes(val,{'numeric'},{'vector','integer'},'signalMask','SelectedCategories')

       if (numel(unique(val))~= numel(val)) || ~all(val == sort(val)) || max(val) > numel(obj.pCategories) || (any(val <=0))
            error(message('signal:internal:segmentation:InvalidSelectedCategoriesContent'));
        end        
    end
    
    %----------------------------------------------------------------------
    function val = validateStringVect(~,val,propName)
        errorFlag = false;
        if ischar(val) || iscellstr(val) || isstring(val)
            val = string(val);
            if ~isvector(val)
                errorFlag = true;
            end
        else
            errorFlag = true;
        end
        
        if errorFlag
            error(message('signal:internal:segmentation:InvalidStrVect',string(propName)));
        end
        val = val(:);
    end
    
end % end methods

%--------------------------------------------------------------------------
% Display control, Copy
%--------------------------------------------------------------------------
methods (Access = protected)
    function propgrp = getPropertyGroups(obj)
        %getPropertyGroups Group properties in order for object display
        
        propList = obj.pPropertyList;
        if ~obj.pSampleRateSpecified
            propList(propList == 'SampleRate') = [];
        end
        if ~obj.SpecifySelectedCategories
            propList(propList == 'SelectedCategories') = [];
        end        
        propgrp = matlab.mixin.util.PropertyGroup(propList);
    end
end % end methods

end % End of class

