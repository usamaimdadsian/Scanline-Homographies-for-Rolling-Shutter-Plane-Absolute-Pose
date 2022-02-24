classdef HorizontalLineExtraction<handle
    
    properties(Access = public)
        c = 0.05; % inlier threshold for m-estimator
        minPoints = 4; % min. no. of feats. that needs to be present on a scanline for it to be considered admissible
        scanline_ids = []; % output
    end    
    
    properties(Access = private)
        ImgPoints = struct ('Data', [], 'Vis', []); % points of 2D features and the associated visibility
        SlidingWindow = cell(0) % empty stub for sliding window
        ImgPointsVis = []; % finally visible image features, derived from the values set in 'ImgPoints'
        numerSlide = 0; % number of fixed-width sliding windows
        m_estimator_option = 0; % 0 -> Huber; Future extension: 1 -> Cauchy; 2 -> Tukey
        prevX; % temporary variable for storing last known parameter (for IRWLS)
        lastWeights; % last known weights from m-estimator
        inliers = {}; % accepted inliers
        rejection_identifier; %  a flag variable indicating the elements/indices to be rejected, arbitrarily chosen very large number
    end    
    
    methods(Access = public)
        function [this] = HorizontalLineExtraction(num_sliding_window, reject_tag, inlier_threshold)
            % HorizontalLineExtraction    this class extracts horizontal
            % scanlines across the height of an image, using a fixed width
            % window. A scanline is intialised from the center of the
            % window (i.e., mean(Y_i)) and minimises the dist. to adjoining
            % feat. pts. using an m-estimator: 
            % $$\min_y \sum_i (y - \rho(Y_i) )^2$$
            %    where \rho() is the robust weight function.
            %
            %This is the constructor of HorizontalLineExtraction class.
            %
            % PARAMETERS:
            %    num_sliding_window: number of sliding windows on image
            %                                           (statically assigned for now)
            %    reject_tag:  a flag variable indicating the elements/indices to be rejected, arbitrarily chosen very large number
            %    inlier_threshold: threshold for inlier from m-estimator,
            %                                   typically 10^-3, increase to accept more points about the
            %                                   scanline, decrease to be more precise
            %
            % RETURNS:
            %    object of this class
            this.numerSlide = num_sliding_window; % number of sliding windows on image, statically assigned for now
            this.rejection_identifier = reject_tag; % a flag variable indicating the elements/indices to be rejected, arbitrarily chosen very large number
            this.c = inlier_threshold; % threshold for inlier from m-estimator
        end
        
        function [scanline_id] = extract_line(this, data, visibility)
            % extract_line    extracting scanlines from data using IRWLS
            % across fixed-width windows along the height of the image
            %
            % PARAMETERS:
            %    data: [Nx2] 3D points
            %    visibility: [Nx1] bi-valued array of 0s and 1s indicating
            %                     if a feature is visible in the current view or not
            % 
            % RETURNS:
            %    <scanline_id> containing the scanline_index (i.e., the y
            %    coordinate) for each feature. unique(scanline_id) will
            %    return the distinct and admissible scanlines
            this.ImgPoints.Data = data; % setting data
            this.ImgPoints.Vis = visibility; % setting visibility
            this.ImgPoints.Data = this.ImgPoints.Data .* this.ImgPoints.Vis; % retaining only visible feat. pts.
            
            [sliding_window] = this.get_window_partition(this.ImgPoints.Data); % partitioning img. height into approrpiate windows (non-overlapping)
            [scanline] = get_scan_line(this); % extracting scanline, IRWLS
            scanline_id = this.extract_scanline_ids(scanline, this.ImgPoints.Data, sliding_window); % re-arranging the scanlines according to feat. pts. dim.
            this.scanline_ids = scanline_id;
        end
    end
    
    methods(Access = private)
        function [sliding_window] = get_window_partition(this, data)
            % get_window_partition    partitioning height along fixed width
            % windows
            %
            % PARAMETERS:
            %    data: [Nx2] 2D data points
            %
            % RETURNS:
            %    sliding_window: window limits, as a cell
            minY = min(data(:,2));
            maxY = max(data(:,2));
            sliding_window = linspace(minY,maxY,this.numerSlide+1);
            this.SlidingWindow = cell(this.numerSlide,2);
            for index = 1:size(data,1)
                y = data(index,2);
                for wIndex = 1:(this.numerSlide)
                    if( (y>=sliding_window(wIndex)) && (y<sliding_window(wIndex+1)) )
                        this.SlidingWindow{wIndex} = vertcat(this.SlidingWindow{wIndex},data(index,:));
                    end
                end
            end
        end
        
        function [scanline_ids] = extract_scanline_ids(this, scanline, data, sliding_window)
            % extract_scanline_ids    post-processing the scanline
            % positions obtained from IRWLS
            % 
            % PARAMETERS:
            %    scanline: estimated scanline positions
            %    data: [Nx2] 2D points
            %    sliding_window: window limits, as a cell
            %
            % RETURNS
            %    scanline_ids: [Nx1] scanline ids, containing the scanline_index (i.e., the y
            %    coordinate) for each feature. unique(scanline_id) will
            %    return the distinct and admissible scanlines
            rIndex = 1;
            N = size(data,1);
            scanline_ids = ones(N,1)*this.rejection_identifier;
            for index = 1:size(data,1)
                y = data(index,2);
                for wIndex = 1:(this.numerSlide)
                    if( (y>=sliding_window(wIndex)) && (y<sliding_window(wIndex+1)) )
                        line_number = scanline(wIndex);
                        if(line_number > this.rejection_identifier)
                            scanline_ids(rIndex,1) = line_number;
                        end
                        rIndex = rIndex + 1;
                    end
                end
            end
            
            for wIndex = 1:(this.numerSlide)
                line_number = scanline(wIndex);
                tIndex = 1;
                for index = 1:size(data,1)
                    if(scanline_ids(index) == line_number)
                        if(~ismember(tIndex,this.inliers{wIndex}))
                            scanline_ids(index) = this.rejection_identifier;
                        end
                        tIndex = tIndex + 1;
                    end
                end
            end
        end
        
        function [scanline] = get_scan_line(this)
            % get_scan_line    extract the actual scanlines by minimising a
            % robust cost
            %
            % RETURNS:
            %    estimated scanline positions
            scanline = ones(this.numerSlide,1)*this.rejection_identifier;
            this.inliers = cell(this.numerSlide,1);
            for wIndex = 1:this.numerSlide
                points = this.SlidingWindow{wIndex};
                if(size(points,1) >= this.minPoints)
                    Y = points(:,2);
                    x0 = mean(Y);
                    this.prevX = x0;
                    fun = @(x)cost_function(this,x, Y);
                    output = @(x,optimValues,state)out_fun(this,x,optimValues,state);
                    options = optimoptions('fminunc','SpecifyObjectiveGradient',false,'OutputFcn',output,...
                        'OptimalityTolerance',10.^-50,'StepTolerance',10^-50,'FiniteDifferenceType', 'central');
                    [x, fval] = fminunc(fun,x0,options);
                    scanline(wIndex) = x;
                    this.inliers{wIndex} = find(this.lastWeights<this.c).';
                end
            end
        end
        
        function [E] = cost_function(this, X, Y)
            % cost_function    optimisation callback function for computing
            % robust cost
            E = 0;
            R = Y-X;
            Rprev = Y - this.prevX;
            this.lastWeights = zeros(size(R,1),1);
            if(this.m_estimator_option == 0)
                % Huber
                k = 1.345*std(Rprev);
                for ii = 1:size(R,1)
                    r = R(ii);
                    rprev = Rprev(ii);
                    if(abs(rprev)<k)
                        w = 1;
                    else
                        w = k/abs(rprev);
                    end
                    this.lastWeights(ii,1) = abs(r);
                    % Weighted residual
                    E = E + (w*r*r);
                end
            end
        end
        
        function [stop] = out_fun(this, x,optimValues,state)
            % out_fun    optimisation callback function for dumping
            % variables at the end of each optim. iteration
            this.prevX = x;
            stop = 0;
        end
        
    end
end