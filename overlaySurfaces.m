function bscan_image = overlaySurfaces(varargin)
    %Overlay surface segmentation onto a bscan image
    %PARAMS:
    %   layer_mat - a matrix of size #Bscans x #Ascans x #Surfaces such as
    %      that generated from processXmlSurfaceFile
    %   bscan_num - (integer) the number of the bscan in the cube set
    %   imagefile - (string) the path to the bscan image (optional)
    %OUTPUTS:
    %   bscan_image - (uint8) an NxMx3 matrix suitable for display with
    %   imshow
    
    %process the inputs
    p=inputParser;
    %setup some default values
    layer_mat = NaN;
    bscan_num = 1;
    imagefile = '';
    colormap = [166,206,227;...
                31,120,180;...
                178,223,138;...
                51,160,44;...
                251,154,153;...
                227,26,28;...
                253,191,111;...
                255,127,0;...
                202,178,214;...
                106,61,154;...
                255,255,153];
    addParamValue(p,'layer_mat',layer_mat,@isnumeric);
    addParamValue(p,'bscan_num',bscan_num,@isnumeric);
    addParamValue(p,'imagefile',imagefile,@ischar);
    
    parse(p,varargin{:});
    
    layer_mat = p.Results.layer_mat;
    bscan_num = p.Results.bscan_num;
    imagefile = p.Results.imagefile;
    
    if ndims(layer_mat) < 2 | isnan(layer_mat)
        error('overlaySurfaces:Params',...
            'Invalid layer matrix')
    end
    
    if size(layer_mat,3) > size(colormap,1)
       wqrning('overlaySurfaces:Params',...
            'A maximum of %i surfaces are supported, colours will be reused',size(colormap,1))
    end 
    
    if bscan_num > size(layer_mat,1) | bscan_num < 1
        error('overlaySurfaces:Params',...
            'Invalid bscan number supplied, number must be between 1 and %i', size(layer_mat,1))
    end
    
    if strcmp(imagefile,'')
        [filename,pathname] = uigetfile({'*.jpg;*.tif;*.png;*.gif','All Image Files';...
          '*.*','All Files' },'Select bscan image.');
        [~,filename,ext] = fileparts(filename);
    else
        [pathname,filename,ext] = fileparts(imagefile);
    end
    
    bscan_image = imread(fullfile(pathname,[filename,ext]));
    
    if size(bscan_image,2) ~= size(layer_mat,2)
        error('overlaySurfaces:Params',...
            'shape of layer matrix does not match bscan image')
    end
    
    if max(layer_mat(:)) > size(bscan_image,1)-1 | min(layer_mat(:)) < 0
        error('overlaySurfaces:Params',...
            'layer matrix contains invalid values')
    end
    
    if size(bscan_image,3) < 3
        %not a color image, convert to RGB
        %assuming bscan_image is type uint8
        %may need to fix this laterz
        bscan_image = repmat(bscan_image,[1 1 3]);
    end
    
    for iLayer=1:size(layer_mat,3)
        coloridx = mod(iLayer-1,size(colormap,1))+1;
        color=colormap(coloridx,:);
        
        ypix = layer_mat(bscan_num,:,iLayer)+1; %I think the output from Iowa reference algorithms is 0 indexed
        for xpix=1:length(ypix)
            bscan_image(ypix(xpix),xpix,:)=color;
        end
    end
    

    
    