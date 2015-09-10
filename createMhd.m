function createMhd(varargin)
% create .raw and .mhd files from a folder of images
%PARAMS:
%   path - full path to input folder
%   outname - output filename (no extension)
%   extension - file extension to search for, defaults to .tif
%   spacing - [x y z] array defining pixel size in m,
%       defaults to 6mm cirrus macular cube
%USAGE:
%   createMhd('path','c:\images\','outname','mydata')

    p=inputParser;
    path='';
    outname='';
    fileext='';
    elementspacing=[0.011719 0.001953 0.046875];
    
    addParamValue(p,'path',path,@ischar);
    addParamValue(p,'outname',outname,@ischar);
    addParamValue(p,'extension',fileext,@ischar);
    addParamValue(p,'spacing',elementspacing,@isnumeric);
    
    parse(p,varargin{:});
    
    path=p.Results.path;
    outname=p.Results.outname;
    
    if isempty(fileext)
        fileext='tif';
    end
    
    if length(elementspacing)~=3
        error('createMhd:InvalidParam','Incorrect element spacing should be [x y z] in m');
    end
    %select the source folder if not specified
    if isempty(path)
        path=uigetdir('/home/tom/Documents/temp/','Select image folder...');
    end
    
    %find the files in the folder
    %try to extract the file index from the filename
    files=dir(fullfile(path,['*.',fileext]));
    file_index=cell(length(files),2);
    for index=1:length(files)
        matches=regexp(files(index).name,['(?<fname>\w*)_(?<index>\d+).',fileext],'names');
        if isempty(outname)
            outname=matches.fname;
        end
        file_index{index,1}=matches.index;
        file_index{index,2}=matches.fname;
    end
    %sort the files into numeric order
    file_index=sortrows(file_index,1);
    
    %load the first file so we can preallocate the output matrix
    image=imread(fullfile(path,[file_index{1,2},'_',file_index{1,1},'.',fileext]));
    image_height = size(image,1);
    image_width = size(image,2);
    
    if ndims(image) > 2
        warning('Image is not bw, using only 1st layer');
    end
    
    im_array = zeros(image_width,image_height,length(file_index),'uint8');
    %loop through the rest of the images
    im_array(:,:,1)=image(:,:,1)';
    for index = 2:length(files)
        image=imread(fullfile(path,[file_index{index,2},'_',file_index{index,1},'.',fileext]));
        im_array(:,:,index)=image(:,:,1)';
    end
    
    %open the output files
    fid_mhd = fopen(fullfile(path,[outname,'.mhd']),'w');
    fid_raw = fopen(fullfile(path,[outname,'.raw']),'w');
    
    %write the mhd file
    fprintf(fid_mhd,'ObjectType = Image\n');
    fprintf(fid_mhd,'NDims = 3\n');
    fprintf(fid_mhd,'BinaryData = True\n');
    fprintf(fid_mhd,'BinaryDataByteOrderMSB = False\n');
    fprintf(fid_mhd,'CompressedData = False\n');
    fprintf(fid_mhd,'TransformMatrix = 1 0 0 0 0 1 0 1 0\n');
    fprintf(fid_mhd,'Offset = 0 0 0\n');
    fprintf(fid_mhd,'CenterOfRotation = 0 0 0\n');
    fprintf(fid_mhd,'AnatomicalOrientation = LAS\n');
    fprintf(fid_mhd,'ElementSpacing = %7.6f %7.6f %7.6f\n',elementspacing);
    fprintf(fid_mhd,'DimSize = %d %d %d\n',...
        image_width,image_height,length(files));
    fprintf(fid_mhd,'ElementType = MET_UCHAR\n');
    fprintf(fid_mhd,'ElementDataFile = %s\n',[outname,'.raw']);
    
    fclose(fid_mhd);
    
    count=fwrite(fid_raw,im_array,'uint8');
    
    fclose(fid_raw);
    