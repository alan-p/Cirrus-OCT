function outmat = processXmlSurfaceFile(filename)
    %Load an XML Surface file generated from the IOWA reference algorithms
    %Return a matrix of size ascans x bscans x layers
    %PARAMS:
    %   filename - optional filename of input xml
    if nargin < 1
        [filename,pathname] = uigetfile({'*.xml','*.*'},'Select xml file.');
        [~,filename,ext] = fileparts(filename);
    else
        [pathname,filename,ext] = fileparts(filename);
    end
    xmlDom = xmlread(fullfile(pathname,[filename,ext]));
    surfaces = xmlDom.getElementsByTagName('surface');
    nSurfaces = surfaces.getLength;
    if nSurfaces == 0
        %no surface data found, throw an error
        error('processXmlSurfaceFile:InvalidFileType',...
            'No surfaces found in file %s.',[filename, ext])
    end
    
    %check the number of bscans and ascans in the first surface
    bscans = surfaces.item(0).getElementsByTagName('bscan');
    nBscans=bscans.getLength;
    if nBscans == 0
       %no bscan data found, throw an error
        error('processXmlSurfaceFile:InvalidFileType',...
            'No bscans found for surface 0 in file %s.',[filename, ext])
    end 
    ascans = bscans.item(0).getElementsByTagName('z');
    nAscans = ascans.getLength;
    if nAscans == 0
       %no ascan data found, throw an error
        error('processXmlSurfaceFile:InvalidFileType',...
            'No ascans found for surface 0, bscan 0 in file %s.',[filename, ext])
    end 
    %preallocate the output matrix
    outmat = zeros(nAscans,nBscans,nSurfaces,'uint16');
    for iSurface=0:nSurfaces-1
        surface=surfaces.item(iSurface);
        surfIdx = str2num(surface.getElementsByTagName('label').item(0).getFirstChild.getData);
        bscans=surface.getElementsByTagName('bscan');
        for iBscan=0:nBscans-1
            bscan=bscans.item(iBscan);
            ascans=bscan.getElementsByTagName('z');
            for iAscan=0:nAscans-1
                outmat(iAscan+1,iBscan+1,iSurface+1)=str2num(ascans.item(iAscan).getFirstChild.getData);
            end
        end
    end
    