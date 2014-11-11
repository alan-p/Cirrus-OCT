function outArray = reshapeArray(inArray,Samples)
    %Reshape an NxM array to [N/S x M/S] x S^2
    [n,m] = size(inArray);
    s=Samples;
    cutsN = round(linspace(1,n+1,Samples+1));
    cutsM = round(linspace(1,m+1,Samples+1));
    %problem with rounding to whole pixels, not all cuts are exactly the
    %same size
    n = max(diff(cutsN));
    m = max(diff(cutsM));
    outArray = NaN((ceil(n/s) * ceil(m/s)), s^2);
    %sprintf('Samples:%i',Samples)
    for iRow = 1:Samples
        for iCol = 1:Samples
            idx = ((iRow - 1) * Samples) + iCol;
            sample = inArray(cutsN(iRow):(cutsN(iRow+1)-1),...
                                      cutsM(iCol):(cutsM(iCol+1)-1));
            outArray(1:numel(sample),idx)=sample(:);
        end
    end
end
        
    