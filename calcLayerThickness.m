function [ data ] = calcLayerThickness( surfaces )
%calcLayerThickness Calulates retinal layer thickness from surface
%coordinates
%   Detailed explanation goes here

s=size(surfaces);

if ~all(s==[128 512 11])
    error('surface data is not expected size');
end


layers=nan(s(1),s(2),s(3)-1);

for iLayer = 2:s(3)
    %each laye thickness is the subtraction from the
    %previous surface
    layers(:,:,iLayer - 1) = surfaces(:,:,iLayer) ...
        - surfaces(:,:,iLayer - 1);
end

data=layers;
end

