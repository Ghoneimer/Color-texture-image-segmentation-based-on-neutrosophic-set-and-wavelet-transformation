function LUV = RGB2ULV(RGB)

% Converting from RGB to Gray level
%Gray = rgb2gray(RGB);

% Converting from RGB to XYZ
XYZ_struct = makecform('srgb2xyz');
XYZ = applycform(RGB,XYZ_struct);

% Converting from XYZ to LUV
LUV_struct = makecform('xyz2uvl');
LUV = applycform(XYZ,LUV_struct);

% figure
% imshow(RGB);
% title('RGB');
% figure
% imshow(Gray);
% title('Gray');
% figure
% imshow(XYZ);
% title('XYZ');
% figure
% imshow(LUV);
% title('LUV');


% Seperating L U V

% figure
% subplot(2,2,1), imagesc(LUV(:,:,1)) , title('L');
% subplot(2,2,2), imagesc(LUV(:,:,2)) , title('U');
% subplot(2,2,3), imagesc(LUV(:,:,3)) , title('V');
end
