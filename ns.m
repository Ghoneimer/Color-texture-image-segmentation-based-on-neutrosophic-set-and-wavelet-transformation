function [Entropy_ns,Alpha, BetaEnhT,BetaEnhI]= ns(G,Enl)
% Calculate local mean value of the image
h = ones(5,5)/25;
g_mean = imfilter(G,h,'same');
% figure;
% subplot(1,2,1),imagesc(G),title('Original');
% subplot(1,2,2),imagesc(g_mean),title('Mean');
% get g_mean max and min value
g_min=min(min(im2double(g_mean)));
g_max=max(max(im2double(g_mean)));
[rows,cols] = size(G);

% En=entropy(G);
% Enmax=-log2(1/(rows*cols));
% Enmin=0;
% AlphaMax=0.1;
% AlphaMin=0.01;
% Alpha=AlphaMin+((AlphaMax-AlphaMin)*(En-Enmin)/(Enmax-Enmin));
% Beta=1-Alpha;

g_min = g_min.*ones(rows,cols);
g_max = g_max.*ones(rows,cols);
%figure,imshow(gmax),title('maxg');
%figure,imshow(gmin),title('ming');
% Absolute Value
o=abs(G-g_mean);
%figure,imshow(o),title('dif');
% get absolute value max and min
mino=min(min(im2double(o)));
maxo=max(max(im2double(o)));
o_min = mino.*ones(rows,cols);
o_max = maxo.*ones(rows,cols);
%figure,imshow(omax),title('maxo');
%figure,imshow(omin),title('mino');
% get True , Flase and Indetermincy
T=(im2double(g_mean)-g_min)./(g_max-g_min);
I=(im2double(o)-o_min)./(o_max-o_min);
F=1-T;

% figure,
% subplot(1,3,1),imshow(T),title('T');
% subplot(1,3,2),imshow(I),title('I');
% subplot(1,3,3),imshow(F),title('F');
% En=entropy(T)+entropy(I)+entropy(F);
% get Entorpy
Enmax=-log2(1/(rows*cols));
Enmin=0;
AlphaMax=0.1;
AlphaMin=0.01;
Alpha=AlphaMin+((AlphaMax-AlphaMin)*(Enl-Enmin)/(Enmax-Enmin));
Beta=1-Alpha;
%Entropy_ns=entropy(I);
fprintf('Alpha= %f , Beta = %f\n', Alpha,Beta);
% T_alpha mean
meanT = imfilter(T,h,'same');
AlphaT=ones(rows,cols);
AlphaT(I<Alpha)=T(I<Alpha);
AlphaT(I>=Alpha)=meanT(I>=Alpha); %required

% F_alpha mean
meanF = imfilter(F,h,'same');
AlphaF=ones(rows,cols);
AlphaF(I<Alpha)=F(I<Alpha);
AlphaF(I>=Alpha)=meanF(I>=Alpha); % required

AlphameanT = imfilter(AlphaT,h,'same');
meanI=abs(AlphaT-AlphameanT);
minmeanI=min(min(im2double(meanI)));
maxmeanI=max(max(im2double(meanI)));
meanImin = minmeanI.*ones(rows,cols);
meanImax = maxmeanI.*ones(rows,cols);
AlphaI=(im2double(meanI)-meanImin)./(meanImax-meanImin); %requierd

% figure,
% subplot(1,3,1),imshow(AlphaT),title('Alpha mean T');
% subplot(1,3,2),imshow(AlphaI),title('Alpha mean I');
% subplot(1,3,3),imshow(AlphaF),title('Alpha mean F');


EnhT=ones(rows,cols);
EnhT(AlphaT<=0.5) = 2*EnhT(AlphaT<=0.5).^2;
EnhT(AlphaT>0.5) = 1-2*(1-EnhT(AlphaT>0.5)).^2;

EnhF=ones(rows,cols);
EnhF(AlphaF<=0.5) = 2*EnhF(AlphaF<=0.5).^2;
EnhF(AlphaF>0.5) = 1-2*(1-EnhF(AlphaF>0.5)).^2;

% figure,
% subplot(1,2,1),imshow(EnhT),title('Enhanced T');
% subplot(1,2,2),imshow(EnhF),title('Enhanced F');

BetaEnhT=ones(rows,cols);
BetaEnhT(AlphaI<Beta)=AlphaT(AlphaI<Beta);
BetaEnhT(AlphaI>=Beta)=EnhT(AlphaI>=Beta); %required

BetaEnhF=ones(rows,cols);
BetaEnhF(AlphaI<Beta)=AlphaF(AlphaI<Beta);
BetaEnhF(AlphaI>=Beta)=EnhF(AlphaI>=Beta); %required

BetaT = imfilter(BetaEnhT,h);
BetaI=abs(BetaEnhT-BetaT);
minbetaI=min(min(im2double(BetaI)));
maxbetaI=max(max(im2double(BetaI)));
betaImin = minbetaI.*ones(rows,cols);
betaImax = maxbetaI.*ones(rows,cols);
BetaEnhI=(im2double(BetaI)-betaImin)./(betaImax-betaImin); %required
Entropy_ns=entropy(BetaEnhI);
% figure,
% subplot(1,3,1),imshow(BetaEnhT),title('Beta Enhanced T');
% subplot(1,3,2),imagesc(BetaEnhI),title('Beta Enhanced I');
% subplot(1,3,3),imshow(BetaEnhF),title('Beta Enhanced F');

% thres = graythresh(T);
% err=1;
% while err>0.0001,
%     mu1 = mean(T(T<=thres));
%     mu2 = mean(T(T>thres));
%     thres2 = (mu1+mu2)/2;
%     err = thres2-thres;
%     thres = thres2;
% end
% tt = thres;
% BinT = im2bw(T,thres);
% thres = graythresh(F);
% err=1;
% while err>0.0001,
%     mu1 = mean(F(F<=thres));
%     mu2 = mean(F(F>thres));
%     thres2 = (mu1+mu2)/2;
%     err = thres2-thres;
%     thres = thres2;
% end
% tf = thres;
% BinF = im2bw(F,thres);
% 
% [Gx, Gy] = gradient(im2double(A));
% eg = sqrt(Gx.*Gx+Gy.*Gy);
% eg=  eg/max(max(eg));
% sd = colfilt(im2double(A),[7 7],'sliding',@std);
% sd = sd/max(max(sd));
% 
% Homogen = 1 - sd.*eg;
% Indet = 1-Homogen;
% 
% lambda = 0.01;
% IO = zeros(size(T));
% IO((T>=tt)&(Indet<lambda)) = 1;
% IE = zeros(size(T));
% IE(((T<tt)|(F<tf))&(Indet>=lambda))=1;
% IB = zeros(size(F));
% IB((F>=tf)&(Indet<lambda))=1;
% 
% BinImage = ones(size(T));
% BinImage((IO==1)|(IB==1)|(IE==0))=0;
% 
% figure,imshow(BinImage),title('Binary Image based on T,I,F');
% 
% BinImage = 1-BinImage;
% cc = bwconncomp(BinImage);
% 
% numr = cc.NumObjects;
% 
% [RA, GA, BA] = deal(Orig);
% SegR = Orig;
% cumlab = zeros(size(T));
% for i=1:numr,
%     
%     labels=zeros(size(T));
%     labels(cc.PixelIdxList{i})=1;
%     labels=bwperim(imerode(labels,strel('disk',3)));
%     if sum(sum(labels))/numel(labels)>0.01,
%         cumlab = cumlab+labels;
%     end
%     
% end
% RA(cumlab>0)=1.0;
% GA(cumlab>0)=0.0;
% BA(cumlab>0)=0.0;
% colorImage = cat(3,RA,GA,BA);
% figure,imshow(colorImage),title('Segmentation Result');
% SegR(cumlab>0)=1.0;
% imwrite(SegR,'segmented.jpg','jpg');
% figure,imshow(SegR),title('Segmentation Result');
% 

%Entns1=entropy(T)+entropy(I)+entropy(F);
%Entns2=entropy(AlphaT)+entropy(AlphaI)+entropy(AlphaF);


end
