close all
clc
clear
img=imread('image.png');
img=rgb2gray(img);
img=im2double(img);
LUV=zeros(size(img));
h = ones(5,5)/25;
MELH=zeros(size(img,1),size(img,2));
MEHL=zeros(size(img,1),size(img,2));
LUV = RGB2ULV(img);
img_gray=rgb2gray(img);
MELH_HL=wavlet_decomposition(img_gray) ;
Enl=0.00001;
err=0.0001;
Entropy_ns=0;
alpha=0.5;
X1=[];
%%%%%%%%%TL%%%%%%%%%%
while 1
[Entropy_ns,alpha1, TL,IL] = ns(LUV(:,:,1),Enl);
fprintf('Enl= %f , Entropy_ns= %f , error= %f\n', Enl,Entropy_ns,(Entropy_ns-Enl)/Enl);
if ((Entropy_ns-Enl)/Enl)< err
    break
else
    Enl=Entropy_ns ;
end
end
disp('---------------------');
X=zeros(size(img,1),size(img,2));
TL_mean = imfilter(TL,h,'same');
X(IL<alpha1)=TL(IL<alpha1);
X(IL>=alpha1)=TL_mean(IL>=alpha1); %required
imshow(X)
X1(:,:,1)=X;
%%%%%%%%%TU%%%%%%%%%%
Enl=0.00001;
while 1
[Entropy_ns,alpha2, TU,IU] = ns(LUV(:,:,2),Enl);
fprintf('Enl= %f , Entropy_ns= %f , error= %f\n', Enl,Entropy_ns,(Entropy_ns-Enl)/Enl);
if ((Entropy_ns-Enl)/Enl)< err
    break
else
    Enl=Entropy_ns ;
end
end
disp('---------------------');
X=zeros(size(img,1),size(img,2));
TU_mean = imfilter(TU,h,'same');
X(IU<alpha2)=TU(IU<alpha2);
X(IU>=alpha2)=TU_mean(IU>=alpha2); %required
imshow(X)
X1(:,:,2)=X;
%%%%%%%%%%TV%%%%%%%%%
Enl=0.00001;
while 1
[Entropy_ns,alpha3, TV,IV] = ns(LUV(:,:,3),Enl);
fprintf('Enl= %f , Entropy_ns= %f , error= %f\n', Enl,Entropy_ns,(Entropy_ns-Enl)/Enl);
if ((Entropy_ns-Enl)/Enl)< err
    break
else
    Enl=Entropy_ns ;
end
end
disp('---------------------');
X=zeros(size(img,1),size(img,2));
TV_mean = imfilter(TV,h,'same');
X(IV<alpha3)=TV(IV<alpha3);
X(IV>=alpha3)=TV_mean(IV>=alpha3); %required

imshow(X)
X1(:,:,3)=X;
%%%%%%%%%TLH%%%%%%%%%%
Enl=0.00001;
while 1
[Entropy_ns,alpha4, TLH,ILH] = ns(MELH_HL(:,:,1),Enl);
fprintf('Enl= %f , Entropy_ns= %f , error= %f\n', Enl,Entropy_ns,(Entropy_ns-Enl)/Enl);
if ((Entropy_ns-Enl)/Enl)< err
    break
else
    Enl=Entropy_ns ;
end
end
disp('---------------------');
X=zeros(size(img,1),size(img,2));
TLH_mean = imfilter(TLH,h,'same');
X(ILH<alpha4)=TLH(ILH<alpha4);
X(ILH>=alpha4)=TLH_mean(ILH>=alpha4); %required

imshow(X)
X1(:,:,4)=X;
%%%%%%%%%THL%%%%%%%%%%
Enl=0.00001;
while 1
[Entropy_ns,alpha5, THL,IHL] = ns(MELH_HL(:,:,2),Enl);
fprintf('Enl= %f , Entropy_ns= %f , error= %f\n', Enl,Entropy_ns,(Entropy_ns-Enl)/Enl);
if ((Entropy_ns-Enl)/Enl)< err
    break
else
    Enl=Entropy_ns ;
end
end
disp('---------------------');
X=zeros(size(img,1),size(img,2));
THL_mean = imfilter(THL,h,'same');
X(IHL<alpha5)=THL(IHL<alpha5);
X(IHL>=alpha5)=THL_mean(IHL>=alpha5); %required
imshow(X)
X1(:,:,5)=X;
%%%%%%%%%%%%%%%%%%%
nrows = size(X1,1);
ncols = size(X1,2);
X1 = reshape(X1,nrows*ncols,5);


nColors = 15;
% repeat the clustering 3 times to avoid local minima
[cluster_idx, cluster_center] = kmeans(X1,nColors,'distance','sqEuclidean', ...
                                        'Replicates',3);
pixel_labels = reshape(cluster_idx,nrows,ncols);
imshow(pixel_labels,[]), title('image labeled by cluster index');
