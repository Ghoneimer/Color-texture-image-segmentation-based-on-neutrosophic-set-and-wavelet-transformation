function MELH_HL = wavlet_decomposition(X)
%[Lo_D,Hi_D,Lo_R,Hi_R] = biorfilt(DF,RF) computes four filters associated 
%with the biorthogonal wavelet specified by decomposition filter DF 
%and reconstruction filter RF. These filters are
[Rf,Df] = biorwavf('bior2.2');
%It is well known in the subband filtering community that if the same FIR 
%filters are used for reconstruction and decomposition, then symmetry and 
%exact reconstruction are incompatible (except with the Haar wavelet). 
%Therefore, with biorthogonal filters, two wavelets are introduced instead 
%of just one
[LoD,HiD,LoR,HiR] = biorfilt(Df,Rf);
% plot
% figure;
% subplot(211); stem(LoD);
% title('Dec. low-pass filter bior2.2');
% subplot(212); stem(HiD);
% title('Dec. high-pass filter bior2.2');
%read image

HA = conv2(X,LoD(:)','same');
HD = conv2(X,HiD(:)','same');
V_LH1 = conv2(HA',HiD(:)','same');
V_LH1 =V_LH1';
H_HL1 = conv2(HD',LoD(:)','same');
H_HL1 = H_HL1';



%wcodemat rescales an input matrix to a specified range for display. 
%If the specified range is the full range of the current colormap
V1img = wcodemat(V_LH1,255,'mat',1);
H1img = wcodemat(H_HL1,255,'mat',1);
% figure;
% imagesc(V1img);
% title('Vertical(LH subband) detail Coef. of Level 1');
% figure;
% imagesc(H1img);
% title('Horizontal(HL subband) detail Coef. of Level 1');

filter_E=1/25*ones(5,5);
MELH=conv2(V_LH1,filter_E,'same');
MEHL=conv2(H_HL1,filter_E,'same');
figure;
imagesc(MELH);
title('MELH');
figure;
imagesc(MEHL);
title('MEHL');
MELH_HL(:,:,1)=MELH;
MELH_HL(:,:,2)=MEHL;
end
