%% Estimating chlorophyll content in the tomato leaves
% Author: Michna Mariusz 
% Date: 01/29/2018

close all
clear all
clc

P = imread('InputImage.jpg');
NumberOfNeurons=7; 


P = im2double(P);
P = imresize(P,[480 640]);
X1 = [];
X2 = [];
height = size(P,1);
width = size(P,2);

% b) HSV

P_hsv = rgb2hsv(P);

% c) HSL

P_hsl =hsv2hsl(P_hsv);

% d) YCbCr

P_YCbCr = rgb2ycbcr(P);

% e) CIELAB

P_lab = rgb2lab(P);

% f) CIELUV

P_luv = luv(P);


%% Leaf segmentation:
% Segmentation based on the created tree model which had been created using 'fitctree' commend
P_seg = zeros(height,width);

 for i = 1:height
 for j = 1:width
     
         if P_luv(i,j,2) < -1.7214;
        if P(i,j,3) < 0.649506
            P_seg(i,j) = 1;
        end
     else
         if P_lab(i,j,2) < -6.38077
            if P(i,j,1) < 0.438254
             P_seg(i,j) = 1;
            end
         else
             if P(i,j,1)< 0.322352
                 if P_lab(i,j,2) < -2.73531
                    if P_YCbCr(i,j,3) < 0.508885
                        P_seg(i,j) = 1; 
                    end
                 else
                     if P_lab(i,j,2) < -1.64771 && P_hsv(i,j,1) < 0.138929
                      P_seg(i,j) = 1;
                     end
                 end
             end
         end
         end

 end
 end

% Additional filtration by morphological erosion

ErosionMask = zeros(25,25); 
ErosionMask(:,:) = 1;
P_seg = imerode(P_seg,ErosionMask);
S = sum(sum(P_seg));

% Average blue (B) coefficient from RGB space
B_sum = sum(sum(P_seg.*P(:,:,3)));
B_avg = B_sum/S; % œrednie B

% Average saturation (S) coefficient from RGB space
S_sum = sum(sum(P_seg.*P_hsv(:,:,2)));
S_avg = S_sum/S;

%% The sigmoid neural network model 

load('TrainingData.mat'); 

Xu = [D(:,3) D(:,5)];
Xt = [B_avg S_avg];
Yu = D(:,(22));
NumberOfTrainingData=size(Yu,1);

Xeu=[ones(NumberOfTrainingData,1) Xu];
Xet = [1 Xt];
e = [];
Hu=[];
H_t=[];

for i=1:NumberOfNeurons
    a=2*(rand(2+1,1)-0.5);
    Gu=tanh(Xeu*a);
    G_t=tanh(Xet*a);
    Hu = [Hu Gu]; 
    H_t = [H_t G_t];
    Weights=pinv(Hu)*Yu;
    F_score = H_t*Weights;
end 

F_score

function Image = hsv2hsl(Image)

   % przejœcie z HSV na HSL, poniewa¿ jest to prostrze w implementacji, ni¿
   
   Max = Image(:,:,3); % gdzie Max = max(R,G,B)
   Min = (1 - Image(:,:,2)).*Max; % gdzie Min = min(R,G,B)
   L = 0.5*(Max + Min);
   
   m = min(L,1-L); 
   
   Image(:,:,2) = 0.5*(Max - Min)./(m + (m == 0)); % S
   % Sprowadza siê to do tego, ¿e jeœli dla danego piksela L = 0 to 
   % S = (Max - Min)/2, a gdy L <= 0.5 wówszas 
   % S = (Max - Min)/(Max + Min), a gdy L > 0.5 
   % S = (Max - Min)/(2 - Max - Min)
   Image(:,:,3) = L;
   return;
end


function Image = luv(Image)
% Convert to CIE L*u*v* (CIELUV)
WhitePoint = [0.950456,1,1.088754];
WhitePointU = (4*WhitePoint(1))./(WhitePoint(1) + 15*WhitePoint(2) + 3*WhitePoint(3));
WhitePointV = (9*WhitePoint(2))./(WhitePoint(1) + 15*WhitePoint(2) + 3*WhitePoint(3));
   % Convert RGB to XYZ
 R = invgammacorrection(Image(:,:,1));
   G = invgammacorrection(Image(:,:,2));
   B = invgammacorrection(Image(:,:,3));

   T = inv([3.2406, -1.5372, -0.4986; -0.9689, 1.8758, 0.0415; 0.0557, -0.2040, 1.057]);
   Image(:,:,1) = T(1)*R + T(4)*G + T(7)*B;  % X 
   Image(:,:,2) = T(2)*R + T(5)*G + T(8)*B;  % Y
   Image(:,:,3) = T(3)*R + T(6)*G + T(9)*B;  % Z%
Denom = Image(:,:,1) + 15*Image(:,:,2) + 3*Image(:,:,3);
U = (4*Image(:,:,1))./(Denom + (Denom == 0));
V = (9*Image(:,:,2))./(Denom + (Denom == 0));
Y = Image(:,:,2)/WhitePoint(2);
L = 116*f(Y) - 16;

Image(:,:,1) = L;                        % L*
Image(:,:,2) = 13*L.*(U - WhitePointU);  % u*
Image(:,:,3) = 13*L.*(V - WhitePointV);  % v*
return;  
end


function R = invgammacorrection(Rp)
R = zeros(size(Rp));
i = (Rp <= 0.0404482362771076);
R(i) = Rp(i)/12.92;
R(~i) = real(((Rp(~i) + 0.055)/1.055).^2.4);
return;
end



function fY = f(Y)
fY = real(Y.^(1/3));
i = (Y < 0.008856);
fY(i) = Y(i)*(841/108) + (4/29);
return;
end