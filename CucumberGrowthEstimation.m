%% Estimating the daily growth of the main stem of cucumber
% Author: Michna Mariusz 
% Date: 01/29/2018

close all
clear all
clc

%% Preparing the data
%
NumberOfFrames = 89;
ProjectNumber = 27;

Q = 'Proj';
ProjectNumber = string(ProjectNumber);
ProjectNumber = char(ProjectNumber);
E1 = '_img0000000';
E2 = '_img000000';
E3 = '_img00000';
s = (1:NumberOfFrames)';
s = string(s);
Y = '.jpg';


 PlantTip_DATA = [];
 Px_p_mm_DATA = [];
 Px_p_mm = 0;
 PlantClip_DATA = [];
 Dc = [];
 Day = 1;
 Data = [];
 


for u = 1:NumberOfFrames; 
    E = E3;
    if u < 10
       E = E1; 
    end
    if u < 100 && u >= 10
       E = E2; 
    end
    
X = char(s(u,1));
Z = [Q ProjectNumber E X Y];
P = imread(Z);
P = im2double(P);

GeneralBrightness = sum(P);
GeneralBrightness = sum(GeneralBrightness);
GeneralBrightness = sum(GeneralBrightness);

if GeneralBrightness > 400000
is_dark = 0;
height = size(P,1);
width = size(P,2);

%% 
% b) HSV

P_hsv = rgb2hsv(P);

%% 
% c) HSL

P_hsl = hsv2hsl(P_hsv);


%% 
% d) YCbCr

P_YCbCr = rgb2ycbcr(P);


%% 
% e) CIELAB

P_lab = rgb2lab(P);


%% 
% f) CIELUV

P_luv = luv(P);
  
%% Determination of the proportionality factor
% Reference object segmentation:
%       Segmentation based on the created tree model which had been created using 'fitctree' commend

P_RefObj = zeros(height,width);
for i = 1:height;
for j = 1:width;
     
     if P_hsv(i,j,1) < 0.836667
         if P_lab(i,j,2) < 21.6311
            if P_hsv(i,j,1) < 0.811625
                if P_hsv(i,j,1) < 0.722222
                   if P_hsv(i,j,1) >= 0.69246 && P_hsv(i,j,2) >= 0.301051
                       P_RefObj(i,j) = 1;
                   end
                else
                    if P_YCbCr(i,j,2) >= 0.523147
                       P_RefObj(i,j) = 1; 
                    end
                end
            else
               if P_YCbCr(i,j,2) >= 0.513092
                   P_RefObj(i,j) = 1;
               end
            end
         else
            if P(i,j,3) >= 0.194118 && P_hsv(i,j,1) < 0.0382555
                P_RefObj(i,j) = 1;
            end
         end
     else
         if P_lab(i,j,2) < 7.6964
            if P_hsv(i,j,1) < 0.852083
                P_RefObj(i,j) = 1;
            end
         else
             P_RefObj(i,j) = 1;
         end
     end
     
end
end


%% 
% Morphological filtration

OpenMask = zeros(20,1);
OpenMask(:,:) = 1;
P_open = imopen(P_RefObj,OpenMask);

i = 1;
j = 1;
Flag = 0;

 while i <= height && Flag == 0
    while j <= width && Flag == 0     
     if P_open(i,j) == 1 
      x_up_RefObj = i;
      y_up_RefObj = j;
      Flag = 1;
     end
     j = j +1;
    end
 j = 1;
 i = i+1;
 end
  
  
i = height;
j = 1;
Flag = 0;

 while i > 0 && Flag == 0
    while j <= width && Flag == 0     
     if P_open(i,j) == 1 
      x_down_RefObj = i;
      y_down_RefObj = j;
      Flag = 1;
     end
     j = j +1;
    end
 j = 1;
 i = i-1;
 end
  
  ReferenceObjectLengthPx = sqrt((x_down_RefObj - x_up_RefObj)^2+(y_down_RefObj - y_up_RefObj)^2); 
  Px_p_mm = 120/ReferenceObjectLengthPx;
  Px_p_mm_DATA = [Px_p_mm_DATA
      Px_p_mm ]; 
  
%% 
% Additional filtering mask determined based on the position of the pattern
  
MaskWidth = 20; % 
k = -2:0.001:10;
k = k';

RefObj_Direction = [(x_down_RefObj - x_up_RefObj) (y_down_RefObj - y_up_RefObj)]; 
  
U = int16(k*RefObj_Direction);
U2 = [];

for i = 1:size(U,1)
RefObj_Direction = int16([x_up_RefObj y_up_RefObj]);
U(i,:) = U(i,:) + RefObj_Direction;
    if U(i,1) >= 1 && U(i,1) <= height
        U2 = [U2
              U(i,:)];
    end
end
U = U2;
FilterMask = zeros(height,width);   
  

for i = 1:size(U,1)
    if U(i,2)-MaskWidth > 0 && U(i,2)+ MaskWidth < width
        FilterMask(U(i,1),(U(i,2)- MaskWidth):(U(i,2)+ MaskWidth)) = 1;
    end 
end
 

%% The tip the plant detection
% Segmentation based on the created tree model which had been created using 'fitctree' commend
P_PlantTip = zeros(height,width); 
 for i = 1:height;
 for j = 1:width;
     
  if P_lab(i,j,2) < -10.3835;
   
    if P_hsv(i,j,1) < 0.35;
       P_PlantTip(i,j) = 1;
    end
    
else
    if P_lab(i,j,2) < -7.48836;
        if P_luv(i,j,2) < -5.04106;
            if P_hsv(i,j,1) < 0.248779;
                P_PlantTip(i,j) = 1;
            end
        else
            if P_hsv(i,j,2) < 0.321148
                if P_lab(i,j,2) < -7.9802;
                    P_PlantTip(i,j) = 1;
                end
            else
                if P_hsl(i,j,2) < 0.262321;
                    if P_luv(i,j,3) >= 10.2297;
                        P_PlantTip(i,j) = 1;
                    end
                else
                    P_PlantTip(i,j) = 1;
                end
            end
        end
    else
        if P_lab(i,j,2) < -6.34943;
            if P(i,j,3) < 0.152941;
                if P_hsv(i,j,1) < 0.354167
                    P_PlantTip(i,j) = 1;
                end
            else
                if P_hsv(i,j,1) < 0.185993;
                    P_PlantTip(i,j) = 1;
                end
            end
        else
            if P(i,j,1) < 0.0215686;
                if P_YCbCr(i,j,2) < 0.491099 
                    P_PlantTip(i,j) = 1;
                else
                    if P_luv(i,j,3) >= 0.942955 && P(i,j,1) < 0.00588235 && P_YCbCr(i,j,2) < 0.493682
                        P_PlantTip(i,j) = 1;
                    end
                end
                
            else
                if P_luv(i,j,3) < 35.6866
                    if P_lab(i,j,2) <-5.2794
                        if P_luv(i,j,2) < -2.3968
                           if P_lab(i,j,2) < -5.28035
                               if P_luv(i,j,2) >= -2.56477 && P_lab(i,j,3) < 8.49256 && P_lab(i,j,2) < -5.61972
                                   P_PlantTip(i,j) = 1;
                               end
                           else
                               P_PlantTip(i,j) = 1;
                           end
                        else
                            if P_hsv(i,j,2) >= 0.293354
                                P_PlantTip(i,j) = 1;
                            end
                        end
                    else
                        if P_luv(i,j,3) < 32.2745
                            if P(i,j,1) < 0.0372549
                                if P_lab(i,j,3) >= 3.64229
                                   P_PlantTip(i,j) = 1; 
                                end
                            else
                                if P_hsl(i,j,3) < 0.230392
                                   if P_YCbCr(i,j,2) < 0.421031
                                        P_PlantTip(i,j) = 1;
                                   else
                                       if P_luv(i,j,3) >= 6.67508
                                          if P_lab(i,j,2) < -3.18724
                                             if P_luv(i,j,3) >= 9.78384
                                                P_PlantTip(i,j) = 1; 
                                             end
                                          else
                                              if P_luv(i,j,3) < 6.68813
                                                  P_PlantTip(i,j) = 1;
                                              else
                                                  if P(i,j,3) < 0.201961 && P_luv(i,j,3) < 6.75471 && P(i,j,1) >= 0.198038
                                                      P_PlantTip(i,j) = 1;
                                                  end
                                              end
                                          end
                                       end
                                   end
                                else
                                    if P_luv(i,j,3) >=25.6809 && P(i,j,1) < 0.445098
                                        P_PlantTip(i,j) = 1;
                                    end
                                end
                            end
                        else
                            if P_hsv(i,j,1) >= 0.122156
                                P_PlantTip(i,j) = 1;
                            end
                        end
                        
                        
                    end
                else
                    if P_hsv(i,j,1) >= 0.124028
                        P_PlantTip(i,j) = 1;
                    end
                end
            end
        end
    end
    
end
  
 end



    end

% Morphological filtration
P_open = imopen(P_PlantTip,OpenMask);
P_open = P_open.*FilterMask;

i = 1;
j = 1;
Flag = 0;

  while i <= height && Flag == 0
    while j <= width && Flag == 0     
     if P_open(i,j) == 1 
      x_PlantTip = i;
      y_PlantTip = j;
      Flag = 1;
     end
     j = j +1;
  end
j = 1;
i = i+1;
  end
 
PlantTip_DATA = [PlantTip_DATA
            x_PlantTip y_PlantTip];
  

%% The lowest clip detection
% Segmentation based on the created tree model which had been created using 'fitctree' commend

P_PlantClip = zeros(height,width);

for i = 1:height;
 for j = 1:width; 
    if P_hsv(i,j,2) < 0.0906158
        if P(i,j,2)>= 0.388235
            P_PlantClip(i,j) = 1;
        end
    else
        if P_hsv(i,j,2) < 0.121132
            if P_luv(i,j,3) < 8.81528
                if P_YCbCr(i,j,3) < 0.485136
                    P_PlantClip(i,j) = 1;
                end
            else
                if P_YCbCr(i,j,2) < 0.479601
                    if P(i,j,1)>= 0.511765 && P_lab(i,j,2) >= -6.65514
                        P_PlantClip(i,j) = 1;
                    end
                else
                    P_PlantClip(i,j) = 1;
                end
            end
        else
            if P(i,j,3)< 0.3
               if P(i,j,3)< 0.256863 && P(i,j,3)>= 0.213725 && P_hsv(i,j,3) < 0.327451 && P_YCbCr(i,j,2) < 0.460956 
                   P_PlantClip(i,j) = 1;
               end
            else
                if P_lab(i,j,2) < -4.90842
                   if P_hsv(i,j,2) < 0.222497
                       if P_YCbCr(i,j,1) >= 0.362426
                           if P_hsv(i,j,1) < 0.234708
                              if P(i,j,1)< 0.65098
                                  if P_lab(i,j,2) >= -5.1376
                                      P_PlantClip(i,j) = 1;
                                  end
                              else
                                  P_PlantClip(i,j) = 1;
                              end
                           else
                               if P_YCbCr(i,j,3) < 0.487475
                                   if P(i,j,1)< 0.337255
                                       P_PlantClip(i,j) = 1;
                                   end
                               else
                                   P_PlantClip(i,j) = 1;
                               end
                           end
                       end
                   else
                      if P(i,j,2) >= 0.413725 && P(i,j,2) < 0.433333 && P(i,j,1)>= 0.382353
                          P_PlantClip(i,j) = 1;
                      end
                   end
                end
            end
        end
    end
 end

end
  
P_PlantClip = P_PlantClip.*FilterMask;
    
OpenMask = zeros(20,1);
OpenMask(:,:) = 1;  
P_sum = P_PlantTip + P_PlantClip;

P_open_sum = imopen(P_sum,OpenMask);

OpenMask = zeros(1,8);
OpenMask(:,:) = 1;
P_open2 = imopen(P_open_sum,OpenMask);

P_open2 = 1 - P_open2;

P_open = P_open2.*P_open_sum.*P_PlantClip;

  i = height;
  j = width;
  Flag = 0;
 while i > 0 && Flag == 0
 while j > 0 && Flag == 0     
     if P_open(i,j) == 1 
      x_PlantClip = i;
      y_PlantClip = j;
      Flag = 1;
     end
     j = j - 1;
          end
j = width;
i = i - 1;
  end
PlantClip_DATA = [PlantClip_DATA   
       x_PlantClip y_PlantClip ];

   
else
    is_dark = is_dark + 1;
end
%%
% Is new day?
if is_dark == 3;
   Day = Day + 1; 
   is_dark = 12;

   Median_PlantTip = median(PlantTip_DATA(:,2),1);
   
   %%
   Median_RefObj = median(Px_p_mm_DATA,1);
   
   %%
   Median_PlantClip = median(PlantClip_DATA(:,2),1);
   
   for t = 1:size(PlantTip_DATA,1)
          if abs(norm(PlantClip_DATA(t,2)-[Median_PlantClip])) < 8 && abs(norm(PlantTip_DATA(t,2)-Median_PlantTip)) < 10 && abs(Px_p_mm_DATA(t,1)-Median_RefObj) < 0.2

            d = sqrt((PlantTip_DATA(t,1) - PlantClip_DATA(t,1))^2 + (PlantTip_DATA(t,2) - PlantClip_DATA(t,2))^2);
            d = Px_p_mm_DATA(t,1)*d; 
            Dc = [Dc
                  d];
          end   
   end
    
   M = median(Dc(:,1),1);
   m = [];
   for t = 1:size(Dc,1)
   
       if abs(Dc(t,1)-M) < 30
       m = [m
            Dc(t,1)];
       end 
   end
   avg = mean(m);
   Data = [Data
            avg]; 
   PlantTip_DATA = []; 
   Dc = [];
   PlantClip_DATA = [];
   Px_p_mm_DATA = [];

end

end

%% 
% Last Day
   if is_dark < 12
   Median_PlantTip = median(PlantTip_DATA,1);
   
   %%
   Median_RefObj = median(Px_p_mm_DATA(:,1),1);
   
   %%
   Median_PlantClip = median(PlantClip_DATA,1);
 
 %% 
 % Rejection of photos which the results differ significantly from the median 
   for t = 1:size(PlantTip_DATA,1) 
       if abs(norm(PlantClip_DATA(t,:)-Median_PlantClip)) < 8 && abs(norm(PlantTip_DATA(t,:)-Median_PlantTip)) < 10 && abs(Px_p_mm_DATA(t,1)-Median_RefObj) < 0.2

            d = sqrt((PlantTip_DATA(t,1) - PlantClip_DATA(t,1))^2 + (PlantTip_DATA(t,2) - PlantClip_DATA(t,2))^2);
            d = Px_p_mm_DATA(t,1)*d; 
            Dc = [Dc
                  d];
       end 
   end
   M = median(Dc(:,1),1);
   m = [];
   for t = 1:size(Dc,1)
       if abs(Dc(t,1)-M) < 20
       m = [m
            Dc(t,1)];
       end
   end
   avg = mean(m);
   Data = [Data
            avg]; 
   end
figure
imshow('Example_Tip.jpg')
figure
imshow('Example_Clip.jpg')
  %% Daily growth of the main stem of cucumber
  ManualMeasurements = [477; 511; 540; 568; 594; 617; 635; 637]; % [mm]
  Delta = [];
  for y = 1:(size(Data,1)-1)
      Delta = [Delta
                (Data(y+1,1) - Data(y,1))
                ];
  end
  
  Delta2 = [];
  for y = 1:(size(ManualMeasurements,1)-1)
      Delta2 = [Delta2
                (ManualMeasurements(y+1,1) - ManualMeasurements(y,1))
                ];
  end

figure, plot(Data,'xb')
title('Distance from the lowest clip to the tip of the plant')
ylabel('a [mm]')
axis([0 size(Data,1)+1 0 max(Data)*1.2])
hold on
grid on
xlabel('day')
plot(ManualMeasurements,'r+')
legend('Estimated','Measured')
  
  
figure, plot(Delta,'xb')
title('Daily growth of the main stem of the plant')
ylabel('delta [mm]')
xlabel('day')
axis([0 size(Delta,1)+1 0 max(Delta)*1.2])
grid on
hold on
plot(Delta2,'r+')
legend('Estimated','Measured')
function Image = hsv2hsl(Image)
   Max = Image(:,:,3); 
   Min = (1 - Image(:,:,2)).*Max; 
   L = 0.5*(Max + Min);
   
   if L <= 0.5
    if L == 0
       S = 0.5*(Max - Min);
    else
       S = (Max - Min)./(Max + Min);
    end
   else
    S = (Max - Min)./(2 - Max + Min);
   end
   Image(:,:,2) = S;
   Image(:,:,3) = L;
   return;
end




function Image = luv(Image)
% Convert to CIE L*u*v* (CIELUV)
B_ref = [0.950456,1,1.088754];
U_ref = (4*B_ref(1))./(B_ref(1) + 15*B_ref(2) + 3*B_ref(3));
V_ref = (9*B_ref(2))./(B_ref(1) + 15*B_ref(2) + 3*B_ref(3));

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
Y = Image(:,:,2)/B_ref(2);
L = 116*f(Y) - 16;

Image(:,:,1) = L;                  % L*
Image(:,:,2) = 13*L.*(U - U_ref);  % u*
Image(:,:,3) = 13*L.*(V - V_ref);  % v*
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