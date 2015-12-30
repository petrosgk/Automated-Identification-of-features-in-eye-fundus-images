%Using MESSIDOR - Base11.

%Read image.
img = imread('img16.tif');
%Get red channel of image, which we'll use to locate the Optical Disk.
imgR = img(:,:,1);
%Get green channel ofimage, which we'll use to locate the Macula.
imgG = img(:,:,2);


%Apply mask.
imgMask = imread('mask.tif');
[m n] = size(imgMask);
for i = 1:m,
    for j = 1:n,
        if imgMask(i,j) == 255
           imgR(i,j) = NaN;
           imgG(i,j) = NaN;
        end
    end
end


%Morphological Closing to make vessels fader.
%Papers: http://www.ijarcsse.com/docs/papers/April2012/Volume_2_issue_4/V2I400117.pdf
%        http://rabida.uhu.es/dspace/bitstream/handle/10272/4446/detecting_the_optic_disc_boundary.pdf
SE = strel('disk',30,8);
imgR = imclose(imgR,SE);
imgG = imclose(imgG,SE);
figure('Name','Morphological Closing - Red Image','NumberTitle','off');
imshow(imgR);
figure('Name','Morphological Closing - Green Image','NumberTitle','off');
imshow(imgG);



%Split the image into 3x4 blocks.
%Paper: http://rabida.uhu.es/dspace/bitstream/handle/10272/4446/detecting_the_optic_disc_boundary.pdf
[mb,nb] = bestblk([m n], 600);
numBlocksVer = m/mb;
numBlocksHor = n/nb;
C_R = mat2cell(imgR, [mb mb mb], [nb nb nb nb]);
C_G = mat2cell(imgG, [mb mb mb], [nb nb nb nb]);
% figure('Name','Divided Red Image','NumberTitle','off');
% imshow(imgR);
% for i = 1:numBlocksVer*numBlocksHor
%     subplot(numBlocksVer,numBlocksHor,i);
%     imshow(C_R{i});
% end
% figure('Name','Divided Green Image','NumberTitle','off');
% imshow(imgG);
% for i = 1:numBlocksVer*numBlocksHor
%     subplot(numBlocksVer,numBlocksHor,i);
%     imshow(C_G{i});
% end


%Isolate the area around the optical disk.
counter = 0;
[m n] = size(C_R);
maxPixels = 0;
%Find the brightest pixel in the image.
brightest = max(max(imgR));
%Find the block containing the optical disk by counting, for each block,
%how many pixels it contains that are equal or greater than 80% of the
%brigthness of the brightest pixel. The block that contains the most of
%those pixels is the block containing the optical disk.
for i = 1:m
    for j = 1:n
        Temp_R = C_R{i,j};
        Temp_G = C_G{i,j};
        counter = 0;
        for k = 1:mb
            for l = 1:nb
                if Temp_R(k,l) >= 0.8*brightest
                    counter = counter + 1;
                end
            end
        end
        if counter > maxPixels
           maxPixels = counter;
           OD_R = Temp_R;
           OD_G = Temp_G;
           ind = [i j];
        end
    end
end
% figure('Name','Optical Disk Area - Red','NumberTitle','off');
% imshow(OD_R);
% figure('Name','Optical Disk Area - Green','NumberTitle','off');
% imshow(OD_G);


%Find the block containing the macula by taking into account the fact that
%the macula is roughly at the same position and distance in relation to the
%optical disk on every image. If the optical disk is on the right of the
%image, the macula will be in the adjacent block on the left of the block
%containing the optical disk and vice versa.
if ((ind(1)-1) == 1) && ((ind(2)-1) == 2)
    MAC = C_G{2,2};
    indMAC = [1 1];
elseif ((ind(1)-1) == 1) && ((ind(2)-1) == 1)
    MAC = C_G{2,3};
    indMAC = [1 2];
end
% figure('Name','Macula','NumberTitle','off');
% imshow(MAC);


%Image thresholding using Otsu's multi-threshold.
%Paper: http://rabida.uhu.es/dspace/bitstream/handle/10272/4446/detecting_the_optic_disc_boundary.pdf
%Segment pixels of OD into 12 classes. IDX is an array containing the class
%each pixel belongs to (1,2,3,4,5,6,7,8,9,10,11,12). Use Red channel.
OD = OD_R;
class = 12;
IDX = otsu(OD,class);
figure('Name','Multi-Thresholding Classes (Optical Disk - Red Channel)','NumberTitle','off');
imagesc(IDX);
[m n] = size(OD);
OD = zeros(m,n);
Temp = zeros(m,n);
for i=1:m
    for j=1:n
        %Make a binary image of the Optical Disk out of the pixels in
        %class 12.
        if IDX(i,j) == class
            OD(i,j) = 1;
        end
    end
end
%If the red channel is too bright and we can't threshold the Optical Disk,
%use the Green channel instead.
numPixels = find(find(OD), 1, 'last' );
if numPixels > 55000 
    OD = OD_G;
    IDX = otsu(OD,class);
    figure('Name','Multi-Thresholding Classes (Optical Disk - Green Channel)','NumberTitle','off');
    imagesc(IDX);
    OD = zeros(m,n);
    for i=1:m
        for j=1:n
            if IDX(i,j) == class
                OD(i,j) = 1;
            end
        end
    end 
end
%If number of pixels in thresholded Optical Disk is less than the specified
%number (meaning initial thresholding didn't provide enough information
%about the optical disk) then use 1 more class to make the binary image.
%Iterate until there's enough information about the optical disk.
class = class - 1;
numPixels = 0;
numPixels_next = 0;
while ((numPixels < 25000) && (numPixels_next < 50000))
     for i=1:m
        for j=1:n
            if IDX(i,j) == class
                OD(i,j) = 1;
            end
        end
     end
    numPixels = find(find(OD), 1, 'last' );
    class = class - 1;
    for i=1:m
        for j=1:n
            if IDX(i,j) == class
                Temp(i,j) = 1;
            end
        end
    end
    numPixels_next = find(find(Temp), 1, 'last' );
end
%Fill holes and discard objects smaller than 50px.
OD = imfill(OD,'holes');
OD = bwareaopen(OD,50);
figure('Name','Thresholded Optical Disk','NumberTitle','off');
imshow(OD);

%Same as above for Cup, using the Green channel, using 16 classes.
OD_CUP = C_G{ind(1),ind(2)};
class = 16;
numPixels_next = 0;
IDX = otsu(OD_CUP,class);
figure('Name','Multi-Thresholding Classes (Cup)','NumberTitle','off');
imagesc(IDX);
OD_CUP = zeros(m,n);
Temp = zeros(m,n);
for i=1:m
    for j=1:n
        if IDX(i,j) == class
            OD_CUP(i,j) = 1;
        end
    end
end
numPixels = find(find(OD_CUP, 1, 'last'));
while ((numPixels < 4000) && (numPixels_next < 8000))
     for i=1:m
        for j=1:n
            if IDX(i,j) == (class - 1)
                OD_CUP(i,j) = 1;
            end
        end
     end
    numPixels = find(find(OD_CUP), 1, 'last');
    class = class - 1;
    for i=1:m
        for j=1:n
            if IDX(i,j) == class
                Temp(i,j) = 1;
            end
        end
    end
numPixels_next = find(find(Temp), 1, 'last' );
end
OD_CUP = imfill(OD_CUP,'holes');
OD_CUP = bwareaopen(OD_CUP,50);
figure('Name','Thresholded Cup','NumberTitle','off');
imshow(OD_CUP);

%Same as above for Macula, using 16 classes.
class = 16;
IDX = otsu(MAC,class);
figure('Name','Multi-Thresholding Classes (Macula)','NumberTitle','off');
imagesc(IDX);
MAC = zeros(m,n);
%We now want to start from the class containing the darker pixels, which is
%class 1.
class = class - 15;
for i=1:m
    for j=1:n
        %Make a binary image of the Macula out of the pixels in class 1
        if IDX(i,j) == class
            MAC(i,j) = 1;
        end
    end
end
numPixels = find(find(MAC), 1, 'last');
while numPixels < 7000
     for i=1:m
        for j=1:n
            if IDX(i,j) == (class + 1)
                MAC(i,j) = 1;
            end
        end
     end
    numPixels = find(find(MAC), 1, 'last');
    class = class + 1;
end
MAC = imfill(MAC,'holes');
% figure('Name','Thresholded Macula','NumberTitle','off');
% imshow(MAC);

%Locate Optical Disk using Connected Components.
%Papers: http://www.ijarcsse.com/docs/papers/April2012/Volume_2_issue_4/V2I400117.pdf       
s = regionprops(OD, 'area','centroid');
%Find the largest CC object and its Center Of Mass.
areas = cat(1, s.Area);
centroids = cat(1, s.Centroid);
[~, i] = max(areas);
%Create and display the label matrix of the CC objects.
LabelCC = labelmatrix(bwconncomp(OD));
figure('Name','Labeled CC Objects with COMs - Optical Disk','NumberTitle','off');
imshow(imcomplement(label2rgb(LabelCC)));
hold on;
plot(centroids(:,1), centroids(:,2), 'black*');
hold off;
%Make a circle of radius R1 around the COM of the largest CC object.
%Every pixel outside this circle is set to be 0.
[m n] = size(OD);
for k = 1:m
    for l = 1:n
        if sqrt((l-centroids(i,1))^2 + (k-centroids(i,2))^2) >= 175
            OD(k,l) = 0;
        end
    end
end
%Merge all the CC objects which have less than a certain distance between
%them so that they form a single object. 
OD_Merged = bwdist(OD) <= 50;
%Find the new COM of this object which will be the COM of the Optical
%Disk.
s = regionprops(OD_Merged,'area','centroid');
areas = cat(1, s.Area);
[~, i] = max(areas);
centroids = cat(1,s.Centroid);
%Display the label matrix.
LabelCC_Merged = labelmatrix(bwconncomp(OD_Merged));
LabelCC_Merged(~OD) = 0;
figure('Name','Merged CC Objects with new COM - Optical Disk','NumberTitle','off');
imshow(imcomplement(label2rgb(LabelCC_Merged)));
hold on;
plot(centroids(i,1), centroids(i,2), 'black*');
hold off;
%Make a smaller circle of R2<R1 around the COM of the optical disk 
%which will be used to locate its borders. Every pixel outside this circle 
%is set to be 0.
for k = 1:m
    for l = 1:n
        if sqrt((l-centroids(i,1))^2 + (k-centroids(i,2))^2) >= 125
            OD(k,l) = 0;
        end
    end
end


%Locate Cup using Connected Components.
s = regionprops(OD_CUP, 'area','centroid');
%Find the largest CC object and its Center Of Mass.
areasCup = cat(1, s.Area);
centroidsCup = cat(1, s.Centroid);
[~, i] = max(areasCup);
% %Create and display the label matrix of the CC objects.
% LabelCC = labelmatrix(bwconncomp(OD_CUP));
% figure('Name','Labeled CC Objects with COMs - Cup','NumberTitle','off');
% imshow(imcomplement(label2rgb(LabelCC)));
% hold on;
% plot(centroidsCup(:,1), centroidsCup(:,2), 'black*');
% hold off;
%Make a circle of radius R1 around the COM of the largest CC object.
%Every pixel outside this circle is set to be 0.
[m n] = size(OD_CUP);
for k = 1:m
    for l = 1:n
        if sqrt((l-centroidsCup(i,1))^2 + (k-centroidsCup(i,2))^2) >= 125
            OD_CUP(k,l) = 0;
        end
    end
end
%Merge all the CC objects so that they form a single object.
OD_Merged = bwdist(OD_CUP) <= 20;
%Find the new COM of this object which will be the center of the Optical
%Disk.
s = regionprops(OD_Merged,'area','centroid');
centroidsCup = cat(1,s.Centroid);
areasCup = cat(1, s.Area);
[~, i] = max(areasCup);
% %Display the label matrix.
% LabelCC_Merged = labelmatrix(bwconncomp(OD_Merged));
% LabelCC_Merged(~OD_CUP) = 0;
% figure('Name','Merged CC Objects with new COM - Cup','NumberTitle','off');
% imshow(imcomplement(label2rgb(LabelCC_Merged)));
% hold on;
% plot(centroidsCup(1,1), centroidsCup(1,2), 'black*');
% hold off;
%Make a smaller circle of R2<R1 around the COM of the Cup which will 
%be used to locate its borders. Every pixel outside this circle is set to 
%be 0.
for k = 1:m
    for l = 1:n
        if sqrt((l-centroidsCup(i,1))^2 + (k-centroidsCup(i,2))^2) >= 100
            OD_CUP(k,l) = 0;
        end
    end
end


%Locate Macula using Connected Components.
s = regionprops(MAC, 'area','centroid');
%Find the largest CC object and its Center Of Mass.
areasMac = cat(1, s.Area);
centroidsMac = cat(1, s.Centroid);
[~, i] = max(areasMac);
% %Create and display the label matrix of the CC objects.
% LabelCC = labelmatrix(bwconncomp(MAC));
% figure('Name','Labeled CC Objects with COMs - Macula','NumberTitle','off');
% imshow(imcomplement(label2rgb(LabelCC)));
% hold on;
% plot(centroidsMac(:,1), centroidsMac(:,2), 'black*');
% hold off;
%Make a circle of radius R1 around the COM of the largest CC object.
%Every pixel outside this circle is set to be 0.
[m n] = size(MAC);
for k = 1:m
    for l = 1:n
        if sqrt((l-centroidsMac(i,1))^2 + (k-centroidsMac(i,2))^2) >= 200
            MAC(k,l) = 0;
        end
    end
end
%Merge all the CC objects so that they form a single object.
MAC_Merged = bwdist(MAC) <= 10;
%Find the new COM of this object which will be the center of the Optical
%Disk.
s = regionprops(MAC_Merged,'area','centroid');
centroidsMac = cat(1,s.Centroid);
areasMac = cat(1, s.Area);
[~, i] = max(areasMac);
% %Display the label matrix.
% LabelCC_Merged = labelmatrix(bwconncomp(MAC_Merged));
% LabelCC_Merged(~MAC) = 0;
% figure('Name','Merged CC Objects with new COM - Macula','NumberTitle','off');
% imshow(imcomplement(label2rgb(LabelCC_Merged)));
% hold on;
% plot(centroidsMac(1,1), centroidsMac(1,2), 'black*');
% hold off;
%Make a smaller circle (R2<R1), now around the COM of the macula we just 
%found. Every pixel outside this circle is set to be 0.
for k = 1:m
    for l = 1:n
        if sqrt((l-centroidsMac(i,1))^2 + (k-centroidsMac(i,2))^2) >= 150
            MAC(k,l) = 0;
        end
    end
end


% Obtain the outlines of the Optical Disc, Cup and Macula.
[m n] = size(imgR);
per_OD = zeros(m,n);
per_OD_CUP = zeros(m,n);
per_MAC = zeros(m,n);
% Put the outlines of the Optical Disk, Cup and Macula at the correct position 
% on the original image dimensions.
% Optical Disk - Cup.
upper = (ind(1)-1)*mb;
lower = (ind(1)-1)*mb + mb -1;
left = (ind(2)-1)*nb;
right = (ind(2)-1)*nb + nb -1;
per_OD(upper:lower, left:right) = bwperim(OD);
imwrite(per_OD, 'per.tif');
smooth = fspecial('disk',15);
OD_CUP = imfilter(OD_CUP,smooth);
per_OD_CUP(upper:lower, left:right) = bwperim(OD_CUP);
% Macula.
upper = indMAC(1)*mb;
lower = indMAC(1)*mb + mb -1;
left = indMAC(2)*nb;
right = indMAC(2)*nb + nb -1;
per_MAC(upper:lower, left:right) = bwperim(MAC);


%Approximate the perimeters with a circle calculated with Kasa's method.
%Paper: AUTOMATIC DETECTION OF OPTIC DISC BASED ON PCA AND STOCHASTIC
%WATERSHED, Sandra Morales, Valery Naranjo, David Perez, Amparo Navea, 
%Mariano Alcaniz.
%2.4: Circular approximation
%Optical Disk.
[X,Y] = find(per_OD);
XY = cat(2,X,Y);
FIT_OD = Kasa(XY);

%Macula.
[X,Y] = find(per_MAC);
XY = cat(2,X,Y);
FIT_MAC = Kasa(XY);

per_OD = zeros(m,n);
for i = 1:m
    for j = 1:n
        if sqrt((i-FIT_OD(1))^2 + (j-FIT_OD(2))^2) <= FIT_OD(3)
            per_OD(i,j) = 1;
        end
    end
end
per_OD = bwperim(per_OD);


%Obtain green channel of the Optical Disk area.
tmp = imgG(:,:);
for i = 1:m
    for j = 1:n
        if sqrt((i-FIT_OD(1))^2 + (j-FIT_OD(2))^2) >= FIT_OD(3)
            tmp(i,j) = 0;
        end
    end
end
figure('Name','Optical Disk area on Green channel','NumberTitle','off');
imshow(tmp);
imwrite(tmp,'Disk.tif');

per_MAC = zeros(m,n);
for i = 1:m
    for j = 1:n
        if sqrt((i-FIT_MAC(1))^2 + (j-FIT_MAC(2))^2) <= FIT_MAC(3)
            per_MAC(i,j) = 1;
        end
    end
end
per_MAC = bwperim(per_MAC);


%Dilate outlines to be more visible.
SE_outl = strel('disk',1,8);
per_OD = imdilate(per_OD,SE_outl);
per_OD_CUP = imdilate(per_OD_CUP,SE_outl);
per_MAC = imdilate(per_MAC,SE_outl);


%Impose outlines on original image.
for i = 1:m,
     for j = 1:n,
        if per_OD(i,j) == 1
            img(i,j,:) = 255;
        end
        if per_OD_CUP(i,j) == 1
            img(i,j,:) = 0;
        end
        if per_MAC(i,j) == 1
            img(i,j,2) = 255;
        end
    end
end
figure('Name','Located Optic Disk and Cup','NumberTitle','off');
imshow(img);
