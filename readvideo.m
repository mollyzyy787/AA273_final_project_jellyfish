%close all; clear all; clc;

threshPercent = 99;
v = VideoReader('jelly_movie_trim.mp4');

location = zeros(2,500);
area_history = zeros(1,500);
height_history = zeros(1,500);
width_history = zeros(1,500);
particle_location = zeros(2,500);
flow_velocity_history = zeros(2,499);

% M2 = VideoWriter('particle.avi');
% open(M2);
for i = 1:500
    img_ori = read(v,i);
    img = double(img_ori);
    img = imcrop(img,[350 122 1000 483]);
    img = 0.299*img(:,:,1) + 0.587*img(:,:,2) + 0.114*img(:,:,3);
    img = uint8(img);
    img_top = imtophat(img,strel('disk',25));
    bi_top = img_top > threshPercent*mean(mean(img))/100;
    CC=bwconncomp(bi_top);
    s=regionprops(CC);
    areas=[s.Area];
    centroids = zeros(length(areas),2);
    boundingbox=zeros(length(areas),4);
    for j=1:length(areas)
        boundingbox(j,:)=s(j).BoundingBox;
        centroids(j,:)=s(j).Centroid;
    end
    [Area,Ind]=sort(areas,'descend'); 
%     imshow(bi_top);
%     hold on;
%     rectangle('Position',[boundingbox(Ind(1),:)],'EdgeColor','r', 'LineWidth', 2);hold on;
%     plot(centroids(Ind(1),1), centroids(Ind(1),2), 'r.','LineWidth',16); 
%     axis off
%     drawnow;
%     currFrame = getframe;
%     writeVideo(M,currFrame);
    location(:,i)=centroids(Ind(1),:);
    area_history(i) = Area(Ind(1));
    height_history(i) = boundingbox(Ind(1),4);
    width_history(i) = boundingbox(Ind(1),3);
    new_im = bi_top;         
    xmin = floor(boundingbox(Ind(1),1));
    ymin = floor(boundingbox(Ind(1),2));
    width = floor(boundingbox(Ind(1),3));
    height = floor(boundingbox(Ind(1),4));
    new_im(ymin:ymin+height-1,xmin:xmin+width-1)=...
        zeros(height,width);
    big_bb = zeros(4,1);
    big_bb(1) = xmin-50;
    big_bb(2) = ymin-50;
    big_bb(3) = width+100;
    big_bb(4) = height+100;
    cropped = imcrop(new_im, big_bb');
    newCC=bwconncomp(cropped);
    new_s=regionprops(newCC);
    new_areas = [new_s.Area];
    new_centroids = zeros(length(new_areas),2);
    for j=1:length(new_areas)
        new_centroids(j,:)=new_s(j).Centroid;
    end
    xc = sum(new_centroids(:,1))/length(new_centroids(:,1))+big_bb(1);
    yc = sum(new_centroids(:,2))/length(new_centroids(:,1))+big_bb(2);
    particle_location(1,i) = xc;
    particle_location(2,i) = yc;
    imshow(new_im)
    rectangle('Position',big_bb,'EdgeColor','g', 'LineWidth', 2);hold on;
    plot(new_centroids(:,1)+big_bb(1), new_centroids(:,2)+big_bb(2), 'r.','LineWidth',16)
    plot(xc,yc,'bo','LineWidth',4)
%     currFrame = getframe;
%     writeVideo(M2,currFrame);
end
%close(M2);


%%
[row_flowV, col_flowV] = size(flow_velocity_history);
for i = 2:col_flowV
    flow_velocity_history(1,i) = particle_location(1,i)-particle_location(1,i-1);
    flow_velocity_history(2,i) = particle_location(2,i)-particle_location(2,i-1);
end

filtered_vpx = flow_velocity_history(1,:);
filtered_vpy = flow_velocity_history(2,:);

for i = 10:length(filtered_vpx)
    filtered_vpx(i) = mean(filtered_vpx(i-9:i));
    filtered_vpy(i) = mean(filtered_vpy(i-9:i));
end
figure;
plot(filtered_vpx)
figure;
plot(filtered_vpy)
 
for i=1:length(location(1,:))
    plot(location(1,i),-location(2,i),'.');hold on;
    title('jellyfish trajectory measurement');
    drawnow;  
end

for i=1:length(particle_location(1,:))
    plot(particle_location(1,i),-particle_location(2,i),'.');hold on;
    drawnow;
end
 
