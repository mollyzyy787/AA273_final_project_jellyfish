banana = ones(50,70);
figure;
imshow(banana)
title('banana 1')

xmin = 20
width = 20
ymin = 30
height = 30

banana(xmin:xmin+width-1,ymin:ymin+height-1)=...
            zeros(width,height);
        
figure;
imshow(banana)
title('banan 2')