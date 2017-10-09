clear; clc;

image_list = dir('./front/*.png');

for i = 1:length(image_list)
    img = imread([image_list(i).folder, '/', image_list(i).name]);
     
end