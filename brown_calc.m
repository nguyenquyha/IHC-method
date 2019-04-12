% Script to compute the ratio of brownish pixels in a set of images
% Author: Ha Q. Nguyen @ Institute of Big Data - Vingroup
% Date Modified: April 02, 2019

%% Prepare the image

close all
clear

% find all subfolders
files = dir(pwd);
dirFlags = [files.isdir];
subdirs = files(dirFlags);
subdirs(1:2) = [];

% lower and upper HSV for brownish pixels
lb = [0.015, 0.15, 0.15];
ub = [0.15, 0.8, 0.8];

brown_table = cell(1000,2);
im_idx = 0;

%% Loop through all imgages in all subfolders

for i = 1:length(subdirs)
    folder_name = subdirs(i).name;
    folder = subdirs(i).folder;
    subfiles = dir(fullfile(folder, folder_name, '*.png'));
    for j = 1:length(subfiles)
        if ~ismember('mask', subfiles(j).name)% exclude the mask files
            
            fname = fullfile(folder, folder_name, subfiles(j).name);
            [r, mask] = brown_ratio(fname, lb, ub);
            
            fname_mask = [fname(1:end-4), '_mask.png'];
            
            % write the mask images into the same folder
            if ~exist(fname_mask,'file')
                imwrite(mask, fname_mask)
            end

            % write the brown ratios into a table
            im_idx = im_idx + 1;
            brown_table{im_idx, 1} = subfiles(j).name;
            brown_table{im_idx, 2} = r;
        end
    end
end

brown_table(im_idx + 1 : end , :) = [];
brown_table = cell2table(brown_table, 'VariableNames',{'Image' , 'Brown_Ratio'});
writetable(brown_table, 'brown_ratio_table.csv');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper funtion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [r, mask] = brown_ratio(fname, lb, ub)
    % Load image
    im = imread(fname);
    % convert from RGB to gray scale
    im_gray = rgb2gray(im);

%   figure, imshow(im)

    %% Find the ball mask
    mask = im_gray < 230 & im_gray > 10;
    % remove the noise outside
    se = strel('disk', 3);
    mask = imopen(mask, se);
    
    % fill the holes
    se2 = strel('disk', 5);    
    mask = imclose(mask, se2);
    
    % imfill(mask,'holes');
    CC = bwconncomp(mask);
    
    % noise removal
    if CC.NumObjects > 1
        object_size = zeros(1,CC.NumObjects);
        for i = 1:CC.NumObjects
            object_size(i) = length(CC.PixelIdxList{i});
        end
        
        [~ , idx] = sort(object_size, 'descend');
        for obj = idx(2:end)
            mask(CC.PixelIdxList{obj}) = 0;
        end
    end
    
   % find the convex hull
    mask = bwconvhull(mask);
    
    %% Count the brownish pixels
    im_hsv = rgb2hsv(mat2gray(im));

    im_hsv_vec = reshape(im_hsv, size(im_hsv,1) * size(im_hsv, 2) , 3);
    
    seg = im_hsv_vec(mask(:), :);
    is_brown = seg > lb & seg < ub;
    num_brown = sum(sum(is_brown, 2) == 3);

    r = num_brown / sum(mask(:));
end