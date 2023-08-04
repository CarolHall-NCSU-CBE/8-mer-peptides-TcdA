%Image Processing of Secondary Screening Beads
%4/16/21
%Adapted from Sahand Saberi and Brandyn Moore

%Updated 6/6/22
%Adapted from previous version to calculate average green fluorescence of
%red halo only. This is more meaningful in the context of peptide:TcdA GTD
%(red) binders, where the exclusion of UDP-Glucose-Fluorescein (green) is
%desired.

%Updated 8/4/23 with comments for publication

clc
clear all
close all
greendirname='Green/';
reddirname='Red/';
%% Image processing
zgreen=dir(fullfile(greendirname,'*.jpg'));
zred=dir(fullfile(reddirname,'*.jpg'));
imred=cell(size(zred,1),1);
imgreen=cell(size(zgreen,1),1);
n=0;
o=0;
for i=1:length(zred)
    filegreen=strcat(greendirname,zgreen(i).name);
    filered=strcat(reddirname,zred(i).name);
    imgreen{i}=imread(filegreen);
    imgreen{i}=im2gray(imgreen{i});
    imred{i}=imread(filered);
    imred{i}=im2gray(imred{i});
    
    T=graythresh(imred{i});
    T2=graythresh(imgreen{i});
    BW=imbinarize(imred{i},T);
    BW_1=imbinarize(imgreen{i},T2);
    
    BW2=bwareaopen(BW,100);
    BW_2=bwareaopen(BW_1,100);
    
    BW3=imfill(BW2,'holes');
    BW_3=imfill(BW_2,'holes');
   
    
    %Visual image of beads with black background
    masked_r=imred{i};
    masked_r(~BW3)=0;
    masked_g=imgreen{i};
    masked_g(~BW3)=0;                   %mask green with same mask as red
    masked_go(~BW_3)=0;
    
    
    %Pixels from beads only (remove background pixels)
    bead_r=double(imred{i}(BW3));
    bead_g=double(imgreen{i}(BW3));     %mask green with same mask as red
    bead_go=double(imgreen{i}(BW_3));
    
    logic=length(find(bead_r));
    logic1=length(find(bead_g));
    

    %% Extract pixels from desired regions of beads & compute means
    k = 0;
    m = 0;
    redpixels = [];
    greenpixels = [];
    greenpixelso = [];
    if logic1~=0&&logic~=0
        pctred = prctile(bead_r(:),90);
        pctgreen = prctile(bead_g(:),50);
        for j = 1:size(bead_r)
            if bead_r(j) > pctred
                k = k+1;
                redpixels(k) = bead_r(j);   %red pixels from 90th prctile red
                greenpixels(k) = bead_g(j); %green pixels from 90th prctile red
            end
        end
        
        n = n+1;
        redmeans(n) = mean(redpixels);      %mean of red pixels from 90th prctile red
        greenmeans(n) = mean(greenpixels);  %mean of green pixels from 90th prctile red

        for l = 1:size(bead_go)
            if bead_go(l) > pctgreen
                m = m+1;
                greenpixelso(m) = bead_go(l);   %green pixels from 50th prctile green
            end
        end
        
        o = o+1;
        greenavg(o) = mean(greenpixelso);       %mean of green pixels from 50th prctile green
    end
end

%% Report statistics
redmean = mean(redmeans)        %mean of means of red pixels from 90th prctile red
redstd = std(redmeans)          %standard deviation of means of red pixels from 90th prctile red
greenmean = mean(greenmeans)    %mean of means of green pixels from 90th prctile red
greenstd = std(greenmeans)      %standard deviation of means of green pixels from 90th prctile red
greenmeannorm = mean(greenmeans)-mean(greenavg) %difference between means of green pixels from 90th prctile red and of green pixels from 50th prctile green 
greenpercentchange = greenmeannorm/mean(greenavg)*100 %percent change between green means