%% Generate thicknessIRL for all participants in pathlist
% Format: all_thicknessIRL(:,:,participant rank in pathlist)
all_thicknessIRL = zeros(512,512,length(pathlist));
for iPath = 1:length(pathlist);
    folder = fullfile(basepath,pathlist{iPath});
    all_thicknessIRL(:,:,iPath) = Cirrus_OCT_thicknessIRL(folder);
end

% Flip right eye participants to left eye orientation
% Eye information should be in Cirrus_OCT_pathlist.m
for participant = [3 4 5 9 10 12 13]
    all_thicknessIRL(:,:,participant)=fliplr(all_thicknessIRL(:,:,participant));
end

%% Generate mean thicknessIRL and standard deviation for each group
% Adjust range in Z-dimension to include appropriate controls/patients.
% Format: con_thicknessIRLMean/con_thicknessIRLStd/pat_thicknessIRLMean/
%         pat_thicknessIRLStd(512,512,1)
con_thicknessIRLMean = mean(all_thicknessIRL(:,:,6:13),3);
pat_thicknessIRLMean = mean(all_thicknessIRL(:,:,1:5),3);
con_thicknessIRLStd = std(all_thicknessIRL(:,:,6:13),0,3);
pat_thicknessIRLStd = std(all_thicknessIRL(:,:,1:5),0,3);

%% Descriptive Statistics Summary (INDIVIDUAL Participants)
% Format: dstat1

dstat1 = zeros(length(pathlist),7,1);
% Max for each participant
for iPath = 1:length(pathlist);
    dstat1(iPath,1,1) = max(max(all_thicknessIRL(:,:,iPath)));
end
% Min for each participant
for iPath = 1:length(pathlist);
    dstat1(iPath,2,1) = min(min(all_thicknessIRL(:,:,iPath)));
end
% Range for each participant
for iPath = 1:length(pathlist);
    dstat1(iPath,3,1) = (dstat1(iPath,1,1)) - (dstat1(iPath,2,1));
end
% Median for each participant
for iPath = 1:length(pathlist);
    dstat1(iPath,4,1) = median(median(all_thicknessIRL(:,:,iPath)));
end
% Mean for each participant
for iPath = 1:length(pathlist);
    dstat1(iPath,5,1) = mean2(all_thicknessIRL(:,:,iPath));
end
% Std for each participant
for iPath = 1:length(pathlist);
    dstat1(iPath,6,1) = std2(all_thicknessIRL(:,:,iPath));
end
% Var for each participant
for iPath = 1:length(pathlist);
    dstat1(iPath,7,1) = var(reshape(all_thicknessIRL(:,:,iPath),[],1));
end

% Formatting table
iPath = 1:length(pathlist);
printmat(dstat1,'Descriptive Statistics Summary (INDIVIDUAL)',...
    ['Participant',num2str(iPath)],'Max Min Range Median Mean Std Var');

%% Descriptive Statistics Summary (GROUPED Participants)
% Format: dstat2
% First, generate mean thicknessIRL for each group.
pat_thicknessIRLMean = mean(all_thicknessIRL(:,:,1:5),3);
con_thicknessIRLMean = mean(all_thicknessIRL(:,:,6:13),3);

dstat2 = zeros(2,5,1);
% Max for each group
dstat2(1,1,1) = max(pat_thicknessIRLMean(:));
dstat2(2,1,1) = max(con_thicknessIRLMean(:));
% Min for each group
dstat2(1,2,1) = min(pat_thicknessIRLMean(:));
dstat2(2,2,1) = min(con_thicknessIRLMean(:));
% Range for each group
dstat2(1,3,1) = (dstat2(1,1,1)) - (dstat2(1,2,1));
dstat2(2,3,1) = (dstat2(2,1,1)) - (dstat2(2,2,1));
% Median for each group
dstat2(1,4,1) = median(pat_thicknessIRLMean(:));
dstat2(2,4,1) = median(con_thicknessIRLMean(:));
% Mean for each group
dstat2(1,5,1) = mean(pat_thicknessIRLMean(:));
dstat2(2,5,1) = mean(con_thicknessIRLMean(:));
% Std for each group
dstat2(1,6,1) = std(pat_thicknessIRLMean(:));
dstat2(2,6,1) = std(con_thicknessIRLMean(:));
% Var for each group
dstat2(1,7,1) = var(pat_thicknessIRLMean(:));
dstat2(2,7,1) = var(con_thicknessIRLMean(:));

% Formatting table
iPath = 1:length(pathlist);
printmat(dstat2,'Descriptive Statistics Summary (GROUPED)',...
    'T1D Control','Max Min Range Median Mean Std Var');

%% Descriptive Statistics Summary (DIFFERENCE between groups)
% Format: dstat3
% First, generate mean thicknessIRL for each group.
pat_thicknessIRLMean = mean(all_thicknessIRL(:,:,1:5),3);
con_thicknessIRLMean = mean(all_thicknessIRL(:,:,6:13),3);

diff_thicknessIRLMean = con_thicknessIRLMean - pat_thicknessIRLMean;

dstat3 = zeros(1,5,1);
% Max
dstat3(1,1,1) = max(diff_thicknessIRLMean(:));
% Min
dstat3(1,2,1) = min(diff_thicknessIRLMean(:));
% Range
dstat3(1,3,1) = range(diff_thicknessIRLMean(:));
% Median
dstat3(1,4,1) = median(diff_thicknessIRLMean(:));
% Mean
dstat3(1,5,1) = mean(diff_thicknessIRLMean(:));
% Std
dstat3(1,6,1) = std(diff_thicknessIRLMean(:));
% Var
dstat3(1,7,1) = var(diff_thicknessIRLMean(:));

% Formatting table
iPath = 1:length(pathlist);
printmat(dstat3,'Descriptive Statistics Summary (DIFFERENCE)',...
    'Difference','Max Min Range Median Mean Std Var');

%% MISCELLANEOUS (for Progress Report)

% PAGE 2-4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IRL thickness of each individual participant:
% Suggest patient / control information gets moved to Cirrus_OCT_pathlist.m

for i=1:5
figure;
imagesc(all_thicknessIRL(:,:,i));
colorbar;
[cmin,cmax] = caxis;
caxis([0,100]);
title(['Participant: ',num2str(i)]);
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
ylabel(colorbar,'Thickness (pixel)');
end

for i=6:13
figure;
imagesc(all_thicknessIRL(:,:,i));
colorbar;
[cmin,cmax] = caxis;
caxis([0,100]);
title(['Participant: ',num2str(i)]);
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
ylabel(colorbar,'Thickness (pixel)');
end

% PAGE 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Boxplot of averaged control minus averaged T1D IRL thickness
all_thicknessIRL_columns = zeros(262144,13,1);
for i = 1:13;
    for j = 2:512;
        all_thicknessIRL_columns(((j-1)*(512)):((j*512)-1),i,1)...
            = all_thicknessIRL(1:512,j,i);
        all_thicknessIRL_columns(1:512,i,1) = all_thicknessIRL(1:512,1,i);
        all_thicknessIRL_columns(262144,i,1) = all_thicknessIRL(512,512,i);
    end
end
boxplot(all_thicknessIRL_columns);
title({'IRL Thickness of Individual Participants'});
xlabel('Participant');
ylabel('IRL Thickness');

% PAGE 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Averaged control and averaged T1D IRL thicknesses
figure;
imagesc(con_thicknessIRLMean(:,:));
colorbar;
[cmin,cmax] = caxis;
caxis([0,100]);
title({'Average IRL Thickness: Control'});
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
ylabel(colorbar,'Thickness (pixel)');
axis square;

figure;
imagesc(pat_thicknessIRLMean(:,:));
colorbar;
[cmin,cmax] = caxis;
caxis([0,100]);
title({'Average IRL Thickness: T1D'});
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
ylabel(colorbar,'Thickness (pixel)');
axis square;

% Averaged control and averaged T1D IRL thicknesses
grouped_thicknessIRL_columns = zeros(262144,2,1);
for j = 2:512;
    grouped_thicknessIRL_columns(((j-1)*(512)):((j*512)-1),1,1)...
        = con_thicknessIRLMean(1:512,j,1);
    grouped_thicknessIRL_columns(1:512,1,1) = con_thicknessIRLMean(1:512,1,1);
    grouped_thicknessIRL_columns(262144,1,1) = con_thicknessIRLMean(512,512,1);
end
for j = 2:512;
    grouped_thicknessIRL_columns(((j-1)*(512)):((j*512)-1),2,1)...
        = pat_thicknessIRLMean(1:512,j,1);
    grouped_thicknessIRL_columns(1:512,2,1) = pat_thicknessIRLMean(1:512,1,1);
    grouped_thicknessIRL_columns(262144,2,1) = pat_thicknessIRLMean(512,512,1);
end
boxplot(grouped_thicknessIRL_columns,'Labels',{'Control','T1D'});
title({'Averaged Control and Averaged T1D IRL Thicknesses'});
xlabel('Group');
ylabel('IRL Thickness');

% PAGE 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histograms of averaged control and averaged T1D IRL thicknesses
figure;
histogram(con_thicknessIRLMean(:,:),'FaceColor',[0 1 0],'EdgeColor','none');
title({'Group: Control'});
xlabel('Thickness (pixel)');
ylabel('Frequency');

figure;
histogram(pat_thicknessIRLMean(:,:),'FaceColor',[1 0 0],'EdgeColor','none');
title({'Group: T1D'});
xlabel('Thickness (pixel)');
ylabel('Frequency');

% Overlay of histograms
figure;
histogram(con_thicknessIRLMean(:,:),'FaceColor',[0 1 0],'EdgeColor','none');
xlabel('Thickness (pixel)');
ylabel('Frequency');
hold on
histogram(pat_thicknessIRLMean(:,:),'FaceColor',[1 0 0],'EdgeColor','none');
title({'Overlay of Averaged Control and Averaged T1D IRL Thickness Frequencies'});
xlabel('Thickness (pixel)');
ylabel('Frequency');
legend('Control','T1D')

% PAGE 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corresponding histogram overlay with location on the retina
bottom_range=57.5;
top_range=70;
overlap=zeros(512,512,3);
con_green=zeros(512,512);
con_green(con_thicknessIRLMean>bottom_range & con_thicknessIRLMean<top_range)=255;
pat_red=zeros(512,512);
pat_red(pat_thicknessIRLMean>bottom_range & pat_thicknessIRLMean<top_range)=255;
overlap(:,:,1)=pat_red;
overlap(:,:,2)=con_green;
image(overlap);
% Primitive version
a=zeros(512,512,3);
b=zeros(512,512);
b(con_thicknessIRLMean>45 & con_thicknessIRLMean<50)=255;
c=zeros(512,512);
c(pat_thicknessIRLMean>45 & pat_thicknessIRLMean<50)=255;
a(:,:,2)=b;
a(:,:,1)=c;
figure;image(a);

% PAGE 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Averaged control minus averaged T1D IRL thickness
figure;
imagesc(diff_thicknessIRLMean(:,:));
colorbar;
[cmin,cmax] = caxis;
caxis([-15,15]);
title({'Difference in IRL Thickness';'(Average T1D Minus Average Control)'});
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
ylabel(colorbar,'Thickness (pixel)');
axis square;

% Boxplot of averaged control minus averaged T1D IRL thicknesses
diff_thicknessIRL_columns = zeros(262144,1,1);
for j = 2:512;
    diff_thicknessIRL_columns(((j-1)*(512)):((j*512)-1),1,1)...
        = diff_thicknessIRLMean(1:512,j,1);
    diff_thicknessIRL_columns(1:512,1,1) = diff_thicknessIRLMean(1:512,1,1);
    diff_thicknessIRL_columns(262144,1,1) = diff_thicknessIRLMean(512,512,1);
end
boxplot(diff_thicknessIRL_columns,'Labels',{'Averaged control minus averaged T1D'});
title({'Averaged Control Minus Averaged T1D IRL Thicknesses'});
ylabel('IRL Thickness');

% PAGE 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Histogram of averaged control minus averaged T1D IRL thickness
histogram(diff_thicknessIRLMean);
title({'Averaged Control Minus Averaged T1D IRL Thicknesses'});
xlabel('Thickness (pixel)');
ylabel('Frequency');

% Quantile-quantile (QQ) plot of averaged control minus averaged T1D IRL
% thicknesses
qqplot(diff_thicknessIRLMean(:));
title({'Averaged Control Minus Averaged T1D IRL Thicknesses'});

% PAGE 11 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QQ plot of p-values comparing control and T1D IRL thicknesses
hvalue_thicknessIRL = zeros(512,512,1);
pvalue_thicknessIRL = zeros(512,512,1);
for iRow = 1:512;
    for iCol = 1:512;
        [h,p] = ttest2(all_thicknessIRL(iRow,iCol,1:5),all_thicknessIRL(iRow,iCol,6:13));
        hvalue_thicknessIRL(iRow,iCol,1) = h;
        pvalue_thicknessIRL(iRow,iCol,1) = p;
    end
end
figure;
qqplot(pvalue_thicknessIRL(:));
title({'P-values Comparing Controls and T1D'});

% Scatter plot of p-values comparing control and T1D IRL thicknesses
figure;
spy(hvalue_thicknessIRL,'k',5);
title({'Locations of Significant IRL Thickness Difference','between Controls and T1Ds'});
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
legend('Significant locations')

%% Generate zscore for each patient
% Adjust range of iPath to include appropriate patients in pathlist.
% Format: zscore(:,:,patient rank in pathlist)
for iPath = 1:5;
    folder = fullfile(basepath,pathlist{iPath});
    zscore(:,:,iPath)=((all_thicknessIRL(:,:,iPath)...
        - con_thicknessIRLMean)./con_thicknessIRLStd);
end

%% Other Miscellaneous 
%Display zscore as mesh plot for each patient
figure;mesh(zscore(:,:,1))
caxis(r)
axis([1 512 1 512 r])

% Playing variance as a function of resolution
figure;
for i = 2:256
image(imresize(thicknessIRLVar{1,i},[256 256]))
pause(0.25)
end

%Creating region of interest (ROI) free-hand
figure;
image(overlap);
axis square;
%Prompts defining the ROI with the current axes of current figure (gca)
ROI=imfreehand(gca);
%Creates logical where 1s are region inside ROI and 0s are region outside
%ROI
BW=createMask(ROI);


%Plotting p-values for patients and controls
pat_matrix=zeros(5,262144,1);
for pat=1:5;
    pat_rows=reshape(all_thicknessIRL(:,:,pat),[1 262144]);
    pat_matrix(pat,:,1)=pat_rows;
end

con_matrix=zeros(8,262144,1);
for con=6:13;
    con_rows=reshape(all_thicknessIRL(:,:,con),[1 262144]);
    con_matrix(con-5,:,1)=con_rows;
end 

pvalue=zeros(512,512,1);
[h,p]=ttest2(pat_matrix,con_matrix);
pvalue=reshape(p,[512 512]);

%Calculate t-test between patients and controls by pixel, and drawing 3
%ROIs
pat_matrix=zeros(5,262144,1);
for pat=1:5;
    pat_rows=reshape(all_thicknessIRL(:,:,pat),[1 262144]);
    pat_matrix(pat,:,1)=pat_rows;
end

con_matrix=zeros(8,262144,1);
for con=6:13;
    con_rows=reshape(all_thicknessIRL(:,:,con),[1 262144]);
    con_matrix(con-5,:,1)=con_rows;
end 

pvalue=zeros(512,512,1);
[h,p]=ttest2(pat_matrix,con_matrix);
pvalue=reshape(p,[512 512]);

figure;
imagesc(pvalue);
axis square;
colormap gray;
title('Significance of Difference in IRL Thickness');
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
ylabel(colorbar,'P-value');

t = 0:pi/20:2*pi;
R0 = 42.66666666666667; x0 = 256; y0 = 256;
xi = R0*cos(t)+x0;
yi = R0*sin(t)+y0;
LineHandler = line(xi,yi,'LineWidth',1,'Color',[0 0 0]);
R1 = 128; x0 = 256; y0 = 256;
xj = R1*cos(t)+x0;
yj = R1*sin(t)+y0;
LineHandler = line(xj,yj,'LineWidth',1,'Color',[0 0 0]);
R3 = 256; x0 = 256; y0 = 256;
xk = R3*cos(t)+x0;
yk = R3*sin(t)+y0;
LineHandler = line(xk,yk,'LineWidth',1,'Color',[0 0 0]);

ROI1=imfreehand(gca);
region1=createMask(ROI1);
ROI2=imfreehand(gca);
region2=createMask(ROI2);
ROI3=imfreehand(gca);
region3=createMask(ROI3);
ROI4=imfreehand(gca);
region4=createMask(ROI3);

%% Settings for automatic figure docking
set(0,'DefaultFigureWindowStyle','docked')
set(0,'DefaultFigureWindowStyle','normal')





% Generating figure to select ROIs, and creating a boxplot and dotplot to
% display individual participant means
pat_matrix=zeros(5,262144,1);
for pat=1:5;
    pat_rows=reshape(all_thicknessIRL(:,:,pat),[1 262144]);
    pat_matrix(pat,:,1)=pat_rows;
end

con_matrix=zeros(8,262144,1);
for con=6:13;
    con_rows=reshape(all_thicknessIRL(:,:,con),[1 262144]);
    con_matrix(con-5,:,1)=con_rows;
end 

pvalue=zeros(512,512,1);
[h,p]=ttest2(pat_matrix,con_matrix);
pvalue=reshape(p,[512 512]);

figure;
imagesc(pvalue);
axis square;
colormap gray;
title('IRL thickness differences between patients and controls by location');
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
ylabel(colorbar,'P-value');

figure;
imagesc(diff_thicknessIRLMean(:,:));
axis square;
colorbar;
[cmin,cmax] = caxis;
caxis([-15,15]);
title({'IRL Thickness (Averaged Control - Averaged T1D)'});
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
ylabel(colorbar,'Thickness (pixel)');

ROI1=imfreehand(gca);
region1=createMask(ROI1);
ROI2=imfreehand(gca);
region2=createMask(ROI2);
ROI3=imfreehand(gca);
region3=createMask(ROI3);
ROI4=imfreehand(gca);
region4=createMask(ROI4);

part1=(all_thicknessIRL(:,:,1));
part2=(all_thicknessIRL(:,:,2));
part3=(all_thicknessIRL(:,:,3));
part4=(all_thicknessIRL(:,:,4));
part5=(all_thicknessIRL(:,:,5));
part6=(all_thicknessIRL(:,:,6));
part7=(all_thicknessIRL(:,:,7));
part8=(all_thicknessIRL(:,:,8));
part9=(all_thicknessIRL(:,:,9));
part10=(all_thicknessIRL(:,:,10));
part11=(all_thicknessIRL(:,:,11));
part12=(all_thicknessIRL(:,:,12));
part13=(all_thicknessIRL(:,:,13));


pat_regions=zeros(5,4,1);
pat_mean=mean(part1(region1));
pat_regions(1,1,1)=pat_mean;
pat_mean=mean(part2(region1));
pat_regions(2,1,1)=pat_mean;
pat_mean=mean(part3(region1));
pat_regions(3,1,1)=pat_mean;
pat_mean=mean(part4(region1));
pat_regions(4,1,1)=pat_mean;
pat_mean=mean(part5(region1));
pat_regions(5,1,1)=pat_mean;

pat_mean=mean(part1(region2));
pat_regions(1,2,1)=pat_mean;
pat_mean=mean(part2(region2));
pat_regions(2,2,1)=pat_mean;
pat_mean=mean(part3(region2));
pat_regions(3,2,1)=pat_mean;
pat_mean=mean(part4(region2));
pat_regions(4,2,1)=pat_mean;
pat_mean=mean(part5(region2));
pat_regions(5,2,1)=pat_mean;

pat_mean=mean(part1(region3));
pat_regions(1,3,1)=pat_mean;
pat_mean=mean(part2(region3));
pat_regions(2,3,1)=pat_mean;
pat_mean=mean(part3(region3));
pat_regions(3,3,1)=pat_mean;
pat_mean=mean(part4(region3));
pat_regions(4,3,1)=pat_mean;
pat_mean=mean(part5(region3));
pat_regions(5,3,1)=pat_mean;

pat_mean=mean(part1(region4));
pat_regions(1,4,1)=pat_mean;
pat_mean=mean(part2(region4));
pat_regions(2,4,1)=pat_mean;
pat_mean=mean(part3(region4));
pat_regions(3,4,1)=pat_mean;
pat_mean=mean(part4(region4));
pat_regions(4,4,1)=pat_mean;
pat_mean=mean(part5(region4));
pat_regions(5,4,1)=pat_mean;


con_regions=zeros(8,4,1);
con_mean=mean(part6(region1));
con_regions(1,1,1)=con_mean;
con_mean=mean(part7(region1));
con_regions(2,1,1)=con_mean;
con_mean=mean(part8(region1));
con_regions(3,1,1)=con_mean;
con_mean=mean(part9(region1));
con_regions(4,1,1)=con_mean;
con_mean=mean(part10(region1));
con_regions(5,1,1)=con_mean;
con_mean=mean(part11(region1));
con_regions(6,1,1)=con_mean;
con_mean=mean(part12(region1));
con_regions(7,1,1)=con_mean;
con_mean=mean(part13(region1));
con_regions(8,1,1)=con_mean;

con_mean=mean(part6(region2));
con_regions(1,2,1)=con_mean;
con_mean=mean(part7(region2));
con_regions(2,2,1)=con_mean;
con_mean=mean(part8(region2));
con_regions(3,2,1)=con_mean;
con_mean=mean(part9(region2));
con_regions(4,2,1)=con_mean;
con_mean=mean(part10(region2));
con_regions(5,2,1)=con_mean;
con_mean=mean(part11(region2));
con_regions(6,2,1)=con_mean;
con_mean=mean(part12(region2));
con_regions(7,2,1)=con_mean;
con_mean=mean(part13(region2));
con_regions(8,2,1)=con_mean;

con_mean=mean(part6(region3));
con_regions(1,3,1)=con_mean;
con_mean=mean(part7(region3));
con_regions(2,3,1)=con_mean;
con_mean=mean(part8(region3));
con_regions(3,3,1)=con_mean;
con_mean=mean(part9(region3));
con_regions(4,3,1)=con_mean;
con_mean=mean(part10(region3));
con_regions(5,3,1)=con_mean;
con_mean=mean(part11(region3));
con_regions(6,3,1)=con_mean;
con_mean=mean(part12(region3));
con_regions(7,3,1)=con_mean;
con_mean=mean(part13(region3));
con_regions(8,3,1)=con_mean;

con_mean=mean(part6(region4));
con_regions(1,4,1)=con_mean;
con_mean=mean(part7(region4));
con_regions(2,4,1)=con_mean;
con_mean=mean(part8(region4));
con_regions(3,4,1)=con_mean;
con_mean=mean(part9(region4));
con_regions(4,4,1)=con_mean;
con_mean=mean(part10(region4));
con_regions(5,4,1)=con_mean;
con_mean=mean(part11(region4));
con_regions(6,4,1)=con_mean;
con_mean=mean(part12(region4));
con_regions(7,4,1)=con_mean;
con_mean=mean(part13(region4));
con_regions(8,4,1)=con_mean;

pat_cols=[1 3 5 7; 1 3 5 7; 1 3 5 7; 1 3 5 7; 1 3 5 7];
con_cols=[2 4 6 8; 2 4 6 8; 2 4 6 8; 2 4 6 8; 2 4 6 8; 2 4 6 8; 2 4 6 8; 2 4 6 8];

figure;
scatter(pat_cols(:,1),pat_regions(:,1),'filled','k');
hold on;
scatter(pat_cols(:,2),pat_regions(:,2),'filled','k');
scatter(pat_cols(:,3),pat_regions(:,3),'filled','k');
scatter(pat_cols(:,4),pat_regions(:,4),'filled','k');
scatter(con_cols(:,1),con_regions(:,1),36,'k');
scatter(con_cols(:,2),con_regions(:,2),36,'k');
scatter(con_cols(:,3),con_regions(:,3),36,'k');
scatter(con_cols(:,4),con_regions(:,4),36,'k');

all_regions=zeros(8,8,1);
where_zeros=find(all_regions==0);
all_regions(where_zeros)=NaN;
all_regions(1:5,1,1)=pat_regions(1:5,1,1);
all_regions(1:5,3,1)=pat_regions(1:5,2,1);
all_regions(1:5,5,1)=pat_regions(1:5,3,1);
all_regions(1:5,7,1)=pat_regions(1:5,4,1);

all_regions(1:8,2,1)=con_regions(1:8,1,1);
all_regions(1:8,4,1)=con_regions(1:8,2,1);
all_regions(1:8,6,1)=con_regions(1:8,3,1);
all_regions(1:8,8,1)=con_regions(1:8,4,1);

hold on;
boxplot(all_regions);
axis square;

title({'IRL Thicknesses by Group and Region'});
xlabel('Region');
ylabel('Thickness (pixel)');

% To compare the regions using t-tests and return the p-values
all_thicknessIRLpatreg1(1,1,1)=mean(part1(region1));
all_thicknessIRLpatreg1(2,1,1)=mean(part2(region1));
all_thicknessIRLpatreg1(3,1,1)=mean(part3(region1));
all_thicknessIRLpatreg1(4,1,1)=mean(part4(region1));
all_thicknessIRLpatreg1(5,1,1)=mean(part5(region1));
all_thicknessIRLconreg1(1,1,1)=mean(part6(region1));
all_thicknessIRLconreg1(2,1,1)=mean(part7(region1));
all_thicknessIRLconreg1(3,1,1)=mean(part8(region1));
all_thicknessIRLconreg1(4,1,1)=mean(part9(region1));
all_thicknessIRLconreg1(5,1,1)=mean(part10(region1));
all_thicknessIRLconreg1(6,1,1)=mean(part11(region1));
all_thicknessIRLconreg1(7,1,1)=mean(part12(region1));
all_thicknessIRLconreg1(8,1,1)=mean(part13(region1));
[h,p]=ttest2(all_thicknessIRLpatreg1,all_thicknessIRLconreg1)

all_thicknessIRLpatreg2(1,1,1)=mean(part1(region2));
all_thicknessIRLpatreg2(2,1,1)=mean(part2(region2));
all_thicknessIRLpatreg2(3,1,1)=mean(part3(region2));
all_thicknessIRLpatreg2(4,1,1)=mean(part4(region2));
all_thicknessIRLpatreg2(5,1,1)=mean(part5(region2));
all_thicknessIRLconreg2(1,1,1)=mean(part6(region2));
all_thicknessIRLconreg2(2,1,1)=mean(part7(region2));
all_thicknessIRLconreg2(3,1,1)=mean(part8(region2));
all_thicknessIRLconreg2(4,1,1)=mean(part9(region2));
all_thicknessIRLconreg2(5,1,1)=mean(part10(region2));
all_thicknessIRLconreg2(6,1,1)=mean(part11(region2));
all_thicknessIRLconreg2(7,1,1)=mean(part12(region2));
all_thicknessIRLconreg2(8,1,1)=mean(part13(region2));
[h,p]=ttest2(all_thicknessIRLpatreg2,all_thicknessIRLconreg2)

all_thicknessIRLpatreg3(1,1,1)=mean(part1(region3));
all_thicknessIRLpatreg3(2,1,1)=mean(part2(region3));
all_thicknessIRLpatreg3(3,1,1)=mean(part3(region3));
all_thicknessIRLpatreg3(4,1,1)=mean(part4(region3));
all_thicknessIRLpatreg3(5,1,1)=mean(part5(region3));
all_thicknessIRLconreg3(1,1,1)=mean(part6(region3));
all_thicknessIRLconreg3(2,1,1)=mean(part7(region3));
all_thicknessIRLconreg3(3,1,1)=mean(part8(region3));
all_thicknessIRLconreg3(4,1,1)=mean(part9(region3));
all_thicknessIRLconreg3(5,1,1)=mean(part10(region3));
all_thicknessIRLconreg3(6,1,1)=mean(part11(region3));
all_thicknessIRLconreg3(7,1,1)=mean(part12(region3));
all_thicknessIRLconreg3(8,1,1)=mean(part13(region3));
[h,p]=ttest2(all_thicknessIRLpatreg3,all_thicknessIRLconreg3)

all_thicknessIRLpatreg4(1,1,1)=mean(part1(region4));
all_thicknessIRLpatreg4(2,1,1)=mean(part2(region4));
all_thicknessIRLpatreg4(3,1,1)=mean(part3(region4));
all_thicknessIRLpatreg4(4,1,1)=mean(part4(region4));
all_thicknessIRLpatreg4(5,1,1)=mean(part5(region4));
all_thicknessIRLconreg4(1,1,1)=mean(part6(region4));
all_thicknessIRLconreg4(2,1,1)=mean(part7(region4));
all_thicknessIRLconreg4(3,1,1)=mean(part8(region4));
all_thicknessIRLconreg4(4,1,1)=mean(part9(region4));
all_thicknessIRLconreg4(5,1,1)=mean(part10(region4));
all_thicknessIRLconreg4(6,1,1)=mean(part11(region4));
all_thicknessIRLconreg4(7,1,1)=mean(part12(region4));
all_thicknessIRLconreg4(8,1,1)=mean(part13(region4));
[h,p]=ttest2(all_thicknessIRLpatreg4,all_thicknessIRLconreg4)

%% Generating all thicknesses for all participants
% First run "Generating all surfaces for all participants"
% Generate all_thickness1

