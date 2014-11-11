%% Cirrus_OCT_runALL
% This script runs all participants listed in "pathlist" to generate the
% following:
%   all_thicknessIRL
%   con_thicknessIRLMean, con_thicknessIRLStd, pat_thicknessIRLMean,
%       pat_thicknessIRLStd
%   zscore
% First, run the following:
Cirrus_OCT_pathlist;

%% Generate thicknessIRL for all participants in pathlist
% Format: all_thicknessIRL(:,:,participant rank in pathlist)
all_thicknessIRL = zeros(512,512,length(pathlist));
for iPath = 1:length(pathlist);
    folder = fullfile(basepath,pathlist{iPath});
    all_thicknessIRL(:,:,iPath) = Cirrus_OCT_thicknessIRL(folder);
end

% Flip right eye participants to left eye orientation
for participant = [3 4 5 9 10 12 13]
    all_thicknessIRL(:,:,participant)=fliplr(all_thicknessIRL(:,:,participant));
end

%% Generate mean thicknessIRL and standard deviation for each group
% Adjust range in Z-dimension to include appropriate controls/patients.
% Format: con_thicknessIRLMean/con_thicknessIRLStd/pat_thicknessIRLMean/
%         pat_thicknessIRLStd(512,512,1)
con_thicknessIRLMean = mean(all_thicknessIRL(:,:,6:13),3);
con_thicknessIRLStd = std(all_thicknessIRL(:,:,6:13),0,3);
pat_thicknessIRLMean = mean(all_thicknessIRL(:,:,1:5),3);
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

%% MISCELLANEOUS

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

figure;
imagesc(con_thicknessIRLMean(:,:));
colorbar;
[cmin,cmax] = caxis;
caxis([0,100]);
title({'Group: Control'});
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
ylabel(colorbar,'Thickness (pixel)');

figure;
imagesc(pat_thicknessIRLMean(:,:));
colorbar;
[cmin,cmax] = caxis;
caxis([0,100]);
title({'Group: T1D'});
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
ylabel(colorbar,'Thickness (pixel)');

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

figure;
imagesc(diff_thicknessIRLMean(:,:));
colorbar;
[cmin,cmax] = caxis;
caxis([-15,15]);
title({'Averaged Control - Averaged T1D'});
xlabel('(Nasal)               Pixel Location               (Temporal)');
ylabel('(Inferior)               Pixel Location               (Superior)');
ylabel(colorbar,'Thickness (pixel)');

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

a=zeros(512,512,3);
b=zeros(512,512);
b(con_thicknessIRLMean>45 & con_thicknessIRLMean<50)=255;
c=zeros(512,512);
c(pat_thicknessIRLMean>45 & pat_thicknessIRLMean<50)=255;
a(:,:,2)=b;
a(:,:,1)=c;
figure;image(a);

histogram(diff_thicknessIRLMean);
title({'Averaged Control Minus Averaged T1D IRL Thicknesses'});
xlabel('Thickness (pixel)');
ylabel('Frequency');

qqplot(diff_thicknessIRLMean(:));
title({'Averaged Control Minus Averaged T1D IRL Thicknesses'});

hvalue_thicknessIRL = zeros(512,512,1);
pvalue_thicknessIRL = zeros(512,512,1);
for iRow = 1:512;
    for iCol = 1:512;
        [h,p] = ttest2(all_thicknessIRL(iRow,iCol,1:5),all_thicknessIRL(iRow,iCol,6:13));
        hvalue_thicknessIRL(iRow,iCol,1) = h;
        pvalue_thicknessIRL(iRow,iCol,1) = p;
    end
end

hvalue_thicknessIRL = zeros(512,512,1);
pvalue_thicknessIRL = zeros(512,512,1);
for iRow = 1:512;
    for iCol = 1:512;
        [h,p] = ttest2(all_thicknessIRL(iRow,iCol,1:5),all_thicknessIRL(iRow,iCol,6:13));
        hvalue_thicknessIRL(iRow,iCol,1) = h;
        pvalue_thicknessIRL(iRow,iCol,1) = p;
    end
end

qqplot(pvalue_thicknessIRL(:));
title({'P-values Comparing Controls and T1D'});

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

%% Miscellaneous 
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