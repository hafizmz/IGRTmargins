clear,clc;close all;
mainfolder = uigetdir;
subfolders = dir(mainfolder);
subfolders = subfolders([subfolders.isdir] & ~startsWith({subfolders.name},'.'));
for i = 1:length(subfolders)
    x = []; y = []; z = []; 
    files = dir(fullfile(mainfolder,subfolders(i).name,'*.xvi'));
    for k = 1:length(files);
        fid = fopen(fullfile(files(k).folder,files(k).name),'r'); 
              while ~feof(fid)
              st = fgetl(fid);           
               if  ~isempty(strfind(st,'CouchShiftLat'))
                 stread = textscan(st,'%s %f','Delimiter',';=');
                  x = [x; stread{2}(1)];
               elseif  ~isempty(strfind(st,'CouchShiftLong'))
                 stread = textscan(st,'%s %f','Delimiter',';=');
                  y = [y; stread{2}(1)];
               elseif ~isempty(strfind(st,'CouchShiftHeight'))
                 stread = textscan(st,'%s %f','Delimiter',';=');
                 z = [z; stread{2}(1)];
                 nonIGRTcm=[x,y,z];            %nonIGRT implies image registration (raw data as extracted from the text file)
                 nonIGRT=nonIGRTcm*10;         %convert the nonIGRT from cm to mm

               end
              end
        fclose(fid);
    end
    %save all the x,y,z error data for each patient as .mat file and save it in the main folder
 save([mainfolder num2str(i) '.mat'],'nonIGRT')
end
for ii = 1:i;
 load([mainfolder num2str(ii) '.mat']) %Load all the .mat files from the main folder
%nonIGRT Protocol
indv_SSE_nonIGRT=mean(nonIGRT);       %Individual systematic setup error  (indv_SSE)
indv_RSE_nonIGRT=std(nonIGRT);        %Individual random setup error  (indv_RSE)
SSE_nonIGRT(ii,:)=indv_SSE_nonIGRT;   %collect indv_SSE for all patients into matrix
RSE_nonIGRT(ii,:)=indv_RSE_nonIGRT;   %collect indv_RSE for all patients into matrix
PoP_SSE_nonIGRT= std(SSE_nonIGRT)     %population sysematic setup error
PoP_RSE_nonIGRT= rms(RSE_nonIGRT)     %population random setup error 
PTV_nonIGRT=2.5*PoP_SSE_nonIGRT+0.7*PoP_RSE_nonIGRT %PTV using van Herk margin recipe

%NAL protocol
First3fxs=nonIGRT(1:3,:);             %Extract error values from the 1st 3 fractions
Correction_NAL=mean(First3fxs);       %cal. average of the error values from the 1st 3 fraction as correction for NAL protocol
Lastfxs=nonIGRT(4:end,:);             %Extract error values from the 4th fraction to the last fraction  
PC_NAL=Lastfxs-Correction_NAL;        %Apply the NAL carrection values
Residual_NAL=[First3fxs;PC_NAL];      %NAL Resisual setup error
indv_SSE_NAL=mean(Residual_NAL);      %NAL individual sysematic setup error
indv_RSE_NAL=std(Residual_NAL);       %NAL individual random setup error
SSE_NAL(ii,:)=indv_SSE_NAL;           %NAL individual sysematic setup error for all patients
RSE_NAL(ii,:)=indv_RSE_NAL;           %NAL individual random setup error for all patients
PoP_SSE_NAL= std(SSE_NAL)             %NAL population sysematic setup error 
PoP_RSE_NAL= rms(RSE_NAL)             %NAL population random setup error
PTV_NAL=2.5*PoP_SSE_NAL+0.7*PoP_RSE_NAL %NAL PTV margin using van Herk margin recipe

%eNAL protocol
% Runing average down the columns.The output is a mtrix in which its 1st row is average of the 1st and 2nd row of original matrix...
%the 2nd row is the average of the 1st-3rd row of the original matrix..until the last row
RunAverage=cell2mat(arrayfun(@(n)mean(nonIGRT(1:n,:)),[2:size(nonIGRT,1)].','uni',0)); 
Correction_eNAL=RunAverage(2:end-1,:);  %correction values for eNAL protocol                                    
PC_eNAL=Lastfxs-Correction_eNAL;        %Apply the correction to the 4th fraction until the last fraction 
Residual_eNAL=[First3fxs;PC_eNAL];      %eNAL residual setup error.
indv_SSE_eNAL=mean(Residual_eNAL);      %eNAL individual sysematic setup error
indv_RSE_eNAL=std(Residual_eNAL) ;      %eNAL individual random setup error
SSE_eNAL(ii,:)=indv_SSE_eNAL;           %eNAL individual sysematic setup error for all patients
RSE_eNAL(ii,:)=indv_RSE_eNAL;           %eNAL individual random setup error for all patients
PoP_SSE_eNAL= std(SSE_eNAL)             %eNAL population sysematic setup error 
PoP_RSE_eNAL= rms(RSE_eNAL)             %eNAL population random setup error
PTV_eNAL=2.5*PoP_SSE_eNAL+0.7*PoP_RSE_eNAL %eNAL PTV margin using van Herk margin recipe

% Online protocol using 3 mm action level
% nonIGRT(any(nonIGRT>3,2),:)=0;        %Replace all the values in a row with zero if any of the element in the row is >3  
% nonIGRT(any(nonIGRT<-3,2),:)=0;       %Replace all the values in a row with zero if any of the element in the row is <-3
% indv_SSE_Online=mean(nonIGRT);        %Online protocol individual sysematic setup error
% indv_RSE_Online=std(nonIGRT);         %Online protocol individual random setup error 
% SSE_Online(ii,:)=indv_SSE_Online;     %Online protocol individual sysematic setup error for all patients
% RSE_Online(ii,:)=indv_RSE_Online;     %Online protocol individual random setup error for all patients
% PoP_SSE_Online= std(SSE_Online);      %Online protocol population systematic setup error 
% PoP_RSE_Online= rms(RSE_Online);      %Online protocol population random setup error
% PTV_Online=2.5*PoP_SSE_Online+0.7*PoP_RSE_Onlin; %Online protocol PTV margin using van Herk margin recipe
end

%plot the results
figure;
subplot(1,2,1) % Bar plot for population  systematic error
PoP_SSE_Error=[PoP_SSE_nonIGRT(1),PoP_SSE_nonIGRT(2),PoP_SSE_nonIGRT(3)];
PoP_SSE_NAL=[PoP_SSE_NAL(1),PoP_SSE_NAL(2),PoP_SSE_NAL(3)];
PoP_SSE_eNAL=[PoP_SSE_eNAL(1),PoP_SSE_eNAL(2),PoP_SSE_eNAL(3)]
Xaxis=[PoP_SSE_Error;PoP_SSE_NAL;PoP_SSE_eNAL];
subplot(1,2,1)
b=bar(Xaxis,'FaceColor','flat');
set(gca, 'XTicklabel',{'Error', 'NAL', 'eNAL'})
for k=1:size(Xaxis,2)
    b(k).CData=k;
end
ylabel('Systematic error (mm)')
legend('x-axis','y-axis','z-axis')
title('PoPulation systematic error')
ylim([0 1.8])

subplot(1,2,2)  % Bar plot for population  random error
PoP_RSE_Error=[PoP_RSE_nonIGRT(1),PoP_RSE_nonIGRT(2),PoP_RSE_nonIGRT(3)];
PoP_RSE_NAL=[PoP_RSE_NAL(1),PoP_RSE_NAL(2),PoP_RSE_NAL(3)];
PoP_RSE_eNAL=[PoP_RSE_eNAL(1),PoP_RSE_eNAL(2),PoP_RSE_eNAL(3)];
Xaxis=[PoP_RSE_Error;PoP_RSE_NAL;PoP_RSE_eNAL];
subplot(1,2,2)
b=bar(Xaxis,'FaceColor','flat');
set(gca, 'XTicklabel',{'Error', 'NAL', 'eNAL'})
for k=1:size(Xaxis,2)
    b(k).CData=k;
end
ylabel('Random error (mm)')
legend('x-axis','y-axis','z-axis')
title('PoPulation random error')
ylim([0 1.8])

% Barplot for PTV margins
PTV_Error=[PTV_nonIGRT(1),PTV_nonIGRT(2),PTV_nonIGRT(3)];
PTV_NAL=[PTV_NAL(1),PTV_NAL(2),PTV_NAL(3)];
PTV_eNAL=[PTV_eNAL(1),PTV_eNAL(2),PTV_eNAL(3)];
Xaxis=[PTV_Error;PTV_NAL;PTV_eNAL];
figure;
subplot(1,2,1)
b=bar(Xaxis,'FaceColor','flat');
set(gca, 'XTicklabel',{'Error', 'NAL', 'eNAL'})
for k=1:size(Xaxis,2)
    b(k).CData=k;
end
ylabel('PTV margin (mm)')
legend('x-axis','y-axis','z-axis')
title('Required PTV margin')

%Calculation for percentage reduction in PTV margin achieved for various protocols relative to nonIGRT protocols
NAL=((PTV_nonIGRT-PTV_NAL)./PTV_nonIGRT)*100;
eNAL=((PTV_nonIGRT-PTV_eNAL)./PTV_nonIGRT)*100;
% Online=((PTV_nonIGRT-PTV_Online)./PTV_nonIGRT)*100; %online protocol is excluded

% Barplots for percentage reduction in PTV margin
NAL=[NAL(1),NAL(2),NAL(3)];
eNAL=[eNAL(1),eNAL(2),eNAL(3)];
Xaxis=[NAL;eNAL];
subplot(1,2,2)
b=bar(Xaxis,'FaceColor','flat');
set(gca, 'XTicklabel',{'NAL', 'eNAL'})
for k=1:size(Xaxis,2)
    b(k).CData=k;
end
ylabel('Percentage reduction in PTV (%)')
legend('x-axis','y-axis','z-axis')
title('Percentage reduction in PTV ')
clear var ans fid files st stread subfolders x y z Xaxis RunAverage Residual_NAL Residual_eNAL...
    b eNAL First3fxs i ii indv_RSE_eNAL indv_RSE_NAL indv_RSE_nonIGRT indv_SSE_eNAL indv_SSE_NAL indv_SSE_eNAL k ...
    Lastfxs NAL nonIGRT nonIGRTcm PC_eNAL PC_NAL indv_SSE_nonIGRT