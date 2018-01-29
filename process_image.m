% clear;
% clc;
% close all;
% chos=0;
% chos=menu('Select an image to process','Exit');
% possibility = 2;
% while chos~=possibility,
%     if chos==1
%         chos=menu('Fingerprint Recognition System','Select image and add to database','Select image for fingerprint recognition','Info','Delete database',...
%             'Fingerprint image: visualization','Gabor Filter: visualization','Exit');
%         %--------------------------------------------------------------------------
%         %--------------------------------------------------------------------------
%         %--------------------------------------------------------------------------
%         % Calculate FingerCode and Add to Database
%         if chos==1
%             clc;
%             close all;
%             selezionato=0;
%             while selezionato==0
%                 [namefile,pathname]=uigetfile({'*.bmp;*.tif;*.tiff;*.jpg;*.jpeg;*.gif','IMAGE Files (*.bmp,*.tif,*.tiff,*.jpg,*.jpeg,*.gif)'},'Chose GrayScale Image');
%                 if namefile~=0
%                     [img,map]=imread(strcat(pathname,namefile));
%                     selezionato=1;
%                 else
%                     disp('Select a grayscale image');
%                 end
%                 if (any(namefile~=0) && (~isgray(img)))
%                     disp('Select a grayscale image');
%                     selezionato=0;
%                 end
%             end
%         end
%     end
% end
% 
%         
%     
%     
% 
% pathname = 'C:\Users\User\Documents\MATLAB\fingercode\test\f0001_01.png';
% img=imread(pathname);
% % figure
% % imshow(img);
% % title('Original');
% resizedimg = imresize(img, 0.5);
% % imshow(resizedimg);
% % title('Resized image');



destinationfolder = 'C:\Users\User\Documents\MATLAB\fingercode\processed_fp\';

% Specify the folder where the files live.
%myFolder = 'C:\Users\User\Desktop\NISTSpecialDatabase4GrayScaleImagesofFIGS\NISTSpecialDatabase4GrayScaleImagesofFIGS\sd04\png_txt\figs_0';
myFolder = 'C:\Users\User\Desktop\NISTSpecialDatabase4GrayScaleImagesofFIGS\NISTSpecialDatabase4GrayScaleImagesofFIGS\sd04\png_txt';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isdir(myFolder)
  errorMessage = sprintf('Error: The following folder does not exist:\n%s', myFolder);
  uiwait(warndlg(errorMessage));
  return;
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '**\*.png'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
  fullFileName = strcat(theFiles(k).folder, '\', baseFileName);
  %fullFileName = fullfile(myFolder, baseFileName);
  % Now do whatever you want with this file name,
  % such as reading it in as an image array with imread()
  img = imread(fullFileName);
  resizedimg = imresize(img, 0.5);
  newFileName = strrep(baseFileName,'.png','.bmp');
  imwrite(resizedimg, strcat(destinationfolder, newFileName));
end
