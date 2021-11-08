%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
%
% Script for converting Amira output files into .txt coordinate
% files for computational toolbox compatibility
%
% Written by Florian Fässler, IST Austria; florian.fAessler@ist.ac.at

% 05/2021
% Please use MATLAB version 2020 or higher
%--------------------------------------------------------------------------

%Input area

%place the files extracted from amira in the same directory as the script

%Type the prefix/basename of your file (eg. TS_01; TS_02, etc.)
name='';
%Your point file should be exported as tab limitted file from amira 
%containing only the header, pointIDs, and point coordinates in X,Y, and Z
%select the suffix/extension name of the points file:
points=horzcat(name,'points.txt');
%Set according to your convention of contour file naming
%Your segment file should be exported as tab limitted file from amira
%containing only the header and pointIDs
%select the suffix/extension name of the segments file: 
contours=horzcat(name,'segments.txt');
%Select the desired name/extension of output file:
output=horzcat(name,'contours_and_points.txt');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

contour_and_points=importdata(points);
contours_formatted=importdata(contours);
out_put_contour_and_points=[];
contours_formatted=contours_formatted(3:length(contours_formatted));
contour_and_points=struct2cell(contour_and_points);
contour_and_points=contour_and_points{1};
for i = 1:length(contours_formatted)
    contours_formatted{i,1}=strcat(',',extractAfter(extractBefore(contours_formatted{i,1}, strlength(contours_formatted{i,1})),1),',');
end
for i = 1:length(contour_and_points(:,1))
    a=contour_and_points(i,1);
    b=num2str(a);
    c=strcat(',',b,',');
    d=find(contains(contours_formatted,c));
    if ~isempty(d)
        contour_and_points(i,1)=d;
        out_put_contour_and_points=vertcat(out_put_contour_and_points, contour_and_points(i,:));
    end
end

writematrix(out_put_contour_and_points,output,'Delimiter','tab');