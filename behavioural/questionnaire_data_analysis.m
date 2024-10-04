%% Analysis of OCC pilot participant questionnaires

clear all;
close all;
clc;

%% Load file
tbl = readtable('/Volumes/methlab_vp/OCC/OCC_ARNEMA.xlsx');
tbl = tbl(2:end, :); % Delete example participant
tbl = tbl(1:20, ["OCCID", "Geschlecht", "Alter", "H_ndigkeit", "H_ndigkeit_Fragebogen_", "Ausbildungsgrad_1_6_Andere_", "OkulareDominanz", "ParadigmVersion"]); % Only keep relevant information

%%



%% Age
% birthdays = rmmissing(tbl(:, "Geburtsdatum"));
% 
% 
% for i = 1:height(birthdays)
%     birth_numdate = datenum(string(birthdays.Geburtsdatum(i)),'DD-mm-YYYY');
%     age=datestr(now-birth_numdate,'DD-mm-YYYY')
% end