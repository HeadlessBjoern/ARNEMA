% Connect to automagic and load project state file

mypath = pwd;
apath = '/Volumes/methlab/4marius_bdf/automagic-master';
addAutomagicPaths
cd(mypath)
load('/Volumes/methlab/Students/Arne/MA/data/automagic_opticat_hp/project_state.mat')