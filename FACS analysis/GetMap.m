function [map,path]=GetMap()

[file_name,path]=uigetfile('.xlsx','Select One Map File');
map=readtable(fullfile(path,file_name));