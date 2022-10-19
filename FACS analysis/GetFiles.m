function file_names=GetFiles(path)

dir_content=dir(path);
dirs=[dir_content.isdir];
number_of_directories=sum(dirs)-2;
dirs=dir_content(dirs);
dir_names=zeros(number_of_directories,1);
dir_names=num2cell(dir_names);
for dd=1:number_of_directories
    dir_names{dd}=dirs(dd+2).name;
end %for dd
dir_content=dir(fullfile(path,dir_names{dd}));
files=~[dir_content.isdir];
number_of_files=sum(files);
file_names=zeros(number_of_files,number_of_directories);
file_names=num2cell(file_names);
for dd=1:number_of_directories
    %[names,path]=uigetfile('.fcs','Select One or More Files','MultiSelect','on');
    dir_content=dir(fullfile(path,dir_names{dd}));
    dir_content=natsortfiles(dir_content);
    files=~[dir_content.isdir];
    files=dir_content(files);
    for ff=1:number_of_files
        file_names{ff,dd}=files(ff).name;
    end %for ff
end %for dd
file_names=cell2table(file_names);
file_names.Properties.VariableNames=dir_names;