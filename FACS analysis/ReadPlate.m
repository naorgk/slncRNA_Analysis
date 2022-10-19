clear all
close all

%make sure to have only one .mat file in the folder!!!
gate_file=dir('*.mat');
gates={};
if(exist(gate_file.name,'file')==2)
    load(gate_file.name,'gates');
else
    ChooseGates();
    gate_file=dir('*.mat');
    load(gate_file.name,'gates');
end

cols={'blue','green','magenta','cyan','red','black','yellow','#FFA500'};

% build ch_pairs for this gate file
ch1s={gates.Gate_channel_X};
ch2s={gates.Gate_channel_Y};
ch_pairs=cellstr([ch1s; ch2s]');
ch_pairs={ch_pairs(1,:),ch_pairs(2,:)};

% set up plots
n_plots=numel(ch_pairs); %number of elements in an array of unknown size.
n_cols=min(n_plots,4);
n_rows=ceil(n_plots/n_cols);

%get plate labels from map
[PlateMap,path]=GetMap();
populations=unique(PlateMap{:,'POPULATION'},'stable');
conditions=unique(PlateMap{:,'CONDITION'},'stable');

%select files and get path
file_names=GetFiles(path);

live=zeros(size(file_names,1),size(file_names,2));
binEdges=round(0:0.05:5,2);
for dd=1:size(file_names,2) %run on all folders
    bins_mCherry=zeros(size(file_names,1),numel(binEdges)-1);
    for ff=1:size(file_names,1) %run on all the files
        fname=file_names{ff,dd}{:};
        plate_sign=regexp(fname,'[^_]+_([A-H]{1,1}\d{1,2})_[^_]+','tokens');
        plate_sign=plate_sign{1}{1};
        fname=(fullfile(path,file_names.Properties.VariableNames{dd},fname));
        % [^_] = not underscore
        % + = once or more
        % _ = looks for underscore
        % ([A-H]{1,1}\d{1,2}) = what we are looking for (looking for A-H once
        % exactly, \d{1,2} = once ore twice of digits)
        
        [fdat, fhdr, ~]=fca_readfcs(fname);
        colnames={fhdr.par.name}; %type of channel
        colnames=strrep(colnames,'-','_'); %matlab doesnt like - in names.
        % convert to table
        ftab=array2table(fdat, 'VariableNames', colnames);
        ghs=[];
        for gg=1:length(gates)
            chX_str=gates(gg).Gate_channel_X; % X,Y are gate channels
            chY_str=gates(gg).Gate_channel_Y;
            nodes = gates(gg).Gate_polygon;
            insiders=inpolygon(ftab{:,chX_str},ftab{:,chY_str},nodes(:,1),nodes(:,2));
            
            %add columns for each gate, denoting insiders: live and single
            tmp=array2table(insiders);
            tmp.Properties.VariableNames={gates(gg).Gate_name};
            ftab=[ftab tmp];
        end
        
        live(ff,dd)=sum(ftab.Gate_1==1)/size(ftab,1);
        
        
        pop=PlateMap{strcmp(PlateMap{:,"WELL"},plate_sign),"POPULATION"}{:};
        cond=PlateMap{strcmp(PlateMap{:,"WELL"},plate_sign),"CONDITION"}{:};
        switch(pop)
            case "E. coli"
                kk=1;
        end %switch
        switch(cond)
            case "PP7-mch C4HSL-"
                row_factor=0;
            case "PP7-mch C4HSL+"
                row_factor=1;
            case "Qb-mch C4HSL-"
                row_factor=2;
            case "Qb-mch C4HSL+"
                row_factor=3;
            case "PP7-4x_Qb-5x+PP7-mch IPTG- C4HSL-"
                row_factor=4;
            case "PP7-4x_Qb-5x+PP7-mch IPTG+ C4HSL-"
                row_factor=5;
            case "PP7-4x_Qb-5x+PP7-mch IPTG- C4HSL+"
                row_factor=6;
            case "PP7-4x_Qb-5x+PP7-mch IPTG+ C4HSL+"
                row_factor=7;
            case "PP7-4x_Qb-5x+Qb-mch IPTG- C4HSL-"
                row_factor=8;
            case "PP7-4x_Qb-5x+Qb-mch IPTG+ C4HSL-"
                row_factor=9;
            case "PP7-4x_Qb-5x+Qb-mch IPTG- C4HSL+"
                row_factor=10;
            case "PP7-4x_Qb-5x+Qb-mch IPTG+ C4HSL+"
                row_factor=11;
            case "Qb-10x+Qb-mch IPTG- C4HSL-"
                row_factor=12;
            case "Qb-10x+Qb-mch IPTG+ C4HSL-"
                row_factor=13;
            case "Qb-10x+Qb-mch IPTG- C4HSL+"
                row_factor=14;
            case "Qb-10x+Qb-mch IPTG+ C4HSL+"
                row_factor=15;
            case "PP7-24x+PP7-mch IPTG- C4HSL-"
                row_factor=16;
            case "PP7-24x+PP7-mch IPTG+ C4HSL-"
                row_factor=17;
            case "PP7-24x+PP7-mch IPTG- C4HSL+"
                row_factor=18;
            case "PP7-24x+PP7-mch IPTG+ C4HSL+"
                row_factor=19;

        end %switch
        kk=kk+size(populations,1)*row_factor;
        
        figure(1)
        subplot(5,4,kk)
        hold on
        h=histogram(log10(ftab{ftab.Gate_1==1 & ftab.Gate_2==1 & ftab.Y2_A>1,"Y2_A"}),binEdges,'Normalization','probability','DisplayStyle','stairs','LineWidth',3);
        bins_mCherry(ff,:)=h.Values;
        xlabel("log10(mCherry)");
        ylabel("Frequency");
        ylim([0,0.1]);
        XT = get(gca, 'XTick');
        title(strcat(cond," in ",pop));

    end %for ff
end %for dd




figure;
ax = axes;
plot(binEdges(1:end-1),mean(bins_mCherry(1:3,:),1),'LineWidth',2); hold on;
plot(binEdges(1:end-1),mean(bins_mCherry(4:6,:),1),'LineWidth',2);
plot(binEdges(1:end-1),mean(bins_mCherry(19:21,:),1),'b+','LineWidth',2); 
plot(binEdges(1:end-1),mean(bins_mCherry(22:24,:),1),'r*','LineWidth',2); 
plot(binEdges(1:end-1),mean(bins_mCherry(55:57,:),1),'k+','LineWidth',2); 
plot(binEdges(1:end-1),mean(bins_mCherry(58:60,:),1),'g*','LineWidth',2); 
legend('PP7 mCherry IPTG+ C4HSL-','PP7 mCherry IPTG+ C4HSL+','PP7-4x-Qb-5x+PP7-mCherry IPTG- C4HSL+',...
    'PP7-4x-Qb-5x+PP7-mCherry IPTG+ C4HSL+','PP7-24x+PP7-mCherry IPTG- C4HSL+',...
    'PP7-24x+PP7-mCherry IPTG+ C4HSL+','Location','NW','FontSize',15);
xlabel("log10(mCherry)");
ylim([0,0.09]);
ylabel("Frequency");

