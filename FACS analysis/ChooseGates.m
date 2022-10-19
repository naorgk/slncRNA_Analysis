function []=ChooseGates()

% classify objects using 2D scatter plot
clear all; close all; fclose all;


% % read/write dir
% wdir=[pwd '\'];%'/Users/sarah/Dropbox/matlab/fcs_read/or_2019/';
% % for defining gates for 200nm beads in all 561nm laser intensities (0, 50%, 100%)
% fname='RC2019-11-27_A1_CoT I 48h.0001.fcs';
[file_name,path]=uigetfile('.fcs','Select One Control File');

% convert fcs data to table so channel names are linked to data
% read .txt of exported features
[fdat, fhdr, ~]=fca_readfcs(fullfile(path,file_name));
colnames={fhdr.par.name}; %type of channel
colnames=strrep(colnames,'-','_'); %matlab doesnt like - in names.
% convert to table
ftab=array2table(fdat, 'VariableNames', colnames);

% colors for ROIs
cols={'b','g','m','c','r'};
gates = struct('Gate_name',{},'Object_type',{},'Gate_offset',{},'Gate_scale',{},'Gate_channel_X',{},'Gate_channel_Y',{},'Gate_polygon',{});

% gate file, new gates are appended so delete/rename if you want to create
% a new set

ls('gate_file*')
gate_file=[pwd '\gate_file_'  datestr(now,'yyyymmdd_HHMM') '.mat'];

% optional offset
offset=0; 
scale_str='log';
ch_pairs={{'FSC_A','SSC_A'}; {'SSC_A','SSC_H'}}; %live and single, respectively.

% set up plots
n_plots=numel(ch_pairs); %number of elements in an array of unknown size.
n_cols=min(n_plots,4);
n_rows=ceil(n_plots/n_cols);

% loop for choosing gates
cur_gate=0;
ah=[]; % axes handles

figure(1);clf;hold on;%set(gcf,'position',[680         403        1035         575]);

% plot the 2D data on requested channel pairs chpairs
for ii=1:length(ch_pairs)
    sh(ii)=subplot(n_rows,n_cols,ii);
    hold on;
    ch1=ch_pairs{ii}{1};
    ch2=ch_pairs{ii}{2};
    xlabel(ch1,'interpreter','none');
    ylabel(ch2,'interpreter','none');
    scatter(ftab{:,ch1},ftab{:,ch2},'.k');
    hold on;
    ah(ii)=gca;
    lh=legend(file_name);set(lh,'interpreter','none');
end

%% choose n gates (n limits loop interations, save is currently outside of loop) - runs on 1 filename
% gate definition
n_gates=2;
if(n_gates>0)
    while(cur_gate<n_gates)
        cur_gate=cur_gate+1;
        if(cur_gate==1)
            subplot(sh(1)); hold on;
            th=text(-0.05,1.05,'Choose subplot for gate definition (click)','units','normalized');
        end
        
        % not needed for first iteration
        subplot(sh(1));hold on;
        set(th,'string','Choose subplot for gate definition (click)','units','normalized');
        
        % choose subplot that we are using for gate
        if(length(ch_pairs)>1)
            n_tries=1;cur_try=0;
            while(cur_try<n_tries)
                cur_try=cur_try+1;
                [xchoice,ychoice]=ginput(1); % makes axes active !!!
                choice=gca;
                ind_choice=find(ah==gca);
                set(th,'string',['ACTIVE subplot = ' num2str(ind_choice(1)) ', choose zoom (press any key when zoom is good)']); 
            end
        end
        
        % zoom with your mouse. when image is okay, press any key
        zoom on;
        pause(); 
        zoom off; % to escape the zoom mode
        set(th,'string',['ACTIVE subplot = ' num2str(ind_choice(1)) ', define polygon with mouse (double click to close)']);
        
        % get ROI from subplot of interest
        subplot(sh(ind_choice));hold on;
        poly = impoly;
        nodes = getPosition(poly);
        
        % adjust nodes
        set(th,'string',['ACTIVE subplot = ' num2str(ind_choice(1)) ', adjust polygon with mouse (click on cmd window and press ENTER when poly is good)']);
        pause();
        nodes = getPosition(poly);
        
        chX=ch_pairs{ind_choice}{1};
        chY=ch_pairs{ind_choice}{2};
        
        % get points from this subplot in this ROI
        insiders=inpolygon(ftab{:,chX}, ftab{:,chY},nodes(:,1),nodes(:,2));
        %ftab_live=ftab(insiders,:); %removing dead cells
        if(sum(insiders)>0)
            % color all points selected on all other subplots
            for ii=1:numel(ch_pairs)
                % color points inside polygon for this plot
                subplot(n_rows,n_cols,ii);hold on;
                ch1=ch_pairs{ii}{1};
                ch2=ch_pairs{ii}{2};
                scatter(ftab{insiders,ch1},ftab{insiders,ch2},[],cols{cur_gate});
            end
        end
        
        % add to gate list: struct('Object_type',{},'Gate_offset',{},'Gate_scale',{},'Gate_channel_X',{},'Gate_channel_Y',{},'Gate_polygon',{});
%         chX=ch_pairs{ind_choice}{1}; 
%         chY=ch_pairs{ind_choice}{2};
        gate_name=['Gate_' num2str(cur_gate)];
        gates(cur_gate)=struct('Gate_name',gate_name,'Object_type',file_name,'Gate_offset',offset,'Gate_scale',scale_str,'Gate_channel_X',chX,'Gate_channel_Y',chY,'Gate_polygon',nodes);
        
    end % gate while
end

%% save all gates - append to existing
if(length(gates)>0)
    new_gates=gates; % new
    if(exist(gate_file,'file')==2)
        load(gate_file,'gates'); % existing
        gates=[gates,new_gates];
    else
        gates=new_gates;
    end
    save(gate_file,'gates');
end

%% load and plot gates - assume all in same scale as gate(1), with same offset
gates={};
if(exist(gate_file,'file')==2)
    load(gate_file,'gates');
else
    return;
end

figure(2); clf; hold on;
offset=gates(1).Gate_offset;
% plot data
for ii=1:length(ch_pairs)
    sh(ii)=subplot(n_rows,n_cols,ii);hold on;
    ch1=ch_pairs{ii}{1};
    ch2=ch_pairs{ii}{2};
    xlabel(ch1,'interpreter','none');
    ylabel(ch2,'interpreter','none');
    scatter(ftab{:,ch1},ftab{:,ch2},'.k');
    %axis([-1e3 1.6e4 -1e3 1.8e4]);
end
subplot(sh(1));hold on;
text(-0.25,1.05,['ALL GATES for ' gate_file ' on ' file_name],'units','normalized','interpreter','none');  
% plot gates
ghs=[];
for gg=1:length(gates)
    chX_str=gates(gg).Gate_channel_X; % X,Y are gate channels
    chY_str=gates(gg).Gate_channel_Y;
    nodes = gates(gg).Gate_polygon;
    insiders=inpolygon(ftab{:,chX_str},ftab{:,chY_str},nodes(:,1),nodes(:,2));
   
    for ii=1:length(ch_pairs)
        subplot(n_rows,n_cols,ii);hold on;
        ch1=ch_pairs{ii}{1};
        ch2=ch_pairs{ii}{2};
       
        % color all points slected on all subplots
        if(sum(insiders)>0)
            sh=scatter(ftab{insiders,ch1},ftab{insiders,ch2},[],cols{mod(gg,length(cols))+1});
        end
       
        % plot gate (gate gg will appear only in one plot)
        if(strcmp(ch1,chX_str) & strcmp(ch2,chY_str))
            gh=plot([nodes(:,1)' nodes(1,1)]',[nodes(:,2)' nodes(1,2)]','linewidth',3,'color',cols{mod(gg,length(cols))+1});
            ghs(gg)=gh;
        end
    end
end
% legend
for ii=1:length(ch_pairs)
    subplot(n_rows,n_cols,1);hold on;
    lh=legend(ghs,gates.Gate_name);
    set(lh,'interpreter','none');
end

%% EOF
