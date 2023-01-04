function [varargout] = plugin_RAMP(fcn,varargin)

if nargin==0
	myfailed('Expects at least one input argument.');
	return;
end;

switch fcn
	case 'getname'
		varargout = cell(1,1);
		
		%Segment with versions >1.636 calls with two input arguments where
		%the second input argument is the handle to the menu item.
		
		uimenu(varargin{1},'Label','Mitral and Tricuspid Respiratory Variation','Callback','plugin_RAMP(''init'')');
		varargout{1} = 'Mitral and Tricuspid Respiratory Variation';
		%    set(varargin{1},'init','');
	case 'getdependencies'
		%Here: List all depending files. This is required if your plugin should
		%be possible to compile to the stand-alone version of Segment.
		varargout = cell(1,4);
		
		%M-files, list as {'hello.m' ...};
		varargout{1} = {};
		
		%Fig-files, list as {'hello.fig' ... };
		varargout{2} = {};
		
		%Mat-files, list as {'hello.mat' ... };
		varargout{3} = {};
		
		%Mex-files, list as {'hello' ...}; %Note i.e no extension!!!
		varargout{4} = {};
		
	otherwise
		macro_helper(fcn,varargin{:}); %Future use to record macros
		[varargout{1:nargout}] = feval(fcn,varargin{:}); % FEVAL switchyard
end;
end

%-----------------------
function calculate_flow()
%-----------------------
%Calculate flow data
global DATA SET RAMP 
gui = DATA.GUI.RAMP;
RAMP=[];

nom = gui.nom;
nop = gui.nop;
rois2take = gui.rois2take;
numrois = gui.numrois;

if (~gui.resamped)
	%gui.velocity = gui.velmean;
	gui.velmean = NaN(SET(nom).TSize,numrois);
	gui.velstd = gui.velmean;
	gui.velmax = gui.velmean;
	gui.velmin = gui.velmean;
	gui.kenergy = gui.velmean;
	gui.area = gui.velmean;
	gui.netflow = gui.velmean;
	gui.posflow = gui.velmean;
	%gui.negflow = gui.velmean;
	gui.phasecorrstring = '';
	if ~isempty(SET(nop).Flow.PhaseCorr)
		if ~isfield(SET(nop).Flow,'PhaseCorrTimeResolved')
			mywarning('Incompatible eddy current correction. Correction reset.',DATA.GUI.Segment);
			SET(nop).Flow.PhaseCorr = [];
			gui.phasecorrstring = '';
		else
			if SET(nop).Flow.PhaseCorrTimeResolved
				gui.phasecorrstring = '[Time-resolved eddy current compensation applied]';
			else
				gui.phasecorrstring = '[Stationary eddy current compensation applied]';
			end;
		end;
	end;
	warnedempty = false;
	a=SET(nom).TSize*2;
    h = waitbar(0,'Calculating flow.');
    set(h,'visible','off');
	myadjust(h,DATA.GUI.Segment);
	set(h,'visible','on');
    counter = 0;
    
	RAMP.ROI(1) = SET(nom).Roi(1);
	RAMP.ROI(2) = SET(nom).Roi(2);
	
	for rloop = 1:numrois
		for tloop = SET(nom).Roi(rois2take(rloop)).T
            counter = counter + 1;
            waitbar(counter/a,h);
			%Create mask
			mask = logical(segment('createmask',...
				gui.outsize,...
				SET(nom).Roi(rois2take(rloop)).Y(:,tloop),...
				SET(nom).Roi(rois2take(rloop)).X(:,tloop)));
			
			%Extract phase image
			temp = SET(nop).IM(:,:,tloop,SET(nom).Roi(rois2take(rloop)).Z);
			
			%If empty phasecorr, the do not add phase correction.
			if isempty(SET(nop).Flow.PhaseCorr)
				veldata = SET(nom).Roi(rois2take(rloop)).Sign*(temp-0.5)*2*SET(nop).VENC;
				
			else
				%Phase correction
				if SET(nop).Flow.PhaseCorrTimeResolved
					%Time resolved phase correction
					veldata = SET(nom).Roi(rois2take(rloop)).Sign*...
						(temp-0.5-SET(nop).Flow.PhaseCorr(:,:,tloop,SET(nom).Roi(rois2take(rloop)).Z))*2*SET(nop).VENC;
				else
					%Stationary phase correction
					veldata = SET(nom).Roi(rois2take(rloop)).Sign*...
						(temp-0.5-SET(nop).Flow.PhaseCorr(:,:,1,SET(nom).Roi(rois2take(rloop)).Z))*2*SET(nop).VENC;
					
				end;
			end;
							
            RAMP.velocity_unmasked(:,:,tloop) = veldata;
			
			veldata = veldata(mask);
			unaltered_vel = temp(mask);
			
			if isempty(veldata)
				if not(warnedempty)
					mywarning('Empty ROI. Please try again',DATA.GUI.Segment);
				end;
				warnedempty = true;
            else
				
				if(rloop == 1) %Mitral
					
                    RAMP.mitral_velocity(tloop,:) = unaltered_vel ;%*SET(nom).Roi(rois2take(rloop)).Sign;
                 
				end
				
				if(rloop == 2) %Tricuspid
					
					RAMP.tricuspid_velocity(tloop,:) = unaltered_vel;%*SET(nom).Roi(rois2take(rloop)).Sign;
                  				
                end;
			end;
           
		end;
		
	
		
	end;
	
	close(h);
	
	timeframes = SET(nom).StartAnalysis:SET(nom).EndAnalysis;
	
	TIncr = SET(nom).TIncr;
	TSize = SET(nom).TSize;
else
	gui.phasecorrstring = '';
	% 	StartAnalysis = 1;
	% 	EndAnalysis = length(gui.t);
	timeframes = 1:gui.tsize;
	TIncr = gui.TIncr;
	TSize = gui.tsize;
end

%DATA.updateaxestables('flow',nom,nop);

end

function ok = init(no)
global RAMP DATA SET NO 
RAMP = [];
ok = false; %inform if unable to launch GUI
if nargin<1
	no=NO;
end;

if isempty(SET(no).Flow)
	myfailed('No flow data in current image stack.',DATA.GUI.Segment);
	return;
end;

if ~isempty(SET(no).Flow)
	nom = SET(no).Flow.MagnitudeNo;
else
	return;
end;

if isempty(SET(nom).Flow.PhaseNo)
	myfailed('No through plane velocity data found.',DATA.GUI.Segment);
	return;
end;

if SET(nom).RoiN==0
	myfailed('No ROIs available.',DATA.GUI.Segment);
	return;
end;

tempnos=[SET(no).Flow.MagnitudeNo SET(no).Flow.PhaseNo];
imissingle=classcheckim(tempnos);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
	return;
end;

if isempty(SET(nom).EndAnalysis)
	SET(nom).EndAnalysis = SET(nom).TSize;
end;

if isempty(SET(nom).StartAnalysis)
	SET(nom).StartAnalysis = 1;
end;

if (SET(nom).EndAnalysis-SET(nom).StartAnalysis)<1
	myfailed('No timespan, adjust Start/End analysis time (under stress menu or drag red bar).',DATA.GUI.Segment);
	return;
end;

nop = SET(nom).Flow.PhaseNo;

if ~isopengui('plugin_RAMP.fig');
	
end;
gui = mygui('plugin_RAMP.fig');
DATA.GUI.RAMP = gui;
gui.TIncr = [];
gui.tsize = [];
gui.resamped = 0;
set(gui.fig,'Name','');
gui.outsize = [SET(nom).XSize SET(nom).YSize];

if not(isfield(SET(nop).Flow,'PhaseCorr'))
	SET(nop).Flow.PhaseCorr = [];
end;

gui.t = SET(nom).TIncr*(-1+(1:SET(nom).TSize));
gui.d = '.';
gui.dt = 9;

gui.nom = nom;
gui.nop = nop;

%Calculate number of ROI's
rois2take = [];
for loop=1:SET(nom).RoiN
	if sum(~isnan(SET(nom).Roi(loop).X(1,:)))>1
		%More than one timeframe
		rois2take = [rois2take loop]; %#ok<AGROW>
	end;
end;
gui.rois2take = rois2take;
gui.numrois = length(rois2take);

%Warningmessage
%mywarning('Experimental feature currently being developed by Karolinska Institutet',DATA.GUI.Segment);

set(gcf,'color',[0.9 0.9 0.9]);
set([gui.handles.text153 gui.handles.text154 gui.handles.uipanel10 gui.handles.uipanel11 gui.handles.mit_samples gui.handles.mit_slider gui.handles.tri_slider gui.handles.tri_samples gui.handles.checkbox_kmax_line gui.handles.checkbox_beat_seg],'BackgroundColor',[0.9 0.9 0.9]);
set([gui.handles.text153 gui.handles.text154 gui.handles.uipanel10 gui.handles.uipanel11 gui.handles.mit_samples gui.handles.mit_slider gui.handles.tri_slider gui.handles.tri_samples gui.handles.checkbox_kmax_line gui.handles.checkbox_beat_seg],'ForegroundColor',[0 0 0]);


calculate_flow();

%Format the table of results
tab = cell(9,2);
tab(1,1) = {'Mitral'};
tab(2,1) = {'Vmax (m/s)'};
tab(3,1) = {'Vmin (m/s)'};
tab(4,1) = {'Variation'};
tab(5,1) = {''};
tab(6,1) = {'Tricuspid'};
tab(7,1) = {'Vmax (m/s)'};
tab(8,1) = {'Vmin (m/s)'};
tab(9,1) = {'Variation'};
set(gui.handles.table, 'Data', tab);
set(gui.handles.table, 'ColumnWidth',{75,65});

colormap(gray);

%Some settings related to ylim and clim

ymin = 40;
ymax = 164;
gui.ymin = ymin;
gui.ymax = ymax;

%Construct flow histogram image to display in mitral axes
mitral_image = flipud(flow_histogram_no_low(RAMP.mitral_velocity));

%Calculate percentiles, there will always be an outlier that will mess up
%the contrast in the histogram. 
percentile = 99; %Set the percentile of values to be max in colormap.
mit_prctl = prctile(mitral_image,percentile,'all');
axes(gui.handles.mit_axes);
clims=[255-mit_prctl 255]; %Black on white, instead of white on black
gui.mit_image = imagesc(255-mitral_image(:,:),clims);

%tricuspid
tricuspid_image = flipud(flow_histogram_no_low(RAMP.tricuspid_velocity));
tri_prctl = prctile(tricuspid_image,percentile,'all');
clims=[255-tri_prctl 255]; %Black on white, instead of white on black
axes(gui.handles.tri_axes);
gui.tri_image = imagesc(255-tricuspid_image(:,:),clims);

set([gui.mit_image gui.tri_image],'PickableParts','none');
set([gui.mit_image gui.tri_image],'HitTest','off');
%Labels
ylabel(gui.handles.mit_axes,'Velocity (m/s)','fontsize',14); 
xlabel(gui.handles.mit_axes,'Time (s)','fontsize',14);

ylabel(gui.handles.tri_axes,'Velocity (m/s)','fontsize',14); 
xlabel(gui.handles.tri_axes,'Time (s)','fontsize',14);

title(gui.handles.mit_axes,'Mitral inflow velocity','FontSize',18,'fontweight','bold');
title(gui.handles.tri_axes,'Tricuspid inflow velocity','FontSize',18,'fontweight','bold');

%tick settings
set([gui.handles.mit_axes gui.handles.tri_axes],'YTick',linspace(0,256,11));
set([gui.handles.mit_axes gui.handles.tri_axes],'YTickLabel',1:-0.2:-1);
set([gui.handles.mit_axes gui.handles.tri_axes],'XTick',linspace(0,625,31));
set([gui.handles.mit_axes gui.handles.tri_axes],'XTickLabel',0:30);

%Callbakcs for visibility panel
set([gui.handles.checkbox_kmax_line gui.handles.checkbox_beat_seg],'Callback',@checkboxes_callback);


%Draw magnitude image and formatting
gui.magnitude = imagesc(SET(nom).IM(:,:,1),'Parent',gui.handles.mag_axes);
title(gui.handles.mag_axes,'Magnitude','FontSize',18,'fontweight','bold');
set([gui.handles.mag_axes],'XTickLabel',[]);
set([gui.handles.mag_axes],'YTickLabel',[]);

%Draw the mmode image which is constructed by simple slicing of
%SET(nom).IM. The value 100 is the x-value of the slice. 
gui.mmd = imagesc(squeeze(SET(nom).IM(:,100,:)),'Parent',gui.handles.mmd_axes);
%Draw the line
RAMP.overview_line = line(gui.handles.mag_axes,[100 100],[0 300],'LineWidth',2,'color',[1 1 0]); 

ylabel(gui.handles.mmd_axes,'Distance (pixels)','fontsize',14);
title(gui.handles.mmd_axes,'Respiration by lung-diaphragm interface M-mode','FontSize',18,'fontweight','bold');
set(gui.magnitude,'ButtonDownFcn',@mag_click_callback);
set(gui.handles.mmd_axes,'XTick',linspace(0,625,31));
set(gui.handles.mmd_axes,'XTickLabel',0:30);

hold(gui.handles.mit_axes,'on');
hold(gui.handles.tri_axes,'on');
hold(gui.handles.mmd_axes,'on');

%Calculate and plot the black velocity-time line in mit_axes and tri_axes
mit_curve = maxk_curve(RAMP.mitral_velocity,10);
tri_curve = maxk_curve(RAMP.tricuspid_velocity,10);

RAMP.mit_curve = mit_curve;
RAMP.tri_curve = tri_curve;

gui.mit_curve = plot(255-mit_curve,'color','black','LineWidth',1,'Parent',gui.handles.mit_axes,'PickableParts','none');
gui.tri_curve = plot(255-tri_curve,'color','black','LineWidth',1,'Parent',gui.handles.tri_axes,'PickableParts','none');

%Find and plot each heartcycle
splits =calculate_negative_flow(RAMP.mitral_velocity,12);

for i=1:length(splits)
    gui.split_line(1,i) = plot([splits(i) splits(i)],[0,256],'Parent',gui.handles.mit_axes,'color',[0.3 0.3 0.3],'PickableParts','none');
    gui.split_line(2,i) = plot([splits(i) splits(i)],[0,256],'Parent',gui.handles.tri_axes,'color',[0.3 0.3 0.3],'PickableParts','none');
                       
end;

RAMP.splits = splits;

%Find peaks and plot
[RAMP.mit_peaks RAMP.mit_inds] = maxk_test(mit_curve,splits,1);

for i=1:length(RAMP.mit_peaks)
    gui.mit_peak(i) = plot(RAMP.mit_inds(i),255-mit_curve(RAMP.mit_inds(i)),'Marker','o','color','r','MarkerSize',5,'MarkerFaceColor','r','LineWidth',1,'Parent',...
        gui.handles.mit_axes,'ButtonDownFcn',@peak_callback,'UserData',i);
end;

[RAMP.tri_peaks RAMP.tri_inds] = maxk_test(tri_curve,splits,1);
for i=1:length(RAMP.tri_peaks)
    gui.tri_peak(i) = plot(RAMP.tri_inds(i),RAMP.tri_peaks(i),'Marker','o','color','b','MarkerSize',5,'MarkerFaceColor','b','LineWidth',1,'Parent',...
    gui.handles.tri_axes,'ButtonDownFcn',@peak_callback,'UserData',i);
end;

set([gui.handles.mit_axes gui.handles.tri_axes],'ylim',[ymin ymax]);
[RAMP.mit_min_val,RAMP.mit_max_val,RAMP.mit_min_ind,RAMP.mit_max_ind]=get_maxmin(RAMP.mit_peaks);
[RAMP.tri_min_val,RAMP.tri_max_val,RAMP.tri_min_ind,RAMP.tri_max_ind]=get_maxmin(RAMP.tri_peaks);

set(gui.mit_peak(RAMP.mit_max_ind),'MarkerSize',10);
set(gui.mit_peak(RAMP.mit_min_ind),'MarkerSize',10);
set(gui.tri_peak(RAMP.tri_max_ind),'MarkerSize',10);
set(gui.tri_peak(RAMP.tri_min_ind),'MarkerSize',10);

gui.l(1) = plot([RAMP.mit_inds(RAMP.mit_max_ind) RAMP.mit_inds(RAMP.mit_max_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','r','LineWidth',2);
gui.l(2) = plot([RAMP.mit_inds(RAMP.mit_min_ind) RAMP.mit_inds(RAMP.mit_min_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','r','LineWidth',2);
gui.l(3) = plot([RAMP.tri_inds(RAMP.tri_max_ind) RAMP.tri_inds(RAMP.tri_max_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','b','LineWidth',2);
gui.l(4) = plot([RAMP.tri_inds(RAMP.tri_min_ind) RAMP.tri_inds(RAMP.tri_min_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','b','LineWidth',2);

set(gui.handles.mit_slider,'Callback',@btn_callback);
set(gui.handles.tri_slider,'Callback',@btn_callback);

%select_maxmin();
update_table();
DATA.GUI.RAMP = gui;
ok = true;
end

function btn_callback(gcbo,event)

global RAMP DATA
gui = DATA.GUI.RAMP;
x=0;
y=0;

% switch gcbo
%     case gui.handles.btn_mit_up
%         x=1;
%     case gui.handles.btn_mit_down
%         x=-1;
%     case gui.handles.btn_tri_up
%         y=1;
%     case gui.handles.btn_tri_down 
%         y=-1;
% end;

set(gui.handles.mit_samples,'String', get(gui.handles.mit_slider,'Value'));
set(gui.handles.tri_samples,'String', get(gui.handles.tri_slider,'Value'));

splits = RAMP.splits;

RAMP.mit_curve = maxk_curve(RAMP.mitral_velocity,round(get(gui.handles.mit_slider,'Value')));
RAMP.tri_curve = maxk_curve(RAMP.tricuspid_velocity,round(get(gui.handles.tri_slider,'Value')));

%mit_curve = maxk_curve(RAMP.mitral_velocity,round(size(RAMP.mitral_velocity,2)*0.1));
%tri_curve = maxk_curve(RAMP.tricuspid_velocity,round(size(RAMP.tricuspid_velocity,2)*0.2));

set(gui.mit_curve,'YData',255-RAMP.mit_curve);
set(gui.tri_curve,'YData',255-RAMP.tri_curve);


%[RAMP.mit_peaks RAMP.mit_inds] = maxk_test(mit_curve,splits,4);
%[RAMP.tri_peaks RAMP.tri_inds] = maxk_test(tri_curve,splits,3);
update_curve();
end

function update_curve()
%Find peaks and plot
global RAMP DATA
gui=DATA.GUI.RAMP;

[RAMP.mit_peaks RAMP.mit_inds] = maxk_test(RAMP.mit_curve,RAMP.splits,2);
[RAMP.tri_peaks RAMP.tri_inds] = maxk_test(RAMP.tri_curve,RAMP.splits,1);
%delete(gui.mit_peak(:));
%delete(gui.tri_peak(:));

for i=1:length(RAMP.mit_peaks)
    %gui.mit_peak(i) = plot(RAMP.mit_inds(i),RAMP.mit_peaks(i),'Marker','o','color','r','MarkerSize',5,'Parent',...
    %    gui.handles.mit_axes,'ButtonDownFcn',@peak_callback,'UserData',i);
    set(gui.mit_peak(i),'YData',RAMP.mit_peaks(i));
    set(gui.mit_peak(i),'XData',RAMP.mit_inds(i));
    
end;


for i=1:length(RAMP.tri_peaks)
    set(gui.tri_peak(i),'YData',RAMP.tri_peaks(i));
    set(gui.tri_peak(i),'XData',RAMP.tri_inds(i));
end;

set([gui.handles.mit_axes gui.handles.tri_axes],'ylim',[gui.ymin gui.ymax]);

[RAMP.mit_min_val,RAMP.mit_max_val,RAMP.mit_min_ind,RAMP.mit_max_ind]=get_maxmin(RAMP.mit_peaks);
[RAMP.tri_min_val,RAMP.tri_max_val,RAMP.tri_min_ind,RAMP.tri_max_ind]=get_maxmin(RAMP.tri_peaks);


set(gui.mit_peak(:),'MarkerSize',5);
set(gui.tri_peak(:),'MarkerSize',5);



set(gui.mit_peak(RAMP.mit_max_ind),'MarkerSize',10);
set(gui.mit_peak(RAMP.mit_min_ind),'MarkerSize',10);
set(gui.tri_peak(RAMP.tri_max_ind),'MarkerSize',10);
set(gui.tri_peak(RAMP.tri_min_ind),'MarkerSize',10);

set(gui.l(1),'XData',[RAMP.mit_inds(RAMP.mit_max_ind) RAMP.mit_inds(RAMP.mit_max_ind)]);
set(gui.l(2),'XData',[RAMP.mit_inds(RAMP.mit_min_ind) RAMP.mit_inds(RAMP.mit_min_ind)]);
set(gui.l(3),'XData',[RAMP.tri_inds(RAMP.tri_max_ind) RAMP.tri_inds(RAMP.tri_max_ind)]);
set(gui.l(4),'XData',[RAMP.tri_inds(RAMP.tri_min_ind) RAMP.tri_inds(RAMP.tri_min_ind)]);

%gui.l(1) = plot([RAMP.mit_inds(RAMP.mit_max_ind) RAMP.mit_inds(RAMP.mit_max_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','white','LineWidth',2);
%gui.l(2) = plot([RAMP.mit_inds(RAMP.mit_min_ind) RAMP.mit_inds(RAMP.mit_min_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','white','LineWidth',2);
%gui.l(3) = plot([RAMP.tri_inds(RAMP.tri_max_ind) RAMP.tri_inds(RAMP.tri_max_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','white','LineWidth',2);
%gui.l(4) = plot([RAMP.tri_inds(RAMP.tri_min_ind) RAMP.tri_inds(RAMP.tri_min_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','white','LineWidth',2);

%select_maxmin();
update_table();
DATA.GUI.RAMP = gui;
end

function update_table()
    global DATA RAMP SET
    gui = DATA.GUI.RAMP;
    table_contents = get(gui.handles.table,'Data');
    gui.clip_array=[];
    %gui.clip_array(1)=str2num(SET(1).OrigFileName);
    mit_max = -(RAMP.mit_max_val-128)/128;
    mit_min = -(RAMP.mit_min_val-128)/128;
    tri_max = -(RAMP.tri_max_val-128)/128;
    tri_min = -(RAMP.tri_min_val-128)/128;
    %mitral
    gui.clip_array(2)=mit_max;
    gui.clip_array(3)=mit_min;
    gui.clip_array(4)=round((mit_max-mit_min)/mit_max*100);
    table_contents(2,2) = {sprintf('%.2f',gui.clip_array(2))};
    table_contents(3,2) = {sprintf('%.2f',gui.clip_array(3))};
    table_contents(4,2) = {sprintf('%d %%',gui.clip_array(4))};
    
    gui.clip_array(5)=tri_max;
    gui.clip_array(6)=tri_min;
    gui.clip_array(7)=round((tri_max-tri_min)/tri_max*100);
    table_contents(7,2) = {sprintf('%.2f',gui.clip_array(5))};
    table_contents(8,2) = {sprintf('%.2f',gui.clip_array(6))};
    table_contents(9,2) = {sprintf('%d %%',gui.clip_array(7))};

    set(gui.handles.table,'Data',table_contents);
    %send to to clipboard
    num2clip(gui.clip_array);

end

function peak_callback(gcbo,event)
    global DATA RAMP
    gui = DATA.GUI.RAMP;
    button = event.Button;
    source = event.Source;
    ind = get(source,'UserData');
    xdata = get(source,'XData');
    switch button
        case 1 %left click set max value
            switch get(source,'Parent')
                case gui.handles.mit_axes
                    set(gui.mit_peak(RAMP.mit_max_ind),'MarkerSize',5);
                    
                    RAMP.mit_max_ind = ind;
                    RAMP.mit_max_val = RAMP.mit_peaks(ind);
                    set(source,'MarkerSize',8);
                    set(gui.l(1),'XData',[xdata xdata]);

                case gui.handles.tri_axes
                    set(gui.tri_peak(RAMP.tri_max_ind),'MarkerSize',5);
                   
                    RAMP.tri_max_ind = ind;
                    RAMP.tri_max_val = RAMP.tri_peaks(ind);
                    set(source,'MarkerSize',8);
                    set(gui.l(3),'XData',[xdata xdata]);
            end;
        case 3  %right click set min value
        switch get(source,'Parent')
                case gui.handles.mit_axes
                    set(gui.mit_peak(RAMP.mit_min_ind),'MarkerSize',5);
                   
                    RAMP.mit_min_ind = ind;
                    RAMP.mit_min_val = RAMP.mit_peaks(ind);
                    set(source,'MarkerSize',8);
                    set(gui.l(2),'XData',[xdata xdata]);
                case gui.handles.tri_axes
                    set(gui.tri_peak(RAMP.tri_min_ind),'MarkerSize',5);
                    
                    RAMP.tri_min_ind = ind;
                    RAMP.tri_min_val = RAMP.tri_peaks(ind);
                    set(source,'MarkerSize',8);
                    set(gui.l(4),'XData',[xdata xdata]);
            end;
    end;
    update_table();
end

function img = flow_histogram(velocities)
    img = zeros([256 size(velocities,1)]);
    for time_index=1:size(velocities,1)
        for pixel_index=1:length(velocities(time_index,:))
            % + 1 is needed here because matlab is awful
            ind = (velocities(time_index,pixel_index))*255+1;
            ind(ind<0)=0;
            ind(ind>255)=255;
            img(round(ind),time_index)=img(round(ind),time_index)+1;
        end;
    end;
end

function img = flow_histogram_no_low(velocities)
    dist = 4;
    img = zeros([256 size(velocities,1)]);
    for time_index=1:size(velocities,1)
        for pixel_index=1:length(velocities(time_index,:))
            % + 1 is needed here because matlab is awful
            
            ind = (velocities(time_index,pixel_index))*255+1;
            ind(ind<0)=0;
            ind(ind>255)=255;
            if ind >128+dist || ind < 128-dist
                img(round(ind),time_index)=img(round(ind),time_index)+1;
            end;
            
        end;
    img(128,time_index)=255;
    end;
    
end

function curve = maxk_curve(velocity,k)
    for i=1:size(velocity,1)
        curve(i)=(sum(maxk(velocity(i,:),k))/k)*255+1;
    end;
end

function [peaks,inds] = maxk_test(velocity,splits,offset)
    ewaves={};
    global DATA
    gui = DATA.GUI.RAMP;
    %B = 1/2*ones(2,1);
    %velocity = 255-velocity;
    %velocity = filter(B,1,velocity);
    
    for i=1:length(velocity)
        if velocity(i)>128
            velocity(i)=velocity(i);
        else
            velocity(i)=128;
        end
    end;
    m = mean(velocity)+offset;
    %plot([0 625],[255-m 255-m],'Parent',gui.handles.tri_axes);
    for i=1:length(velocity)
        if velocity(i)>m
            velocity(i)=velocity(i);
        else
            velocity(i)=m;
        end
    end;
    for i=1:length(splits)-1
        ewaves{i}=velocity(splits(i):splits(i+1));
    end;

    for i=1:size(ewaves,2)
        for ii=1:size(ewaves{i})
            t = splits(i)+ii+1;
        
            [pks,locs,w,h] = findpeaks(ewaves{i});
            w
            if(isempty(pks))
               peaks(i)=-1
               inds(i)=t-9;
            else
            [maxv,maxind] = max(pks)    
            peak_e_ind = locs(maxind)%min(locs);
            peaks(i) =255-ewaves{i}(peak_e_ind);     
            inds(i) = peak_e_ind+t-3;
            
            end;
        end;
    end;
end

function splits = calculate_negative_flow(velocity,min_peak_distance)
    neg_velocity = min(velocity,[],2);
    neg_velocity = neg_velocity(neg_velocity<0.5);
    B = 1/5*ones(5,1);
    neg_velocity = filter(B,1,neg_velocity);
    [peak,loc]=findpeaks(-neg_velocity,'MinPeakDistance',min_peak_distance);
    splits = loc-3;
end 
        
function draw_mmode(x)
global DATA SET
gui=DATA.GUI.RAMP;
nom = gui.nom;
nop = gui.nop;

%gui.phase = imagesc(squeeze(SET(nom).IM(:,x,:)),'Parent',gui.handles.mmd_axes);
set(gui.mmd,'CData',squeeze(SET(nom).IM(:,x,:)));
end

function checkboxes_callback(varargin)
    global DATA
    gui = DATA.GUI.RAMP;
    switch varargin{1}.Tag 
        case 'checkbox_kmax_line'
            set([gui.mit_curve gui.tri_curve],'Visible',get(varargin{1},'Value'));

        case 'checkbox_beat_seg'
            set(gui.split_line(:),'Visible',get(varargin{1},'Value'));
    end;

end

function mag_click_callback(gcbo,event)
    global DATA RAMP;
    gui = DATA.GUI.RAMP;
    pos = get(get(gco,'Parent'),'CurrentPoint');
    h = gco;

    t = pos(1,1);
    draw_mmode(round(t));

    set(RAMP.overview_line,'Xdata',[t t]);
    %uistack(RAMP.overview_line,'top');
end

function [max_val,min_val,max_ind,min_ind] = get_maxmin(peaks)
[max_val,max_ind] = max(peaks);
[min_val,min_ind] = min(peaks);
end

function select_maxmin()
global DATA SET RAMP
gui = DATA.GUI.RAMP;
table_contents = get(gui.handles.table,'Data');

mit_peaks = RAMP.mit_peaks;
tri_peaks = RAMP.tri_peaks;

mit_inds = RAMP.mit_inds;
tri_inds = RAMP.tri_inds;

h = findobj('YData',max(mit_peaks));
set(h,'MarkerSize',8);

h = findobj('YData',min(mit_peaks));
set(h,'MarkerSize',8);

h = findobj('YData',max(tri_peaks));
set(h,'MarkerSize',8);

h = findobj('YData',min(tri_peaks));
set(h,'MarkerSize',8);

v1 = -(mit_peaks-128)/128;
v2 = -(tri_peaks-128)/128;

gui.clip_array=[];
gui.clip_array(1)=str2num(SET(1).OrigFileName);

%mitral
gui.clip_array(2)=max(v1);
gui.clip_array(3)=min(v1);
gui.clip_array(4)=round((max(v1)-min(v1))/max(v1)*100);
table_contents(2,2) = {sprintf('%.2f',gui.clip_array(2))};
table_contents(3,2) = {sprintf('%.2f',gui.clip_array(3))};
table_contents(4,2) = {sprintf('%d %%',gui.clip_array(4))};

%tricuspid
gui.clip_array(5)=max(v2);
gui.clip_array(6)=min(v2);
gui.clip_array(7)=round((max(v2)-min(v2))/max(v2)*100);
table_contents(7,2) = {sprintf('%.2f',gui.clip_array(5))};
table_contents(8,2) = {sprintf('%.2f',gui.clip_array(6))};
table_contents(9,2) = {sprintf('%d %%',gui.clip_array(7))};

set(gui.handles.table,'Data',table_contents);
%send to to clipboard
num2clip(gui.clip_array);
end

function mit_tri_click_callback(gcbo,event)
global DATA SET
gui = DATA.GUI.RAMP;

button = event.Button;
source = event.Source;
gui.clip_array=[];
gui.clip_array(1)=str2num(SET(1).OrigFileName);
if source == gui.mit_image
    pos = get(gui.handles.mit_axes,'CurrentPoint');
    offset=1;
end;
if source == gui.tri_image
    pos = get(gui.handles.tri_axes,'CurrentPoint');
    offset=3;
end;

switch button
    case 1
        set(gui.p(offset),'XData',pos(1,1));
        set(gui.p(offset),'YData',pos(1,2)); 
        set(gui.l(offset),'Xdata',[pos(1,1) pos(1,1)]);
    case 3
        set(gui.p(offset+1),'XData',pos(1,1));
        set(gui.p(offset+1),'YData',pos(1,2));
        set(gui.l(offset+1),'Xdata',[pos(1,1) pos(1,1)]);
end;
table_contents = get(gui.handles.table,'Data');


    if get(gui.p(1),'XData') ~=-1 && get(gui.p(2),'XData') ~= -1
        v1 = [-(get(gui.p(1),'YData')-128)/128 -(get(gui.p(2),'YData')-128)/128]; 
        gui.clip_array(2)=max(v1);
        gui.clip_array(3)=min(v1);
        gui.clip_array(4)=round((max(v1)-min(v1))/max(v1)*100);
        table_contents(2,2) = {sprintf('%.2f',gui.clip_array(2))};
        table_contents(3,2) = {sprintf('%.2f',gui.clip_array(3))};
        table_contents(4,2) = {sprintf('%d %%',gui.clip_array(4))};
    end;


    if get(gui.p(3),'XData') ~=-1 && get(gui.p(4),'XData') ~= -1
        v2 = [-(get(gui.p(3),'YData')-128)/128 -(get(gui.p(4),'YData')-128)/128]; 
        gui.clip_array(5)=max(v2);
        gui.clip_array(6)=min(v2);
        gui.clip_array(7)=round((max(v2)-min(v2))/max(v2)*100);
        table_contents(7,2) = {sprintf('%.2f',gui.clip_array(5))};
        table_contents(8,2) = {sprintf('%.2f',gui.clip_array(6))};
        table_contents(9,2) = {sprintf('%d %%',gui.clip_array(7))};
    end;
    
    %update the table
    set(gui.handles.table,'Data',table_contents);
    %send data to clipboard
    num2clip(gui.clip_array);
end

function arraystring = num2clip(array)


arraystring = num2str(array);
arraystring(:,end+1) = char(10); 
arraystring = reshape(arraystring',1,prod(size(arraystring))); 
arraystringshift = [' ',arraystring]; 
arraystring = [arraystring,' ']; 


arraystring = arraystring((double(arraystring)~=32 | double(arraystringshift)~=32) & ~(double(arraystringshift==10) & double(arraystring)==32) );
arraystring(double(arraystring)==32) = char(9); %convert the space characters to tab characters
arraystring(double(arraystring)==46) = char(44);
arraystring(double(arraystring)==10) = char(0);
arraystring(double(arraystring)==12) = char(0);
arraystring(double(arraystring)==13) = char(0);

clipboard('copy',arraystring); %copy the result to the clipboard ready for pasting

end

function draw_overview(x)
global DATA SET NO RAMP  
gui = DATA.GUI.RAMP;
nom = gui.nom;
nop = gui.nop;

set(gui.magnitude,'CData',SET(nom).IM(:,:,x));
%set(gui.phase,'CData',SET(nop).IM(:,:,x));


%ymin=RAMP.ymin;
%xmin=RAMP.xmin;
% offset=1;
% hold(gui.handles.pha_axes,'on');
% 
%gui.mitralroi = plot(gui.handles.pha_axes,RAMP.ROI(1).Y(:,1)-ymin+offset,RAMP.ROI(1).X(:,1)-xmin+offset);
%gui.tricuspidroi = plot(gui.handles.pha_axes,RAMP.ROI(2).Y(:,1)-ymin+offset,RAMP.ROI(2).X(:,1)-xmin+offset);
%gui.mitral_marker = plot(gui.handles.pha_axes,RAMP.mitral.constrained_y(x)-ymin+offset,RAMP.mitral.constrained_x(x)-xmin+offset,'o');
%gui.tricuspid_marker = plot(gui.handles.pha_axes,RAMP.tricuspid.constrained_y(x)-ymin+offset,RAMP.tricuspid.constrained_x(x)-xmin+offset,'o');
% 

%plot(gui.handles.pha_axes,RAMP.mitral.index_y(x)-RAMP.ymin+offset,RAMP.mitral.index_x(x)-RAMP.xmin+offset,'o','color',[1 0 0]);

%imshow(uint8(DATA.ViewIM{1}(:,:,round(x))),'Parent',gui.handles.pha_axes,'InitialMagnification',300)

%colormap(gui.handles.mag_axes,gray)
%set(RAMP.overview,'ButtonDownFcn',@image_mouseclick_callback);

%hold(gui.handles.mag_axes,'on');

end