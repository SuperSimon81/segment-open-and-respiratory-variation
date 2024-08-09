function [varargout] = plugin_RAMP(fcn,varargin)

if nargin==0
    myfailed('Expects at least one input argument.');
    return;
end

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
end
end

function ok = init(no)
global RAMP DATA SET NO
ok = false; %inform if unable to launch GUI

if nargin<1
    no=NO;
end

if ~isempty(SET(no).Flow)
    nom = SET(no).Flow.MagnitudeNo;
    nop = SET(nom).Flow.PhaseNo;
else
    myfailed('No flow data in current image stack.',DATA.GUI.Segment);
    return;
end

if isempty(SET(nop).Flow.PhaseCorr)
    myfailed('No eddy current compensation applied.',DATA.GUI.Segment);
    return;
end

if isempty(SET(nom).Flow.PhaseNo)
    myfailed('No through plane velocity data found.',DATA.GUI.Segment);
    return;
end

if SET(nom).RoiN==0
    myfailed('No ROIs available. make sure the first ROI is the mitral orifice and the second ROI the tricuspid orifice',DATA.GUI.Segment);
    return;
end

%tempnos=[SET(no).Flow.MagnitudeNo SET(no).Flow.PhaseNo];
imissingle=classcheckim([SET(no).Flow.MagnitudeNo SET(no).Flow.PhaseNo]);%checks so that SET(tempnos).IM is single and can also convert from int16 to singel if user wants
if not(imissingle)
    return;
end

if isempty(SET(nom).EndAnalysis)
    SET(nom).EndAnalysis = SET(nom).TSize;
end

if isempty(SET(nom).StartAnalysis)
    SET(nom).StartAnalysis = 1;
end

if (SET(nom).EndAnalysis-SET(nom).StartAnalysis)<1
    myfailed('No timespan, adjust Start/End analysis time (under stress menu or drag red bar).',DATA.GUI.Segment);
    return;
end

RAMP.nom = nom;
RAMP.nop = nop;
RAMP.numrois = 2;
RAMP.outsize = [SET(nom).XSize SET(nom).YSize];
calculate_flow();

RAMP.gui = mygui('plugin_RAMP.fig');
gui = RAMP.gui;
%set(gui.fig,'Name','');

%Set some gui colors
set(gcf,'color',[0.9 0.9 0.9]);
set([gui.handles.uipanel11 gui.handles.checkbox_kmax_line gui.handles.checkbox_beat_seg],'BackgroundColor',[0.9 0.9 0.9]);
set([gui.handles.uipanel11 gui.handles.checkbox_kmax_line gui.handles.checkbox_beat_seg],'ForegroundColor',[0 0 0]);

%Warningmessage, might be needed if compiled
%mywarning('Experimental feature being developed by Karolinska Institutet',DATA.GUI.Segment);

%Format the table of results
tab = cell(9,2);
tab(:,1) = {'Mitral'; 'Vmax (m/s)'; 'Vmin (m/s)'; 'Variation'; ''; 'Tricuspid'; 'Vmax (m/s)'; 'Vmin (m/s)'; 'Variation'};
set(gui.handles.table, 'Data', tab, 'ColumnWidth', {75, 65});
colormap(gray);

%Some settings related to ylim and clim
ymin = 40;
ymax = 164;
gui.ymin = ymin;
gui.ymax = ymax;

%Construct flow histogram image to display in mitral axes
mitral_image = flipud(flow_histogram_no_low(RAMP.mitral_velocity));

%There will always be an outlier that will mess up
%the contrast in the histogram.
percentile = 99; %Set the percentile of values to be max in colormap.

%mitral
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
set(gui.handles.mag_axes, 'XTickLabel', [], 'YTickLabel', []);

%Draw the mmode image which is constructed by simple slicing of
%SET(nom).IM. The value 100 is the x-value of the slice.
gui.mmd = imagesc(squeeze(SET(nom).IM(:,100,:)),'Parent',gui.handles.mmd_axes);
%Draw the line
RAMP.overview_line = line(gui.handles.mag_axes,[100 100],[0 300],'LineWidth',2,'color',[1 1 0]);

ylabel(gui.handles.mmd_axes,'Distance (pixels)','fontsize',14);
title(gui.handles.mmd_axes,'Respiration by lung-diaphragm interface M-mode','FontSize',18,'fontweight','bold');
set(gui.magnitude,'ButtonDownFcn',@mag_click_callback);
set(gui.handles.mmd_axes, 'XTick', linspace(0,625,31), 'XTickLabel', 0:30);

hold(gui.handles.mit_axes,'on');
hold(gui.handles.tri_axes,'on');
hold(gui.handles.mmd_axes,'on');

%Calculate and plot the black velocity-time line in mit_axes and tri_axes
mit_curve = max_percentile_curve(RAMP.mitral_velocity,95);
tri_curve = max_percentile_curve(RAMP.tricuspid_velocity,95);

RAMP.mit_curve = mit_curve;
RAMP.tri_curve = tri_curve;

gui.mit_curve = plot(255-mit_curve,'color','black','LineWidth',1,'Parent',gui.handles.mit_axes,'PickableParts','none');
gui.tri_curve = plot(255-tri_curve,'color','black','LineWidth',1,'Parent',gui.handles.tri_axes,'PickableParts','none');

%Find and plot each heartcycle
splits =calculate_negative_flow(RAMP.mitral_velocity,12);

for i=1:length(splits)
    gui.split_line(1,i) = plot([splits(i) splits(i)],[0,256],'Parent',gui.handles.mit_axes,'color',[0.3 0.3 0.3], "ButtonDownFcn",@splits_callback,'UserData',i);
    gui.split_line(2,i) = plot([splits(i) splits(i)],[0,256],'Parent',gui.handles.tri_axes,'color',[0.3 0.3 0.3],"ButtonDownFcn",'','UserData',i);
end

RAMP.splits = splits;

%Find peaks and plot
[RAMP.mit_peaks, RAMP.mit_inds] = find_ewave_peaks(mit_curve,splits,1);

for i=1:length(RAMP.mit_peaks)
    gui.mit_peak(i) = plot(RAMP.mit_inds(i),255-mit_curve(RAMP.mit_inds(i)),'Marker','o','color','r','MarkerSize',5,'MarkerFaceColor','r','LineWidth',1,'Parent',...
        gui.handles.mit_axes,'UserData',i);
    set(gui.mit_peak(i),'ButtonDownFcn', @peak_callback)
    %optional line in mmode for every peak
    gui.mmode_lines(i) = plot([RAMP.mit_inds(i) RAMP.mit_inds(i)],[0,256],'Parent',gui.handles.mmd_axes,'color','w','LineWidth',1,'LineStyle',':');

end

[RAMP.tri_peaks, RAMP.tri_inds] = find_ewave_peaks(tri_curve,splits,1);
for i=1:length(RAMP.tri_peaks)
    gui.tri_peak(i) = plot(RAMP.tri_inds(i),RAMP.tri_peaks(i),'Marker','o','color','b','MarkerSize',5,'MarkerFaceColor','b','LineWidth',1,'Parent',...
        gui.handles.tri_axes,'UserData',i);

end

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

set(gui.handles.mit_axes, 'ButtonDownFcn', @curve_click_callback);
set(gui.handles.tri_axes, 'ButtonDownFcn', @curve_click_callback);

update_table();
RAMP.gui = gui;
ok = true;
end

%Core functions
function calculate_flow()
global DATA SET RAMP

nom = RAMP.nom;
nop = RAMP.nop;
numrois = RAMP.numrois;

a = SET(nom).TSize * 2;
h = waitbar(0, 'Calculating flow.');
set(h, 'visible', 'off');
myadjust(h, DATA.GUI.Segment);
set(h, 'visible', 'on');

RAMP.Roi = SET(nom).Roi(1:2);

counter = 0;

for rloop = 1:numrois

    for tloop = RAMP.Roi(rloop).T
        counter = counter + 1;
        waitbar(counter/a, h);

        mask = logical(segment('createmask', RAMP.outsize, RAMP.Roi(rloop).Y(:,tloop), RAMP.Roi(rloop).X(:,tloop)));
        temp = SET(nop).IM(:,:,tloop,RAMP.Roi(rloop).Z);

        veldata = RAMP.Roi(rloop).Sign * (temp - 0.5) * 2 * SET(nop).VENC;

        if ~isempty(SET(nop).Flow.PhaseCorr)
            if SET(nop).Flow.PhaseCorrTimeResolved
                correction = SET(nop).Flow.PhaseCorr(:,:,tloop,RAMP.Roi(rloop).Z);
            else
                correction = SET(nop).Flow.PhaseCorr(:,:,1,RAMP.Roi(rloop).Z);
            end
            veldata = RAMP.Roi(rloop).Sign * (temp - 0.5 - correction) * 2 * SET(nop).VENC;
        end

        RAMP.velocity_unmasked(:,:,tloop) = veldata;
        veldata = veldata(mask);

        if isempty(veldata) && not(warnedempty)
            mywarning('Empty ROI. Please try again', DATA.GUI.Segment);
            warnedempty = true;
        else
            if rloop == 1
                RAMP.mitral_velocity(tloop, :) = temp(mask);
            elseif rloop == 2
                RAMP.tricuspid_velocity(tloop, :) = temp(mask);
            end
        end
    end
end
close(h);
end

function img = flow_histogram(velocities)
% Construct a flow histogram or spectral plot where each voxel contributes one intensity value on the y-level of its velocity.

img = zeros([256 size(velocities,1)]);
for time_index=1:size(velocities,1)
    for pixel_index=1:length(velocities(time_index,:))
        % + 1 is needed here because matlab is awful
        ind = (velocities(time_index,pixel_index))*255+1;
        ind(ind<0)=0;
        ind(ind>255)=255;
        img(round(ind),time_index)=img(round(ind),time_index)+1;
    end
end
end

function img = flow_histogram_no_low(velocities)
% Construct a flow histogram or spectral plot where each voxel contributes one intensity value on the y-level of its velocity.
% Exclude -4 to 4 v velocites like in Doppler echo

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
        end

    end
    img(128,time_index)=255;
end

end

function curve= max_percentile_curve(velocity,percentile)
%constructs the vt curve from the of voxels with above or equal to the nth
%percentile in the timeframe
for i=1:size(velocity,1)
    value = prctile(velocity(i,:),percentile);

    all_values = velocity(i,:);
    values_above_percentile = all_values(all_values>=value);
    %length(values_above_percentile)
    curve(i) = sum(values_above_percentile)/length(values_above_percentile)*255+1;

end
end

function curve = maxk_curve(velocity,k)
%constructs the vt curve from the k nr of voxels with the highest velocity
%in the timeframe

for i=1:size(velocity,1)
    curve(i)=(sum(maxk(velocity(i,:),k))/k)*255+1;
end
end

function [peaks,inds] = find_ewave_peaks(velocity,splits,offset)
%This function is used to identify the peaks of the ewave.
ewaves={};
global RAMP
gui = RAMP.gui;

for i=1:length(velocity)
    if velocity(i)>128
        velocity(i)=velocity(i);
    else
        velocity(i)=128;
    end
end
m = mean(velocity)+offset;


for i=1:length(velocity)
    if velocity(i)>m
        velocity(i)=velocity(i);
    else
        velocity(i)=m;
    end
end
for i=1:length(splits)-1
    ewaves{i}=velocity(splits(i):splits(i+1));
end




for i=1:size(ewaves,2)
    %for ii=1:size(ewaves{i})
    t = splits(i);

    %bias = linspace(1,0.9,length(ewaves{i}));
    %test = ewaves{i}.*bias;

    [pks,locs,w,h] = findpeaks(ewaves{i});

    if(isempty(pks))
        peaks(i)=-1
        inds(i)=t-9;
    else
        [maxv,maxind] = max(w);
        peak_e_ind = locs(maxind);%min(locs);
        peaks(i) =255-ewaves{i}(peak_e_ind);
        inds(i) = peak_e_ind+t-1;

    end
    %end
end

end

function splits = calculate_negative_flow(velocity,min_peak_distance)
%this function is used to determine the negative peaks of mitral flow,
%which corresponds to aortic outflow during systole. This is used to
%separate each heartbeat

global DATA SET RAMP
gui=RAMP.gui;


neg_velocity = sum(mink(velocity,25,2),2)/25;
neg_velocity = neg_velocity(neg_velocity<0.5);
%Plot for debugging
%gui.test = plot(255-neg_velocity*255,'color','red','LineWidth',1,'Parent',gui.handles.mit_axes,'PickableParts','none');
B = 1/5*ones(5,1);
neg_velocity = filter(B,1,neg_velocity);
[peak,loc]=findpeaks(-neg_velocity,'MinPeakDistance',min_peak_distance);
splits = loc-3;
end

function draw_mmode(x)
global SET RAMP
gui=RAMP.gui;
nom = RAMP.nom;
nop = RAMP.nop;

%gui.phase = imagesc(squeeze(SET(nom).IM(:,x,:)),'Parent',gui.handles.mmd_axes);
set(gui.mmd,'CData',squeeze(SET(nom).IM(:,x,:)));
end

function [max_val,min_val,max_ind,min_ind] = get_maxmin(peaks)
[max_val,max_ind] = max(peaks);
[min_val,min_ind] = min(peaks);
end

%Callback functions
function splits_callback(h,~)
global RAMP
gui = RAMP.gui;

set (gcf, 'WindowButtonMotionFcn', @mouseMove);
set(gcf,"WindowButtonUpFcn", @click)
    function mouseMove (object, eventdata)
        C = get (gca, 'CurrentPoint');
        set(h,'XData',[C(1) C(1)]);
        i = get(h,'UserData');
        set(gui.split_line(2,i),'XData',[C(1) C(1)]);
    end

    function click(object, eventdata)
        set (gcf, 'WindowButtonMotionFcn','' );
        set (gcf, 'WindowButtonUpFcn', '');
        i = get(h,'UserData');
        x = get(h,'XData');
        RAMP.splits(i) = x(1);


        %Find peaks and plot
        [RAMP.mit_peaks RAMP.mit_inds] = find_ewave_peaks(RAMP.mit_curve,round(RAMP.splits),1);
        
        delete(gui.mit_peak(:));
        delete(gui.tri_peak(:));
        
        
        for i=1:length(RAMP.mit_peaks)
            
            gui.mit_peak(i) = plot(RAMP.mit_inds(i),RAMP.mit_peaks(i),'Marker','o','color','r','MarkerSize',5,'MarkerFaceColor','r','LineWidth',1,'Parent',...
            gui.handles.mit_axes,'UserData',i);

        end

        [RAMP.tri_peaks RAMP.tri_inds] = find_ewave_peaks(RAMP.tri_curve,round(RAMP.splits),1);
        for i=1:length(RAMP.tri_peaks)
            gui.tri_peak(i) = plot(RAMP.tri_inds(i),RAMP.tri_peaks(i),'Marker','o','color','b','MarkerSize',5,'MarkerFaceColor','b','LineWidth',1,'Parent',...
                gui.handles.tri_axes,'UserData',i);
        end

        %    set([gui.handles.mit_axes gui.handles.tri_axes],'ylim',[ymin ymax]);
        [RAMP.mit_min_val,RAMP.mit_max_val,RAMP.mit_min_ind,RAMP.mit_max_ind]=get_maxmin(RAMP.mit_peaks);
        [RAMP.tri_min_val,RAMP.tri_max_val,RAMP.tri_min_ind,RAMP.tri_max_ind]=get_maxmin(RAMP.tri_peaks);

        set(gui.mit_peak(RAMP.mit_max_ind),'MarkerSize',10);
        set(gui.mit_peak(RAMP.mit_min_ind),'MarkerSize',10);
        set(gui.tri_peak(RAMP.tri_max_ind),'MarkerSize',10);
        set(gui.tri_peak(RAMP.tri_min_ind),'MarkerSize',10);
        %delete(gui.l(:));
        
        %set(gui.l(1),'XData',[RAMP.mit_inds(RAMP.mit_max_ind) RAMP.mit_inds(RAMP.mit_max_ind)])
        %gui.l(1) = plot([RAMP.mit_inds(RAMP.mit_max_ind) RAMP.mit_inds(RAMP.mit_max_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','r','LineWidth',2);
        %gui.l(2) = plot([RAMP.mit_inds(RAMP.mit_min_ind) RAMP.mit_inds(RAMP.mit_min_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','r','LineWidth',2);
        %gui.l(3) = plot([RAMP.tri_inds(RAMP.tri_max_ind) RAMP.tri_inds(RAMP.tri_max_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','b','LineWidth',2);
        %gui.l(4) = plot([RAMP.tri_inds(RAMP.tri_min_ind) RAMP.tri_inds(RAMP.tri_min_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','b','LineWidth',2);
        update_mmode_lines();
        

        
        update_table();
    end

end

function peak_callback(~,event)
%To mnanually switch peaks if something went wrong
global RAMP
gui = RAMP.gui;
button = event.Button;
source = event.Source;
ind = get(source,'UserData');
xdata = get(source,'XData')
ydata = get(source,'YData')
switch button
    case 1 %left click set max value
        asdasd
        text(gui.handles.mit_axes,xdata(1),ydata(1) - 10,0,num2str(-(RAMP.mit_peaks(ind)-128)/128));



end

end

function curve_click_callback(gcbo,event)

global DATA RAMP
gui = RAMP.gui;
x = round(event.IntersectionPoint(1,1));

for i=1:length(RAMP.splits)
    dist(i)= RAMP.splits(i) - x;

end
dist(dist<0)=nan;

[val, ind] = min(dist);
ind = ind -1;


if x>RAMP.splits(1) && x<RAMP.splits(end)
    if gcbo == gui.handles.mit_axes
        RAMP.mit_inds(ind) = x;
        RAMP.mit_peaks(ind) = 255-RAMP.mit_curve(x);
        set(gui.mit_peak(ind),'XData',[x x]);
        set(gui.mit_peak(ind),'YData',[RAMP.mit_peaks(ind) RAMP.mit_peaks(ind)]);
    end

    if gcbo == gui.handles.tri_axes
        RAMP.tri_inds(ind) = x;
        RAMP.tri_peaks(ind) = 255-RAMP.tri_curve(x);
        set(gui.tri_peak(ind),'XData',[x x]);
        set(gui.tri_peak(ind),'YData',[RAMP.tri_peaks(ind) RAMP.tri_peaks(ind)]);
    end
end

[RAMP.mit_min_val,RAMP.mit_max_val,RAMP.mit_min_ind,RAMP.mit_max_ind]=get_maxmin(RAMP.mit_peaks);
[RAMP.tri_min_val,RAMP.tri_max_val,RAMP.tri_min_ind,RAMP.tri_max_ind]=get_maxmin(RAMP.tri_peaks);

set([gui.mit_peak(:); gui.tri_peak(:)], 'MarkerSize', 5);
set([gui.mit_peak(RAMP.mit_max_ind), gui.mit_peak(RAMP.mit_min_ind), gui.tri_peak(RAMP.tri_max_ind), gui.tri_peak(RAMP.tri_min_ind)], 'MarkerSize', 10);

update_mmode_lines();
update_table();

end

function checkboxes_callback(varargin)
global DATA RAMP
gui = RAMP.gui;
switch varargin{1}.Tag
    case 'checkbox_kmax_line'
        set([gui.mit_curve gui.tri_curve],'Visible',get(varargin{1},'Value'));

    case 'checkbox_beat_seg'
        set(gui.split_line(:),'Visible',get(varargin{1},'Value'));
end

end

function mag_click_callback(gcbo,event)
%Select what slice to use as mmode
global DATA RAMP;
gui = RAMP.gui;
pos = get(get(gco,'Parent'),'CurrentPoint');
h = gco;

t = pos(1,1);
draw_mmode(round(t));

set(RAMP.overview_line,'Xdata',[t t]);
%uistack(RAMP.overview_line,'top');
end

%Update functions that update the gui after callbacks have been triggered 
function update_mmode_lines()
global RAMP
gui = RAMP.gui;

set(gui.l(1),'XData',[RAMP.mit_inds(RAMP.mit_max_ind) RAMP.mit_inds(RAMP.mit_max_ind)]);
set(gui.l(2),'XData',[RAMP.mit_inds(RAMP.mit_min_ind) RAMP.mit_inds(RAMP.mit_min_ind)]);
set(gui.l(3),'XData',[RAMP.tri_inds(RAMP.tri_max_ind) RAMP.tri_inds(RAMP.tri_max_ind)]);
set(gui.l(4),'XData',[RAMP.tri_inds(RAMP.tri_min_ind) RAMP.tri_inds(RAMP.tri_min_ind)]);
end

function update_curve()
% Update graphics and plots after changing something
%Find peaks and plot
global RAMP DATA
gui=RAMP.gui;

[RAMP.mit_peaks RAMP.mit_inds] = find_ewave_peaks(RAMP.mit_curve,RAMP.splits,2);
[RAMP.tri_peaks RAMP.tri_inds] = find_ewave_peaks(RAMP.tri_curve,RAMP.splits,1);
%delete(gui.mit_peak(:));
%delete(gui.tri_peak(:));

for i=1:length(RAMP.mit_peaks)
    set(gui.mit_peak(i),'YData',RAMP.mit_peaks(i));
    set(gui.mit_peak(i),'XData',RAMP.mit_inds(i));
    %Optional mmode_lines
    set(gui.mmode_lines(i),'XData',RAMP.mit_inds(i));


    %%gui.mit_peak(i) = plot(RAMP.mit_inds(i),RAMP.mit_peaks(i),'Marker','o','color','r','MarkerSize',5,'Parent',...
    %%    gui.handles.mit_axes,'ButtonDownFcn',@peak_callback,'UserData',i);


end
%
%
% for i=1:length(RAMP.tri_peaks)
%     set(gui.tri_peak(i),'YData',RAMP.tri_peaks(i));
%     set(gui.tri_peak(i),'XData',RAMP.tri_inds(i));
% end
%
% set([gui.handles.mit_axes gui.handles.tri_axes],'ylim',[gui.ymin gui.ymax]);
%
% [RAMP.mit_min_val,RAMP.mit_max_val,RAMP.mit_min_ind,RAMP.mit_max_ind]=get_maxmin(RAMP.mit_peaks);
% [RAMP.tri_min_val,RAMP.tri_max_val,RAMP.tri_min_ind,RAMP.tri_max_ind]=get_maxmin(RAMP.tri_peaks);
%
%
% set(gui.mit_peak(:),'MarkerSize',5);
% set(gui.tri_peak(:),'MarkerSize',5);
%
% set(gui.mit_peak(RAMP.mit_max_ind),'MarkerSize',10);
% set(gui.mit_peak(RAMP.mit_min_ind),'MarkerSize',10);
% set(gui.tri_peak(RAMP.tri_max_ind),'MarkerSize',10);
% set(gui.tri_peak(RAMP.tri_min_ind),'MarkerSize',10);
%
% set(gui.l(1),'XData',[RAMP.mit_inds(RAMP.mit_max_ind) RAMP.mit_inds(RAMP.mit_max_ind)]);
% set(gui.l(2),'XData',[RAMP.mit_inds(RAMP.mit_min_ind) RAMP.mit_inds(RAMP.mit_min_ind)]);
% set(gui.l(3),'XData',[RAMP.tri_inds(RAMP.tri_max_ind) RAMP.tri_inds(RAMP.tri_max_ind)]);
% set(gui.l(4),'XData',[RAMP.tri_inds(RAMP.tri_min_ind) RAMP.tri_inds(RAMP.tri_min_ind)]);

%gui.l(1) = plot([RAMP.mit_inds(RAMP.mit_max_ind) RAMP.mit_inds(RAMP.mit_max_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','white','LineWidth',2);
%gui.l(2) = plot([RAMP.mit_inds(RAMP.mit_min_ind) RAMP.mit_inds(RAMP.mit_min_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','white','LineWidth',2);
%gui.l(3) = plot([RAMP.tri_inds(RAMP.tri_max_ind) RAMP.tri_inds(RAMP.tri_max_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','white','LineWidth',2);
%gui.l(4) = plot([RAMP.tri_inds(RAMP.tri_min_ind) RAMP.tri_inds(RAMP.tri_min_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color','white','LineWidth',2);

%select_maxmin();
update_table();
DATA.GUI.RAMP = gui;
end

function update_table()
global RAMP
gui = RAMP.gui;
table_contents = get(gui.handles.table,'Data');


%Clipboard array to be exported
gui.clip_array=[];

%First column is the OrigFileName from Segment file format as a number
% gui.clip_array(1)=str2num(SET(1).OrigFileName);
%SET(1).OrigFileName

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

for i=1:length(RAMP.mit_peaks)

    gui.clip_array(7+i) = -(RAMP.mit_peaks(i)-128)/128;
end

for i=1:length(RAMP.tri_peaks)

    gui.clip_array(8+length(RAMP.mit_peaks)+i) = -(RAMP.tri_peaks(i)-128)/128;

end
num2clip(gui.clip_array);
end

%Utility funct
function arraystring = num2clip(array)
%Copies the result to the clipboard for pasting
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

