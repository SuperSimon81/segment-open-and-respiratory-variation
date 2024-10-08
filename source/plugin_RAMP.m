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
RAMP = [];
RAMP.debug = true;
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
    myfailed('No through plane velocity data found.',DATA.GUI.Segment);segment
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
RAMP.filename= SET(1).OrigFileName;
set(gui.fig,'Name',RAMP.filename);

%Enable debug
if RAMP.debug == true
    set([gui.handles.panel_debug],'Visible',1); 
end

%Set some gui colors
set(findall(gcf, 'Style', 'pushbutton'),'BackgroundColor',[0.9 0.9 0.9],'ForegroundColor',[0 0 0]); 
set(findall(gcf, 'Style', 'checkbox'),'BackgroundColor',[0.9 0.9 0.9],'ForegroundColor',[0 0 0]); 
set(findall(gcf, 'Style', 'text'),'BackgroundColor',[0.9 0.9 0.9],'ForegroundColor',[0 0 0]);
set(findall(gcf, 'Style', 'edit'),'BackgroundColor',[1 1 1],'ForegroundColor',[0 0 0]);
set([gui.handles.panel_visibility, gui.handles.panel_results, gui.handles.panel_debug],'BackgroundColor',[0.9 0.9 0.9],'ForegroundColor',[0 0 0]); 
set(gcf,'color',[0.9 0.9 0.9]);
set(gui.handles.result_mitral_variation,'String','asd');
%Warningmessage, might be needed if compiled
%mywarning('Experimental feature being developed by Karolinska Institutet',DATA.GUI.Segment);


colormap(gray);

%Some settings related to ylim and clim
ymin = -0.2;
ymax = 1;
RAMP.ymin = ((ymin+2)*512)/4;
RAMP.ymax = ((ymax+2)*512)/4;

%Construct spectral plot image to display in mitral axes
mitral_spectral_plot = spectral_plot(RAMP.mitral.velocity);

%There will always be an outlier that will mess up
%the contrast in the histogram. This sets some windowing
hpercentile = 98; 
lpercentile = 1;

hmit_prctl = prctile(mitral_spectral_plot(:),hpercentile);
lmit_prctl = prctile(mitral_spectral_plot(:),lpercentile);

axes(gui.handles.mit_axes);
clims=[0 hmit_prctl+8]; %Black on white, instead of white on black

gui.mit_image = imagesc(mitral_spectral_plot(:,:),clims);
set(gui.handles.mit_axes,'YDir','normal');
colormap(gui.handles.mit_axes,flipud(gray));
%tricuspid
tricuspid_spectral_plot = spectral_plot(RAMP.tricuspid.velocity);
htri_prctl = prctile(tricuspid_spectral_plot(:),hpercentile);
ltri_prctl = prctile(tricuspid_spectral_plot(:),lpercentile);

clims=[0 htri_prctl+8]; %Black on white, instead of white on black
axes(gui.handles.tri_axes);

gui.tri_image = imagesc(tricuspid_spectral_plot(:,:),clims);
set(gui.handles.tri_axes,'YDir','normal');
colormap(gui.handles.tri_axes,flipud(gray));

set([gui.mit_image gui.tri_image],'PickableParts','none');
set([gui.mit_image gui.tri_image],'HitTest','off');

%Labels
ylabel(gui.handles.mit_axes,'Velocity (m/s)','fontsize',14);
xlabel(gui.handles.mit_axes,'Time (s)','fontsize',14);

ylabel(gui.handles.tri_axes,'Velocity (m/s)','fontsize',14);
xlabel(gui.handles.tri_axes,'Time (s)','fontsize',14);

title(gui.handles.mit_axes,'Mitral inflow velocity','FontSize',18,'fontweight','bold');
title(gui.handles.tri_axes,'Tricuspid inflow velocity','FontSize',18,'fontweight','bold');

%Callbakcs for visibility panel
set(findall(gcf, 'Style', 'checkbox'),'Callback',@checkboxes_callback);

%debug
set(gui.handles.button_debug,'Callback',@button_debug_callback);
set(gui.handles.button_clipboard,'Callback',@button_clipboard_callback);


%Draw magnitude image and formatting
gui.magnitude = imagesc(SET(nom).IM(:,:,1),'Parent',gui.handles.mag_axes);
title(gui.handles.mag_axes,'Magnitude','FontSize',18,'fontweight','bold');
set(gui.handles.mag_axes, 'XTickLabel', [], 'YTickLabel', []);

%Draw the mmode image which is constructed by simple slicing of
%SET(nom).IM. The default value 100 is the x-value of the slice, empirically works for most images.
gui.mmd = imagesc(squeeze(SET(nom).IM(:,100,:)),'Parent',gui.handles.mmd_axes);
%Draw the line
RAMP.overview_line = line(gui.handles.mag_axes,[100 100],[0 300],'LineWidth',2,'color',[1 1 1]);

ylabel(gui.handles.mmd_axes,'Distance (pixels)','fontsize',14);
title(gui.handles.mmd_axes,'Respiration by lung-diaphragm interface M-mode','FontSize',18,'fontweight','bold');
set(gui.magnitude,'ButtonDownFcn',@mag_click_callback);
set(gui.handles.mmd_axes, 'XTick', linspace(0,625,31), 'XTickLabel', 0:30);

hold(gui.handles.mit_axes,'on');
hold(gui.handles.tri_axes,'on');
hold(gui.handles.mmd_axes,'on');

%Calculate and plot the black velocity-time line in mit_axes and tri_axes
%The percentile is set to 95 empirically proved to work well. It can be
%discussed. Observed differences are around 1-2 %-units in the derived variation.
RAMP.mitral.curve = max_percentile_curve(RAMP.mitral.velocity,95);
RAMP.tricuspid.curve = max_percentile_curve(RAMP.tricuspid.velocity,95);

gui.mit_curve_plot = plot(256+128*RAMP.mitral.curve,'color','black','LineWidth',1,'Parent',gui.handles.mit_axes,'PickableParts','none','Visible',get(gui.handles.checkbox_vt_curve,'Value'));
gui.tri_curve_plot = plot(256+128*RAMP.tricuspid.curve,'color','black','LineWidth',1,'Parent',gui.handles.tri_axes,'PickableParts','none','Visible',get(gui.handles.checkbox_vt_curve,'Value'));

%tick settings
set([gui.handles.mit_axes gui.handles.tri_axes],'YTick',linspace(RAMP.ymin,RAMP.ymax,7));
set([gui.handles.mit_axes gui.handles.tri_axes],'YTickLabel',ymin:0.2:ymax);
set([gui.handles.mit_axes gui.handles.tri_axes],'XTick',linspace(0,625,31));
set([gui.handles.mit_axes gui.handles.tri_axes],'XTickLabel',0:30);
%Find and plot each heartcycle
splits =segment_heart_cycles(-RAMP.mitral.velocity,12);

for i=1:length(splits)
    gui.split_line(1,i) = plot([splits(i) splits(i)],[0,512],'Parent',gui.handles.mit_axes,'color',[0.3 0.3 0.3],'LineWidth',0.5, "ButtonDownFcn",@splits_callback,'UserData',i,'Visible',get(gui.handles.checkbox_beat_seg,'Value'));
    gui.split_line(2,i) = plot([splits(i) splits(i)],[0,512],'Parent',gui.handles.tri_axes,'color',[0.3 0.3 0.3],'LineWidth',0.5,"ButtonDownFcn",@splits_callback,'UserData',i,'Visible',get(gui.handles.checkbox_beat_seg,'Value'));
end

RAMP.splits = splits;

%Find peaks and plot
[RAMP.mitral.peaks, RAMP.mitral.indices,debug_curve] = find_ewave_peaks(RAMP.mitral.curve,splits);

gui.mitral_debug = plot(256+128*debug_curve,'color','r','LineWidth',0.5,'Parent',gui.handles.mit_axes,'PickableParts','none','Visible',get(gui.handles.checkbox_vt_curve,'Value'));

%[RAMP.tricuspid.peaks, RAMP.tricuspid.indices,debug_curve] = find_ewave_peaks(RAMP.tricuspid.curve,splits);
[RAMP.tricuspid.peaks, RAMP.tricuspid.indices] = copy_peaks(RAMP.mitral.indices,-RAMP.tricuspid.curve);
gui.tricuspid_debug = plot(256+128*debug_curve,'color','r','LineWidth',0.5,'Parent',gui.handles.tri_axes,'PickableParts','none','Visible',get(gui.handles.checkbox_vt_curve,'Value'));


plot_peaks();



set([gui.handles.mit_axes gui.handles.tri_axes],'ylim',[RAMP.ymin RAMP.ymax]);
[RAMP.mitral.max_val,RAMP.mitral.min_val,RAMP.mitral.max_ind,RAMP.mitral.min_ind]=get_maxmin(RAMP.mitral.peaks);
[RAMP.tricuspid.max_val,RAMP.tricuspid.min_val,RAMP.tricuspid.max_ind,RAMP.tricuspid.min_ind]=get_maxmin(RAMP.tricuspid.peaks);

update_marker_size();

set([gui.handles.mit_axes gui.handles.tri_axes], 'ButtonDownFcn', @curve_click_callback);
set([gui.handles.mit_axes gui.handles.tri_axes], 'ButtonDownFcn', @curve_click_callback);
set([gui.handles.mit_axes gui.handles.tri_axes], 'ButtonDownFcn', @curve_click_callback);

%[RAMP.exploc RAMP.insploc RAMP.inflections RAMP.respiration] = calculate_respiratory_curve(112);

color = 'w';    
gui.l(1) = plot([RAMP.mitral.indices(RAMP.mitral.max_ind) RAMP.mitral.indices(RAMP.mitral.max_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color',color,'LineWidth',1);
 
gui.l(2) = plot([RAMP.mitral.indices(RAMP.mitral.min_ind) RAMP.mitral.indices(RAMP.mitral.min_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color',color,'LineWidth',1);
  
gui.l(3) = plot([RAMP.tricuspid.indices(RAMP.tricuspid.max_ind)+1 RAMP.tricuspid.indices(RAMP.tricuspid.max_ind)+1],[0,256],'Parent',gui.handles.mmd_axes,'color',color,'LineWidth',1);

gui.l(4) = plot([RAMP.tricuspid.indices(RAMP.tricuspid.min_ind)+1 RAMP.tricuspid.indices(RAMP.tricuspid.min_ind)+1],[0,256],'Parent',gui.handles.mmd_axes,'color',color,'LineWidth',1);


update_table();
%RAMP.gui = gui;



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

        %veldata = RAMP.Roi(rloop).Sign * (temp - 0.5) * 2 * SET(nop).VENC;

        if ~isempty(SET(nop).Flow.PhaseCorr)
            if SET(nop).Flow.PhaseCorrTimeResolved
                mywarning('Time resolved eddy current correction is not supported.', DATA.GUI.Segment);
            warnedempty = true;
            else
                correction = SET(nop).Flow.PhaseCorr(:,:,1,RAMP.Roi(rloop).Z);
            end
            veldata = RAMP.Roi(rloop).Sign * (temp - 0.5 - correction) *2* SET(nop).VENC;
            %old_veldata = (temp-0.5-SET(nop).Flow.PhaseCorr(:,:,1,SET(nom).Roi(rois2take(rloop)).Z))*2*SET(nop).VENC
        end

        RAMP.velocity_unmasked(:,:,tloop) = veldata;
        veldata = veldata(mask);
        try
        if isempty(veldata) && not(warnedempty)
            mywarning('Empty ROI. Please try again', DATA.GUI.Segment);
            warnedempty = true;
        else
            if rloop == 1
                RAMP.mitral.velocity(tloop, :) = veldata/100; %m/s not cm/s
            elseif rloop == 2
                RAMP.tricuspid.velocity(tloop, :) = veldata/100; %m/s not cm/s
            end
        end
        catch
        close(h);
        
        
        return
        end
    end
end
close(h);
end

function img = spectral_plot(velocities)
% Construct a flow histogram or spectral plot where each voxel contributes one intensity value on the y-level of its velocity.
% Exclude -4 to 4 v velocites like in Doppler 




dist = 4;
img = zeros([512 size(velocities,1)]);

for time_index=1:size(velocities,1)
    for pixel_index=1:length(velocities(time_index,:))
        %linear transformation to bring [-2 2] to [0 to 512]
        ind = (((velocities(time_index,pixel_index))+2)*512)/4;
        ind(ind<1)=1;
        %if ind>512
        %    disp(ind);
        %    disp('ding');
        %end
        
        %ind(ind>255)=255;
        
        if ind >256+dist || ind < 256-dist
            img(round(ind),time_index)=img(round(ind),time_index)+1;
        end

    end
end


maximum = max(max(img));

%img = (img/maximum)*255


for time_index=1:size(velocities,1)    
    img(256,time_index)=maximum;
end
end

function curve = max_percentile_curve(velocity,percentile)
%constructs the vt curve from the of voxels with above the nth
%percentile in each the timeframe
for i=1:size(velocity,1)
    value = prctile(velocity(i,:),percentile);
    
    all_values = velocity(i,:);
    values_above_percentile = all_values(all_values>value);
    
    curve(i) = sum(values_above_percentile)/length(values_above_percentile);
    
end
end

function [peaks, inds,debug_curve] = find_ewave_peaks(velocity, splits)
% This function is used to identify the peaks of the ewave from a velocity signal.
% Inputs:
%   velocity - The input signal array from which peaks need to be identified.
%   splits - Indices in the velocity array that demarcate different segments.

ewaves = {}; % Initialize a cell array to store segments of the velocity signal.
global RAMP 
gui = RAMP.gui; 

% Only keep positive values 
for i = 1:length(velocity)
    if velocity(i) < 0
        velocity(i) = 0;
    end
end

% Calculate the threshold as the mean of the modified velocity plus an offset.




m = mean(velocity) + std(velocity);


% Adjust the velocity such that all values below the threshold m are set to m.
for i = 1:length(velocity)
    if velocity(i) < m
        velocity(i) = m;
    end
end

debug_curve= velocity;


% Split the velocity into segments according to the indices provided in splits.
for i = 1:length(splits) - 1
    ewaves{i} = velocity(splits(i):splits(i + 1));
    
end

% Initialize arrays to store the peaks and their corresponding indices.
peaks = zeros(1, size(ewaves, 2));
inds = zeros(1, size(ewaves, 2));


% Find peaks within each segment.
for i = 1:size(ewaves, 2)
    t = splits(i); % Starting index of the current segment.

    % Find peaks in the current segment using the findpeaks function.
    [pks, locs, w, p] = findpeaks(ewaves{i});

    if isempty(pks) % If no peaks are found,
        peaks(i) = -1; % assign -1 to indicate no peak.
        inds(i) = t; % Assign an index which is adjusted from the start by -9. 
        
    else
        %[maxv, maxind] = max(p); % Find the peak with the maximum prominence as defined in findpeaks().
        [maxv,maxind] = min(locs); %find the peak that is leftmost.
        peak_e_ind = locs(maxind); % Find the index of this peak within the segment.
        peaks(i) = ewaves{i}(peak_e_ind); % Calculate the peak value 
        inds(i) = peak_e_ind + t - 1; % Adjust the index relative to the whole velocity array.
    end
    
    
end

    %Remove peaks with value 256, not an ideal solution
    %k = find(inds==256)
    %peaks = [peaks(1:k-1), peaks(k+1:end)];
    %inds = [inds(1:k-1), inds(k+1:end)];
    %splits = [splits(1:k-1), splits(k+1:end)];
   
end

function [peaks, inds] = copy_peaks(input_inds,input)
inds = input_inds;
input=-input; %because tricuspid
steps_to_check = 3;
for i=1:length(inds)
    
        [peak, ind] = findpeaks(input(inds(i)-steps_to_check:inds(i)+steps_to_check));
        if isempty(peak)
        peaks(i)=input(inds(i));
        else
        [max_peak, max_ind] = max(peak);
        peaks(i)=max_peak;
        inds(i)=inds(i)+ind(max_ind)-steps_to_check-1;
        end
    
end


end

function splits = segment_heart_cycles(velocity, min_peak_distance)
% This function calculates the locations of negative peaks in mitral flow velocity,
% which correspond to the aortic outflow during systole. This information is used to
% delineate individual heartbeats.
% Inputs:
%   velocity - An array of velocity measurements.
%   min_peak_distance - Minimum number of timeframes between successive peaks to avoid
%                       detection of closely spaced peaks.

global RAMP 
gui = RAMP.gui; % Access the GUI component from the global RAMP structure.

neg_velocity = sum(maxk(velocity, 5, 2), 2) / 5; % Smoothes the lowest 25 velocity values by averaging them.
%neg_velocity = neg_velocity(neg_velocity < 0.5); % Retain only values of the smoothed array that are less than 0.5.

% Optionally, uncomment the following line to plot the negative velocity for debugging purposes.
%gui.test = plot(256+256*neg_velocity, 'color', 'red', 'LineWidth', 1, 'Parent', gui.handles.mit_axes, 'PickableParts', 'none');

% Define a moving average filter with a window size of 5 and a coefficient of 1/5.
B = 1/5 * ones(5, 1);
neg_velocity = filter(B, 1, neg_velocity); % Apply the moving average filter to smooth the signal further.

% Find peaks in the negative of the smoothed negative velocity signal.
[peak, loc,w,p] = findpeaks(neg_velocity, 'MinPeakProminence',0.01,'MinPeakDistance', min_peak_distance); % Detect peaks, treating negative troughs as peaks.
% Adjust peak locations to account for potential lead-in effect of filtering and other processing steps.
splits = loc - 2; % Subtract 3 from the locations to compensate for any delays introduced by processing.
%splits = splits(splits>3);
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

function [exploc,insploc,inflections,respiration] = calculate_respiratory_curve(y)
global DATA RAMP;
gui = RAMP.gui;
ywindowsize = 75;
pos = RAMP.overview_line.XData(1);
posy = round(y);
pos = round(pos);
count=0;

for posloop=pos-5:pos+5
safe_posy = posy(posy>0+round(ywindowsize/2) & posy<round(size(DATA.ViewIM{1},1)-round(ywindowsize/2)));
safe_posloop = posloop(posloop>0 & posloop<round(size(DATA.ViewIM{1},2)));
disp(safe_posloop);
disp(safe_posy);

image = squeeze(uint8(DATA.ViewIM{1}(safe_posy-round(ywindowsize/2):safe_posy+round(ywindowsize/2),round(safe_posloop),:)));
y=y-5;
count=count+1;

exploc = [];
insploc = [];
inflections = [];
respiration = [];
BW = grayconnected(image,70,500,8);

for i = 1:size(BW,2)
    for ii = 1:size(BW,1)
            if BW(ii,i) == 1
                c(count,i)=ii;
                
                break;
        end
        
    end
end



end


d=round(mean(c,1));
respiration = d;

filt = 8;
B = 1/filt*ones(filt,1);
e = filter(B,1,d);

%The mean peak distance here is arbitrary but is constrained by the physiologically
%plausible breathing frequency which odc affects the duration of expiration and
%inspiration.

[exp,exploc]= findpeaks(e,'MinPeakDistance',40,'MinPeakProminence',0.5);

if ~isempty(exp)  
    for i =1:length(exploc)-1
    insploc(i) = exploc(i)+ round((exploc(i+1)-exploc(i))/2);
    end
else
return
end
lines = findobj(gui.handles.mmd_axes, 'Type', 'line','-and','LineWidth',1); % Find all line objects with width 1.5
delete(lines); % Delete all of them
gui.resp_curve = plot(gui.handles.mmd_axes,[1:exploc(1)],d(1:exploc(1))+y,'Color','w','LineWidth',1,'Visible',get(gui.handles.checkbox_resp_curve,'Value'));

for i = 1:length(exploc)
    
    try
    gui.resp_curve = [gui.resp_curve plot(gui.handles.mmd_axes,[exploc(i):insploc(i)],d(exploc(i):insploc(i))+y,'Color','r','LineWidth',1,'Visible',get(gui.handles.checkbox_resp_curve,'Value'))];
catch
    
end

try
    gui.resp_curve = [gui.resp_curve plot(gui.handles.mmd_axes,[insploc(i):exploc(i+1)],d(insploc(i):exploc(i+1))+y,'Color',[1 0 0],'LineWidth',1,'Visible',get(gui.handles.checkbox_resp_curve,'Value'))];
catch
    
end

end
if (max(exploc)>max(insploc))
    gui.resp_curve = [gui.resp_curve plot(gui.handles.mmd_axes,[exploc(end):size(d,2)],d(exploc(end):size(d,2))+y,'Color',[0 0 1],'LineWidth',1,'Visible',get(gui.handles.checkbox_resp_curve,'Value'))];
 
else
    gui.resp_curve = [gui.resp_curve plot(gui.handles.mmd_axes,[insploc(end):size(d,2)],d(insploc(end):size(d,2))+y,'Color',[1 0 0],'LineWidth',1,'Visible',get(gui.handles.checkbox_resp_curve,'Value'))];
 
    
end

%plot(gui.handles.mmd_axes,exploc,exp+y,'bo','MarkerSize',10);
%plot(gui.handles.mmd_axes,insploc,d(insploc)+y,'go','MarkerSize',10);

inflections = sort([exploc , insploc , 0 , size(BW,2)]);

%Determine which is the first inflection and add a 0 before it to deal with
%the first half of the very first resp cycle
[minv,~] = min([exploc,insploc]);
if ismember(minv,exploc) 
    insploc = [0,insploc];
else
    exploc = [0,exploc];

end



%mitral

gui.l = gobjects(1,4); 
%Check if mitral_max occurs during expiration
if ~isempty(insploc) && ~isempty(exploc)
    
a = insploc - RAMP.mitral.indices(RAMP.mitral.max_ind);
b = exploc - RAMP.mitral.indices(RAMP.mitral.max_ind);
[~,ind] = min([min(a(a>0)),min(b(b>0))]);
disp(ind)
if ind == 1 
    color = 'b';
else
    color = 'r';
end
gui.l(1) = plot([RAMP.mitral.indices(RAMP.mitral.max_ind) RAMP.mitral.indices(RAMP.mitral.max_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color',color,'LineWidth',1);

%Check if mitral_min occurs during inspiration
a = insploc - RAMP.mitral.indices(RAMP.mitral.min_ind);
b = exploc - RAMP.mitral.indices(RAMP.mitral.min_ind);
[~,ind] = min([min(a(a>0)),min(b(b>0))]);
disp(ind)
if ind == 1 
    color = 'r';
else
    color = 'b';
end
gui.l(2) = plot([RAMP.mitral.indices(RAMP.mitral.min_ind) RAMP.mitral.indices(RAMP.mitral.min_ind)],[0,256],'Parent',gui.handles.mmd_axes,'color',color,'LineWidth',1);

%Check if tricuspid_max occurs during inspiration
a = insploc - RAMP.tricuspid.indices(RAMP.tricuspid.max_ind);
b = exploc - RAMP.tricuspid.indices(RAMP.tricuspid.max_ind);
[~,ind] = min([min(a(a>0)),min(b(b>0))])
%disp(ind)
if ind == 1 
    color = 'r';
else
    color = 'b';
end
gui.l(3) = plot([RAMP.tricuspid.indices(RAMP.tricuspid.max_ind)+1 RAMP.tricuspid.indices(RAMP.tricuspid.max_ind)+1],[0,256],'Parent',gui.handles.mmd_axes,'color',color,'LineWidth',1);

%Check if tricuspid_min occurs during expiration
a = insploc - RAMP.tricuspid.indices(RAMP.tricuspid.min_ind);
b = exploc - RAMP.tricuspid.indices(RAMP.tricuspid.min_ind);
[~,ind] = min([min(a(a>0)),min(b(b>0))]);
disp(ind)
if ind == 1 
    color = 'b';
else
    color = 'r';
end
gui.l(4) = plot([RAMP.tricuspid.indices(RAMP.tricuspid.min_ind)+1 RAMP.tricuspid.indices(RAMP.tricuspid.min_ind)+1],[0,256],'Parent',gui.handles.mmd_axes,'color',color,'LineWidth',1);
else
    

end

end






function [max_type, min_type, nr_max_insp, nr_max_exp, nr_min_insp, nr_min_exp] = consolidate_results(peaks,indices,curve_min_index,curve_max_index, inflections,respiration)

mean_value = mean(peaks);
for i=1:length(indices)
    distances = [];
    for ii=1:length(inflections)
        distances(ii) = inflections(ii)-indices(i);
    end
   [min_value min_index] = min(distances(distances>0));
   value = respiration(min_value) - respiration(inflections(min_index+1));
 
    if peaks(i)<mean_value
        if value<0 %negative means inspiration
            peaks_type(i) = 0; %below mean during inspiration
        else
            peaks_type(i) = 1; %below mean during expiration
        end    
    else
        if value<0 %negative means inspiration
            peaks_type(i) = 3; %above mean during inspiration
        else
            peaks_type(i) = 4; %above mean during expiration
        end  
    end

if i == curve_min_index
    min_type = peaks_type(i);
end


if i == curve_max_index
    max_type = peaks_type(i);
end
    
AD(i) = abs(mean_value - indices(i));  
 
end
nr_min_insp = numel(peaks_type(peaks_type==0));
nr_min_exp= numel(peaks_type(peaks_type==1));
nr_max_insp= numel(peaks_type(peaks_type==2));
nr_max_exp= numel(peaks_type(peaks_type==3));

end

%Callback functions
function button_debug_callback(varargin)
global RAMP

%hFig = findobj('type', 'figure', 'name', 'YourFigureName');  % Replace 'YourFigureName' with the actual name if known

% Capture the figure with getframe
frame = getframe(RAMP.gui.fig);

% Convert the frame to an image
image = frame2im(frame);

target = '/Users/simon/Documents/SoftwareArticleRAMP/Analys20240810/screenshots/';
fullFilePath = [target, RAMP.filename,'.png']
% Save the image to a file
imwrite(image, fullFilePath);

end

function peak_callback(h, ~)
    global RAMP
    gui = RAMP.gui;
    f = ancestor(h, 'figure');
    switch f.SelectionType
        case 'normal'   % Left mouse button
            disp('Left click');
            
        case 'alt'      % Right mouse button
            disp('Right click');
            k = find(RAMP.mitral.indices == round(h.XData)); % Find the index of the peak to be deleted
            if ~isempty(k)
                % Delete the peak and its index for mitral and tricuspid
                
                  RAMP.mitral.peaks = [RAMP.mitral.peaks(1:k-1), RAMP.mitral.peaks(k+1:end)];
            RAMP.mitral.indices = [RAMP.mitral.indices(1:k-1), RAMP.mitral.indices(k+1:end)];
            RAMP.tricuspid.peaks = [RAMP.tricuspid.peaks(1:k-1), RAMP.tricuspid.peaks(k+1:end)];
            RAMP.tricuspid.indices = [RAMP.tricuspid.indices(1:k-1), RAMP.tricuspid.indices(k+1:end)];
         
                
                
                %RAMP.mitral.peaks(k) = [];
                %RAMP.mitral.indices(k) = [];
                %RAMP.tricuspid.peaks(k) = [];
                %RAMP.tricuspid.indices(k) = [];
                
                % Re-index the peaks and update the plot
                delete(gui.mit_peak(k));
                delete(gui.tri_peak(k));
                
                % Re-draw the peaks with updated indices
                plot_peaks();

                % Recalculate max/min values after deletion
                [RAMP.mitral.max_val, RAMP.mitral.min_val, RAMP.mitral.max_ind, RAMP.mitral.min_ind] = get_maxmin(RAMP.mitral.peaks);
                [RAMP.tricuspid.max_val, RAMP.tricuspid.min_val, RAMP.tricuspid.max_ind, RAMP.tricuspid.min_ind] = get_maxmin(RAMP.tricuspid.peaks);

                % Update the marker sizes to reflect changes
                update_marker_size();
                
                % Update the GUI table and M-mode lines
                update_table();
                update_mmode_lines();
            end
    end
end


function peak_callback2(h,~)
global RAMP
gui = RAMP.gui;
f = ancestor(h, 'figure');
    switch f.SelectionType
        case 'normal'   % Left mouse button
            disp('Left click');
            
        case 'alt'      % Right mouse button
            disp('Right click');
            k= find(RAMP.mitral.indices==round(h.XData));
            RAMP.mitral.peaks = [RAMP.mitral.peaks(1:k-1), RAMP.mitral.peaks(k+1:end)];
            RAMP.mitral.indices = [RAMP.mitral.indices(1:k-1), RAMP.mitral.indices(k+1:end)];
            RAMP.tricuspid.peaks = [RAMP.tricuspid.peaks(1:k-1), RAMP.tricuspid.peaks(k+1:end)];
            RAMP.tricuspid.indices = [RAMP.tricuspid.indices(1:k-1), RAMP.tricuspid.indices(k+1:end)];
           % RAMP.splits = [RAMP.splits(1:k-1), RAMP.splits(k+1:end)];
            
            
            
             disp(length(gui.mit_peak));
             delete(gui.mit_peak);
             delete(gui.tri_peak);
            plot_peaks();

        
   
                
         [RAMP.mitral.max_val,RAMP.mitral.min_val,RAMP.mitral.max_ind,RAMP.mitral.min_ind]=get_maxmin(RAMP.mitral.peaks);
        [RAMP.tricuspid.max_val,RAMP.tricuspid.min_val,RAMP.tricuspid.max_ind,RAMP.tricuspid.min_ind]=get_maxmin(RAMP.tricuspid.peaks);

        set(gui.mit_peak(RAMP.mitral.max_ind),'MarkerSize',12);
        set(gui.mit_peak(RAMP.mitral.min_ind),'MarkerSize',12);
        set(gui.tri_peak(RAMP.tricuspid.max_ind),'MarkerSize',12);
        set(gui.tri_peak(RAMP.tricuspid.min_ind),'MarkerSize',12);
       
        update_table();
        disp(length(gui.mit_peak));
    end
end

function splits_callback(h,~)
global RAMP
gui = RAMP.gui;

    f = ancestor(h, 'figure');
    switch f.SelectionType
        case 'normal'   % Left mouse button
            disp('Left click');
            set (gcf, 'WindowButtonMotionFcn', @mouseMove);
            set(gcf,"WindowButtonUpFcn", @click)
        case 'alt'      % Right mouse button
            %debug
            %disp('Right click');
    end
    
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
        [RAMP.mitral.peaks RAMP.mitral.indices,debug_curve] = find_ewave_peaks(RAMP.mitral.curve,round(RAMP.splits));
        [RAMP.tricuspid.peaks RAMP.tricuspid.indices,debug_curve] = find_ewave_peaks(RAMP.tricuspid.curve,round(RAMP.splits));
       
        delete(gui.mit_peak);
        delete(gui.tri_peak);
        
        

        plot_peaks();

        %    set([gui.handles.mit_axes gui.handles.tri_axes],'ylim',[ymin ymax]);
        [RAMP.mitral.max_val,RAMP.mitral.min_val,RAMP.mitral.max_ind,RAMP.mitral.min_ind]=get_maxmin(RAMP.mitral.peaks);
        [RAMP.tricuspid.max_val,RAMP.tricuspid.min_val,RAMP.tricuspid.max_ind,RAMP.tricuspid.min_ind]=get_maxmin(RAMP.tricuspid.peaks);

      
        update_marker_size();
        update_mmode_lines();
        
        %calculate_respiratory_curve(64)

        update_table();
    end

end


function plot_peaks()
global RAMP

gui = RAMP.gui;

gui.mit_peak = gobjects(1, length(RAMP.mitral.peaks));

for i=1:length(RAMP.mitral.peaks)
    gui.mit_peak(i) = plot(RAMP.mitral.indices(i),256+128*RAMP.mitral.curve(RAMP.mitral.indices(i)),'Marker','o','color','r','MarkerSize',8,'LineWidth',1.5,'Parent',...
        gui.handles.mit_axes,'UserData',i,'Visible',get(gui.handles.checkbox_peaks,'Value'));
    set(gui.mit_peak(i),'ButtonDownFcn', @peak_callback)
end

gui.tri_peak = gobjects(1, length(RAMP.tricuspid.peaks));

for i=1:length(RAMP.tricuspid.peaks)
    gui.tri_peak(i) = plot(RAMP.tricuspid.indices(i),256+128*RAMP.tricuspid.curve(RAMP.tricuspid.indices(i)),'Marker','o','color','r','MarkerSize',8,'LineWidth',1.5,'Parent',...
        gui.handles.tri_axes,'UserData',i,'Visible',get(gui.handles.checkbox_peaks,'Value'));
    set(gui.tri_peak(i),'ButtonDownFcn', @peak_callback)
end


end
    


function update_marker_size()
global RAMP

gui = RAMP.gui;
    


for i=1:length(gui.mit_peak)
    set(gui.mit_peak(i),'MarkerSize', 8);
    set(gui.tri_peak(i),'MarkerSize', 8); 
    set(gui.mit_peak(i),'LineWidth', 1);
    set(gui.tri_peak(i),'LineWidth', 1); 
end

set([gui.mit_peak(RAMP.mitral.max_ind), gui.mit_peak(RAMP.mitral.min_ind), gui.tri_peak(RAMP.tricuspid.max_ind), gui.tri_peak(RAMP.tricuspid.min_ind)], 'MarkerSize', 15);
set([gui.mit_peak(RAMP.mitral.max_ind), gui.mit_peak(RAMP.mitral.min_ind), gui.tri_peak(RAMP.tricuspid.max_ind), gui.tri_peak(RAMP.tricuspid.min_ind)], 'LineWidth', 2);


end


function curve3_click_callback(gcbo, event)
    global RAMP
    gui = RAMP.gui;
    x = round(event.IntersectionPoint(1, 1));  % The x-coordinate where the click occurred

    % Calculate the distances between the click location and each mitral peak index
    dist_to_mitral = abs(RAMP.mitral.indices - x);

    % Find the index of the closest peak
    [~, ind] = min(dist_to_mitral);
    
    % Ensure the click is within reasonable proximity of the nearest peak (optional)
    if min(dist_to_mitral) < 5  % Set an appropriate threshold for selecting a peak
        % Update the mitral peak's position
        if gcbo == gui.handles.mit_axes
            RAMP.mitral.indices(ind) = x;
            RAMP.mitral.peaks(ind) = RAMP.mitral.curve(x);
            set(gui.mit_peak(ind), 'XData', [x x]);
            set(gui.mit_peak(ind), 'YData', [256 + 128 * RAMP.mitral.peaks(ind) 256 + 128 * RAMP.mitral.peaks(ind)]);
        end

        % Update the tricuspid peak's position if the click occurred in the tricuspid axes
        if gcbo == gui.handles.tri_axes
            RAMP.tricuspid.indices(ind) = x;
            RAMP.tricuspid.peaks(ind) = RAMP.tricuspid.curve(x);
            set(gui.tri_peak(ind), 'XData', [x x]);
            set(gui.tri_peak(ind), 'YData', [256 + 128 * RAMP.tricuspid.peaks(ind) 256 + 128 * RAMP.tricuspid.peaks(ind)]);
        end
    end

    % Recalculate max/min values after the update
    [RAMP.mitral.max_val, RAMP.mitral.min_val, RAMP.mitral.max_ind, RAMP.mitral.min_ind] = get_maxmin(RAMP.mitral.peaks);
    [RAMP.tricuspid.max_val, RAMP.tricuspid.min_val, RAMP.tricuspid.max_ind, RAMP.tricuspid.min_ind] = get_maxmin(RAMP.tricuspid.peaks);

    % Update the GUI to reflect the new peak sizes and positions
    update_marker_size();
    update_mmode_lines();
    update_table();
end





function curve_click_callback(gcbo, event)
    global RAMP
    gui = RAMP.gui;
    x = round(event.IntersectionPoint(1, 1));

    % Calculate distances to find the closest peak segment
    dist = RAMP.splits - x;
    dist(dist < 0) = NaN;
    [~, ind] = min(dist);
    ind = ind - 1; % Adjust for the correct segment
    
    if x > RAMP.splits(1) && x < RAMP.splits(end)
        if gcbo == gui.handles.mit_axes
            RAMP.mitral.indices(ind) = x;
            RAMP.mitral.peaks(ind) = RAMP.mitral.curve(x);
            set(gui.mit_peak(ind), 'XData', [x x]);
            set(gui.mit_peak(ind), 'YData', [256 + 128 * RAMP.mitral.peaks(ind) 256 + 128 * RAMP.mitral.peaks(ind)]);
        elseif gcbo == gui.handles.tri_axes
            RAMP.tricuspid.indices(ind) = x;
            RAMP.tricuspid.peaks(ind) = RAMP.tricuspid.curve(x);
            set(gui.tri_peak(ind), 'XData', [x x]);
            set(gui.tri_peak(ind), 'YData', [256 + 128 * RAMP.tricuspid.peaks(ind) 256 + 128 * RAMP.tricuspid.peaks(ind)]);
        end
    end

    % Recalculate max/min values and update accordingly
    [RAMP.mitral.max_val, RAMP.mitral.min_val, RAMP.mitral.max_ind, RAMP.mitral.min_ind] = get_maxmin(RAMP.mitral.peaks);
    [RAMP.tricuspid.max_val, RAMP.tricuspid.min_val, RAMP.tricuspid.max_ind, RAMP.tricuspid.min_ind] = get_maxmin(RAMP.tricuspid.peaks);

    % Update the marker sizes, M-mode lines, and table
    update_marker_size();
    update_mmode_lines();
    update_table();
end




function curve_click_callback2(gcbo,event)

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
        RAMP.mitral.indices(ind) = x;
        RAMP.mitral.peaks(ind) = RAMP.mitral.curve(x);
        set(gui.mit_peak(ind),'XData',[x x]);
        set(gui.mit_peak(ind),'YData',[256+128*RAMP.mitral.peaks(ind) 256+128*RAMP.mitral.peaks(ind)]);
    end

    if gcbo == gui.handles.tri_axes
        RAMP.tricuspid.indices(ind) = x;
        RAMP.tricuspid.peaks(ind) = RAMP.tricuspid.curve(x);
        set(gui.tri_peak(ind),'XData',[x x]);
        set(gui.tri_peak(ind),'YData',[256+128*RAMP.tricuspid.peaks(ind) 256+128*RAMP.tricuspid.peaks(ind)]);
    end
end

[RAMP.mitral.max_val,RAMP.mitral.min_val,RAMP.mitral.max_ind,RAMP.mitral.min_ind]=get_maxmin(RAMP.mitral.peaks);
[RAMP.tricuspid.max_val,RAMP.tricuspid.min_val,RAMP.tricuspid.max_ind,RAMP.tricuspid.min_ind]=get_maxmin(RAMP.tricuspid.peaks);

%set([gui.mit_peak(:); gui.tri_peak(:)], 'MarkerSize', 8);

update_marker_size();
update_mmode_lines();
update_table();

end


function button_clipboard_callback(varargin)

global RAMP


if RAMP.debug == true
    set([gui.handles.panel_debug],'Visible',0); 
else
    set([gui.handles.panel_debug],'Visible',1); 
end




end

function checkboxes_callback(varargin)
global RAMP
gui = RAMP.gui;
switch varargin{1}.Tag
    case 'checkbox_vt_curve'
        set([gui.mit_curve_plot gui.tri_curve_plot],'Visible',get(varargin{1},'Value'));
    case 'checkbox_beat_seg'
        set(gui.split_line(:),'Visible',get(varargin{1},'Value'));
    case 'checkbox_mmode_helplines'
        %set(gui.mmode_lines(:),'Visible',get(varargin{1},'Value'));
        set(gui.l(:),'Visible',get(varargin{1},'Value'));
    
       
    case 'checkbox_resp_curve'
        
        set(gui.tricuspid_debug,'Visible',get(varargin{1},'Value'));
        set(gui.mitral_debug,'Visible',get(varargin{1},'Value'));
    
        %lines = findobj(gui.handles.mmd_axes, 'Type', 'line','-and','LineWidth',1.5);
        %set(lines,'Visible',get(varargin{1},'Value'));
    
    case 'checkbox_peaks'
    
        set( gui.mit_peak ,'Visible',get(varargin{1},'Value'));
        set( gui.tri_peak ,'Visible',get(varargin{1},'Value'));
end

end

function mag_click_callback(gcbo,event)
%Select what slice to use as mmode
global DATA RAMP;
gui = RAMP.gui;
pos = get(get(gco,'Parent'),'CurrentPoint');
h = gco;

t = pos(1,1);
y = round(pos(1,2));
draw_mmode(round(t(t>0)));

set(RAMP.overview_line,'Xdata',[t t]);

%[RAMP.exploc RAMP.insploc RAMP.inflections RAMP.respiration] = calculate_respiratory_curve(y);


%uistack(RAMP.overview_line,'top');
end

%Update functions that update the gui after callbacks have been triggered 
function update_mmode_lines()
global RAMP
gui = RAMP.gui;


set(gui.l(1),'XData',[RAMP.mitral.indices(RAMP.mitral.max_ind) RAMP.mitral.indices(RAMP.mitral.max_ind)]);
set(gui.l(2),'XData',[RAMP.mitral.indices(RAMP.mitral.min_ind) RAMP.mitral.indices(RAMP.mitral.min_ind)]);
set(gui.l(3),'XData',[RAMP.tricuspid.indices(RAMP.tricuspid.max_ind) RAMP.tricuspid.indices(RAMP.tricuspid.max_ind)]);
set(gui.l(4),'XData',[RAMP.tricuspid.indices(RAMP.tricuspid.min_ind) RAMP.tricuspid.indices(RAMP.tricuspid.min_ind)]);


end

function update_table()
global RAMP SET
gui = RAMP.gui;

RAMP.clip_array=[];
%[RAMP.mitral.max_type RAMP.mitral.min_type RAMP.mitral.nr_max_insp RAMP.mitral.nr_max_exp RAMP.mitral.nr_min_insp RAMP.mitral.nr_min_exp]=consolidate_results(RAMP.mitral.peaks,RAMP.mitral.indices,RAMP.mitral.min_ind,RAMP.mitral.max_ind,RAMP.inflections,RAMP.respiration);


% RAMP.clip_array(5) = RAMP.mitral.max_type;
% RAMP.clip_array(6) = RAMP.mitral.min_type;
% RAMP.clip_array(7) = RAMP.mitral.nr_max_insp;
% RAMP.clip_array(8) = RAMP.mitral.nr_max_exp;
% RAMP.clip_array(9) = RAMP.mitral.nr_min_insp;
% RAMP.clip_array(10) = RAMP.mitral.nr_min_exp;
% 
% [RAMP.tricuspid.max_type RAMP.tricuspid.min_type RAMP.tricuspid.nr_max_insp RAMP.tricuspid.nr_max_exp RAMP.tricuspid.nr_min_insp RAMP.tricuspid.nr_min_exp]=consolidate_results(RAMP.tricuspid.peaks,RAMP.tricuspid.indices,RAMP.tricuspid.min_ind,RAMP.tricuspid.max_ind,RAMP.inflections,RAMP.respiration);
% RAMP.clip_array(14) = RAMP.tricuspid.max_type;
% RAMP.clip_array(15) = RAMP.tricuspid.min_type;
% RAMP.clip_array(16) = RAMP.tricuspid.nr_max_insp;
% RAMP.clip_array(17) = RAMP.tricuspid.nr_max_exp;
% RAMP.clip_array(18) = RAMP.tricuspid.nr_min_insp;
% RAMP.clip_array(19) = RAMP.tricuspid.nr_min_exp;

%table_contents = get(gui.handles.table,'Data');

%Clipboard array to be exported


%First column is the OrigFileName from Segment file format as a number
%RAMP.clip_array(1)=str2num(RAMP.filename);

mit_max = RAMP.mitral.max_val;
mit_min = RAMP.mitral.min_val;
tri_max = RAMP.tricuspid.max_val;
tri_min = RAMP.tricuspid.min_val;

%mitral
RAMP.clip_array(2)=mit_max;
RAMP.clip_array(3)=mit_min;
RAMP.clip_array(4)=(mit_max-mit_min)/mit_max;
RAMP.clip_array(5)=tri_max;
RAMP.clip_array(6)=tri_min;
RAMP.clip_array(7)=(tri_max-tri_min)/tri_max;


set(gui.handles.result_mitral_vmax,'String',sprintf('%.2f',mit_max));
set(gui.handles.result_mitral_vmin,'String',sprintf('%.2f',mit_min));
set(gui.handles.result_mitral_variation,'String',sprintf('%d%%',round(RAMP.clip_array(4)*100)));

set(gui.handles.result_tricuspid_vmax,'String',sprintf('%.2f',tri_max));
set(gui.handles.result_tricuspid_vmin,'String',sprintf('%.2f',tri_min));
set(gui.handles.result_tricuspid_variation,'String',sprintf('%d%%',round(RAMP.clip_array(7)*100)));


%gui.handles.result_mitral_variation = sprintf('%.2f',RAMP.clip_array(4));
%table_contents(2,2) = {sprintf('%.2f',RAMP.clip_array(2))};
%table_contents(3,2) = {sprintf('%.2f',RAMP.clip_array(3))};
%table_contents(4,2) = {sprintf('%d %%',RAMP.clip_array(4))};


%table_contents(7,2) = {sprintf('%.2f',RAMP.clip_array(5))};
%table_contents(8,2) = {sprintf('%.2f',RAMP.clip_array(6))};
%table_contents(9,2) = {sprintf('%d %%',RAMP.clip_array(7))};

% %The individual peaks
% 
% %for i=1:length(RAMP.mitral.peaks)
% 
% %    RAMP.clip_array(19%+i) = RAMP.mitral.peaks(i);
% %end
% 
% %for i=1:length(RAMP.tricuspid.peaks)
% 
%     RAMP.clip_array(19+length(RAMP.mitral.peaks)+i+1) = RAMP.tricuspid.peaks(i);
% 
% end
number2clipboard(RAMP.clip_array);
%filename	mit_max	mit_min	mit_var	mit_max_type	mit_min_type	mit_nr_max_insp	mit_nr_max_exp	mit_nr_min_insp	mit_nr_min_exp	tri_min	tri_max	tri_var	tri_max_type	tri_min_type	tri_nr_max_insp	tri_nr_max_exp	tri_nr_min_insp	tri_nr_min_exp	individual peaks separated by 0																																																																												
end

%Utility funct
function arraystring = number2clipboard(array)
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
%arraystring(double(arraystring)==0) = '';
clipboard('copy',arraystring); %copy the result to the clipboard ready for pasting

end

