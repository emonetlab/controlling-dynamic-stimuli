
addpath('../src/')
dm = dataManager;

load(dm.getPath('63facb76b4a9543c8ced04d4030ebf54'));

close all

% Create figure
main_fig = figure('Position',[1  1  1239 1293]);
clear ax
ax.cartoon = axes('Parent',main_fig,'Position',[0.06 0.82 0.25 0.14]); hold on


ax.dist(1) = axes('Parent',main_fig,'Position',[0.42 0.83 0.14 0.14]); hold on
ax.dist(2) = axes('Parent',main_fig,'Position',[0.62 0.83 0.14 0.14]); hold on

ax.flash(1) = axes('Parent',main_fig,'Position',[0.10 0.62 0.13 0.14]); hold on
ax.flash(2) = axes('Parent',main_fig,'Position',[0.365 0.62 0.13 0.14]);hold on
ax.flash(3) = axes('Parent',main_fig,'Position',[0.62 0.62 0.13 0.14]);hold on


ax.step1_mfc = axes('Parent',main_fig,'Position',[0.09 0.49 0.5 0.08]);hold on
ax.step1_pid = axes('Parent',main_fig,'Position',[0.09 0.42 0.5 0.08]);hold on
ax.step2_mfc = axes('Parent',main_fig,'Position',[0.09 0.31 0.5 0.08]);hold on
ax.step2_pid = axes('Parent',main_fig,'Position',[0.09 0.23 0.5 0.08]);hold on
ax.final_pid = axes('Parent',main_fig,'Position',[0.09 0.06 0.5 0.08]);hold on
ax.inset = axes('Parent',main_fig,'Position',[.36 .1 .08 .08]);hold on
ax.step1_pid_r2 = axes('Parent',main_fig,'Position',[0.66 0.44 0.1 0.1]);hold on
ax.step2_pid_r2 = axes('Parent',main_fig,'Position',[0.66 0.26 0.1 0.1]);hold on
ax.final_cv = axes('Parent',main_fig,'Position',[0.66 0.06 0.1 0.1]);hold on


I = imread('../images/nat_stim.png');
figlib.showImageInAxes(ax.cartoon,I)


% make a base axes to draw stuff on 
ax.base = axes;
ax.base.Position = [0 0  1 1];
axis(ax.base,'off')
uistack(ax.base,'bottom');




% show dilution -> PID map
dil_levels = logspace(log10(nsb.min_dilution),log10(nsb.max_dilution),nsb.n_levels);
x = []; y = [];
for i = 1:nsb.n_levels
	x = [x dil_levels(i)*ones(1,nsb.n_reps)];
	y = [y max(nsb.dose_response_data(i).PID')];
end
plot(ax.flash(1),x,y,'k+')
ff_peak_PID = fit([zeros(10,1); x(:)],[zeros(10,1); y(:)],'smoothingspline','SmoothingParam',.999);
plot(ax.flash(1),x,ff_peak_PID(x),'r')
xlabel(ax.flash(1),'Diluton (%)')
ylabel(ax.flash(1),'PID (V)')

% show dil -> peak time
x = []; y = [];
for i = 1:nsb.n_levels
	[~,loc] = max(nsb.dose_response_data(i).PID');
	loc = loc - floor(nsb.pre_pulse_buffer*nsb.sampling_rate);
	loc = loc/nsb.sampling_rate;
	x = [x dil_levels(i)*ones(1,nsb.n_reps)];
	y = [y loc];
end
plot(ax.flash(2),x,y,'k+')
y(y<0) = NaN;
y(y>3*nanmean(y)) = NaN;
x(isnan(y)) = [];
y(isnan(y)) = [];
ff_time_to_peak = fit(x(:),y(:),'poly4');
plot(ax.flash(2),x,ff_time_to_peak(x),'r')

xlabel(ax.flash(2),'Dilution (%)')
ylabel(ax.flash(2),'Time to peak stim (s)')



% compute whiff durations based on threshold crfilesepings 
x = []; 
t_whiff_on = [];
t_whiff_off = [];
t_whiff = [];
for i = 1:nsb.n_levels
	S = nsb.dose_response_data(i).PID';
	for j = 1:size(S,2)
		x = [x dil_levels(i)];
		temp1 = find(nsb.dose_response_data(i).PID(j,:)>nsb.c0,1,'first');
		if isempty(temp1)
			temp1 = NaN;
		end
		temp1 = temp1 - floor(nsb.pre_pulse_buffer*nsb.sampling_rate);
		temp1 = temp1/nsb.sampling_rate;
		t_whiff_on = [t_whiff_on temp1];

		temp2 = find(nsb.dose_response_data(i).PID(j,:)>nsb.c0,1,'last');
		if isempty(temp2)
			temp2 = NaN;
		end
		temp2 = temp2 - floor(nsb.pre_pulse_buffer*nsb.sampling_rate);
		temp2 = temp2/nsb.sampling_rate;
		t_whiff_off = [t_whiff_off temp2];

		t_whiff = [t_whiff (temp2 - temp1)];
	end
end

x = x(:);
t_whiff_on = t_whiff_on(:);
t_whiff_off = t_whiff_off(:);
t_whiff = t_whiff(:);

plot(ax.flash(3),x,t_whiff,'k+')
xlabel(ax.flash(3),'Dilution (%)')
ylabel(ax.flash(3),'Whiff duration (s)')
ff_time_whiff = fit(x(~isnan(t_whiff)),t_whiff(~isnan(t_whiff)),'smoothingspline','SmoothingParam',.9);
plot(ax.flash(3),x,ff_time_whiff(x),'r')

% show distributions
[S,pS, t_blanks, pt_blanks] = genPDF(nsb);
plot(ax.dist(1),pS,S,'k')
set(ax.dist(1),'XScale','log','YScale','log')
xlabel(ax.dist(1),'Stimulus Amplitude (V)')
ylabel(ax.dist(1),'p.d.f')

plot(ax.dist(2),pt_blanks,t_blanks,'k')
set(ax.dist(2),'XScale','log','YScale','log')
xlabel(ax.dist(2),'Blank duration (ms)')
ylabel(ax.dist(2),'p.d.f')




% time = (1:length(nsb.dose_response_data(5).MFC(3,:)))*1e-4;
% plot(ax.mfc_ex,time,nsb.dose_response_data(5).MFC(3,:),'k')
% plot(ax.mfc_ex,time,nsb.dose_response_data(2).MFC(3,:),'b')
% set(ax.mfc_ex,'XLim',[.9 1.4])
% axis(ax.mfc_ex,'off')
% ylabel(ax.mfc_ex,'MFC')

% valve_signal = time*0;
% valve_signal(1e4:1e4+500) = 1;
% plot(ax.valve_ex,time,valve_signal,'k')
% set(ax.valve_ex,'XLim',[.9 1.4],'YLim',[-.1 1.1])
% axis(ax.valve_ex,'off')
% ylabel(ax.valve_ex,'Valve')

% plot(ax.pid_ex,time,nsb.dose_response_data(5).PID(3,:),'k')
% plot(ax.pid_ex,time,nsb.dose_response_data(2).PID(3,:),'b')
% set(ax.pid_ex,'XLim',[.9 1.4])
% axis(ax.pid_ex,'off')
% ylabel(ax.pid_ex,'PID')


% load final data
load(dm.getPath('046b595b034b7afdd709fc2194479987'),'-mat')
PID = data.PID(:,1:10:end)';
dt = 1e-3;
time = dt:dt:100;

plot(ax.final_pid,time,PID,'r')
xlabel(ax.final_pid,'Time (s)')
ylabel(ax.final_pid,'PID (V)')
set(ax.final_pid,'YLim',[-.5 7],'XLim',[0 100])

% also plot in the inset
plot(ax.inset,time,PID,'r')
set(ax.inset,'XLim',[43.52 43.66])
ax.inset.YLim = [-.1 3.5];

% show variability
valve_ons = floor(veclib.computeOnsOffs(nsb.valve_ansatz)/10);
v = NaN*valve_ons;
mu = NaN*valve_ons;


for i = 1:length(valve_ons)
	whiff_peaks = max(PID(valve_ons(i)+30:valve_ons(i)+110,:));
	v(i) = std(whiff_peaks)./mean(whiff_peaks);
	mu(i) = mean(whiff_peaks);
end

scatter(ax.final_cv,mu,v,64,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.5)
xlabel(ax.final_cv,'Whiff peak (V)')
ylabel(ax.final_cv,'CV')
set(ax.final_cv,'XScale','log','box','off','XLim',[1e-2 10],'YLim',[0.01 1],'YScale','log')


;;;;;;;; ;;;; ;;;;;;;;   ;;;;;;  ;;;;;;;; 
;;        ;;  ;;     ;; ;;    ;;    ;;    
;;        ;;  ;;     ;; ;;          ;;    
;;;;;;    ;;  ;;;;;;;;   ;;;;;;     ;;    
;;        ;;  ;;   ;;         ;;    ;;    
;;        ;;  ;;    ;;  ;;    ;;    ;;    
;;       ;;;; ;;     ;;  ;;;;;;     ;;    

 ;;;;;;  ;;;;;;;; ;;;;;;;; ;;;;;;;;  
;;    ;;    ;;    ;;       ;;     ;; 
;;          ;;    ;;       ;;     ;; 
 ;;;;;;     ;;    ;;;;;;   ;;;;;;;;  
      ;;    ;;    ;;       ;;        
;;    ;;    ;;    ;;       ;;        
 ;;;;;;     ;;    ;;;;;;;; ;;        

%  -- PID
dt = nsb.T/length(nsb.optimization_data(1).PID);
time = dt:dt:nsb.T;
plot(ax.step1_pid,time,nsb.optimization_data(1).PID,'r')
set(ax.step1_pid,'YLim',[-.5 7],'XLim',[0 100])

% MFC
dt = nsb.T/length(nsb.optimization_data(1).PID);
time = dt:dt:nsb.T;
plot(ax.step1_mfc,time,nsb.optimization_data(1).mfc1_flow*nsb.mfc1_units,'b')
ylabel(ax.step1_mfc,{'MFC Flow','(mL/min)'})


% how good is this?
[desired_whiff_intensity, actual_whiff_intensity] = nsb.measureWhiffAmplitudes(nsb.optimization_data(1).PID);
plot(ax.step1_pid_r2,[1e-5 1e5],[1e-5 1e5],'k--')
scatter(ax.step1_pid_r2,desired_whiff_intensity,actual_whiff_intensity,64,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.5)
xlabel(ax.step1_pid_r2,'Target whiff amplitude (V)')
ylabel(ax.step1_pid_r2,'Actual whiff amplitude (V)')
set(ax.step1_pid_r2,'XLim',[1e-2 10],'YLim',[1e-2 10],'YScale','log','XScale','log')
r2 = corr(desired_whiff_intensity, actual_whiff_intensity)^2;
title(ax.step1_pid_r2,['r^2 = ' strlib.oval(r2)])
ylabel(ax.step1_pid,'PID (V)')


;;          ;;;     ;;;;;;  ;;;;;;;; 
;;         ;; ;;   ;;    ;;    ;;    
;;        ;;   ;;  ;;          ;;    
;;       ;;     ;;  ;;;;;;     ;;    
;;       ;;;;;;;;;       ;;    ;;    
;;       ;;     ;; ;;    ;;    ;;    
;;;;;;;; ;;     ;;  ;;;;;;     ;;    


 ;;;;;;  ;;;;;;;; ;;;;;;;; ;;;;;;;;  
;;    ;;    ;;    ;;       ;;     ;; 
;;          ;;    ;;       ;;     ;; 
 ;;;;;;     ;;    ;;;;;;   ;;;;;;;;  
      ;;    ;;    ;;       ;;        
;;    ;;    ;;    ;;       ;;        
 ;;;;;;     ;;    ;;;;;;;; ;;        

%  -- PID
dt = nsb.T/length(nsb.optimization_data(end).PID);
time = dt:dt:nsb.T;
plot(ax.step2_pid,time,nsb.optimization_data(end).PID,'r')
set(ax.step2_pid,'YLim',[-.5 7],'XLim',[0 100])

% MFC
dt = nsb.T/length(nsb.optimization_data(end).PID);
time = dt:dt:nsb.T;
plot(ax.step2_mfc,time,nsb.optimization_data(end).mfc1_flow*nsb.mfc1_units,'b')
ylabel(ax.step2_mfc,{'MFC Flow','(mL/min)'})
% how good is this?
[desired_whiff_intensity, actual_whiff_intensity] = nsb.measureWhiffAmplitudes(nsb.optimization_data(end).PID);
plot(ax.step2_pid_r2,[1e-5 1e5],[1e-5 1e5],'k--')
scatter(ax.step2_pid_r2,desired_whiff_intensity,actual_whiff_intensity,64,'MarkerFaceColor',[1 0 0],'MarkerEdgeColor',[1 0 0],'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.5)
xlabel(ax.step2_pid_r2,'Target whiff amplitude (V)')
ylabel(ax.step2_pid_r2,'Actual whiff amplitude (V)')
set(ax.step2_pid_r2,'XLim',[1e-2 10],'YLim',[1e-2 10],'YScale','log','XScale','log')
r2 = corr(desired_whiff_intensity, actual_whiff_intensity)^2;
title(ax.step2_pid_r2,['r^2 = ' strlib.oval(r2)])
ylabel(ax.step2_pid,'PID (V)')




 ;;;;;;   ;;;;;;;   ;;;;;;  ;;     ;; ;;;;;;;; ;;;;;;;; ;;;;  ;;;;;;  
;;    ;; ;;     ;; ;;    ;; ;;;   ;;; ;;          ;;     ;;  ;;    ;; 
;;       ;;     ;; ;;       ;;;; ;;;; ;;          ;;     ;;  ;;       
;;       ;;     ;;  ;;;;;;  ;; ;;; ;; ;;;;;;      ;;     ;;  ;;       
;;       ;;     ;;       ;; ;;     ;; ;;          ;;     ;;  ;;       
;;    ;; ;;     ;; ;;    ;; ;;     ;; ;;          ;;     ;;  ;;    ;; 
 ;;;;;;   ;;;;;;;   ;;;;;;  ;;     ;; ;;;;;;;;    ;;    ;;;;  ;;;;;;  


figlib.pretty('PlotLineWidth',1.5,'LineWidth',1.5,'FontSize',13);
axlib.separate(ax.dist);




% draw some boxes
ax.base.YLim = [0 1];
ax.base.XLim = [0 1];
clear r

% for target distributions
curvature = .15;
r(1) = rectangle(ax.base,'Position',[0 0 .5 .5],'Curvature',curvature);
r(1).Position = [.35 .79 .44 .2];
r(1).LineWidth = 3;
c = lines;
r(1).EdgeColor = [1 1 1];
r(1).FaceColor = [c(1,:) .3];

% for flash reponse
r(2) = rectangle(ax.base,'Position',[0 0 .5 .5],'Curvature',curvature);
r(2).Position = [.03 .58 .76 .2];
r(2).LineWidth = 3;
r(2).EdgeColor = [1 1 1];
r(2).FaceColor = [c(3,:) .3];


% for step 1
r(3) = rectangle(ax.base,'Position',[0 0 .5 .5],'Curvature',curvature);
r(3).Position = [.03 .39 .76 .18];
r(3).LineWidth = 3;
r(3).EdgeColor = [1 1 1];
r(3).FaceColor = [c(4,:) .3];



% for step 10
r(4) = rectangle(ax.base,'Position',[0 0 .5 .5],'Curvature',curvature);
r(4).Position = [.03 .2 .76 .18];
r(4).LineWidth = 3;
r(4).EdgeColor = [1 1 1];
r(4).FaceColor = [c(4,:) .3];


set(ax.flash(1),'XScale','log','YScale','log','XLim',[nsb.min_dilution 105])
set(ax.flash(2),'XScale','log','YLim',[0 .3],'XLim',[nsb.min_dilution 105])
set(ax.flash(3),'XScale','log','XLim',[nsb.min_dilution 105])

set(ax.step1_mfc,'XColor','k','YLim',[0 260])
set(ax.step2_mfc,'XColor','k','YLim',[0 260])

ax.step1_pid.XColor = 'w';
ax.step1_mfc.XTick = [];
ax.step1_pid.XTick = [];
ax.step2_pid.XColor = 'w';
ax.step2_pid.XTick = [];
ax.step2_mfc.XColor = 'k';
ax.step2_mfc.XTick = [];

ax.step1_pid.Position = [.09 .395 .5 .08];
ax.step1_mfc.Position = [0.0900 0.4800 0.5000 0.0800];
ax.step2_pid.Position(2) = .205;
ax.step2_mfc.Position(2) = .29;
ax.step2_pid_r2.Position(2) = .25;

ax.inset.Box = 'on';

% axlib.label(ax.cartoon,'a','y_offset',-.01,'x_offset',-.01);