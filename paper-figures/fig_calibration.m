

clear l

%%
% This makes a figure showing calibration curves for different odorants 

dm = dataManager;



c = lines(4);

 ;;;;;;;     ;;;     ;;;;;;  
;;     ;;   ;; ;;   ;;    ;; 
       ;;  ;;   ;;  ;;       
 ;;;;;;;  ;;     ;; ;;       
;;        ;;;;;;;;; ;;       
;;        ;;     ;; ;;    ;; 
;;;;;;;;; ;;     ;;  ;;;;;;  

addpath('../src')
load(dm.getPath('55091fea105ec2c5f802d71a06059478'))

pc.calibrate;
handles = pc.handles;

handles.ax.r2 = subplot(4,3,12); hold on
handles.ax.r2.Position(2) = .08;
l(1) = plot(handles.ax.r2,1/pc.PID_to_odorant_flux_fit.p1,pc.PID_to_odorant_flux_fit_r2,'r+');
xlabel(handles.ax.r2,'PID sensitivity (V/\mu mol/s)')
ylabel(handles.ax.r2,'r^2 of fit')



;;;;;;;; ;;;;;;;; ;;     ;;    ;;;    ;;    ;;  ;;;;;;;  ;;       
;;          ;;    ;;     ;;   ;; ;;   ;;;   ;; ;;     ;; ;;       
;;          ;;    ;;     ;;  ;;   ;;  ;;;;  ;; ;;     ;; ;;       
;;;;;;      ;;    ;;;;;;;;; ;;     ;; ;; ;; ;; ;;     ;; ;;       
;;          ;;    ;;     ;; ;;;;;;;;; ;;  ;;;; ;;     ;; ;;       
;;          ;;    ;;     ;; ;;     ;; ;;   ;;; ;;     ;; ;;       
;;;;;;;;    ;;    ;;     ;; ;;     ;; ;;    ;;  ;;;;;;;  ;;;;;;;; 


load(dm.getPath('b8ac3ac9f5a97e24fe4bf6be9680552e'))
pc.calibrate(false);
x = [pc.calibration_data.mean_PID];
y = [pc.calibration_data.mean_flux]*1e6;

plot(handles.ax.PID_to_flux,x,y,'bo')
plot(handles.ax.PID_to_flux,x,pc.PID_to_odorant_flux_fit(x),'b')

l(2) = plot(handles.ax.r2,1/pc.PID_to_odorant_flux_fit.p1,pc.PID_to_odorant_flux_fit_r2,'bo');


 ;;;;;;;  ;;;;;;;;  ;;     ;; ;;;;;;;; 
;;     ;; ;;     ;; ;;     ;;    ;;    
       ;; ;;     ;; ;;     ;;    ;;    
 ;;;;;;;  ;;;;;;;;  ;;     ;;    ;;    
;;        ;;     ;; ;;     ;;    ;;    
;;        ;;     ;; ;;     ;;    ;;    
;;;;;;;;; ;;;;;;;;   ;;;;;;;     ;;    


load(dm.getPath('f59fd9547e84279a927e615a397cebe3'))

pc.calibrate(false)
x = [pc.calibration_data.mean_PID];
y = [pc.calibration_data.mean_flux]*1e6;

plot(handles.ax.PID_to_flux,x,y,'gd')
plot(handles.ax.PID_to_flux,x,pc.PID_to_odorant_flux_fit(x),'g')

l(3) = plot(handles.ax.r2,1/pc.PID_to_odorant_flux_fit.p1,pc.PID_to_odorant_flux_fit_r2,'gd');

handles.ax.PID_to_flux.XLim = [.01 14];
handles.ax.r2.YLim = [0.5 1];
handles.ax.r2.XLim = [0 5];

;;;;;;;;    ;;;     ;;;;;;  
;;         ;; ;;   ;;    ;; 
;;        ;;   ;;  ;;       
;;;;;;;  ;;     ;; ;;       
      ;; ;;;;;;;;; ;;       
;;    ;; ;;     ;; ;;    ;; 
 ;;;;;;  ;;     ;;  ;;;;;;  

load(dm.getPath('b64f9de37719ba493c380983b7035859'))
pc.calibrate(false)
x = [pc.calibration_data.mean_PID];
y = [pc.calibration_data.mean_flux]*1e6;

xx = linspace(0,10,100);
plot(handles.ax.PID_to_flux,x,y,'ms')
plot(handles.ax.PID_to_flux,xx,pc.PID_to_odorant_flux_fit(xx),'m')

l(4) = plot(handles.ax.r2,1/pc.PID_to_odorant_flux_fit.p1,pc.PID_to_odorant_flux_fit_r2,'ms');

handles.ax.PID_to_flux.XLim = [.01 14];
handles.ax.r2.YLim = [0.5 1];
handles.ax.r2.XLim = [0 5];

legend(l,{'ethyl acetate','ethanol','2-butanone','pentyl acetate'},'Location','southeast')

axlib.label(handles.ax.r2,'i','x_offset',0.01,'y_offset',-.01)


figlib.pretty('FontSize',12);

handles.ax.r2.Box = 'off';
