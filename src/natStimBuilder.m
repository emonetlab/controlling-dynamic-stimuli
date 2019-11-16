% class that helps you construct naturalistic stimuli
% natStimBuilder
% class to help you build "naturalistic" odorant stimuli
% with sparse whiffs of varying heights
% natStimBuilder assumes that you have a 2MFC set up to dilute odorants


classdef natStimBuilder < handle

	properties
		% output config -- should match kontroller config and physical reality
		main_air_channel = 1;
		mfc1_channel = 2;  % goes to odorant
		mfc2_channel = 3;  % dilutes odorant with air
		valve_channel = 4;

		% data stores: stores all the data during running this 
		dose_response_data
		optimization_data
		ControlParadigm
		mfc1_ansatz
		valve_ansatz
		target_PID
		tweak_these
        n_successful_hits

		% ~~~~~~~~~~~~~~~ parameters for MFCs ~~~~~~~~~~~~~`
		min_dilution = 1; % in percent 
		max_dilution = 99; % in percent
		n_levels = 10; % how many pulses do you want to use to build the dose-response?
		main_air_flow = 2000; % mL/min
		secondary_air_flow = 200; % mL/min
		n_reps = 3; 

		mfc1_units = 100; % mL/min/V
		mfc2_units = 100; % mL/min/V
		mfc_main_air_units = 1000; % mL/min/V

		% ~~~~~~~~~~~~~~~~ timing parameters ~~~~~~~~~~~
		pulse_duration = 50e-3; % seconds 
		pre_pulse_buffer = 1; % how long do you want the MFC to keep a setpoint before the valve turns on? (seconds) applies only to the initial dose-response estimation 
		sampling_rate = 10000; % let's keep it the same as ephys for simplicty 
		T = 100; % seconds, which is how long our "natural" stimulus will be
		min_blank_duration = .1; % blanks below this duration are not allowed

		

	
		c0 = .1; % signal threshold, used in computing statistics 
		c_max = 7; % maximum signal we expect, used in computing the normalisation factors for power law distributions.

		handles

		random_seed = 1441025
		optimization_step = 0

		% function fits
		ff_peak_PID
		ff_time_to_peak
		ff_time_whiff_on
		ff_time_whiff_off
		ff_time_whiff

	end % end props 

	methods

		function nsb = natStimBuilder()
			% check if there is a saved object on disk
			if exist('natural_stimulus_builder.mat','file') == 2
				load('natural_stimulus_builder.mat')
			end
		end % init function 

		function estimateDoseResponse(nsb)
			handles.dose_response_fig = figure('outerposition',[0 0 1200 800],'PaperUnits','points','PaperSize',[1200 800]); hold on

			for i = 1:6
				handles.ax(i) = subplot(2,3,i); hold on
			end
			xlabel(handles.ax(1),'Time')
			ylabel(handles.ax(1),'PID (V)')
			handles.dose_response_data = plot(handles.ax(2),NaN,NaN,'k+');
			dil_levels = logspace(log10(nsb.min_dilution),log10(nsb.max_dilution),nsb.n_levels);
			c = parula(nsb.n_levels);

			if isempty(nsb.dose_response_data)

				ControlParadigm = buildControlParadigm(nsb);

				% first turn the MFC off and measure the "background"
				temp = ControlParadigm(1);
				temp.Outputs(nsb.valve_channel,:) = 0;
				temp.Outputs(nsb.mfc1_channel,:) = 0;
				temp.Outputs(nsb.mfc2_channel,:) = 0;
				zero_data = kontroller('ControlParadigm',temp,'RunTheseParadigms',1,'w',nsb.sampling_rate);
				min_PID = min(zero_data.PID);
				std_PID = std(zero_data.PID);

				for i = 1:nsb.n_levels
					data(i) = kontroller('ControlParadigm',ControlParadigm(i),'RunTheseParadigms',ones(nsb.n_reps,1),'w',nsb.sampling_rate);
					plot(handles.ax(1),data(i).PID','Color',c(i,:));
					drawnow
					data(i).PID = data(i).PID - min_PID;
					handles.dose_response_data.XData = [handles.dose_response_data.XData dil_levels(i)*ones(1,nsb.n_reps)];
					handles.dose_response_data.YData = [handles.dose_response_data.YData max(data(i).PID')];
					drawnow 
				end

				nsb.dose_response_data = data;
			end
			cla(handles.ax(1))
			time = (1:length(nsb.dose_response_data(1).PID))/nsb.sampling_rate;
			for i = 1:nsb.n_levels
				plot(handles.ax(1),time,nsb.dose_response_data(i).PID','Color',c(i,:));
			end

			% plot peak PID for each case
			x = []; y = [];
			for i = 1:nsb.n_levels
				x = [x dil_levels(i)*ones(1,nsb.n_reps)];
				y = [y max(nsb.dose_response_data(i).PID')];
			end
			plot(handles.ax(2),x,y,'k+')
			ff_peak_PID = fit([zeros(10,1); x(:)],[zeros(10,1); y(:)],'smoothingspline','SmoothingParam',.999);
			plot(handles.ax(2),x,ff_peak_PID(x),'r')
			xlabel(handles.ax(2),'Diluton (%)')
			ylabel(handles.ax(2),'PID (V)')
			set(handles.ax(2),'XScale','log','YScale','log','XLim',[nsb.min_dilution 105])

			% compute the time to peak, time of threshold crfileseping,etc. 
			x = []; y = [];
			for i = 1:nsb.n_levels
				[~,loc] = max(nsb.dose_response_data(i).PID');
				loc = loc - floor(nsb.pre_pulse_buffer*nsb.sampling_rate);
				loc = loc/nsb.sampling_rate;
				x = [x dil_levels(i)*ones(1,nsb.n_reps)];
				y = [y loc];
			end
			plot(handles.ax(3),x,y,'k+')
			y(y<0) = NaN;
			y(y>3*nanmean(y)) = NaN;
			x(isnan(y)) = [];
			y(isnan(y)) = [];
			ff_time_to_peak = fit(x(:),y(:),'poly4');
			plot(handles.ax(3),x,ff_time_to_peak(x),'r')
			set(handles.ax(3),'XScale','log','YLim',[0 .3],'XLim',[nsb.min_dilution 105])
			xlabel(handles.ax(3),'Dilution (%)')
			ylabel(handles.ax(3),'Time to peak stim (s)')

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

			plot(handles.ax(4),x,t_whiff_on,'k+')
			set(handles.ax(4),'XScale','log','XLim',[nsb.min_dilution 105])
			xlabel(handles.ax(4),'Dilution (%)')
			ylabel(handles.ax(4),'Time to crfilesep threshold (s)')
			ff_time_whiff_on = fit(x(~isnan(t_whiff_on)),t_whiff_on(~isnan(t_whiff_on)),'smoothingspline','SmoothingParam',.99);
			plot(handles.ax(4),x,ff_time_whiff_on(x),'r')

			plot(handles.ax(5),x,t_whiff_off,'k+')
			set(handles.ax(5),'XScale','log','XLim',[nsb.min_dilution 105])
			xlabel(handles.ax(5),'Dilution (%)')
			ylabel(handles.ax(5),'Time to return below threshold (s)')
			ff_time_whiff_off = fit(x(~isnan(t_whiff_off)),t_whiff_off(~isnan(t_whiff_off)),'smoothingspline','SmoothingParam',.9);
			plot(handles.ax(5),x,ff_time_whiff_off(x),'r')

			plot(handles.ax(6),x,t_whiff,'k+')
			set(handles.ax(6),'XScale','log','XLim',[nsb.min_dilution 105])
			xlabel(handles.ax(6),'Dilution (%)')
			ylabel(handles.ax(6),'Whiff duration (s)')
			ff_time_whiff = fit(x(~isnan(t_whiff)),t_whiff(~isnan(t_whiff)),'smoothingspline','SmoothingParam',.9);
			plot(handles.ax(6),x,ff_time_whiff(x),'r')

			figlib.pretty('fs',16);
			
			% save the fits for later use
			nsb.ff_peak_PID = ff_peak_PID;
			nsb.ff_time_whiff_on = ff_time_whiff_on;
			nsb.ff_time_whiff = ff_time_whiff;
			nsb.ff_time_whiff_on = ff_time_whiff_on;
			nsb.ff_time_whiff_off = ff_time_whiff_off;
			nsb.ff_time_to_peak = ff_time_to_peak;

			nsb.handles = handles;
            
            nsb.c_max = max(nsb.dose_response_data(end).PID(:));

            % minor cosmetic fixes

            nsb.handles.ax(1).Position(1) = .1;
            nsb.handles.ax(4).Position(1) = .1;
            nsb.handles.ax(3).Position(1) = .75;
            nsb.handles.ax(6).Position(1) = .75;

            nsb.handles.ax(1).XLim = [0.5 2];
            nsb.handles.ax(4).YLim(1) = 0;
            nsb.handles.ax(5).YLim(1) = 0;
            nsb.handles.ax(6).YLim(1) = 0;

            figlib.label('x_offset',-.01,'font_size',20);

            deintersectAxes(nsb.handles.ax(2))
            deintersectAxes(nsb.handles.ax(3))
            deintersectAxes(nsb.handles.ax(4))
            deintersectAxes(nsb.handles.ax(5))
            deintersectAxes(nsb.handles.ax(6))

		end % end estimateDoseResponse


		function ControlParadigm = buildControlParadigm(nsb)
			all_dilutions = logspace(log10(nsb.min_dilution),log10(nsb.max_dilution),nsb.n_levels);
			all_dilutions = all_dilutions/100;
			mfc1_setpoints = ((all_dilutions)*nsb.secondary_air_flow)/nsb.mfc1_units;
			mfc2_setpoints = ((1 - all_dilutions)*nsb.secondary_air_flow)/nsb.mfc2_units;
			for i = 1:nsb.n_levels
				ControlParadigm(i).Name = ['mfc:' strlib.oval(100*all_dilutions(i))];
				ControlParadigm(i).Outputs = zeros(4,floor(5*nsb.sampling_rate));
				ControlParadigm(i).Outputs(nsb.main_air_channel,:) = nsb.main_air_flow/nsb.mfc_main_air_units;
				ControlParadigm(i).Outputs(nsb.mfc1_channel,:) = mfc1_setpoints(i);
				ControlParadigm(i).Outputs(nsb.mfc2_channel,:) = mfc2_setpoints(i);
				a = floor(nsb.pre_pulse_buffer*nsb.sampling_rate);
				z = a + floor(nsb.pulse_duration*nsb.sampling_rate);
				ControlParadigm(i).Outputs(nsb.valve_channel,a:z) = 1;
			end
		end % end buildControlParadigm


		function buildNatStim(nsb)

			% online function that uses a mixture of numerical calculations and 
			% live experiments to build and iteratively tune control paradigms to 
			% generate natural stimuli with desired stimuli.

			if ~isempty(nsb.handles) && any(strcmp(fieldnames(nsb.handles),'reroll_button'))
				% clear all the axes
				try
					nsb.handles.mfc_plot1.XData = NaN;
					nsb.handles.mfc_plot2.XData = NaN;
					cla(nsb.handles.ax(4))
					cla(nsb.handles.ax(5))
				catch
				end
			else
				nsb.buildNatStimUI;
			end

		
			% generate the pdfs for whiff intensity and blank duration 
			[S,pS, t_blanks, pt_blanks] = genPDF(nsb);

			% sample from these distributions
			whiff_amplitudes = sampleFromDist(nsb,S,pS);
			blank_durations = sampleFromDist(nsb,t_blanks,pt_blanks);

			plot(nsb.handles.ax(4),S,pS,'k')
			xlabel(nsb.handles.ax(4),'Whiff intensity (V)')
			ylabel(nsb.handles.ax(4),'p.d.f')
			set(nsb.handles.ax(4),'YScale','log','XScale','log')

			plot(nsb.handles.ax(5),t_blanks,pt_blanks,'k')
			xlabel(nsb.handles.ax(5),'Blank durations (s)')
			ylabel(nsb.handles.ax(5),'p.d.f')
			set(nsb.handles.ax(5),'YScale','log','XScale','log')

			% incrementally build the simulated PID signal, the MFC control signal, and the valve signal
			PID = 0*(1:(nsb.T*nsb.sampling_rate));
			mfc1 = 0*PID;
			mfc2 = 0*PID;
			valve = 0*PID;

			goon = 1;
			a = 5*nsb.sampling_rate;
			previous_valve_off = 1;
			n_whiffs = 0;

			nsb.handles.mfc_plot2.XData = (1:length(PID))/nsb.sampling_rate;
			nsb.handles.mfc_plot1.XData = (1:length(PID))/nsb.sampling_rate;
            nsb.handles.mfc_plot3.XData = (1:length(PID))/nsb.sampling_rate;
			nsb.handles.mfc_plot2.YData = NaN*(1:length(PID));
			nsb.handles.mfc_plot1.YData = NaN*(1:length(PID));
			y2 = NaN*(1:length(PID));

			while goon
				% pick a whiff_amplitude
				this_whiff_amplitude = whiff_amplitudes(1);
				whiff_amplitudes(1) = [];

				% find the dilution to achieve this PID target
				desired_dilution = cfinvert(nsb.ff_peak_PID,this_whiff_amplitude,[0 100],1e3);

				if isnan(desired_dilution)
					warning('Cannot invert MFC -> PID function. Your c_max is probably too high')
					desired_dilution = 100;
				end

				% estimate time to threshold crossing 
				estimated_time_to_threshold_crfileseping = nsb.ff_time_whiff_on(desired_dilution);

				% subtract this from the the valve opening time
				a_PID = a;
				z_PID = a + floor(nsb.pulse_duration*nsb.sampling_rate);
				a = a - floor(nsb.sampling_rate*estimated_time_to_threshold_crfileseping);
				z = a + floor(nsb.pulse_duration*nsb.sampling_rate);

				PID(a_PID:z_PID) = this_whiff_amplitude;

				% move the MFCs to the desired setpoints
				mfc1_setpoint = ((desired_dilution)/100*nsb.secondary_air_flow)/nsb.mfc1_units;
				mfc2_setpoint = (((100 - desired_dilution)/100)*nsb.secondary_air_flow)/nsb.mfc2_units;
				mfc1(previous_valve_off:z) = mfc1_setpoint;
				mfc2(previous_valve_off:z) = mfc2_setpoint;

				% turn the valve on
				valve(a:z) = 1;

				% construct markers for graph
				y2(a) = 0;
				y2(a+1) = 5;

				% update previous valve off
				previous_valve_off = z;

				% estimate whiff duration
				estimated_whiff_duration = floor(nsb.sampling_rate*nsb.ff_time_whiff_off(desired_dilution));

				% pick a new blank duration 
				this_blank_duration = blank_durations(1);
				blank_durations(1) = [];

				% advance index 
				a = floor(this_blank_duration*nsb.sampling_rate + estimated_whiff_duration + a);

				if a > length(PID)
					goon = 0;
				end

				n_whiffs = n_whiffs + 1;

			end

			% update the graph
			nsb.handles.mfc_plot2.YData = y2;
			nsb.handles.mfc_plot1.YData = mfc1;

			nsb.handles.PID_plot1.XData = (1:length(PID))/nsb.sampling_rate;
			nsb.handles.PID_plot1.YData = PID;
			nsb.target_PID = PID;
            nsb.handles.ax(2).YLim = [0 1.1*(max(nsb.target_PID))];
            
			set(nsb.handles.fig,'Name',['# of whiffs = ' strlib.oval(n_whiffs)])

			% clear the data and the ControlParadigms, because we will build this anew when we run the optimization 

			nsb.ControlParadigm = [];
			nsb.mfc1_ansatz = mfc1;
			nsb.valve_ansatz = valve;
            
            nsb.n_successful_hits = zeros(n_whiffs,1);
           

		end % end buildNatStim


		function buildNatStimUI(nsb)
			nsb.handles = [];
			nsb.handles.fig = figure('outerposition',[100 100 700 700],'CloseRequestFcn',@nsb.cleanup_handles,'NumberTitle','off'); hold on

			% make a re-roll button
			nsb.handles.reroll_button = uicontrol(nsb.handles.fig,'Units','normalized','Position',[.01 .95 .1 .05],'String','Re-roll','Style','pushbutton','Callback',@nsb.reroll);

			% make some other controls
			nsb.handles.optim_step_picker = uicontrol(nsb.handles.fig,'Units','normalized','Position',[.2 .95 .2 .05],'String',{'Step 0'},'Style','popupmenu','Callback',@nsb.showStepData);
			nsb.updateStepPicker;

			nsb.handles.run_next_step = uicontrol(nsb.handles.fig,'Units','normalized','Position',[.6 .95 .2 .05],'String','Run next Step...','Style','pushbutton','Callback',@nsb.runNextStep);

			nsb.handles.ax(1) = subplot(3,1,1); hold on
			xlabel(nsb.handles.ax(1),'Time')
			ylabel(nsb.handles.ax(1),'MFC setpoint')
			nsb.handles.mfc_plot1 = plot(nsb.handles.ax(1),NaN,NaN,'k-');
			nsb.handles.mfc_plot2 = plot(nsb.handles.ax(1),NaN,NaN,'b-','LineWidth',3);
            nsb.handles.mfc_plot3 = plot(nsb.handles.ax(1),NaN,NaN,'r-');
			nsb.handles.ax(2) = subplot(3,1,2); hold on
			xlabel(nsb.handles.ax(2),'Time')
			ylabel(nsb.handles.ax(2),'PID')
			nsb.handles.PID_plot1 = plot(nsb.handles.ax(2),NaN,NaN,'k-');
			nsb.handles.PID_plot2 = plot(nsb.handles.ax(2),NaN,NaN,'r-');
			nsb.handles.ax(3) = subplot(3,3,7); hold on
			nsb.handles.ax(4) = subplot(3,3,8); hold on
			nsb.handles.ax(5) = subplot(3,3,9); hold on

			% link the top two
			linkaxes(nsb.handles.ax(1:2),'x')
			nsb.handles.ax(1).XLim = [0 nsb.T];
			nsb.handles.ax(1).YLim = [0 5];
			

			figlib.pretty('fs',12)

        end
        
        function [] = showStepData(nsb,src,~)

        	if src == nsb.handles.optim_step_picker
        		nsb.optimization_step = src.Value;
        	end

            this_step = nsb.optimization_step;
            nsb.handles.optim_step_picker.Value = this_step;
            if this_step > length(nsb.optimization_data)
            	return
            else
                PID = nsb.optimization_data(this_step).PID;
                this_mfc = nsb.optimization_data(this_step).mfc1_ansatz;
                this_flow = nsb.optimization_data(this_step).mfc1_flow;
            end
            
            nsb.handles.PID_plot2.XData = (1:length(PID))/nsb.sampling_rate;
			nsb.handles.PID_plot2.YData = PID;
            
            nsb.handles.mfc_plot3.XData = (1:length(PID))/nsb.sampling_rate;
            nsb.handles.mfc_plot3.YData = this_flow;
            
            nsb.handles.mfc_plot1.XData = (1:length(PID))/nsb.sampling_rate;
            nsb.handles.mfc_plot1.YData = this_mfc;
            
            % estimate how close each whiff is to target
            [desired_whiff_intensity, actual_whiff_intensity] = nsb.measureWhiffAmplitudes(PID);

			% show the results of that 
			cla(nsb.handles.ax(3))
			hold(nsb.handles.ax(3),'on')
			plot(nsb.handles.ax(3),desired_whiff_intensity,actual_whiff_intensity,'ko','MarkerSize',5)
			% only tweak the worst half of all whiffs
			xlabel(nsb.handles.ax(3),'Desired Whiff intensity (V)')
			ylabel(nsb.handles.ax(3),'Actual Whiff intensity (V)')
            try
            	r2 = corr(desired_whiff_intensity,actual_whiff_intensity)^2;
                title(nsb.handles.ax(3),['r^2 = ' strlib.oval(r2)])
            catch
            end
			set(nsb.handles.ax(3),'XScale','log','YScale','log')
			plot(nsb.handles.ax(3),[nsb.c0 nsb.c_max],[nsb.c0 nsb.c_max],'k--')
			axis(nsb.handles.ax(3),'square')

            
        end

        function [desired_whiff_intensity, actual_whiff_intensity] = measureWhiffAmplitudes(nsb,PID)
        	target_PID = nsb.target_PID;
			target_PID(target_PID>0) = 1;
			[ons,offs] = veclib.computeOnsOffs(target_PID);
			target_PID = nsb.target_PID;

			desired_whiff_intensity = NaN(length(ons),1);
			actual_whiff_intensity = NaN(length(ons),1);
			for i = 1:length(ons)
				desired_whiff_intensity(i) = target_PID(ons(i)+1);
				actual_whiff_intensity(i) = max(PID(ons(i):offs(i)));
			end
        end
        
        function [] = updateStepPicker(nsb,~,~)
            update_string = {''};
            for i = 1:nsb.optimization_step
                update_string{i} = ['Step ' strlib.oval(i)];
            end
            nsb.handles.optim_step_picker.String = update_string;
        end

		function [] = save(nsb)
			% shut down all figures
			nsb.cleanup_handles;

			save('natural_stimulus_builder.mat','nsb')

		end

		function [] = runNextStep(nsb,~,~)

			nsb.optimization_step = nsb.optimization_step + 1;

			% construct the control paradigm
			mfc1 = nsb.mfc1_ansatz;
			mfc1_flow = mfc1*nsb.mfc1_units;
			mfc2_flow = nsb.secondary_air_flow - mfc1_flow;
			mfc2 = mfc2_flow/nsb.mfc2_units;
            mfc2(mfc2<5/200) = 5/200; 
            mfc2(end) = 0;
            mfc1(end) = 0;

			ControlParadigm.Name = 'dummy';
			ControlParadigm.Outputs = zeros(4,floor(nsb.T*nsb.sampling_rate));
			ControlParadigm.Outputs(nsb.main_air_channel,:) = nsb.main_air_flow/nsb.mfc_main_air_units;
			ControlParadigm.Outputs(nsb.mfc1_channel,:) = mfc1;
			ControlParadigm.Outputs(nsb.mfc2_channel,:) = mfc2;
			ControlParadigm.Outputs(nsb.valve_channel,:) = nsb.valve_ansatz;

			% run it
			if ispc
				optim_data = kontroller('ControlParadigm',ControlParadigm,'RunTheseParadigms',1,'w',nsb.sampling_rate);
			else
				% this bullshit hack is only for development
				clear optim_data
				optim_data.MFC = nsb.optimization_data(nsb.optimization_step).mfc1_flow;
				optim_data.PID = nsb.optimization_data(nsb.optimization_step).PID;

			end

			% extract the data, stash it in the nsb object
			nsb.optimization_data(nsb.optimization_step).mfc1_flow = optim_data.MFC;
			nsb.optimization_data(nsb.optimization_step).PID = optim_data.PID;
			nsb.optimization_data(nsb.optimization_step).mfc1_ansatz = mfc1;

			% tweak the ansatz in an attempt to improve the stimulus
			nsb.tweakAnsatz;
            
            % update the step picker UI
            nsb.updateStepPicker;
            
			% update the graph
			nsb.showStepData([]);

		end

		function tweakAnsatz(nsb)
			% estimate how close each whiff is to target
			target_PID = nsb.target_PID;
			target_PID(target_PID > 0) = 1;
			[ons,offs] = computeOnsOffs(target_PID);
			target_PID = nsb.target_PID;

           
			% make a matrix of all PIDs and all mfc1_ansatzes
			all_PID = vertcat(nsb.optimization_data(1:nsb.optimization_step).PID);
			PID = nsb.optimization_data(nsb.optimization_step).PID;
			all_mfc1_ansatz = vertcat(nsb.optimization_data(1:nsb.optimization_step).mfc1_ansatz);
			new_mfc1_ansatz = nsb.optimization_data(nsb.optimization_step).mfc1_ansatz;
            new_mfc1_ansatz(new_mfc1_ansatz > 5) = 5;
            new_mfc1_ansatz(new_mfc1_ansatz < 0) = 0;
            [valve_ons,valve_offs] = computeOnsOffs(nsb.valve_ansatz);
            
            [desired_whiff_intensity, actual_whiff_intensity] = measureWhiffAmplitudes(nsb)

			% only tweak whiffs with an error of more than 10%
			nsb.tweak_these = max([desired_whiff_intensity./actual_whiff_intensity actual_whiff_intensity./desired_whiff_intensity]') > 1.1;
            nsb.n_successful_hits(~nsb.tweak_these) = nsb.n_successful_hits(~nsb.tweak_these) + 1;
            
            % if something has more than 2 successful hits, freeze it
            nsb.tweak_these(nsb.n_successful_hits > 1) = 0;
            
            
            
            
			% f = figure; hold on
			% clear ax
			% for i = 1:sum(nsb.tweak_these)
			% 	ax(i) = autoPlot(sum(nsb.tweak_these),i); hold on
			% 	set(ax(i),'YLim',[0 5])
			% end

			c = 0;

            
            
			for i = 1:length(ons)
				if ~nsb.tweak_these(i)
					continue
				end
				c = c + 1;
				desired_whiff_intensity = target_PID(ons(i)+1);
				actual_whiff_intensity = max(all_PID(:,ons(i):offs(i))');
                if i == 1
                    applied_mfc_setpoints = max(all_mfc1_ansatz(:,1:valve_ons(i)-1)');
                else
                    applied_mfc_setpoints = max(all_mfc1_ansatz(:,valve_offs(i-1)+1:valve_ons(i)-1)');
                end
                % disp(['tweaking whiff #' strlib.oval(i)])
                % disp(['Error fold change:' strlib.oval(desired_whiff_intensity/actual_whiff_intensity)])
                
                % show all the data so far
                % plot(ax(c),actual_whiff_intensity,applied_mfc_setpoints,'k+')
                % plot(ax(c),[desired_whiff_intensity desired_whiff_intensity],[0 5],'r--')

                % title(ax(c),['t = ' strlib.oval(ons(i)/nsb.sampling_rate,3)])

                if length(actual_whiff_intensity) > 1
                	ff = fit(actual_whiff_intensity(:),applied_mfc_setpoints(:),'poly1');
                	% show the curve fit
                    x = linspace(nsb.c0,nsb.c_max,100);
                    %plot(ax(c),x,ff(x),'r')
                    new_setpoint = ff(desired_whiff_intensity);  
                else
                	if min(actual_whiff_intensity) > desired_whiff_intensity
						% actual whiffs always greater than desired whiff,
						% halve the control setpoint going in here
                        fudge_factor = .7;
						
					elseif max(actual_whiff_intensity) < desired_whiff_intensity
						% actual whiffs always smaller than desired whiffs,
						% increase the control setpoint going in here 
                        fudge_factor = 1.5;
                    end
                    if i == 1
                        new_setpoint = mean(new_mfc1_ansatz(1:valve_offs(i)-1)*fudge_factor);
                    else
                        new_setpoint = mean(new_mfc1_ansatz(valve_offs(i-1):valve_offs(i)-1)*fudge_factor);
                    end
                    
                end
                if new_setpoint < 0
                    new_setpoint = 0;
                end
                if new_setpoint > 5
                    new_setpoint = 5;
                end
                %plot(ax(c),[nsb.c0 nsb.c_max],[new_setpoint new_setpoint],'b--')

                if ~isnan(new_setpoint)
                    if i == 1
                        new_mfc1_ansatz(1:valve_offs(i)-1) = new_setpoint;
                    else
                        new_mfc1_ansatz(valve_offs(i-1):valve_offs(i)-1) = new_setpoint;
                    end
                    
                end
			end


			% update the ansatz 
			assert(max(new_mfc1_ansatz) <= 5, 'New MFC ansatz has values > 5')
			assert(min(new_mfc1_ansatz) >= 0, 'New MFC ansatz has values < 0')
            nsb.mfc1_ansatz = new_mfc1_ansatz;
		end

		function cleanup_handles(nsb,~,~)
			if isempty(nsb.handles)
				return
			end
			fn = fieldnames(nsb.handles);
			for i = 1:length(fn)
				try
					delete(nsb.handles.(fn{i}));
				catch
				end
				try
					rmfield(nsb.handles,fn{i});
				catch
				end
			end

			try
				nsb.handles = [];
			end
			try
				rmfield(nsb,'handles');
			end

		end

		function [S,pS, t_blanks, pt_blanks] = genPDF(nsb)
			S  = linspace(nsb.c0,nsb.c_max,1e4);
			pS = (1/(log(nsb.c_max/nsb.c0))).*(1./S);

			t_blanks = linspace(nsb.min_blank_duration,nsb.T,1e4);
			pt_blanks = (sqrt(nsb.min_blank_duration)/2).*(t_blanks.^-1.5);
		end

		function rand_samples = sampleFromDist(nsb,x,px)
			% sample from these PDFs 
			RandStream.setGlobalStream(RandStream('mt19937ar','Seed',nsb.random_seed)); 
			rand_samples = pdfrnd(x,px,1e3);
		end

		function reroll(nsb,~,~)
			nsb.random_seed = randi(1e7);
			nsb.buildNatStim;

		end

	end % end methods 
	

end % end classdef 