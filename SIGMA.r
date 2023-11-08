##########################################################################
# Description: 	This is a complete version of the SIGMA model
#       	with the application to the Jinghe River Watershed in China
# Author: 	Dr. Xander Wang, P.Eng., CCCCA at UPEI
# Date:   	February 6, 2023
# Location: 	Charlottetown, PE, Canada
# Website:	http://climatesmartlab.ca
# Contact:	xiuquan.wang@gmail.com or xxwang@upei.ca
##########################################################################

#: create a new environment to host all parameters
FM.env <- new.env()

###########################
# Step 1: Data preparation and processing
# -------------------------
# Note that for our case study: 
#    We have hourly peak flow data at the watershed outlet
#    but we only have daily precipitation data for the weather stations
#    we have to convert daily precipitation to hourly/sub-hourly precipitation
# -------------------------

####################################
# Step 2: Prepare the model settings
# ----------------------------------
# 2.1: 4 zones including Huanxian + Xifengzhen + Pingliang + Changwu
# Their sequences in the vector: c(1, 2, 3, 4)

# Total number of zones = > 4
p_nZones = 4

# Scale factor fa ==> to be calibrated
p_fa = c(0.001, 0.001, 0.001, 0.001)

# Total area for each zone ==> to be measured in GIS, units: km^2
p_totalArea = c(10588.75, 8445.09, 14204.97, 9961.89) #: Sz
p_totalArea = p_totalArea * 1000000 #: convert from km^2 to m^2

# Two distances for each zone ==> to be measured in GIS
p_a = c(162.53, 157.58, 151.95, 190.28) #: longest distance within the zone, units: km
p_b = c(64.52, 39.34, 104.71, 53.56) #: the distance perpendicular to the longest one, units: km
p_fs = p_a / p_b #: calculate the shape factor, units: km/km = m/m

# Shape info for each water tank
p_waterTank_Length = sqrt(p_totalArea * p_fa * p_fs) #: L
p_waterTank_Width = sqrt((p_totalArea * p_fa) / p_fs)  #: W  
p_waterTank_Area = p_waterTank_Length * p_waterTank_Width #: Sw

# Info to be measured in GIS
p_D = c(265.98, 183.82, 218.47, 96.05) * 1000  #: Distance from weather station along river channel to the watershed outlet, units: m
p_A = c(216.01, 216.01, 216.01, 216.01)  #: Cross-sectional area of the river channel in each zone, units: m^2
p_Pw = c(117.46, 117.46, 117.46, 117.46)  #: Wetted perimeter of the river channel in each zone, units: m
p_S = c(0.038, 0.066, 0.052, 0.116)  #: Average slope for the river channel between the weather station and the watershed outlet (units: m/m)
p_R = p_A / p_Pw #: calculate the average hydraulic radius
p_n = c(0.04, 0.04, 0.04, 0.04)  #: Manning coefficient of roughness, range -> [0.02, 0.06]

# Average water velocity
p_v = (p_R ^ (2/3) * p_S ^ (1/2)) / p_n

# Travel time of the water tank to the watershed outlet
p_t = p_D / p_v

# Duration of the single-pulse signal
p_deltaT = p_waterTank_Length / p_v

# Apply uniform or triangular pattern to daily precipitation
p_dailyP_pattern = "uniform" 	#: you can switch this between "uniform" and "triangular"
# ----------------------------------


#########################################################
## Function to generate small step rainfall time series 
## Assume that there is only one peak
## Apply simply linear interpolation for both sides
f_generate_rainfall_linear = function(total_rainfall, total_steps)
{
	#: to store the rainfall values
	o_output_rainfall = rep(0, total_steps)
	
	#: right in the middle
	n_peak_step = round(total_steps / 2, digits = 0)
	n_peak_rainfall = (total_rainfall / total_steps) * 2
	o_output_rainfall[n_peak_step] = n_peak_rainfall
	
	#######################
	#: left side
	for (iLeft in 1:(n_peak_step-1))
	{
		#: use simply liner interpolation: y = ax + b
		n_left_slope = (n_peak_rainfall - 0) / (n_peak_step - 1)
		n_left_intercept = (-1) * n_left_slope
		o_output_rainfall[iLeft] = n_left_slope * iLeft + n_left_intercept
	}
	
	############################
	#: right side
	for (iRight in (n_peak_step+1):total_steps)
	{
		#: use simply liner interpolation: y = ax + b
		n_right_slope = (0 - n_peak_rainfall) / (total_steps - n_peak_step)
		n_right_intercept = (-1) * n_right_slope * total_steps
		o_output_rainfall[iRight] = n_right_slope * iRight + n_right_intercept 
		
	}

	return(o_output_rainfall)
}
###########################################


o_years = c("1970", "1973", "1975", "1977", "1992", "1994", "1995", "1996", "1998")
s_labels = c("a", "b", "c", "d", "e", "f", "g", "h", "i")
s_folder = "events_data/"
png(filename = paste("SIGMA_9years_", p_dailyP_pattern, "_comparison.png", sep = ""), units = "in", width = 20, height = 13, pointsize = 24, res = 300)
m <- matrix(c(1:9), nrow = 3, ncol = 3, byrow = TRUE)
layout(mat = m, heights = c(1/3, 1/3, 1/3), widths = c(1/3, 1/3, 1/3))

#: to store the simulated peak flow for each year: datetime, seconds, flowrate
o_sim_peakflows = matrix(0, length(o_years), 3)
	
for (iyear in 1:length(o_years))
{
	#####################################
	# Step 3: Event specific data
	# ----------------------------------
	# Start and end datetime -> to be consistent with the precipitation records
	e_startDatetime = ""
	e_endDatetime = ""

	# Precipitation data, units: mm
	e_precipDataFile = paste(s_folder, o_years[iyear], "_precip_daily.csv", sep = "") 	# format: Datetime,Station1,Station2,Station3,Station4
	e_precipData = read.csv(file = e_precipDataFile, header = T)
	e_eventInitialDateTime = paste(e_precipData[1, 1], " 00:00:00", sep = "")	# initial datatime for this event
	e_precipDataInterval = 24 * 60 * 60											# record internal represented in seconds

	# Extract information about this event
	e_totalPEventDuration = nrow(e_precipData) * e_precipDataInterval  # total precipitation event duration in seconds

	# ----------------------------------

	#####################################
	# Step 4: Flow data output matrix, units: m3/s
	# ----------------------------------
	# Model time step in seconds
	p_timeStep = 60		# 60 seconds by default ==> output step for flow rate
	p_currentRunTime = 0  # in seconds from the beginning
	p_totalRunTime = round(e_totalPEventDuration * 1.5, digits = 0) # increase the duration by 50% to account for the delay in outflow
	p_totalSteps = round(p_totalRunTime / p_timeStep, digits = 0)

	# Columes for the matrix
	# NOTE: no baseflow will be generated in this model because: 
	#       the baseflow (estimated as the initial flow right before the precipitation event) has been subtracted from the observed outflow
	o_colNames = c("Datetime", "SecondsFromBeginning", "HX_flowsignal", "XF_flowsignal", "PL_flowsignal", "CW_flowsignal", "Totalflow") 
	o_flowRate = matrix(0, p_totalSteps, length(o_colNames))
	colnames(o_flowRate) = o_colNames

	# To store the reproduced precipitation data for small steps
	o_smallsteps_pdata = rep(0, p_totalSteps)
	o_smallsteps_pdata_uniform = rep(0, p_totalSteps)

	# Assign the output datetime for the matrix
	o_initialDatetime = as.POSIXct(e_eventInitialDateTime, format = "%Y-%m-%d %H:%M:%S")
	for (iStep in 1:p_totalSteps)
	{
		p_currentRunTime = p_timeStep * iStep
		o_currentDatetime = o_initialDatetime + p_currentRunTime
		o_flowRate[iStep, 1] = format(o_currentDatetime, "%Y-%m-%d %H:%M:%S")
		o_flowRate[iStep, 2] = p_currentRunTime
	}

	# ----------------------------------


	#####################################
	# Step 5: Generate signals for each zone and each time step
	# ----------------------------------
	for (iZone in 1:p_nZones)
	{
		# process each zone separately and store the signal into the corresponding column of the output matrix
		for (iPrecip in 1:nrow(e_precipData))
		{
			# for each line of precipitation record, divide it into small duration which is the same as the model time step
			if (e_precipDataInterval == p_timeStep)
			{
				#: same interval, no need to divide
				#: generate the single-pulse signal (travel time = p_t[iZone], duration = p_deltaT[iZone], flow rate = Q)
				n_precip = e_precipData[iPrecip, 1 + iZone] # as the 1st colume is Datetime
				if (n_precip > 0)
				{
					#: if it didn't rain, no signal to generate
					n_flowrate = n_precip * p_waterTank_Width[iZone] * p_v[iZone]
					n_traveltime = p_t[iZone]
					n_duration = p_deltaT[iZone]
					
					#: convert the travel time + precipitation datatime to model step
					n_currentRunTime = iPrecip * e_precipDataInterval + n_traveltime
					n_currentStep = round(n_currentRunTime / p_timeStep, digits = 0)
					n_signalTotalSteps = round(n_duration / p_timeStep, digits = 0)
					
					#: store signal (plus 2 --> as the first two columns stores datetime and totalseconds)
					for (mStep in 1:n_signalTotalSteps)
					{
						o_flowRate[n_currentStep + mStep - 1, 2 + iZone] = n_flowrate
					}
				}
			}
			else
			{
				#: precipitation internal is great than model time step, need to divide
				n_totalInternalSteps = round(e_precipDataInterval / p_timeStep, digits = 0)
				n_precip = e_precipData[iPrecip, 1 + iZone] # as the 1st colume is Datetime
				if (n_precip > 0)
				{
					#: if it didn't rain, no signal to generate
					n_internalprecip = n_precip / n_totalInternalSteps # uniform distribution for precipitation during the entire interval
					
					#################
					#: convert daily rainfall to hourly rainfall, central peak patern + left/right peak pattern + random pattern
					o_smallstepprecips = f_generate_rainfall_linear(n_precip, n_totalInternalSteps)
					n_smallPcurrentRunTime = iPrecip * e_precipDataInterval
					n_smallPcurrentStep = round(n_smallPcurrentRunTime / p_timeStep, digits = 0)
					o_smallsteps_pdata[n_smallPcurrentStep:(n_smallPcurrentStep+n_totalInternalSteps-1)] = o_smallstepprecips
					o_smallsteps_pdata_uniform[n_smallPcurrentStep:(n_smallPcurrentStep+n_totalInternalSteps-1)] = n_internalprecip
					#################
					
					for (iSmall in 1:n_totalInternalSteps)
					{
						#: handle each small event
						n_smallflowrate = 0
						if (p_dailyP_pattern == "uniform")
							n_smallflowrate = n_internalprecip * p_waterTank_Width[iZone] * p_v[iZone]  #==> for uniform pattern
						if (p_dailyP_pattern == "triangular")
							n_smallflowrate = o_smallstepprecips[iSmall] * p_waterTank_Width[iZone] * p_v[iZone]  #==> for triangular pattern
							
						n_smalltraveltime = p_t[iZone]
						n_smallduration = p_deltaT[iZone]
						
						#: convert the travel time + precipitation datatime to model step
						n_small_to_currentRunTime = iPrecip * e_precipDataInterval + p_timeStep * iSmall + n_smalltraveltime
						n_small_to_currentStep = round(n_small_to_currentRunTime / p_timeStep, digits = 0)
						n_small_to_signalTotalSteps = round(n_smallduration / p_timeStep, digits = 0)
						
						#: store signal (plus 2 --> as the first two columns stores datetime and totalseconds)
						for (smStep in 1:n_small_to_signalTotalSteps)
						{
							o_flowRate[n_small_to_currentStep + smStep - 1, 2 + iZone] = n_smallflowrate
						}
					}
				}
			}
		}
	}

	# add all signals together
	o_flowRate[ , 7] = as.numeric(o_flowRate[ , 3]) + as.numeric(o_flowRate[ , 4]) + as.numeric(o_flowRate[ , 5]) + as.numeric(o_flowRate[ , 6])

	# ----------------------------------

	# output
	write.csv(o_flowRate, file = paste(s_folder, o_years[iyear], "_flowrate_sigma.csv", sep = ""), row.names = FALSE, quote = FALSE)

	# comparison
	o_obs_flowrate = read.csv(file = paste(s_folder, o_years[iyear], "_flowrate.csv", sep = ""), header = T)
	
	#: flow rate plot
	n_max_flowrate = max(as.numeric(o_obs_flowrate[ , 2]), as.numeric(o_flowRate[ , 7]))
	par(mar=c(3.5, 3.5, 2, 1.2), mgp=c(1.8, 0.5, 0))
	plot(x = as.numeric(o_flowRate[ , 2]), y = as.numeric(o_flowRate[ , 7]) + as.numeric(o_obs_flowrate[1, 2]), 
		main = paste("(", s_labels[iyear], ") Peak flow in ", o_years[iyear], sep = ""), 
		ylab = expression(paste("Q (", m^3, "/s)", sep = "")), xlab = paste("Time in seconds (since ", e_eventInitialDateTime, ")", sep = ""), 
		type = "o", col = "#0072BD", ylim = c(0, n_max_flowrate * 1.2))
	lines(x = as.numeric(o_obs_flowrate[ , 1]), y = as.numeric(o_obs_flowrate[ , 2]), type = "o", col = "#D95319")
	if (iyear == 1)
	{
		legend("topright", legend = c("Simulated", "Observed"), col = c("#0072BD", "#D95319"), lty = c(1, 1), pch = 1, cex = 0.95)
	}
	
	#: store the peak flow info
	n_maxindex = which.max(as.numeric(o_flowRate[ , 7]))
	o_sim_peakflows[iyear, 1] = o_flowRate[n_maxindex, 1]
	o_sim_peakflows[iyear, 2] = o_flowRate[n_maxindex, 2]
	o_sim_peakflows[iyear, 3] = o_flowRate[n_maxindex, 7]
}

write.csv(o_sim_peakflows, file = paste(s_folder, "peakflow_sim_", p_dailyP_pattern, ".csv", sep = ""), row.names = FALSE, quote = FALSE)
dev.off()
