"""
GenX: An Configurable Capacity Expansion Model
Copyright (C) 2021,  Massachusetts Institute of Technology
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
A complete copy of the GNU General Public License v2 (GPLv2) is available
in LICENSE.txt.  Users uncompressing this from an archive may not have
received this license file.  If not, see <http://www.gnu.org/licenses/>.
"""



function write_core_behaviors(EP::Model, inputs::Dict, symbol::Symbol, SET::Vector{Int}, filename::AbstractString)
	dfTS = inputs["dfTS"]
	T = inputs["T"]

	resources = by_rid_df(SET, :Resource, dfTS)
	zones = by_rid_df(SET, :Zone, dfTS)

	df = DataFrame(Resource = resources, Zone = zones)
	event = value.(EP[symbol][SET,:]).data
	df.Sum = vec(sum(event, dims=2))

	df = hcat(df, DataFrame(event, :auto))
	auxNew_Names=[:Resource; :Zone; :Sum; [Symbol("t$t") for t in 1:T]]
	rename!(df,auxNew_Names)
	total = DataFrame(["Total" 0 sum(df[!,:Sum]) zeros(1,T)], :auto)
	total[:, 4:T+3] .= sum(event, dims=1)
	rename!(total,auxNew_Names)
	df = vcat(df, total)

	CSV.write(filename, dftranspose(df, false), writeheader=false)
	return df
end

function write_scaled_values(EP::Model, inputs::Dict, symbol::Symbol, SET::Vector{Int}, filename::AbstractString, msf)
	dfTS = inputs["dfTS"]
	T = inputs["T"]

	resources = by_rid_df(SET, :Resource, dfTS)
	zones = by_rid_df(SET, :Zone, dfTS)

	df = DataFrame(Resource = resources, Zone=zones)
	quantity = value.(EP[symbol][SET, :]).data * msf
	df.AnnualSum = quantity * inputs["omega"]

	df = hcat(df, DataFrame(quantity, :auto))
	auxNew_Names=[:Resource; :Zone; :AnnualSum; [Symbol("t$t") for t in 1:T]]
	rename!(df,auxNew_Names)
	total = DataFrame(["Total" 0 sum(df.AnnualSum) zeros(1,T)], :auto)
	total[:, 4:T+3] .= sum(quantity, dims=1)
	rename!(total,auxNew_Names)
	df = vcat(df, total)
	CSV.write(filename, dftranspose(df, false), writeheader=false)

	return df
end

function write_thermal_storage_system_max_dual(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
	dfTS = inputs["dfTS"]

	FIRST_ROW = 1
	if dfTS[FIRST_ROW, :System_Max_Cap_MWe_net] >= 0
		val = -1*dual.(EP[:cCSystemTot])
		val *= setup["ParameterScale"] == 1 ? ModelScalingFactor : 1.0
		df = DataFrame(:System_Max_Cap_MW_th_dual => val)
		filename = joinpath(path, "System_Max_TS_Cap_dual.csv")
		CSV.write(filename, dftranspose(df, false), writeheader=false)
	end
end


@doc raw"""
	write_capacity(path::AbstractString, inputs::Dict, setup::Dict, EP::Model))

Function for writing the diferent capacities for the different generation technologies (starting capacities or, existing capacities, retired capacities, and new-built capacities).
"""
function write_thermal_storage(path::AbstractString, inputs::Dict, setup::Dict, EP::Model)
	# Capacity decisions
	dfGen = inputs["dfGen"]
	dfTS = inputs["dfTS"]
	T = inputs["T"]


	TSResources = dfTS[!,:Resource]
	TSG = length(TSResources)
	corecappower = zeros(TSG)
	for i in 1:TSG
		corecappower[i] = first(value.(EP[:vCCAP][dfTS[i,:R_ID]]))
	end

	corecapenergy = zeros(TSG)
	for i in 1:TSG
		corecapenergy[i] = first(value.(EP[:vTSCAP][dfTS[i,:R_ID]]))
	end

	rhcapacity = zeros(TSG)
	for i in 1:TSG
		if dfTS[i, :RH] == 1
			rhcapacity[i] = first(value.(EP[:vRHCAP][dfTS[i,:R_ID]]))
		else
			rhcapacity[i] = 0
		end
	end

	dfCoreCap = DataFrame(
		Resource = TSResources, Zone = dfTS[!,:Zone],
		CorePowerCap = corecappower[:],
		TSEnergyCap = corecapenergy[:],
		RHPowerCap = rhcapacity[:]
	)

	# set a single scalar to avoid future branching
	msf = setup["ParameterScale"] == 1 ? ModelScalingFactor : 1
	dfCoreCap.CorePowerCap = dfCoreCap.CorePowerCap * msf
	dfCoreCap.TSEnergyCap = dfCoreCap.TSEnergyCap * msf
	dfCoreCap.RHPowerCap = dfCoreCap.RHPowerCap * msf
	CSV.write(joinpath(path,"TS_capacity.csv"), dfCoreCap)

	THERMAL_STORAGE = dfTS.R_ID
	RH = dfTS[dfTS.RH .==1, :R_ID]
	FUS = dfTS[dfTS.FUS .== 1, :R_ID]
	NONFUS = dfTS[dfTS.FUS .== 0, :R_ID]
	### CORE POWER TIME SERIES ###
	dfCorePwr = write_scaled_values(EP, inputs, :vCP, THERMAL_STORAGE, joinpath(path, "TS_CorePwr.csv"), msf)

	### THERMAL SOC TIME SERIES ###
	dfTSOC = write_scaled_values(EP, inputs, :vTS, THERMAL_STORAGE, joinpath(path, "TS_SOC.csv"), msf)

	### RESISTIVE HEATING TIME SERIES ### 
	if !isempty(RH)
		dfRH = write_scaled_values(EP, inputs, :vRH, RH, joinpath(path, "TS_RH.csv"), msf)
	end

	### FUSION SPECIFIC OUTPUTS ###
	if !isempty(FUS)
		### RECIRCULATING POWER TIME SERIES ###
		dfRecirc = write_scaled_values(EP, inputs, :eTotalRecircFus, FUS, joinpath(path, "TS_Recirc.csv"), msf)

		### CORE STARTS, SHUTS, COMMITS, and MAINTENANCE TIMESERIES ###
		dfFStart = write_core_behaviors(EP, inputs, :vFSTART, FUS, joinpath(path, "TS_FUS_start.csv"))
		dfFShut = write_core_behaviors(EP, inputs, :vFSHUT, FUS, joinpath(path, "TS_FUS_shut.csv"))
		dfFCommit = write_core_behaviors(EP, inputs, :vFCOMMIT, FUS, joinpath(path, "TS_FUS_commit.csv"))

		if setup["OperationWrapping"] == 0 && !isempty(get_maintenance(inputs))
			dfMaint = write_core_behaviors(EP, inputs, :vFMDOWN, FUS,  joinpath(path, "TS_FUS_maint.csv"))
			dfMShut = write_core_behaviors(EP, inputs, :vFMSHUT, FUS, joinpath(path, "TS_FUS_maintshut.csv"))
		end
	end

	### NON FUS CORE STARTS, SHUTS, COMMITS ###
	if (!isempty(NONFUS) && setup["UCommit"] > 0)
		dfNStart = write_core_behaviors(EP, inputs, :vCSTART, NONFUS, joinpath(path, "TS_NONFUS_start.csv"))
		dfNShut = write_core_behaviors(EP, inputs, :vCSHUT, NONFUS, joinpath(path, "TS_NONFUS_shut.csv"))
		dfNCommit = write_core_behaviors(EP, inputs, :vCCOMMIT, NONFUS, joinpath(path, "TS_NONFUS_commit.csv"))
	end

	# Write dual values of certain constraints
	write_thermal_storage_system_max_dual(path, inputs, setup, EP)

end
