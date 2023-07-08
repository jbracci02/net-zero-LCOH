# H2 Techno-economic Analysis

# Justin Bracci, Evan Sherwin, Adam Brandt
# Department of Energy Resources Engineering
# Stanford University
# Updated 7/2023

# In this file we generate and compare levelized electricity-based and fossil-based H2 
# production costs for a 250 metric ton/day facility placed in Sacramento, CA.
# We also calculate life-cycle emissions for the different production pathways and can
# choose whether there is a cost for emissions

# For the electricity-based pathways, we use the Gurobi linear optimizer 
# to minimize the cost of H2 production by correctly sizing production components.
# These components include solar PV, battery and hydrogen storage, an electrolyzer, 
# and a grid connection. We look at different levels of reliability of electricity-based
# H2 production as well as consider life-cycle emissions for the project.

# Second, we generate fossil-based levelized H2 production costs for a 250 ton/day
# facility. We examine LCOH for both grey and blue hydrogen and include full and
# partial CO2 capture schemes.

# Finally, we compare the LCOH of all net-zero, hourly reliable hydrogen production pathways.

# In this version of the model, we conduct 50 Monte Carlo runs and use a 96% Confidence Interval
# for construction of the error bars in each pathway LCOH.


####### ADMIN AND PACKAGES ######

# Importing necessary packages. Uncomment all packages the first time running code.
import Pkg;
# Pkg.add("JuMP");

# Add solver for nonlinear optimmization
# Pkg.add("Ipopt")

# Add MATLAB
# Pkg.add("MATLAB")

# Pkg.add("Distributions")

# Add Statistics
# Pkg.add("Statistics")

# Add StatsPlot
Pkg.add("StatsPlots")

# Use packages as needed
using JuMP
using Gurobi
using MATLAB
using Plots
using Statistics
using StatsPlots
using Distributions

##### Start of Large While Loop to Generate Error Bounds

# Initialize matrix that will store the LCOH for each pathway in each iteration
LCOHallEB = []

# Set Number of Model Runs [this will increase by 1 with each iteration]
m=1

# Keep running the model until 50 model runs are performed
while m <= 50

    # Loading capacity factor data from MATLAB
    file = MatFile("CF_Data_SAM.mat")
    CapFactorHR = get_variable(file, "CA_CF_2012") # Solar Ratio is fraction of Output each hour compared to capacity
    close(file)

    # Loading grid data from MATLAB
    file = MatFile("CA_grid_data_cambium.mat")
    CIgrid1 = get_variable(file, "emissions_MC")/1000   # CA Grid Carbon Intensity 2035 (kg CO2/kWh)
    CIgrid_adjusted1 = get_variable(file, "emissions_adjusted_MC")/1000   # CA Grid Carbon Intensity 2035 with embodied (kg CO2/kWh)
    LMPgrid1 = get_variable(file, "total_cost_enduse_MC")/1000   # CA Grid Hourly Marginal Price ($/kWh)
    CIgrid4 = get_variable(file, "emissions_MC_LRE")/1000   # CA Grid Carbon Intensity 2035 (kg CO2/kWh)
    CIgrid_adjusted4 = get_variable(file, "emissions_adjusted_MC_LRE")/1000   # CA Grid Carbon Intensity 2035 with embodied (kg CO2/kWh)
    LMPgrid4 = get_variable(file, "total_cost_enduse_MC_LRE")/1000   # CA Grid Hourly Marginal Price ($/kWh)
    CIgrid5 = get_variable(file, "emissions_MC_HRE")/1000   # CA Grid Carbon Intensity 2035 (kg CO2/kWh)
    CIgrid_adjusted5 = get_variable(file, "emissions_adjusted_MC_HRE")/1000   # CA Grid Carbon Intensity 2035 with embodied (kg CO2/kWh)
    LMPgrid5 = get_variable(file, "total_cost_enduse_MC_HRE")/1000   # CA Grid Hourly Marginal Price ($/kWh)
    close(file)

    ###### GLOBAL CONVERSIONS AND DATA #######################
    HRSperDAY = 24

    KWperMW = 1000

    KgPerTon = 1000

    DAYSperYR = 365

    MONTHSperYR = 12

    DAYSperMONTH = [31,28,31,30,31,30,31,31,30,31,30,31]
    DAYSperMONTHcum = [0,31,59,90,120,151,181,212,243,273,304,334,365]

    MJperKWH = 3.6

    MJperTHERM = 105.5


    # Reliability Type
    types = ["yearly" "monthly" "daily" "hourly" "hourly_GS" "hourly_GS2"]

    # Hours in a Year
    nHours = length(CapFactorHR)

    # Hydrogen Supply Required (ton/day)
    H2SupplyDaily = 250

    # H2 Energy Content (kWh/kg)
    H2Energy = 33.3


    # PART I: GREEN H2 LCOH ANALYSIS

    ###### PARAMETERS and DATA ###############

    # Carbon Intensity and electricity pricing of the Grid
    c = rand(1:3)
    if c == 1
        LMPgrid = LMPgrid1
        CIgrid = CIgrid1
        CIgrid_adjusted = CIgrid_adjusted1
    elseif c == 2
        LMPgrid = LMPgrid4
        CIgrid = CIgrid4
        CIgrid_adjusted = CIgrid_adjusted4
    else
        LMPgrid = LMPgrid5
        CIgrid = CIgrid5
        CIgrid_adjusted = CIgrid_adjusted5   
    end

    # Electrolyzer Efficiency
    ElectEff = Base.rand(TriangularDist(0.6,0.7,0.65))

    # Electrolyzer Energy Input (kWh/kg)
    ElectInput = H2Energy ./ ElectEff

    # Electrolyzer Facility Yearly Energy Demand (kWh/yr)
    YearlyEnergyDemand = H2SupplyDaily * KgPerTon .*(ElectInput) .* DAYSperYR 

    # Capital Cost for Electrolyzers ($/kWdc)
    CapElect = Base.rand(TriangularDist(340,915,340))

    # Capital Cost for Solar Farm ($/kWdc)
    CapSolar = Base.rand(TriangularDist(400,950,400))

    # Project Life (yrs)
    YrTot = Base.rand(TriangularDist(20,30,25))

    # Lifetime of Electrolyzers (yrs)
    YrElect = Base.rand(TriangularDist(20,30,25))

    # Capital Recovery Factor
    WACC = Base.rand(TriangularDist(0.05,0.10,0.08))
    CapRF1 = WACC/(1-(1+WACC)^(-YrTot))

    # Fixed O&M Solar PV (fraction/yr)
    SolarOM = 0.02

    # Fixed O&M Electrolyzer (fraction/yr)
    ElectOM = 0.075

    # Energy Requirements to Compress Hydrogen (kWh/kg)
    CompEnergy = 1.2

    # Hydrogen Compressors Capital Cost ($/kW)
    # H2CompCapex = 2400

    # Fixed O&M H2 Storage (%/yr)
    H2StorageOM = 0.01

    # Battery Storage Capex ($/kWh)
    BattCapex = Base.rand(TriangularDist(100,350,250))

    # Battery Storage Efficiency
    BattEff = Base.rand(TriangularDist(0.8,0.9,0.85))

    # H2 Tank Storage Efficiency
    TankEff = 1

    # Fixed O&M Battery Storage (%/yr)
    BattStorageOM = 0.025

    # Hydrogen Storage Capital Cost ($/kg)
    H2StoreCapex = Base.rand(TriangularDist(200,830,500))

    # Land Area Requirement for Solar Farm (acres/MW installed)
    SolarArea = 7.5

    # Land Area Cost for Solar Farm ($/acre/year)
    SolarLandCost = Base.rand(TriangularDist(500,1250,750))

    # Water Cost ($/tonne)
    WaterCost = 1

    # Water Usage (kg H2O / kg H2)
    WaterUse = 15

    # Grid Interconnection Fee ($/kW electricity)
    GridConnectFee = Base.rand(TriangularDist(45,340,180))

    # Fixed O&M Grid Connection (fraction/yr)
    GridOM = 0.01

    # DAC Cost ($/kg CO2)
    DACCost = Base.rand(TriangularDist(0.1,0.6,0.1))

    # Life Cycle Emissions Solar PV (kg CO2/kWh electricity)
    LCESolarPV = Base.rand(TriangularDist(0,0.05,0.04))

    # Ramping Cost ($ per change in KWe- [1] and kg H2 [2] each hour)
    RampCost1 = 0.005
    RampCost2 = 0.05

    # SMUD Grid Infrastructure Fixed Charge  ($/month)
    GridFixedCharge = 203.15

    # SMUD Site Infrastructure Charge ($/max KW over year/month)
    # [not used in final model]
    SiteCharge = 0

    # Demand Charge ($/max kW each month/month)
    DemandCharge = Base.rand(TriangularDist(0,30,10))

    # Net Surplus Compensation Rate ($/kWh)
    NSCR = LMPgrid.*0.99

    # Turn CO2 Mitigation On (CO2Lever = 1 -> ON | 0 -> OFF)
    CO2Lever = 1

    # Fraction Electricity Allowed from the Grid in PV/Storage/Grid** pathway
    GridAllow = 0.1

    # Initialize LCOH breakdown vectors for Plotting
    Capex_Solar = zeros(6,1)
    Solar_Cost = zeros(6,1)
    Opex_Solar = zeros(6,1)
    Capex_Electrolyzer = zeros(6,1)
    Opex_Electrolyzer = zeros(6,1)
    Capex_Storage = zeros(6,1)
    Opex_Storage = zeros(6,1)
    Capex_Batt = zeros(6,1)
    Opex_Batt = zeros(6,1)
    Capex_Comp = zeros(6,1)
    Opex_Comp = zeros(6,1)
    Opex_Water = zeros(6,1)
    Land = zeros(6,1)
    DAC = zeros(6,1)
    Capex_Grid = zeros(6,1)
    Opex_Grid = zeros(6,1)
    Grid_Compensation = zeros(6,1)
    LCOHg = zeros(6,1)

    # Initializing Electricity-Based Operations Data Vectors for Plotting
    Fraction_Elect_Generated_from_Solar = zeros(6,1)
    Fraction_Elect_Curtailed = zeros(6,1)
    Fraction_Elect_Used = zeros(6,1)
    Fraction_Elect_Generated_from_Grid = zeros(6,1)
    Fraction_Elect_Generated_from_Batt = zeros(6,1)
    Fraction_H2_Used_Directly = zeros(6,1)
    Fraction_H2_Used_from_Storage = zeros(6,1)

    # Initializing Electricity-Based Emissions Data Vectors for Plotting
    CO2Gridg = zeros(6,1)
    CO2Solarg = zeros(6,1)

    # Initializing Data Vectors used to Compare Electricity- and Fossil-Based Pathways
    LCOHall = zeros(7,1)
    CO2Solarall = zeros(7,1)
    CO2Gridall = zeros(7,1)
    CO2NatGasLeakall = zeros(7,1)
    CO2NatGasOtherall = zeros(7,1)
    CO2Directall = zeros(7,1)

    # Initializing Data Vectors for Calculating Grid Use Prices and Carbon Intensity
    AvgGridPrice = zeros(2,1)
    AvgGridCI = zeros(2,1)

    ###### Electricity-Based H2 Production Optimization Model Under Various Reliability Levels #############

    for type in types
        if type == "yearly"
            runs = 1
        elseif type == "monthly"
            runs = 2
        elseif type == "daily"
            runs = 3
        elseif type == "hourly"
            runs = 4
        elseif type == "hourly_GS"
            runs = 5
        elseif type == "hourly_GS2"
            runs = 6
        end
    
        # using the Gurobi solver
        m2 = Model(Gurobi.Optimizer)
    
        ###### DECISION VARIABLES #####
        # Installed Solar (kW)
        @variable(m2, SolarBuild >= 0)
        # Installed Electrolyzer Capacity (kW)
        @variable(m2, ElectBuild >= 0)
        # Storage Tank Size (kg)
        @variable(m2, H2TankSize >= 0)
        # Battery Size (kWh)
        @variable(m2, BattSize >= 0)
        # Land Area Required for PV Facility (acres)
        @variable(m2, SolarAreaTotal >= 0)
        # Electricity Used Each Hour to make H2 (kWh e)
        @variable(m2, ElectUsed[1:nHours] >= 0)
        # Electricity Curtailed Each Hour (kWh e)
        @variable(m2, ElectCurtailed[1:nHours] >= 0)
        # Hydrogen Used Each Hour (kg)
        @variable(m2, H2UsedDirectly[1:nHours] >= 0)
        # Hydrogen Stored Each Hour (kg H2)
        @variable(m2, H2Stored[1:nHours] >= 0)
        # Hydrogen Out of Storage (kg H2)
        @variable(m2, H2UsedfromStorage[1:nHours] >= 0)
        # Hydrogen in Storage (kg)
        @variable(m2, H2InStorage[1:nHours + 1] >= 0)
        # Electricity Stored Each Hour in Battery (kWh)
        @variable(m2, ElectStored[1:nHours] >= 0)
        # Electricity Out of Battery (kWh)
        @variable(m2, ElectUsedfromStorage[1:nHours] >= 0)
        # Electricity in Storage (kWh)
        @variable(m2, ElectInStorage[1:nHours + 1] >= 0)
        
        # Electricity Pulled from the Grid each Hour of the Year (kW/kWh)
        # Grid Connection Size (kW)
        if runs == 5 || runs == 6
            @variable(m2, GridUse[1:nHours] >= 0)
            @variable(m2, GridConnection >= 0)
        end
        # Annual Life Cycle CO2 Emissions
        @variable(m2, AnnualCO2 >= 0)
        # Electrolyzer Ramping (Change in kW per hour)
        @variable(m2, ElectrolyzerRamp[1:nHours - 1] >= -10E9)
        # Battery Ramping (Change in kWh stored per hour)
        @variable(m2, BatteryRamp[1:nHours - 1] >= -10E9)
        # H2 Storage Out Ramping (Change in H2 Used from Storage per hour)
        @variable(m2, H2OutRamp[1:nHours - 1] >= -10E8)
        # Annual Water Costs
        @variable(m2, H2OAnnualCost >= 0)
    
        ###### LINEAR EXPRESSIONS USING DECISON VARIABLES #######
    
        ##### CONSTRAINTS: GENERAL ################
    
        # Solar Electricity must be equal to energy used plus energy curtailed
        if runs != 5 && runs != 6
            for i = 1:nHours
                    @constraint(m2, CapFactorHR[i]*SolarBuild + ElectUsedfromStorage[i] == ElectUsed[i] + ElectCurtailed[i])
            end
        else
            for i = 1:nHours
                @constraint(m2, CapFactorHR[i]*SolarBuild + ElectUsedfromStorage[i] + GridUse[i] == ElectUsed[i] + ElectCurtailed[i])
            end
        end
    
        # Electricity Used each hour must be less than Electrolyzer size
        if runs != 5 && runs != 6
            for i = 1:nHours
                    @constraint(m2, ElectUsed[i] - H2Stored[i]*CompEnergy - ElectStored[i] <= ElectBuild)
            end
        else
            for i = 1:nHours
                    @constraint(m2, ElectUsed[i] - H2Stored[i]*CompEnergy - ElectStored[i] <= ElectBuild)
            end
        end
    
        # Electricity used each hour to make H2 either is stored or used directly 
        for i = 1:nHours
            @constraint(m2, ElectUsed[i] == H2UsedDirectly[i]*ElectInput + H2Stored[i]*(ElectInput + CompEnergy) + ElectStored[i])
        end
    
        if runs == 6
            @constraint(m2, sum(GridUse[i] for i= 1:nHours) <= GridAllow*sum(ElectUsed[i] + ElectCurtailed[i] for i=1:nHours))
        end
    
        if type == "yearly"
            # Must meet yearly hydrogen demand constraint
            @constraint(m2, sum(H2UsedDirectly[i] + H2UsedfromStorage[i] for i = 1:nHours) == H2SupplyDaily*KgPerTon*DAYSperYR)
        elseif type == "monthly"
            # Must meet monthly hydrogen demand constraint
            for i = 1:MONTHSperYR
                first = (DAYSperMONTHcum[i]*HRSperDAY + 1)
                second = DAYSperMONTHcum[i + 1]*HRSperDAY
                @constraint(m2, sum(H2UsedDirectly[j] + H2UsedfromStorage[j] for j = first:second) == H2SupplyDaily*KgPerTon*DAYSperMONTH[i])
            end
        elseif type == "daily"
            # Must meet daily hydrogen demand constraint
            for i = 1:DAYSperYR
                first = (i - 1)*HRSperDAY + 1
                second = i*HRSperDAY
                @constraint(m2, sum(H2UsedDirectly[j] + H2UsedfromStorage[j] for j = first:second) == H2SupplyDaily*KgPerTon)
            end
        elseif runs == 4 || runs == 5 || runs == 6
            # Must meet hourly hydrogen demand constraint
            for i = 1:nHours
                @constraint(m2, H2UsedDirectly[i] + H2UsedfromStorage[i] == H2SupplyDaily*KgPerTon/HRSperDAY)
            end
        end
    
        if runs == 5 || runs == 6
        # Can not pull from the grid more than Connection Capacity
            for i = 1:nHours
                @constraint(m2, GridUse[i] <= GridConnection)
            end
            # Can sell from grid more than Connection Capacity
            for i = 1:nHours
                @constraint(m2, ElectCurtailed[i] <= GridConnection)
            end
        end
    
        # H2 and Battery Storage Balance
        for i = 1:nHours
                @constraint(m2, H2InStorage[i + 1] - H2InStorage[i] == H2Stored[i] - H2UsedfromStorage[i] - H2Stored[i]*(1 - TankEff))
                @constraint(m2, ElectInStorage[i + 1] - ElectInStorage[i] == ElectStored[i] - ElectUsedfromStorage[i] - ElectStored[i]*(1 - BattEff))
        end
    
        # H2 and Electricity storage at start and end of run must be equal
        @constraint(m2, H2InStorage[1]==H2InStorage[end])
        @constraint(m2, ElectInStorage[1]==ElectInStorage[end])
    
        # H2/Battery Storage Can not Exceed Tank Size
        for i = 1:nHours
                @constraint(m2, H2InStorage[i] <= H2TankSize)
                @constraint(m2, H2Stored[i] <= H2TankSize)
                @constraint(m2, H2UsedfromStorage[i] <= H2TankSize)
                @constraint(m2, ElectInStorage[i] <= BattSize)
                @constraint(m2, ElectStored[i] <= BattSize)
                @constraint(m2, ElectUsedfromStorage[i] <= BattSize)
        end
    
        # Land Area Calculation
        @constraint(m2, SolarAreaTotal == SolarBuild/KWperMW*SolarArea)
    
        # Total Annual CO2 Emissions
        if runs != 5 && runs != 6
            @constraint(m2, AnnualCO2 == LCESolarPV*sum(CapFactorHR[i].*SolarBuild for i = 1:nHours))
        else
            @constraint(m2, AnnualCO2 == LCESolarPV*sum(CapFactorHR[i].*SolarBuild for i = 1:nHours) + sum(GridUse[i]*CIgrid[i] for i = 1:nHours))
        end
    
        # Total Annual Water Cost
        @constraint(m2, H2OAnnualCost == WaterCost/KgPerTon*WaterUse*sum(H2UsedfromStorage[i] + H2UsedDirectly[i] for i = 1:nHours))
    
    
        # Electrolyzer Ramping
        for i = 1:nHours - 1
                @constraint(m2, ElectrolyzerRamp[i] >= (ElectUsed[i] - ElectUsed[i + 1]))
                @constraint(m2, ElectrolyzerRamp[i] >= -(ElectUsed[i] - ElectUsed[i + 1]))
        end
    
        # Battery Ramping
        for i = 1:nHours - 2
            @constraint(m2, BatteryRamp[i] >= (ElectInStorage[i] - ElectInStorage[i + 1]))
            @constraint(m2, BatteryRamp[i] >= -(ElectInStorage[i] - ElectInStorage[i + 1]))
        end
    
        # H2 Out Ramping
        for i = 1:nHours - 1
                @constraint(m2, H2OutRamp[i] >= (H2UsedfromStorage[i] - H2UsedfromStorage[i + 1]))
                @constraint(m2, H2OutRamp[i] >= -(H2UsedfromStorage[i] - H2UsedfromStorage[i + 1]))
        end
    
        #### Objective ####
    
        # Minimize the Levelized Cost of Produced Hydrogen
        if runs != 5 && runs != 6 # without grid connection
            @objective(m2, Min, ((CapRF1 + SolarOM)*CapSolar*SolarBuild 
                            + (CapRF1 + ElectOM)*CapElect*ElectBuild
                                            + (CapRF1 + H2StorageOM)*H2TankSize*H2StoreCapex
                                            + (CapRF1 + BattStorageOM)*BattSize*BattCapex
                                                    + SolarAreaTotal*SolarLandCost
                                                            + AnnualCO2*DACCost*CO2Lever
                                                            + H2OAnnualCost
                                                                    + sum(ElectrolyzerRamp[i] for i = 1:nHours-1)*RampCost1
                                                                    + sum(BatteryRamp[i] for i = 1:nHours-2)*RampCost1
                                                                    + sum(H2OutRamp[i] for i = 1:nHours-1)*RampCost2) 
                                                                            /(YearlyEnergyDemand/ElectInput))
        else # with grid connection
            @objective(m2, Min, ((CapRF1 + SolarOM)*CapSolar*SolarBuild 
                            + (CapRF1 + ElectOM)*CapElect*ElectBuild
                                            + (CapRF1 + H2StorageOM)*H2TankSize*H2StoreCapex
                                            + (CapRF1 + BattStorageOM)*BattSize*BattCapex
                                                    + sum(GridUse[i]*LMPgrid[i] for i = 1:nHours)
                                                    + GridConnection*(DemandCharge+SiteCharge)*MONTHSperYR
                                                    - sum(ElectCurtailed[i]*NSCR[i] for i = 1:nHours)
                                                    + (CapRF1 + GridOM)*GridConnectFee*GridConnection
                                                    + SolarAreaTotal*SolarLandCost
                                                            + AnnualCO2*DACCost*CO2Lever
                                                            + H2OAnnualCost
                                                                    + sum(ElectrolyzerRamp[i] for i = 1:nHours-1)*RampCost1
                                                                    + sum(BatteryRamp[i] for i = 1:nHours-2)*RampCost1
                                                                    + sum(H2OutRamp[i] for i = 1:nHours-1)*RampCost2) 
                                                                            /(YearlyEnergyDemand/ElectInput))
        end
        
        # Solve the model or redo the model run if the model does not converge
        optimize!(m2)
        try
            ObjValue = objective_value(m2)
        catch
            break
        end

        # Store Optimal Values for Decision Variables
        OptimalSolarBuild = value.(SolarBuild)
        OptimalElectBuild = value.(ElectBuild)
        OptimalLandBuild = value.(SolarAreaTotal)
        OptimalElectUsed = value.(ElectUsed)
        OptimalElectCurtailed = value.(ElectCurtailed)
        CO2perYR = value.(AnnualCO2)
        OptimalWaterCost = value.(H2OAnnualCost)
        OptimalH2UsedDirectly = value.(H2UsedDirectly)
        OptimalH2Stored = value.(H2Stored)
        OptimalH2UsedfromStorage = value.(H2UsedfromStorage)
        OptimalH2TankSize = value.(H2TankSize)
        OptimalH2InStorage = value.(H2InStorage)
        OptimalElectStored = value.(ElectStored)
        OptimalElectUsedfromStorage = value.(ElectUsedfromStorage)
        OptimalBattSize = value.(BattSize)
        OptimalElectInStorage = value.(ElectInStorage)

        if runs == 5 || runs == 6
            OptimalGridBuild = value.(GridConnection)
            OptimalGridUse = value.(GridUse)
            CO2perYR = CO2perYR + sum(OptimalGridUse.*(CIgrid_adjusted.-CIgrid))
            AvgGridPrice[runs - 4] = sum(OptimalGridUse.*LMPgrid)/sum(OptimalGridUse)
            AvgGridCI[runs - 4] = sum(OptimalGridUse.*CIgrid_adjusted)/sum(OptimalGridUse)
        end
        
        # Storing Infrastructure Cost Values for each electricity-based pathway
        Capex_Solar[runs] = CapRF1*CapSolar*OptimalSolarBuild/(YearlyEnergyDemand/ElectInput)
        Solar_Cost[runs] = (CapRF1+SolarOM)*CapSolar*OptimalSolarBuild/(YearlyEnergyDemand)
        Opex_Solar[runs] = SolarOM*CapSolar*OptimalSolarBuild/(YearlyEnergyDemand/ElectInput)
        Capex_Electrolyzer[runs] = CapRF1*CapElect*OptimalElectBuild/(YearlyEnergyDemand/ElectInput)
        Opex_Electrolyzer[runs] = ElectOM*CapElect*OptimalElectBuild/(YearlyEnergyDemand/ElectInput)
        Capex_Storage[runs] = CapRF1*OptimalH2TankSize*H2StoreCapex/(YearlyEnergyDemand/ElectInput)
        Opex_Storage[runs] = H2StorageOM*OptimalH2TankSize*H2StoreCapex/(YearlyEnergyDemand/ElectInput)
        Capex_Batt[runs] = CapRF1*OptimalBattSize*BattCapex/(YearlyEnergyDemand/ElectInput)
        Opex_Batt[runs] = BattStorageOM*OptimalBattSize*BattCapex/(YearlyEnergyDemand/ElectInput)
        
        if runs == 5 || runs == 6
            Capex_Grid[runs] = CapRF1*OptimalGridBuild*GridConnectFee/(YearlyEnergyDemand/ElectInput)
            Opex_Grid[runs] = (GridOM*OptimalGridBuild*GridConnectFee + sum(OptimalGridUse.*LMPgrid) + OptimalGridBuild*(DemandCharge+SiteCharge)*MONTHSperYR + GridFixedCharge*MONTHSperYR)/(YearlyEnergyDemand/ElectInput)
            Grid_Compensation[runs] = (-1)*sum(OptimalElectCurtailed.*NSCR)/(YearlyEnergyDemand/ElectInput)
        end
        Land[runs] = OptimalLandBuild*SolarLandCost/(YearlyEnergyDemand/ElectInput)
        Opex_Water[runs] = OptimalWaterCost/(YearlyEnergyDemand/ElectInput)
        DAC[runs] = CO2Lever*CO2perYR*DACCost/(YearlyEnergyDemand/ElectInput)

        # Storing emissions data for each pathway
        CO2Solarg[runs] = LCESolarPV*sum(CapFactorHR[i].*OptimalSolarBuild for i = 1:nHours)/(H2SupplyDaily*KgPerTon*DAYSperYR)
        if runs > 4
            CO2Gridg[runs] = sum(OptimalGridUse[i]*CIgrid_adjusted[i] for i = 1:nHours)/(H2SupplyDaily*KgPerTon*DAYSperYR)
            CO2Gridall[runs+1] = CO2Gridg[runs]
        end

        # Storing high-level operations data for each pathway
        Fraction_Elect_Curtailed[runs] = sum(OptimalElectCurtailed)/sum(OptimalElectUsed + OptimalElectCurtailed)
        Fraction_Elect_Used[runs] = sum(OptimalElectUsed)/sum(OptimalElectUsed + OptimalElectCurtailed)

        if runs != 5 && runs !=6
            Fraction_Elect_Generated_from_Solar[runs] = sum(OptimalSolarBuild.*CapFactorHR)/sum(OptimalElectUsed + OptimalElectCurtailed)
            Fraction_Elect_Generated_from_Batt[runs] = sum(OptimalElectUsedfromStorage)/sum(OptimalElectUsed + OptimalElectCurtailed)
            Fraction_H2_Used_Directly[runs] = sum(OptimalH2UsedDirectly)/sum(OptimalH2UsedDirectly + OptimalH2UsedfromStorage)
            Fraction_H2_Used_from_Storage[runs] = sum(OptimalH2UsedfromStorage)/sum(OptimalH2UsedDirectly + OptimalH2UsedfromStorage)
        else
            Fraction_Elect_Generated_from_Solar[runs] = sum(OptimalSolarBuild.*CapFactorHR)/sum(OptimalElectUsed + OptimalElectCurtailed)
            Fraction_Elect_Generated_from_Batt[runs] = sum(OptimalElectUsedfromStorage)/sum(OptimalElectUsed + OptimalElectCurtailed)
            Fraction_Elect_Generated_from_Grid[runs] = sum(OptimalGridUse)/sum(OptimalElectUsed + OptimalElectCurtailed)
            Fraction_H2_Used_Directly[runs] = sum(OptimalH2UsedDirectly)/sum(OptimalH2UsedDirectly + OptimalH2UsedfromStorage)
            Fraction_H2_Used_from_Storage[runs] = sum(OptimalH2UsedfromStorage)/sum(OptimalH2UsedDirectly + OptimalH2UsedfromStorage)
        end
        
        # Daily operations data processing for each electricity-based pathway
        DailyECurtailedFraction = zeros(365,6)
        DailyESolarFraction = zeros(365,6)
        DailyEUsedFraction = zeros(365,6)
        DailyEGridFraction = zeros(365,6)
        DailyH2UsedDirectlyFraction = zeros(365,6)
        DailyH2UsedfromStorageFraction = zeros(365,6)
        DailyH2Stored = zeros(365,6)
        DailyH2InStorage = zeros(365,6)
        DailyH2UsedfromStorage = zeros(365,6)

        for i = 1:DAYSperYR
            first = (i - 1)*HRSperDAY + 1
            second = i*HRSperDAY
            DailyECurtailedFraction[i,runs] = -sum(OptimalElectCurtailed[first:second])/sum(OptimalElectUsed[first:second] + OptimalElectCurtailed[first:second])
            DailyEUsedFraction[i,runs]  = -sum(OptimalElectUsed[first:second])/sum(OptimalElectUsed[first:second] + OptimalElectCurtailed[first:second])
            DailyH2UsedDirectlyFraction[i,runs] = sum(OptimalH2UsedDirectly[first:second])/sum(OptimalH2UsedDirectly[first:second] + OptimalH2UsedfromStorage[first:second])
            DailyH2UsedfromStorageFraction[i,runs] = sum(OptimalH2UsedfromStorage[first:second])/sum(OptimalH2UsedDirectly[first:second] + OptimalH2UsedfromStorage[first:second])
            DailyH2Stored[i,runs] = sum(OptimalH2Stored[first:second])
            DailyH2UsedfromStorage[i,runs] = sum(OptimalH2UsedfromStorage[first:second])
            
            if i == 365
                DailyH2InStorage[i,runs] = sum(OptimalH2InStorage[first:second-1])
            else
                DailyH2InStorage[i,runs] = sum(OptimalH2InStorage[first:second])
            end

            if runs != 5 && runs !=6
                DailyESolarFraction[i,runs] = sum(OptimalSolarBuild.* CapFactorHR[first:second])/sum(OptimalElectUsed[first:second] + OptimalElectCurtailed[first:second])
            else
                DailyESolarFraction[i,runs] = sum(OptimalSolarBuild.* CapFactorHR[first:second])/sum(OptimalElectUsed[first:second] + OptimalElectCurtailed[first:second])
                DailyEGridFraction[i,runs]  = sum(OptimalGridUse[first:second])/sum(OptimalElectUsed[first:second] + OptimalElectCurtailed[first:second])
            end
        end

        # Hourly operations data processing for each electricity-based pathway
        HourlyECurtailed = zeros(7*HRSperDAY,6)
        HourlyESolar = zeros(7*HRSperDAY,6)
        HourlyEUsed = zeros(7*HRSperDAY,6)
        HourlyEGrid = zeros(7*HRSperDAY,6)
        HourlyH2Stored = zeros(7*HRSperDAY,6)
        HourlyH2InStorage = zeros(7*HRSperDAY,6)
        HourlyH2UsedfromStorage = zeros(7*HRSperDAY,6)
        HourlyH2UsedDirectly = zeros(7*HRSperDAY,6)

        # Assume a week in May (130-137 days of year)
        first = 130*(HRSperDAY) + 1
        second = first + HRSperDAY*7 - 1
        HourlyECurtailed[:,runs] = -OptimalElectCurtailed[first:second]
        HourlyEUsed[:,runs]  = -OptimalElectUsed[first:second]
        HourlyH2UsedDirectly[:,runs] = OptimalH2UsedDirectly[first:second]
        HourlyH2UsedfromStorage[:,runs] = OptimalH2UsedfromStorage[first:second]
        HourlyH2Stored[:,runs] = OptimalH2Stored[first:second]
        HourlyH2InStorage[:,runs] = OptimalH2InStorage[first:second]
        HourlyH2UsedfromStorage[:,runs] = OptimalH2UsedfromStorage[first:second]

        if runs != 5 && runs !=6
            HourlyESolar[:,runs] = OptimalSolarBuild.* CapFactorHR[first:second]
        else
            HourlyESolar[:,runs] = OptimalSolarBuild.* CapFactorHR[first:second]
            HourlyEGrid[:,runs]  = OptimalGridUse[first:second]
        end

    end

    # Calculation LCOH for each pathway
    LCOHg = Capex_Solar + Opex_Solar + Capex_Electrolyzer + Opex_Electrolyzer + Capex_Storage + Opex_Storage + Capex_Batt + Opex_Batt + Capex_Grid + Opex_Grid + Grid_Compensation + Land + Opex_Water + DAC
    
    # Skip run if model does not converge correctly
    if 0 in LCOHg
        continue
    end

    LCOHall[1:3] = Capex_Solar[4:6] + Opex_Solar[4:6] + Capex_Electrolyzer[4:6] + Opex_Electrolyzer[4:6] + Capex_Storage[4:6] + Opex_Storage[4:6] + Capex_Batt[4:6] + Opex_Batt[4:6] + Capex_Grid[4:6] + Opex_Grid[4:6] + Grid_Compensation[4:6] + Land[4:6] + Opex_Water[4:6] + DAC[4:6]

    ############ START OF FOSSIL-BASED  TECHNO-ECONOMIC ANALYSIS #############


    ###### GLOBAL CONVERSIONS #######################
    HRSperDAY = 24

    KWperMW = 1000

    KgPerTon = 1000

    DAYSperYR = 365

    MONTHSperYR = 12

    MJperKWH = 3.6

    MJperTHERM = 105.5

    DAYSperMONTHav = 30.4


    ###### PARAMETERS and DATA ###############

    # Hydrogen Plant Production  (ton/day)
    H2Capacity = H2SupplyDaily

    # Hydrogen Plant Capacity Factor
    # H2CapFactor = 0.9

    # Hydrogen Supply (ton/day)
    H2SupplyDaily = H2SupplyDaily

    # H2 Energy Content (kWh/kg)
    H2Energy = 33.3

    # CI Grid (Consider same input as green H2 model)
    # CIgrid = CIgrid

    # SMR Plant Life (yrs)
    YrSMR = Base.rand(TriangularDist(20,40,30))

    # Methane Leakage Rate
    CH4Leakage = Base.rand(Uniform(0.0,0.04))

    # CO2 Removal Cost ($/kg CO2)
    # DACCost = 0.200  [same as used in electricity-based model]

    # Methane GWP Potency
    GWPPotency = rand((30,85))

    # Methane GWP 100
    GWP20= 85
    # Methane GWP 100
    GWP100= 30

    # Natural Gas Energy Content (MJ/kg)
    NatGasEnergy = 52

    # Other Life Cycle Natural Gas Emissions (kg CO2/ kg NG)
    LCEmissionsNG = Base.rand(TriangularDist(0.1,0.5,0.3))

    # Weighted Average Cost of Capital (WACC)
    # SMR Equipment Capital Recovery Factor (%)
    WACC = Base.rand(TriangularDist(0.04,0.1,0.04))
    CapRF2 = WACC/(1-(1+WACC)^(-YrSMR))

    # Natural Gas Prices CA ($/therm gas)
    GasPrice = Base.rand(TriangularDist(0.3,1.5,0.3))

    # Max Grid use (kWh/ kg H2)
    MaxGridUse = 4.0

    # SMUD Electricity Price CA (/kWh)
    ElectricityPrice = (mean(LMPgrid) + (CapRF2 + GridOM)*GridConnectFee/DAYSperYR/HRSperDAY + (SiteCharge+DemandCharge)/(DAYSperYR/MONTHSperYR)/HRSperDAY)
    
    # SMR/CO2 Capture Scale Exponent (N/A)
    ScaleExp = Base.rand(TriangularDist(0.55,0.8,0.55))

    #### SMR w/o CCS PARAMETERS ##########

    # SMR Baseline Size (tonne/day H2)
    SMRBaselineSize = 483

    # SMR Baseline Cost ($/(tonne/day))
    SMRBaselineCost = 763E3

    # SMR Capital Cost ($)
    CapExSMR = (H2Capacity/SMRBaselineSize)^ScaleExp * SMRBaselineCost * SMRBaselineSize

    # SMR Fixed O&M Charges ($/(kg/day) each year)
    SMRFixedOM = 22 * (H2Capacity/SMRBaselineSize)^ScaleExp / (H2Capacity/SMRBaselineSize)

    # SMR Electricity Usage (kWh electricity per kg H2 produced)
    SMRElectricityInput = 0.65

    # SMR Natural Gas Usage (kg gas / kg H2 produced)
    SMRInputNatGas = 3.53

    # Non-Energy Variable O&M ($/(kg/day) each yr)
    SMRVarOM = 13.7

    # Direct Plant Emissions (kg CO2/kg H2)
    SMRDirectEmissions = 9.3

    # SMR CO2e Emissions Calculation (kg CO2/kg H2)
    CO2Grey = (CH4Leakage*SMRInputNatGas*GWPPotency + mean(CIgrid_adjusted)*SMRElectricityInput + SMRDirectEmissions + SMRInputNatGas*LCEmissionsNG)
    CO2Gridall[1] = mean(CIgrid_adjusted)*SMRElectricityInput
    CO2NatGasLeakall[1] = CH4Leakage*SMRInputNatGas*GWPPotency
    CO2Directall[1] = SMRDirectEmissions
    CO2NatGasOtherall[1] = SMRInputNatGas*LCEmissionsNG


    # SMR H2 LCOH Calculation
    CapEx_SMR = CapRF2*CapExSMR/H2Capacity/DAYSperYR/KgPerTon
    FixedOM_SMR = SMRFixedOM/DAYSperYR
    OpEx_NatGas_SMR = SMRInputNatGas*NatGasEnergy/MJperTHERM*GasPrice
    OpEx_Electricity_SMR = SMRElectricityInput*ElectricityPrice
    OpEx_Other_Var_SMR = SMRVarOM/DAYSperYR
    DAC_SMR = CO2Grey*DACCost

    #### SMR with 55% PROCESS CO2 CAPTURE PARAMETERS #####

    # SMR55 Baseline Size (tonne/day H2)
    SMR55BaselineSize = 483

    # SMR55 Baseline Cost ($/(tonne/day))
    SMR55BaselineCost = 800E3

    # SMR55 Capital Cost ($)
    CapExSMR55 = (H2Capacity/SMR55BaselineSize)^ScaleExp * SMR55BaselineCost * SMR55BaselineSize


    # SMR Fixed O&M Charges (as fraction CAPEX per yr)
    #SMR55FixedOM = 0.03

    # SMR Fixed O&M Charges ($/(kg/day) each year)
    SMR55FixedOM = 32 * (H2Capacity/SMR55BaselineSize)^ScaleExp / (H2Capacity/SMR55BaselineSize)


    # CO2 Captured (kg CO2 per kg H2)
    CO2CapturedSMR55 = 5.2

    # CO2 Emissions Released (kg CO2 per kg H2)
    DirectEmissionsSMR55 = 4.1

    # Process Capture Baseline Size (tonne/day CO2 captured)
    CC55BaselineSize = 1090

    # Process Capture Baseline Capital Cost ($)
    CC55BaselineCost = 75E6

    # Process CO2 capture Capital Cost ($)
    CapExCC55 = (H2Capacity*CO2CapturedSMR55/CC55BaselineSize)^ScaleExp * CC55BaselineCost

    # Process CO2 capture Fixed O&M Charges (as fraction CC CAPEX per yr)
    CC55FixedOM = 0.07

    # SMR with Process CO2 capture Electricity input (kWh electricity per kg H2 capacity)
    SMR55ElectricityInput = 1.5

    # SMR with Process CO2 capture Natural Gas input (kg NG / kg H2 capacity)
    SMR55InputNatGas = 3.56

    # Non-Energy Variable O&M ($/(kg/day) each yr)
    SMR55VarOM = 20

    # Cost for CO2 Transport and Storage ($/ton CO2 captured)
    CO2TandSCost = Base.rand(TriangularDist(10,20,10))

    # SMR55 CO2e Emissions Calculation (kg CO2/kg H2)
    CO2SMR55 = (CH4Leakage*SMR55InputNatGas*GWPPotency + mean(CIgrid_adjusted)*SMR55ElectricityInput + DirectEmissionsSMR55 + SMR55InputNatGas*LCEmissionsNG)
    CO2Gridall[2] = mean(CIgrid_adjusted)*SMR55ElectricityInput
    CO2NatGasLeakall[2] = CH4Leakage*SMR55InputNatGas*GWPPotency
    CO2Directall[2] = DirectEmissionsSMR55
    CO2NatGasOtherall[2] = SMR55InputNatGas*LCEmissionsNG

    # SMR55 H2 LCOH Calculation
    CapEx_SMR55 = CapRF2*CapExSMR55/H2Capacity/DAYSperYR/KgPerTon
    CapEx_CC55 = CapRF2*CapExCC55/H2Capacity/DAYSperYR/KgPerTon
    FixedOM_SMR55 = SMR55FixedOM/DAYSperYR
    OpEx_NatGas_SMR55 = SMR55InputNatGas*NatGasEnergy/MJperTHERM*GasPrice
    OpEx_Electricity_SMR55 = SMR55ElectricityInput*ElectricityPrice
    OpEx_Other_Var_SMR55 = SMR55VarOM/DAYSperYR
    CO2TandSSMR55 = CO2TandSCost/KgPerTon*CO2CapturedSMR55
    DAC_SMR55 = DACCost*CO2SMR55

    ### SMR with 96% POST-COMBUSTION and PROCESS CO2 CAPTURE PARAMETERS ######

    # SMRCCS Baseline Size (tonne/day H2)
    SMRCCSBaselineSize = 483

    # SMR Baseline Cost ($/(tonne/day))
    SMRCCSBaselineCost = 1124E3

    # Process Capture Baseline Cost ($/(tonne/day))
    ProCCBaselineCost = 85E3

    # Post-Combustion Capture Baseline Cost ($/(tonne/day))
    PostCCBaselineCost = 647E3

    # SMR Capital Cost ($)
    CapExSMRCCS = (H2Capacity/SMRCCSBaselineSize)^ScaleExp * SMRCCSBaselineCost * SMRCCSBaselineSize

    # Process Capture Capital Cost ($)
    CapExProCCS = (H2Capacity/SMRCCSBaselineSize)^ScaleExp * ProCCBaselineCost * SMRCCSBaselineSize

    # Post-Combustion Capture Capital Cost ($)
    CapExPostCCS = (H2Capacity/SMRCCSBaselineSize)^ScaleExp * PostCCBaselineCost * SMRCCSBaselineSize

    # SMRCCS Fixed O&M Charges ($/(kg/day) each year)
    SMRCCSFixedOM = 48.3 * (H2Capacity/SMRCCSBaselineSize)^ScaleExp / (H2Capacity/SMRCCSBaselineSize)

    # SMRCCS Electricity Usage (kWh electricity per kg H2 produced)
    SMRCCSElectricityInput = 2.04

    # SMRCCS Natural Gas Usage (kg gas / kg H2 produced)
    SMRCCSInputNatGas = 3.75

    # Non-Energy Variable O&M ($/(kg/day) each yr)
    SMRCCSVarOM = 31.8

    # CO2 Captured (kg CO2 per kg H2)
    CO2CapturedSMRCCS = 10.1

    # CO2 Emissions Released (kg CO2 per kg H2)
    SMRCCSDirectEmissions = 0.4

    # SMRCCS CO2e Emissions Calculation (kg CO2/kg H2)
    CO2SMRCCS = (CH4Leakage*SMRCCSInputNatGas*GWPPotency + mean(CIgrid_adjusted)*SMRCCSElectricityInput + SMRCCSDirectEmissions + SMRCCSInputNatGas*LCEmissionsNG)
    CO2Gridall[3] = mean(CIgrid_adjusted)*SMRCCSElectricityInput
    CO2NatGasLeakall[3] = CH4Leakage*SMRCCSInputNatGas*GWPPotency
    CO2Directall[3] = SMRCCSDirectEmissions
    CO2NatGasOtherall[3] = SMRCCSInputNatGas*LCEmissionsNG

    # SMR-CCS H2 LCOH Calculation
    CapEx_SMRCCS = CapRF2*CapExSMRCCS/H2Capacity/DAYSperYR/KgPerTon
    CapEx_ProCCS = CapRF2*CapExProCCS/H2Capacity/DAYSperYR/KgPerTon
    CapEx_PostCCS = CapRF2*CapExPostCCS/H2Capacity/DAYSperYR/KgPerTon
    FixedOM_SMRCCS = SMRCCSFixedOM/DAYSperYR
    OpEx_NatGas_SMRCCS = SMRCCSInputNatGas*NatGasEnergy/MJperTHERM*GasPrice
    OpEx_Electricity_SMRCCS = SMRCCSElectricityInput*ElectricityPrice
    OpEx_Other_Var_SMRCCS = SMRCCSVarOM/DAYSperYR
    CO2TandSSMR = CO2TandSCost/KgPerTon*CO2CapturedSMRCCS
    DAC_SMRCCS = DACCost*CO2SMRCCS

    ## ATR with 94.5% CO2 CAPTURE PARAMETERS #####

    # ATRCCS Baseline Size (tonne/day H2)
    ATRBaselineSize = 660

    # ATR Baseline Cost ($/(tonne/day))
    ATRBaselineCost = 839E3

    # ATR Capture Baseline Cost ($/(tonne/day))
    ATRCCBaselineCost = 220E3

    # ATR Air Separation Unit Baseline Cost ($/(tonne/day))
    ASUBaselineCost = 408E3

    # ATR Capital Cost ($)
    CapExATR = (H2Capacity/ATRBaselineSize)^ScaleExp * ATRBaselineCost * ATRBaselineSize

    # Process ATR Capture Capital Cost ($)
    CapExATRCCS = (H2Capacity/ATRBaselineSize)^ScaleExp * ATRCCBaselineCost * ATRBaselineSize

    # ATR Air Separation Unit Capital Cost ($)
    CapExASU = (H2Capacity/ATRBaselineSize)^ScaleExp * ASUBaselineCost * ATRBaselineSize

    # ATR Fixed O&M Charges ($/(kg/day) each year)
    ATRCCSFixedOM = 37.47 * (H2Capacity/ATRBaselineSize)^ScaleExp / (H2Capacity/ATRBaselineSize)

    # ATRCCS Electricity Usage (kWh electricity per kg H2 produced)
    ATRCCSElectricityInput = 4.0

    # ATR Natural Gas Usage (kg gas / kg H2 produced)
    ATRCCSInputNatGas = 3.52

    # Non-Energy Variable O&M ($/(kg/day) each yr)
    ATRCCSVarOM = 23.3

    # CO2 Captured (kg CO2 per kg H2)
    CO2CapturedATRCCS = 8.6

    # CO2 Emissions Released (kg CO2 per kg H2)
    ATRCCSDirectEmissions = 0.5

    # ATRCCS CO2e Emissions Calculation (kg CO2/kg H2)
    CO2ATRCCS = (CH4Leakage*ATRCCSInputNatGas*GWPPotency + mean(CIgrid_adjusted)*ATRCCSElectricityInput + ATRCCSDirectEmissions + ATRCCSInputNatGas*LCEmissionsNG)
    CO2Gridall[4] = mean(CIgrid_adjusted)*ATRCCSElectricityInput
    CO2NatGasLeakall[4] = CH4Leakage*ATRCCSInputNatGas*GWPPotency
    CO2Directall[4] = ATRCCSDirectEmissions
    CO2NatGasOtherall[4] = ATRCCSInputNatGas*LCEmissionsNG

    # ATR-CCS H2 LCOH Calculation
    CapEx_ATR = CapRF2*CapExATR/H2Capacity/DAYSperYR/KgPerTon
    CapEx_ATRCCS = CapRF2*CapExATRCCS/H2Capacity/DAYSperYR/KgPerTon
    CapEx_ASU = CapRF2*CapExASU/H2Capacity/DAYSperYR/KgPerTon
    FixedOM_ATR = ATRCCSFixedOM/DAYSperYR
    OpEx_NatGas_ATR = ATRCCSInputNatGas*NatGasEnergy/MJperTHERM*GasPrice
    OpEx_Electricity_ATR = ATRCCSElectricityInput*ElectricityPrice
    OpEx_Other_Var_ATR = ATRCCSVarOM/DAYSperYR
    CO2TandSATR = CO2TandSCost/KgPerTon*CO2CapturedATRCCS
    DAC_ATRCCS = DACCost*CO2ATRCCS


    # LCOH Breakdowns
    #[Grey H2, Blue H2 (SMR1), Blue H2 (SMR2), Blue H2 (ATR)]
    CapExReformerv = [CapEx_SMR, CapEx_SMR55, CapEx_SMRCCS, CapEx_ATR]
    CapExProCCv = [0, CapEx_CC55, CapEx_ProCCS, CapEx_ATRCCS]
    CapExPostCCv = [0, 0, CapEx_PostCCS, 0]
    CapExASUv = [0, 0, 0, CapEx_ASU]
    FixedOMv = [FixedOM_SMR, FixedOM_SMR55, FixedOM_SMRCCS, FixedOM_ATR]
    VarOMv = [OpEx_Other_Var_SMR, OpEx_Other_Var_SMR55, OpEx_Other_Var_SMRCCS, OpEx_Other_Var_ATR]
    OpEx_NatGasv = [OpEx_NatGas_SMR, OpEx_NatGas_SMR55, OpEx_NatGas_SMRCCS, OpEx_NatGas_ATR]
    OpEx_Electricityv = [OpEx_Electricity_SMR, OpEx_Electricity_SMR55, OpEx_Electricity_SMRCCS, OpEx_Electricity_ATR]
    CO2TandSv = [0, CO2TandSSMR55, CO2TandSSMR, CO2TandSATR]
    DACv = [DAC_SMR, DAC_SMR55, DAC_SMRCCS, DAC_ATRCCS]

    GreySMRLCOH = CapEx_SMR+FixedOM_SMR+OpEx_Other_Var_SMR+OpEx_NatGas_SMR+OpEx_Electricity_SMR
    BlueSMR55LCOH = CapEx_SMR55 + CapEx_CC55 + FixedOM_SMR55 + OpEx_Electricity_SMR55 + OpEx_Other_Var_SMR55 + OpEx_NatGas_SMR55 + CO2TandSSMR55
    BlueSMRLCOH = CapEx_SMRCCS+CapEx_ProCCS+CapEx_PostCCS+FixedOM_SMRCCS+OpEx_Other_Var_SMRCCS+OpEx_NatGas_SMRCCS+OpEx_Electricity_SMRCCS+CO2TandSSMR
    BlueATRLCOH = CapEx_ATR+CapEx_ATRCCS+CapEx_ASU+FixedOM_ATR+OpEx_Other_Var_ATR+OpEx_NatGas_ATR+OpEx_Electricity_ATR+CO2TandSATR
    LCOHsFossil = [GreySMRLCOH, BlueSMR55LCOH, BlueSMRLCOH, BlueATRLCOH]
    
    # Below LCOH vectors needed to compare results to electricity-based pathways
    LCOHall[4] = GreySMRLCOH + DACv[1]
    LCOHall[5] = BlueSMR55LCOH + DACv[2]
    LCOHall[6] = BlueSMRLCOH + DACv[3]
    LCOHall[7] = BlueATRLCOH + DACv[4]


    # Push all LCOH results into a matrix to be stored for error bar creation
    push!(LCOHallEB,vcat(LCOHg,LCOHall[4:7]))
    print(vcat(LCOHg,LCOHall[4:7]))

    # Move to next iteration of monte carlo simulation
    global m = m + 1
end

# Reordering LCOH matrix
LCOHallEBmatrix = zeros(50,10)
for i in 1:50
    for j in 1:10
        LCOHallEBmatrix[i,j] = LCOHallEB[i][j]
    end
end

# Choosing low and high end LCOH results from Monte Carlo simulation to be error bounds for LCOH Figure 4.
# We use a 96% Confidence Interval for error bar construction.
using Base.Sort

lowEB = zeros(10,1)
highEB = zeros(10,1)
for k = 1:10
    lowEB[k] = sort(LCOHallEBmatrix[:,k]; alg=Sort.PartialQuickSort(5))[2]
    highEB[k] = sort(LCOHallEBmatrix[:,k]; alg=Sort.PartialQuickSort(5),rev=true)[2]
end

# Uncomment lines below if first time running code.
# import Pkg;
# Pkg.add("FileIO")
# Pkg.add("JLD2")

# Saving error bounds so they can be loaded into primary model.
using FileIO, JLD2
FileIO.save("LCOHallEBmatrix_CA.jld2", "LCOHallEBmatrix", LCOHallEBmatrix)
FileIO.save("lowEB_CA.jld2", "lowEB", lowEB)
FileIO.save("highEB_CA.jld2", "highEB", highEB)


