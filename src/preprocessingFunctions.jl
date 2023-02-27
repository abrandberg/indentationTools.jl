### Packages used
using Statistics
using NaNStatistics
using Optim
using DelimitedFiles
using DataFrames
using CSV
using Distributions
using Plots

"""
ReadIBW(filename) 

Reads a .IBW binary file, extracting the position and force signal of the AFM indenter. 

filename is the string to a file which either needs to be absolute or on the path.
"""
function readIBW(filename::String)

    A = open( filename, "r")   
    firstByte = read(A, Int8)
    close(A)
    
    if firstByte == 0
        machineFormat = 'b'
    else
        machineFormat = 'l'
    end
    
    A = open( filename, "r")
    version         = read(A, Int16)  
    checksum        = ltoh(read(A,UInt16)) # Does not work           
    wfmSize         = ltoh(read(A, Int32))
    formulaSize     = ltoh(read(A, Int32))
    noteSize        = ltoh(read(A, Int32)) # Does not work
    dataEUnitsSize  = ltoh(read(A, Int32))
    dimEUnitsSize   = Array{Int32, 1}(undef, 4); read!(ltoh(A),dimEUnitsSize)
    dimLabelsSize   = Array{Int32, 1}(undef, 4); read!(ltoh(A),dimLabelsSize)
    sIndicesSize    = read(A, Int32)
    optionSize1     = read(A, Int32)
    optionSize2     = read(A, Int32)   
    ignore          = read(A, Int32)
    CreationDate    = read(A, UInt32)
    modData         = read(A, UInt32)
    npnts           = read(A, Int32)
    
    type = read(A, Int16)
    if type == 2
        datatype = "single"
    elseif type == 4
        datatype = "double"
    else
        println("ERROR")
    end

    read(A, Int16)
    for aLoop = 1:6
        read(A, Char)
    end
    read(A, Int16)

    for aLoop = 1:32
        read(A, Char)
    end
    for aLoop = 1:6
        read(A, Int32)
    end

    sfA = Array{Float64, 1}(undef, 4)
    read!(A, sfA)
    sfB = Array{Float64, 1}(undef, 4)
    read!(A, sfB)

    dUnits = Array{Char, 1}(undef,4)
    for aLoop = 1:4
        dUnits[aLoop] = read(A, Char)
    end
    dUnits = join(dUnits)

    xUnits  = Array{Char, 1}(undef, 16)
    for aLoop = 1:16
        xUnits[aLoop] = read(A, Char)
    end
    xUnits = join(xUnits)

    fsValid = read(A, Int16)
    whpad3 = read(A, Int16)
    read(A, Float64)
    read(A, Float64)
    for aLoop = 1:26
        read(A, Int32)
    end
    for aLoop = 1:3
        read(A, Int16)
    end
    read(A, Char)
    read(A, Char)
    read(A, Int32)
    read(A, Int16)
    read(A, Int16)
    read(A, Int32)
    modDate = read(A, Int32)
    creationDate = read(A, Int32)
    wdata = Array{Float32, 1}(undef, npnts); read!(A,wdata)

    close(A)

    x0 = Float64(sfB[1]*fsValid)
    dx = Float64(sfA[1]*fsValid)
    return wdata#, npnts# , dUnits, x0, dx, xUnits

end



#data, dUnits, x0, dx, xUnits = ReadIBW(filename)
"""
IBWtoTXT(filename)

Takes a filename, calls ReadIBW and formats the returned data into x and y signals where
x is the position of the sensor and y is the deflection of the sensor.
"""
function IBWtoTXT(filename::String)
    data = readIBW(filename)
    lengthOfData = length(data) ÷ 3

    y = data[lengthOfData+1:2*lengthOfData]
    x = data[2*lengthOfData+1:3*lengthOfData]
    return [x y]
end


function offsetAndDriftCompensation(xy::Matrix{Float32}, endRangeBaseFit=25, noiseMultiplier=4.0)
    #endRangeBaseFit = 25
    breakForContact = true
    aLoop = 0
    basefit = 0

    while breakForContact
        aLoop += 1

        sIdx = min(10+endRangeBaseFit*(aLoop-1), size(xy,1)-endRangeBaseFit);
        eIdx = min(10+endRangeBaseFit*(aLoop), size(xy,1)-1);

        dPdt = 100*(xy[eIdx,2] - xy[sIdx,2]) / endRangeBaseFit;

        if dPdt > 3.3  && aLoop > 3
            basefit = endRangeBaseFit*(aLoop-1);
            breakForContact = false;
        elseif endRangeBaseFit*aLoop > 1000
            basefit = endRangeBaseFit*7;
            breakForContact = false;
        end
    end
   
    zeroLineM = ones(length(xy[10:basefit,1])) \ xy[10:basefit,2]
    xy[:,2] .-= zeroLineM

    # Locate the point where the ramp begins.
    noise = max(0.1, std(xy[10:basefit,2]) );
    rampStartIdx = max(1, findfirst(x -> x > noiseMultiplier*noise, xy[:,2]) - 1);
    xZero = xy[rampStartIdx,1];
    xy[:,1] .-= xZero;                                       # Shift the curve by the value at the start of the ramp.
    xy[:,2] .-= xy[rampStartIdx,2];
    
    zeroLineD = xy[10:basefit,1] \ xy[10:basefit,2];

    zeroLine = [zeroLineD zeroLineM];

    return xy , basefit , zeroLine , rampStartIdx, xZero   
end


"""
subdirImport(placeToLook::String,stringToMatch::String) imports the file names matching some criterion in a 
target directory. 
"""
function subdirImport(placeToLook::String,stringToMatch::String)
    return filter(x-> contains(x, stringToMatch) , readdir(placeToLook))
end

"""
findStartOfHold(xy::Matrix{Float32}, directionOfSearch::String) takes as input a 2xN vector
and finds the first/last entry larger than a specific value.
"""
function findStartOfHold(xy::Matrix{Float32}, directionOfSearch::String)
    # 1. Determine range of deflection values.
    # 2. Bin the values (heuristic bin size at the moment)
    # 3. Determine the most common deflection value at the highest load levels under the assumption that this bin will contain 
    #    the hold sequence.
    # 4. Determine the mean value of all values larger than this bin value.
    # 5. Find the first time the vector exceeds this value.
    # 6. This is taken as the first value in the hold sequence.
    sensorRange = ( maximum(xy[:,2]) - minimum(xy[:,2]) ) .*10e3./maximum(xy[:,2])
    vecLengthTemp = Int64(round(4.0*sensorRange))
    edgesOfHist = range(minimum(xy[:,2]), maximum(xy[:,2]), length = vecLengthTemp)
    histTemp = histcounts(xy[:,2] , edgesOfHist)
    
    tailVecStart = Int64(round(0.9*vecLengthTemp))

    peakIdx = argmax(histTemp[tailVecStart:end])
    peakIdx += tailVecStart-1
    
    idxTemp = map(x -> x > edgesOfHist[peakIdx], xy[:,2])
    
    meanOfPlateau = mean(xy[idxTemp,2])
    stdOfPlateau = std(xy[idxTemp,2])

    if cmp(directionOfSearch,"first") == 0
        returnIdx = findfirst(x -> x ≥ meanOfPlateau-stdOfPlateau, xy[:,2])
    elseif cmp(directionOfSearch,"last") == 0
        returnIdx = findlast(x -> x ≥ meanOfPlateau-stdOfPlateau, xy[:,2])
    end

    return returnIdx
end


function determineThermalCreep(xy::Matrix{Float32},sampleRate::Int64,thermalHoldTime::Int64,ctrl::control,noiseMultiplier::Float64)
    sensorRange = maximum(xy[:,2]) - minimum(xy[:,2])
    vecLengthTemp = Int64(max(100,round(0.1*sensorRange)))
    edgesOfHist = range(minimum(xy[:,2]), maximum(xy[:,2]), length = vecLengthTemp)
    peakIdx = argmax(histcounts(xy[:,2] , edgesOfHist))
    
    idxTemp = (xy[:,2] .< edgesOfHist[peakIdx+1]) .& (xy[:,2] .> edgesOfHist[max(1,peakIdx-1)])
    meanOfPlateau = mean(xy[idxTemp,2])
    stdOfPlateau = std(xy[idxTemp,2])
    
    thermalHoldStartIdx = findlast(x -> x > meanOfPlateau+noiseMultiplier*stdOfPlateau, xy[:,2])
    thermalHoldEndIdx = findlast(x -> x > meanOfPlateau-noiseMultiplier*stdOfPlateau, xy[:,2])

    thermalHoldStartIdx += Int64(round(sampleRate*0.5))
    thermalHoldEndIdx -= Int64(round(sampleRate*0.5))
    
    # Ensure that the thermal hold sequence contains at least 25 seconds (out of the 30 secounds
    # specified).
    if thermalHoldStartIdx > thermalHoldEndIdx
        ctrl.verboseMode && println("Missed thermal hold. Increasing search range.")
        
        while thermalHoldTime+thermalHoldStartIdx > thermalHoldEndIdx
            noiseMultiplier += 1.0
            thermalHoldStartIdx = findlast(x -> x > meanOfPlateau+noiseMultiplier*stdOfPlateau, xy[:,2])
            thermalHoldEndIdx = findlast(x -> x > meanOfPlateau-noiseMultiplier*stdOfPlateau, xy[:,2])

            try
                thermalHoldStartIdx += sampleRate
                thermalHoldEndIdx -= sampleRate
            catch
                println("Failure in the termal hold calculation")
                return 0.0 , length(xy[:,2])
            end
        end
        ctrl.verboseMode && println("Thermal hold found using multiplier $noiseMultiplier")
    end 
    
    # Fit a function of displacement (due to thermal fluctuation)
    # h_thermal(time) = A1 + A2*time^A3 
    thermalHoldDisplacement = xy[thermalHoldStartIdx:thermalHoldEndIdx,1];
    thermalHoldTime = collect(1:length(thermalHoldDisplacement))./sampleRate;
    
    deltaDisp = (xy[thermalHoldEndIdx,1] - xy[thermalHoldStartIdx,1])

    thermalCreepFun(x) = xy[thermalHoldStartIdx,1] .+ x[1].*(thermalHoldTime)#.^x[2]
    thermalHoldMinFun(x) = sqrt(sum(( (thermalCreepFun(x) .- thermalHoldDisplacement)./thermalHoldDisplacement ).^2))

    result = optimize(thermalHoldMinFun, [deltaDisp], BFGS(),Optim.Options(time_limit = 5.0))
    thermal_p = result.minimizer

    # Estimate the thermal drift rate by taking the median of the differentiated h_thermal
    # dhtdt = d(h_thermal(time))/d(time)
    # The functional form accounts for any viscous effects lingering from the unloading at the start
    # of the thermal hold, while the median provides a roboust average of the thermal drift rate.
    dhtdt = thermal_p[1]
    return dhtdt , thermalHoldStartIdx
end


function determineCreepDuringHold(xy_hold,sampleRate::Int64)
    # Fit the creep during the hold time at maximum load to a linear spring-dashpot model
    #
    # Equation (17c) - not labeled in [1]
    # 
    # h(t) = h_i + \beta * t.^(1/3)
    #
    #   h(t)    - Displacement as a function of time during hold sequence.
    #   h_i     - Fitting constant
    #   \beta   - Fitting constant
    #   t       - Time

    holdTimeVals = collect(1:length(xy_hold[:,1])) ./ sampleRate
    # Generate time signal

    deltaDisp = (xy_hold[end,1] - xy_hold[1,1])

    # Define fitting functions
    hOftFun(x) = xy_hold[1,1] .+ x[1].*holdTimeVals.^x[2]
    minFcn(x) = sqrt( sum( ( (hOftFun(x) - xy_hold[:,1])./ xy_hold[:,1] ).^2 ) )
    holdDrift = optimize(minFcn, [deltaDisp 0.33], BFGS(), Optim.Options(time_limit = 3.0))
    crp_p = holdDrift.minimizer

    h_dot_tot = crp_p[1] * crp_p[2]*holdTimeVals[end]^(crp_p[2] - 1.0)

    return max(h_dot_tot, 0.0)
end



function calculateMachineCompliance(indentationSet::metaInfoExperimentalSeries,hyperParameters,ctrl::control)
    resultNames = subdirImport(indentationSet.targetDir,".ibw")
    resultNames = resultNames[6:end]
    collectCompliances = []
    collectAreas = []
    for file in resultNames
        ctrl.verboseMode && println(file)

        try
            currentCompliance , currentArea = extractSingleComplianceExperiment(indentationSet,hyperParameters,ctrl,file)
            if !isnan(currentCompliance) && !isnan(currentArea)
                push!(collectCompliances,currentCompliance)
                push!(collectAreas,currentArea)
            end
        catch
            currentCompliance = NaN
            currentArea = NaN
        end
    end

    collectAreas = Float32.(collectAreas)
    collectCompliances = Float32.(collectCompliances)
    #println(collectAreas)
    #println(collectCompliances)
    squaredInverseArea = 1.0 ./sqrt.(collectAreas)
    
    if length(squaredInverseArea) > 0
        effectiveCompliance = [squaredInverseArea[:] ones(length(squaredInverseArea))] \ collectCompliances
        println(effectiveCompliance) 
    else
        println("----> No accepted measurements")
    end

    return [squaredInverseArea[:]  collectCompliances[:]] # effectiveCompliance[2]
end


function extractSingleComplianceExperiment(indentationSet::metaInfoExperimentalSeries,hyperParameters,ctrl::control,resultFile::String)
    
    xy = IBWtoTXT(indentationSet.targetDir*resultFile)
    # Import signal
    xy .*= 1e9     
    # Convert to nano-meters

    xy , basefit , zeroLine , rampStartIdx, xZero  = offsetAndDriftCompensation(xy)
    # Find initial contact

    xy[:,1] .-= xy[:,2]
    xy[:,2] .*= indentationSet.springConstant
    # Convert displacement-deflection matrix to indentation-force matrix

    # Determine the start of the hold time at circa max force.
    if cmp( lowercase(hyperParameters.controlLoop) , "force") == 0
        holdStartIdx = findStartOfHold(xy,"first")
    elseif cmp( lowercase(hyperParameters.controlLoop) , "displacement") == 0
        holdStartIdx = argmax(xy[:,2])
    else
        throw(DomainError(hyperParameters.controlLoop, "controlLoop setting not defined."))
    end

    # Split into loading and unloading.
    xy_unld1 = xy[holdStartIdx:end,:];

    #Determine the end of the hold time.
    if cmp( lowercase(hyperParameters.controlLoop) , "force") == 0
        unloadStartIdx = findStartOfHold(xy_unld1,"last")
    elseif cmp( lowercase(hyperParameters.controlLoop) , "displacement") == 0
        unloadStartIdx = max(argmax(xy_unld1[:,1]) , findStartOfHold(xy_unld1,"last"))
    else
        throw(DomainError(hyperParameters.controlLoop, "controlLoop setting not defined."))
    end

    if ctrl.plotMode
        display(plot!([xy_unld1[unloadStartIdx,1]], [xy_unld1[unloadStartIdx,2]],seriestype = :scatter, label = :none))
        plotd = plot(xy[:,1], xy[:,2], xlims = (0.0, maximum(xy[:,1])), xlab = "Indentation [nm]", ylab = "Force [uN]", label = "Signal")
        plot!([xy[holdStartIdx,1]], [xy[holdStartIdx,2]], 
                             seriestype = :scatter, lab = "Start of hold", legend = :topleft)
        plot!([xy_unld1[unloadStartIdx,1]], [xy_unld1[unloadStartIdx,2]],seriestype = :scatter, label = "Start of unload")
        println("$(indentationSet.targetDir)$(resultFile[1:end-4]).png")
        plot!(size=(500,500))
        savefig(plotd,"$(indentationSet.targetDir)$(resultFile[1:end-4]).png")
    end

    # Split into two new pieces
    xy_hold = xy_unld1[1:unloadStartIdx-1,:]
    xy_unld = xy_unld1[unloadStartIdx:end,:]

    dhtdt , thermalHoldStartIdx = determineThermalCreep(xy,hyperParameters.sampleRate,indentationSet.thermalHoldTime,ctrl,hyperParameters.noiseMultiplier)
    if thermalHoldStartIdx > length(xy_unld[:,1])
        println("Could not find a thermal hold period. Assuming no thermal hold.")
        thermalHoldStartIdx = unloadStartIdx + hyperParameters.sampleRate*2
    end
    xy_unld5 = xy_unld[1:thermalHoldStartIdx,:]

    #xy_unld5 = xy_unld[1:min(Int64(round(2000*0.90)),size(xy_unld,1)),:];  # OBS 2000 is hard coded!

    # Fitting of the unloading curve.
    stiffness_fit = Array{Float64}(undef,1)    
    tempLen = minimum([hyperParameters.unloadingFitRange, length(xy_unld5[:,1])])

    dispVals = xy_unld5[1:tempLen ,1]
    forceVals = xy_unld5[1:tempLen,2]
    Fmax = xy_unld5[1,2]               # Maximum force during unloading

    if cmp(hyperParameters.unloadingFitFunction,"Oliver-Pharr") == 0
        Dmax = xy_unld5[1,1]               
        # Maximum indentation depth during unloading
        
        function unloadFitFun(fitCoefs)
            return fitCoefs[1].*(dispVals .- fitCoefs[2]).^fitCoefs[3] .- forceVals
        end
        function unloadFitMinFun(fitCoefs)
            sqrt( sum( (unloadFitFun(fitCoefs) ./ forceVals).^2 ) )
        end

        lx = [-Inf, -Inf , 0.0]; ux = [Inf, minimum(dispVals)-1e-2 , Inf];
        dfc = TwiceDifferentiableConstraints(lx, ux)
        resultFit = optimize(unloadFitMinFun, dfc, [1.0, 1.0, 1.0], IPNewton())
        uld_p = resultFit.minimizer
        println(uld_p)
        stiffness_fit = uld_p[1]*uld_p[3]*(Dmax - uld_p[2]).^(uld_p[3] - 1.0)
    
        if ctrl.plotMode
            plot(xlabel = "Indentation [nm]", ylabel = "Force [uN]", size = (500 , 500) , dpi = 600, legend = :topleft)
            plot!(dispVals, forceVals, label = "Signal")
            plot!(dispVals, unloadFitFun(uld_p).+forceVals , label = "Fit")# P(h) = $(round(uld_p[1],digits=2)(h - $(round(uld_p[2], digits = 2)) )^{$(round(uld_p[3],digits = 2))} )")
            savefig("$(indentationSet.targetDir)$(resultFile[1:end-4])_unloadFit.png")
        end

    elseif cmp(hyperParameters.unloadingFitFunction, "AP-OP") == 0
        Dmax = xy_unld5[1,1]               
        # Maximum indentation depth during unloading
        
        function unloadFitFunAP(fitCoefs)
            return Fmax.*((dispVals .- fitCoefs[1])./(Dmax .- fitCoefs[1]) ).^( fitCoefs[2] ) .- forceVals
        end
        function unloadFitMinFunAP(fitCoefs)
            sqrt( sum( (unloadFitFunAP(fitCoefs) ./ forceVals).^2 ) )
        end

        resultFit = optimize(unloadFitMinFunAP, [ Dmax.*0.5 , 2.0], NewtonTrustRegion())
        uld_p = Optim.minimizer(resultFit)
        stiffness_fit = Fmax.*uld_p[2]*(Dmax .- uld_p[1]).^(-1.0)


        if ctrl.plotMode && uld_p[1] > 0.0 && uld_p[2] > 0.0
            plotd = plot(dispVals, forceVals, xlabel = "Indentation [nm]" , ylabel = "Force [uN]" , label = "Signal")
            plot!(dispVals, unloadFitFunAP(uld_p).+forceVals , label = "Fit")# \$F(z)=F_{max}((z - $(round(uld_p[1],digits = 1)))/(D_{max} - $(round(uld_p[1],digits = 1))) )^{$(round(uld_p[2],digits = 1))} \$", legend = :topleft)
            plot!(size=(500,500))
            println("$(indentationSet.targetDir)$(resultFile[1:end-4])_unloadFit.png")
            savefig(plotd,"$(indentationSet.targetDir)$(resultFile[1:end-4])_unloadFit.png")
        end
    end
    stiffness = stiffness_fit
    dhtdt = 0;

    maxIndentation = median(xy_unld5[1,1]) - dhtdt*(length(xy[rampStartIdx:holdStartIdx,1])+length(xy_hold[:,1]))/hyperParameters.sampleRate;  #%OBS OBS OBS

    if cmp(indentationSet.indenterType,"pyramid") == 0
        x0 = maxIndentation - 0.72*Fmax/stiffness;
    elseif cmp(indentationSet.indenterType,"hemisphere") == 0
        x0 = maxIndentation - 0.75*Fmax/stiffness;
    end

    area_xy = readdlm(indentationSet.areaFile, ' ', Float64, '\n')
    # % Determine the area by loading the calibration data and fitting a polynom to the data.        

    if (x0 > 100.0)
        area_fit_end = length(area_xy[:,1])
    elseif (x0 < 100.0)
        area_fit_end = findfirst( x -> x > 100,area_xy[:,1])
    end
    
    tempVec = area_xy[1:area_fit_end,1]
    p_area = [tempVec.^2 tempVec tempVec.^0.5 tempVec.^0.25 tempVec.^0.125] \ area_xy[1:area_fit_end,2]
    
    unloadArea = [x0^2 x0 x0^0.5 x0^0.25 x0^0.125] * p_area
    unloadArea = unloadArea[1]
    

    if stiffness < 0.0
        stop
    end

    if uld_p[1] < 0.0 || uld_p[2] < 1.3 || uld_p[2] > 2.5
        return NaN, NaN
    else
        return 1/stiffness , unloadArea
    end
    # Assign outputs
end

function importNI_forceDisplacementData(filename::String)   
    df = CSV.File(filename, decimal = ',', skipto = 72, delim = "\t", header = ["Time", "Pd_mm", "Fn_mN","FnRef_mN","SegmentID"]) |> DataFrame
    dx = Matrix(df)

    return dx[:,2:3]
end
function importNI_forceDisplacementData_v2(filename::String)   
    df = CSV.File(filename, decimal = '.', skipto = 7, delim = "\t", header = ["Depth_nm", "Load_uN", "Time_s","Depth_V","Load_V"]) |> DataFrame
    dx = Matrix(df)

    return dx[:,1:2]
end

function areaCheck(indentationSet , ctrl)
    area_xy = readdlm(indentationSet.areaFile, ' ', Float64, '\n')
    # % Determine the area by loading the calibration data and fitting a polynom to the data.        
  
    tempVec = area_xy[:,1]
    p_area = [tempVec.^2 tempVec tempVec.^0.5 tempVec.^0.25 tempVec.^0.125] \ area_xy[:,2]

    area_coneIndenter(indentationDepth) = 24.5.*indentationDepth.^2.0
    area_halfSphere(indentationDepth) = π.*(2.0.*indentationDepth.*300.0 .- indentationDepth.^2.0)
    area_halfSphere2(indentationDepth) = π.*(300.0 .*indentationDepth)

    if ctrl.plotMode
        plot(xlabel = "\$z\$ [nm]", ylabel = "Area \$A(z)\$ [nm]^2", size = (500,500), dpi = 600)
        plot!(area_xy[:,1] , area_xy[:,2], label = split(indentationSet.areaFile, '/')[end])
        plot!(area_xy[:,1] , area_coneIndenter(area_xy[:,1]), label = "Vickers/Berkovich")
        plot!(area_xy[:,1] , area_halfSphere(area_xy[:,1]), label = "Halfsphere, R = 300 nm", legend = :bottomright)
        plot!(area_xy[:,1] , area_halfSphere2(area_xy[:,1]), label = "Halfsphere contact, R = 300 nm", legend = :bottomright)
        ylims!(0.0 , maximum(area_xy[:,2]))
        savefig("$(indentationSet.areaFile[1:end-4]).png")
    end
end

function plotVersusTime(xy, sampleRate, saveName)
    plot(xlabel = "Time [s]", ylabel = "Force [nN]", size = (500,500) , dpi = 600)
    plot!( (collect(1:length(xy[:,2])).-1)./sampleRate , xy[:,2] , label = "Force signal")
    savefig("$(saveName).png")
end



function signalImporter( indentationSet, ctrl, resultFile, hyperParameters)
    if cmp(lowercase(indentationSet.indentationDataType), "afm") == 0
        xy = IBWtoTXT(indentationSet.targetDir*resultFile)
        # Import signal
        xy .*= 1e9     
        # Convert to nano-meters
        ctrl.plotMode && display(plot([xy[1:100:end,1]],[xy[1:100:end,2]]))

        xy , ~ , ~ , rampStartIdx, ~  = offsetAndDriftCompensation(xy)
        # Find initial contact
        ctrl.plotMode && display(plot!(xy[1:100:end,1],xy[1:100:end,2]))

        xy[:,1] .-= xy[:,2]
        xy[:,2] .*= indentationSet.springConstant
        xy[:,1] .-= hyperParameters.machineCompliance.*xy[:,2];
        # Convert displacement-deflection matrix to indentation-force matrix

    elseif cmp( lowercase( indentationSet.indentationDataType), "ni") == 0
        xy = importNI_forceDisplacementData(indentationSet.targetDir*resultFile)   
        # Import data
        xy[:,2] *= 1.0e6
        # Convert force to nano-Newtons
        
        xy = Float32.(xy)
        # Convert to Float32
        
        rampStartIdx = 1
        # Software handles rampStart, so set to 1.

        ctrl.plotMode && display(plot([xy[:,1]],[xy[:,2]]))

    elseif cmp( lowercase( indentationSet.indentationDataType), "ni_v2") == 0
        xy = importNI_forceDisplacementData_v2(indentationSet.targetDir*resultFile)   
        # Import data
        xy[:,2] *= 1.0e3
        # Convert force to nano-Newtons
        
        xy = Float32.(xy)
        # Convert to Float32

        xy , ~ , ~ , rampStartIdx, ~  = offsetAndDriftCompensation(xy)
        # Find initial contact
        
        ctrl.plotMode && display(plot([xy[:,1]],[xy[:,2]]))

    else 
        
    end
    return xy , rampStartIdx

end