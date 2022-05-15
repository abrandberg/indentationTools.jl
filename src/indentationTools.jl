module indentationTools

################################################################################
## Structs used
struct hyperParameters
    sampleRate              ::Int64          # [Hz] Sample rate.
    unloadingFitRange       ::Int64          # [-] Vector samples lengths (measured from start of unloading) to be included in unload fit.
    unloadingFitFunction    ::String         # [string] Function to use when fitting.
    compensateCreep         ::Bool           # [Bool] Compensate creep using Feng's method, y/n
    constrainHead           ::Int            # Experimental, not implemented
    constrainTail           ::Int            # Experimental, not implemented
    machineCompliance       ::Float64        # Machine compliance
    noiseMultiplier         ::Float64        # Number of standard deviations to add when thermal hold is missed.
end

struct control
    plotMode                ::Bool           # Activates plotting of intermediate results
    verboseMode             ::Bool           # Verbose output
end

struct metaInfoExperimentalSeries
    designatedName          ::String         # Designated name
    relativeHumidity        ::Float64        # Relative humidity
    indenterType            ::String
    indentationNormal       ::String
    springConstant          ::Float64
    areaFile                ::String
    targetDir               ::String
    thermalHoldTime         ::Int64          # Should be changed to something of unit "second".
    indentationDataType     ::String         # "afm" for classical CC or "ni" for force-displacement
end

################################################################################
## Import and preprocess data functions
include("preprocessingFunctions.jl")

function modulusfitter(indentationSet::metaInfoExperimentalSeries,hyperParameters,ctrl::control,resultFile::String)


    if cmp(indentationSet.indentationDataType, "afm") == 0
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

    elseif cmp(indentationSet.indentationDataType, "ni") == 0
        xy = importNI_forceDisplacementData(indentationSet.targetDir*resultFile)   
        # Import data
        xy = Float32.(xy)
        # Convert to Float32
        xy[:,2] *= 1.0e6
        # Convert force to nano-Newtons
        
        rampStartIdx = 1
        # Software handles rampStart, so set to 1.

        ctrl.plotMode && display(plot([xy[:,1]],[xy[:,2]]))
    end

    
        
    #################################################################
    # Determine the start of the hold time at circa max force.
    # 
    # 1. Determine range of deflection values.
    # 2. Bin the values (heuristic bin size at the moment)
    # 3. Determine the most common deflection value at the highest load levels under the assumption that this bin will contain 
    #    the hold sequence.
    # 4. Determine the mean value of all values larger than this bin value.
    # 5. Find the first time the vector exceeds this value.
    # 6. This is taken as the first value in the hold sequence.
    holdStartIdx = findStartOfHold(xy,"first")
    ctrl.plotMode && display(plot(xy[:,1], xy[:,2], xlims = (0.0, maximum(xy[:,1])), xlab = "Indentation [nm]", ylab = "Force [uN]", legend = false))
    ctrl.plotMode && display(plot!([xy[holdStartIdx,1]], [xy[holdStartIdx,2]], 
                             seriestype = :scatter, lab = "Start of hold", legend = :topleft))

    # Split into loading and unloading.
    xy_unld1 = xy[holdStartIdx:end,:];

    #Determine the end of the hold time.
    unloadStartIdx = findStartOfHold(xy_unld1,"last")
    ctrl.plotMode && display(plot!([xy_unld1[unloadStartIdx,1]], [xy_unld1[unloadStartIdx,2]],seriestype = :scatter))

    # Split into two new pieces
    xy_hold = xy_unld1[1:unloadStartIdx-1,:];
    xy_unld = xy_unld1[unloadStartIdx:end,:];

    # Accept only indentations that had positive creep. "Negative" creep (indenter moves outwards 
    # during hold sequence) can occur if the thermal drift is substantial, but this typically
    # indicates that the system was not in equilibrium (since the thermal drift dominates the creep)
    # and furthermore it messes up the mathematical framework if you accept such indentations (see
    # Cheng & Cheng articles.)
    condition1 = xy[holdStartIdx,1] < xy_unld1[unloadStartIdx,1] 
    ctrl.verboseMode && println("condition1 = $condition1")
    
    # Accept only monotonously increasing load-displacement curves. A curve may show weird behaviour
    # and our solution is to simply drop the curve in that case. 
    condition2 = minimum(xy[rampStartIdx:holdStartIdx,1]) > (0.0-eps())
    condition3 = maximum(xy_unld[:,1]) < (xy_unld[1,1]+0.5)


    #ctrl.plotMode && display(title!(string(condition3)))
    
    if condition1 #&& (length(xy_unld[:,1]) > 285) #&& condition2 && condition3

        ctrl.plotMode && display(plot(xy_unld[:,1], xy_unld[:,2]))
        ctrl.verboseMode && println(length(xy_unld[:,1]))
        # Make sure that the thermal hold sequence is not included in the unloading curve.

        if cmp(indentationSet.indentationDataType, "afm") == 0
            dhtdt , thermalHoldStartIdx = determineThermalCreep(xy,hyperParameters.sampleRate,indentationSet.thermalHoldTime,ctrl,hyperParameters.noiseMultiplier)

            if thermalHoldStartIdx > length(xy_unld[:,1])
                println("Could not find a thermal hold period. Assuming no thermal hold.")
                #println("Check results carefully.")
                #thermalHoldStartIdx = length(xy_unld[:,1])
                thermalHoldStartIdx = unloadStartIdx + hyperParameters.sampleRate*2
            end
        else
            dhtdt = 0.0
            thermalHoldStartIdx = length(xy_unld[:,1])
        end
        #println(thermalHoldStartIdx)
        xy_unld5 = xy_unld[1:thermalHoldStartIdx,:]
        ctrl.plotMode && display(plot!(xy_unld5[:,1], xy_unld5[:,2]))

        # Fitting of the unloading curve.
        stiffness_fit = Array{Float64}(undef,1)    
        tempLen = minimum([hyperParameters.unloadingFitRange, length(xy_unld5[:,1])])
        dispVals = xy_unld5[1:tempLen ,1]
        ctrl.verboseMode && println(length(dispVals))
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

            lx = [0.0, 0.0 , 0.0]; ux = [Inf, minimum(dispVals)-1e-2 , Inf];
            dfc = TwiceDifferentiableConstraints(lx, ux)
            resultFit = optimize(unloadFitMinFun, dfc, [1.0, 1.0, 1.0], IPNewton())
            uld_p = resultFit.minimizer
            stiffness_fit = uld_p[1]*uld_p[3]*(Dmax - uld_p[2]).^(uld_p[3] - 1.0)


        elseif cmp(hyperParameters.unloadingFitFunction, "AP-OP") == 0
            Dmax = xy_unld5[1,1]               
            # Maximum indentation depth during unloading
            
            function unloadFitFunAP(fitCoefs)
                return Fmax.*((dispVals .- fitCoefs[1])./(Dmax .- fitCoefs[1]) ).^( fitCoefs[2] ) .- forceVals
            end
            function unloadFitMinFunAP(fitCoefs)
                sqrt( sum( (unloadFitFunAP(fitCoefs) ./ forceVals).^2.0 ) )
            end

            ctrl.plotMode && display(plot(dispVals, forceVals, label = :none))

            resultFit = optimize(unloadFitMinFunAP, [ Dmax.*0.5 , 2.0], NewtonTrustRegion())
            uld_p = Optim.minimizer(resultFit)
            stiffness_fit = Fmax.*uld_p[2].*(Dmax .- uld_p[1]).^(-1.0)


            if ctrl.plotMode && uld_p[1] > 0.0 && uld_p[2] > 0.0
                plotd = plot(dispVals, forceVals, xlabel = "Indentation [nm]" , ylabel = "Force [uN]" , label = "Signal")
                plot!(dispVals, unloadFitFunAP(uld_p).+forceVals , label = "Fit \$F(z)=F_{max}((z - $(round(uld_p[1],digits = 1)))/(D_{max} - $(round(uld_p[1],digits = 1))) )^{$(round(uld_p[2],digits = 1))} \$", legend = :topleft)
                plot!(size=(500,500))
                println("$(indentationSet.targetDir)$(resultFile[1:end-4])_unloadFit.png")
                savefig(plotd,"$(indentationSet.targetDir)$(resultFile[1:end-4])_unloadFit.png")
            end


        elseif cmp(hyperParameters.unloadingFitFunction, "Feng") == 0
            
            unloadFitFun2(fitCoefs) = fitCoefs[1] .+ fitCoefs[2].*forceVals.^0.5 + fitCoefs[3].*forceVals.^fitCoefs[4] .- dispVals
            unloadFitMinFun2(fitCoefs) = sqrt(sum( (unloadFitFun2(fitCoefs) ./ dispVals).^2) )
            resultFit = optimize(unloadFitMinFun2, [1.0 1.0 1.0 1.0], BFGS())
            uld_p = resultFit.minimizer
            stiffness_fit = inv(( 0.5*uld_p[2].*Fmax.^-0.5 + uld_p[4]*uld_p[3].*Fmax.^(uld_p[4] - 1.0) ))
            
        else
            println("Not implemented.")
            stop
        end

        h_dot_tot = determineCreepDuringHold(xy_hold,hyperParameters.sampleRate)       
        dPdt = [1/hyperParameters.sampleRate .* collect(0:(length(xy_unld5[:,1])-1)) ones(length(xy_unld5[:,1]))] \ xy_unld5[:,2]

        
        if hyperParameters.compensateCreep
            stiffness = inv(1/stiffness_fit + h_dot_tot/(abs(dPdt[1]))); 
            #stiffness = inv(-hyperParameters.machineCompliance + 1/stiffness_fit + h_dot_tot/(abs(dPdt[1]))); 
        else
            stiffness = stiffness_fit;
            #println(-min(1e6,1/hyperParameters.machineCompliance))
            #println(1/stiffness_fit)
            
            #stiffness = inv(-hyperParameters.machineCompliance + 1/stiffness_fit ); 
            #println(stiffness)
        end

        # Equation (2) in [1]
        maxIndentation = xy_unld5[1,1] - dhtdt*(length(xy[rampStartIdx:holdStartIdx,1])+length(xy_hold[:,1]))/hyperParameters.sampleRate;  #%OBS OBS OBS

        if cmp(indentationSet.indenterType,"pyramid") == 0
            x0 = maxIndentation - 0.72*Fmax/stiffness;
        elseif cmp(indentationSet.indenterType,"hemisphere") == 0
            x0 = maxIndentation - 0.75*Fmax/stiffness;
        end
        x0 < 0.0 && return 0.0 , maxIndentation , x0 , 0.0 , stiffness , uld_p[3]


        if cmp(indentationSet.areaFile, "vickers") == 0 || cmp(indentationSet.areaFile, "berkovich") == 0
            area_xy(indentationDepth) = 24.5.*indentationDepth.^2
            # Define function
            unloadArea = area_xy(x0)
            # Extract value

            # N.A. Sakharova, J.V. Fernandes, J.M. Antunes, M.C. Oliveira,
            # Comparison between Berkovich, Vickers and conical indentation tests: 
            # A three-dimensional numerical simulation study,
            # International Journal of Solids and Structures,
            # Volume 46, Issue 5, 2009, Pages 1095-1104,
            # https://doi.org/10.1016/j.ijsolstr.2008.10.032

        else
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
        end
        unloadArea < 0.0 && return 0.0 , maxIndentation , x0 , unloadArea , stiffness , uld_p[3]
        
        # % Equation (1) in [1]
        Er = sqrt(pi)/(2.0)/sqrt(unloadArea) / ( 1.0/stiffness )

        if cmp(indentationSet.indenterType,"pyramid") == 0
            Er = Er/1.05;
        end
        return Er , maxIndentation , x0 ,  unloadArea , stiffness  , uld_p[2]
    else
        #println(condition1)
        #println(condition2)
        #println(condition3)
        return 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0
    end

    ctrl.verboseMode && println(Er)
end

################################################################################
# Export functions
    # Data import
    export readIBW
    export IBWtoTXT
    export importNI_forceDisplacementData
    export subdirImport

    # Preprocessing of signal
    export offsetAndDriftCompensation
    export findStartOfHold
    export determineCreepDuringHold
    export determineThermalCreep

    # Principal functions/functionality
    export modulusfitter
    export control
    export hyperParameters
    export metaInfoExperimentalSeries
    
    # Control
    export calculateMachineCompliance
    export areaCheck
end
