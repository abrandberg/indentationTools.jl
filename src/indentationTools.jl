################################################################################
#
# indentationTools.jl
# 
# Created by : August Brandberg
# Date       : 2022-05-22
#
# Contents of this file:
#
# SECTION       CONTENTS
#       1       Struct definitions.
#       2       Import and preprocess data functions
#       3       Main functions
module indentationTools

################################################################################
## SECTION 1
## Struct definitions
"""
hyperParameters(sampleRate::Int, unloadingFitRange::Int, unloadingFitFunction::String, constrainHead::Int, constrainTail::Int, machineCompliance::Float64, noiseMultiplier::Float64)

Sets instructions for how the files should be treated. The assumption is that all files in an indentation set share sample rate and that all files should be treated equally.

- sampleRate : [Hz] {> 0} Force-indentation data sample rate.
- unloadingFitRange : [-] {> 1} Number of samples to be included in the fit.
- unloadingFitFunction : [] Should be "Oliver-Pharr", "AP-OP" or "Feng".
- compensateCreep : [] Whether to attempt to use Ngan's method to compensate for creep during the indentation.
- constrainHead : Experimental, not implemented.
- constrainTail : Experimental, not implemented.
- machineCompliance : [m/N] {≥ 0.0} Machine compliance.
- noiseMultiplier : [-] {> 0} Number of signal standard deviations to add per loop when the thermal hold is missed.

"""
struct hyperParameters
    sampleRate              ::Int            # [Hz] Sample rate.
    unloadingFitRange       ::Int            # [-] Vector samples lengths (measured from start of unloading) to be included in unload fit.
    unloadingFitFunction    ::String         # [String] Function to use when fitting.
    compensateCreep         ::Bool           # [Bool] Compensate creep using Feng's method, y/n
    constrainHead           ::Bool           # [Bool] Experimental, not implemented
    constrainTail           ::Bool           # [Bool] Experimental, not implemented
    machineCompliance       ::Float64        # [m/N] Machine compliance
    noiseMultiplier         ::Float64        # [-] Number of standard deviations to add when thermal hold is missed.
    controlLoop             ::String         # Control for the indentation. Can be "force" or "nanoindentation"

    # Input checks
    hyperParameters(sampleRate, unloadingFitRange, unloadingFitFunction, 
                    compensateCreep, constrainHead, constrainTail, 
                    machineCompliance, noiseMultiplier , controlLoop) = 
                                                        sampleRate ≤ 0        ? throw(DomainError(sampleRate, "sampleRate must be a positive integer."))           : 
                                                        unloadingFitRange ≤ 1 ? throw(DomainError(unloadingFitRange , "unloadingFitRange must be larger than 1.")) :  
                                                        unloadingFitFunction ∉ ["Oliver-Pharr" , "AP-OP" , "Feng"] ? throw(DomainError(unloadingFitFunction , "unloadingFitFunction $unloadingFitFunction is not implemented.")) :                
                                                        machineCompliance < 0.0 ?  throw(DomainError(machineCompliance , "machineCompliance must be ≥ 0.0."))      :
                                                        noiseMultiplier ≤ 0.0 ?  throw(DomainError(noiseMultiplier , "noiseMultiplier must be > 0.0."))            :
                                                        controlLoop ∉ ["force" , "displacement"] ?  throw(DomainError(controlLoop , "controlLoop $controlLoop is not implemented.")) :
                    new(sampleRate, unloadingFitRange, unloadingFitFunction, compensateCreep, constrainHead, constrainTail, machineCompliance, noiseMultiplier, controlLoop)   
end


"""
control(plotMode::Bool,verboseMode::Bool) 

controls what outputs are generated when code from the package indentationTools is run.

- plotMode    - If true, generate and save plots during the processing of files.
- verboseMode - If true, report back information to the terminal.

Example:
```
control(true, false)        # Generate plots, but no verbose output.
```
"""
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
## SECTION 2
## Import and preprocess data functions
include("preprocessingFunctions.jl")

################################################################################
## SECTION 3
## Main functions
function modulusfitter(indentationSet::metaInfoExperimentalSeries,hyperParameters,ctrl::control,resultFile::String)


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
        xy = Float32.(xy)
        # Convert to Float32
        xy[:,2] *= 1.0e6
        # Convert force to nano-Newtons
        
        rampStartIdx = 1
        # Software handles rampStart, so set to 1.

        ctrl.plotMode && display(plot([xy[:,1]],[xy[:,2]]))
    end

    # Determine the start of the hold time at circa max force.
    if cmp( lowercase(hyperParameters.controlLoop) , "force") == 0
        holdStartIdx = findStartOfHold(xy,"first")
    elseif cmp( lowercase(hyperParameters.controlLoop) , "displacement") == 0
        holdStartIdx = argmax(xy[:,2])
    else
        throw(DomainError(hyperParameters.controlLoop, "controlLoop setting not defined."))
    end
    
    
    # Split into loading and unloading.
    xy_unld1 = xy[holdStartIdx:end,:]

    #Determine the end of the hold time.
    if cmp( lowercase(hyperParameters.controlLoop) , "force") == 0
        unloadStartIdx = findStartOfHold(xy_unld1,"last")
    elseif cmp( lowercase(hyperParameters.controlLoop) , "displacement") == 0
        unloadStartIdx = max(argmax(xy_unld1[:,1]) , findStartOfHold(xy_unld1,"last"))
    else
        throw(DomainError(hyperParameters.controlLoop, "controlLoop setting not defined."))
    end
    
    # Split into two new pieces
    xy_hold = xy_unld1[1:unloadStartIdx-1,:]
    xy_unld = xy_unld1[unloadStartIdx:end,:]

    if ctrl.plotMode
        plotd = plot(xy[:,1], xy[:,2], xlims = (0.0, maximum(xy[:,1])), xlab = "Indentation [nm]", ylab = "Force [uN]", label = "Signal")
        plot!([xy[holdStartIdx,1]], [xy[holdStartIdx,2]], seriestype = :scatter, lab = "Start of hold")
        plot!([xy_unld1[unloadStartIdx,1]], [xy_unld1[unloadStartIdx,2]], seriestype = :scatter, lab = "Start of unload", legend = :topleft)
        plot!(size=(500,500))
        savefig(plotd,"$(indentationSet.targetDir)$(resultFile[1:end-4])_signal.png")
    end



    # Accept only indentations that had positive creep. "Negative" creep (indenter moves outwards 
    # during hold sequence) can occur if the thermal drift is substantial, but this typically
    # indicates that the system was not in equilibrium (since the thermal drift dominates the creep)
    # and furthermore it messes up the mathematical framework if you accept such indentations (see
    # Cheng & Cheng articles.)
    condition1 = xy[holdStartIdx,1] < xy_unld1[unloadStartIdx,1] 
        
    if condition1

        ctrl.plotMode && display(plot(xy_unld[:,1], xy_unld[:,2]))
        ctrl.verboseMode && println(length(xy_unld[:,1]))
        # Make sure that the thermal hold sequence is not included in the unloading curve.

        if cmp(indentationSet.indentationDataType, "afm") == 0
            dhtdt , thermalHoldStartIdx = determineThermalCreep(xy,hyperParameters.sampleRate,indentationSet.thermalHoldTime,ctrl,hyperParameters.noiseMultiplier)

            if thermalHoldStartIdx > length(xy_unld[:,1])
                println("Could not find a thermal hold period. Assuming no thermal hold.")
                thermalHoldStartIdx = unloadStartIdx + hyperParameters.sampleRate*2
            end
        else
            dhtdt = 0.0
            thermalHoldStartIdx = length(xy_unld[:,1])
        end
        xy_unld5 = xy_unld[1:thermalHoldStartIdx,:]
        

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

            lx = [0.0, 0.0 , 0.0]; ux = [Inf, minimum(dispVals)-1e-2 , Inf];
            dfc = TwiceDifferentiableConstraints(lx, ux)
            resultFit = optimize(unloadFitMinFun, dfc, [1.0, 1.0, 1.0], IPNewton())
            uld_p = resultFit.minimizer
            stiffness_fit = uld_p[1]*uld_p[3]*(Dmax - uld_p[2]).^(uld_p[3] - 1.0)

            if ctrl.plotMode #&& uld_p[1] > 0.0 && uld_p[2] > 0.0
                plotd = plot(dispVals, forceVals, xlabel = "Indentation [nm]" , ylabel = "Force [uN]" , label = "Signal")
                plot!(dispVals, unloadFitFun(uld_p).+forceVals , label = "Fit \$F(z)= $(round(uld_p[1],digits = 1))(z - $(round(uld_p[2],digits = 1)))^{$(round(uld_p[3],digits = 1))} \$", legend = :topleft)
                plot!(size=(500,500))
                println("$(indentationSet.targetDir)$(resultFile[1:end-4])_unloadFit.png")
                savefig(plotd,"$(indentationSet.targetDir)$(resultFile[1:end-4])_unloadFit.png")
            end

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


        elseif cmp(lowercase(hyperParameters.unloadingFitFunction), "feng") == 0
            
            unloadFitFun2(fitCoefs) = fitCoefs[1] .+ fitCoefs[2].*forceVals.^0.5 + fitCoefs[3].*forceVals.^fitCoefs[4] .- dispVals
            unloadFitMinFun2(fitCoefs) = sqrt(sum( (unloadFitFun2(fitCoefs) ./ dispVals).^2) )
            resultFit = optimize(unloadFitMinFun2, [1.0 1.0 1.0 1.0], BFGS())
            uld_p = resultFit.minimizer
            stiffness_fit = inv(( 0.5*uld_p[2].*Fmax.^-0.5 + uld_p[4]*uld_p[3].*Fmax.^(uld_p[4] - 1.0) ))
           
        end

        h_dot_tot = determineCreepDuringHold(xy_hold,hyperParameters.sampleRate)       
        dPdt = [1/hyperParameters.sampleRate .* collect(0:(length(xy_unld5[:,1])-1)) ones(length(xy_unld5[:,1]))] \ xy_unld5[:,2]

        # Apply creep compensation
        hyperParameters.compensateCreep ? stiffness = inv(1/stiffness_fit + h_dot_tot/(abs(dPdt[1]))) : stiffness = stiffness_fit

        # Equation (2) in [1]
        maxIndentation = xy_unld5[1,1] - dhtdt*(length(xy[rampStartIdx:holdStartIdx,1])+length(xy_hold[:,1]))/hyperParameters.sampleRate

        if cmp(lowercase(indentationSet.indenterType),"pyramid") == 0
            x0 = maxIndentation - 0.72*Fmax/stiffness;
        elseif cmp(lowercase(indentationSet.indenterType),"hemisphere") == 0
            x0 = maxIndentation - 0.75*Fmax/stiffness;
        end
        x0 < 0.0 && return 0.0 , maxIndentation , x0 , 0.0 , stiffness , uld_p[3]


        if cmp(lowercase(indentationSet.areaFile), "vickers") == 0 || cmp(lowercase(indentationSet.areaFile), "berkovich") == 0
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

        #elseif cmp(lowercase(indentationSet.areaFile), "halfsphere") == 0 || cmp(lowercase(indentationSet.areaFile), "sphere") == 0
            # Needs radius of the indenter to work!

            # area_xy(indentationDepth) = π.*(2.0.*indentationDepth.*300.0 .- indentationDepth.^2.0)
            # # Define function
            # unloadArea = area_xy(x0)
            # # Extract value

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

        if cmp( lowercase( indentationSet.indenterType),"pyramid") == 0
            Er = Er/1.05;
        end
        return Er , maxIndentation , x0 ,  unloadArea , stiffness  , uld_p[2]
    else
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
