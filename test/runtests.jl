using indentationTools
using Test


include("testFunctions.jl")

#importNI_forceDisplacementData("test/testData/anonNI.TXT")
@testset "importNI_forceDisplacementData " begin
     1000 == length(importNI_forceDisplacementData("test/testData/anonNI.TXT"))
end

@testset "Test readIBW                   " begin
    # TO DO:
    # ERROR HANDLING
    @test  begin # "Read back written data: small array"
        filename = "Test1.ibw"
        inputData = Float32.(rand(100,3))
        writeIBW(filename,inputData)
        all(reshape(inputData, : , 1) .== readIBW(filename))
    end
    rm("Test1.ibw")
    @test  begin # "Read back written data: big array"
        filename = "Test2.ibw"
        inputData = Float32.(rand(10_000,3))
        writeIBW(filename,inputData)
        all(reshape(inputData, : , 1) .== readIBW(filename))
    end
    rm("Test2.ibw")
    @test begin # "Data supplied as a Float64 matrix"
        filename = "Test3.ibw"
        inputData = Float64.(rand(100,3))
        writeIBW(filename,inputData)
        all(Float32.(reshape(inputData, : , 1)) .== readIBW(filename))
    end
    rm("Test3.ibw")

    ### Here comes a series of regression tests designed on real data files.
    ### Obviously any failing test means only that behavior has changed, not
    ### necessarily that the previous implementation was correct. This is
    ### because the real .IBW file format is "black-box" to me.
    @test begin
        # Regression test 1
        inputVec = readIBW("test/testData/anonMG.ibw")
        length(inputVec) == 254541
    end
    @test begin
        # Regression test 2
        inputVec = readIBW("test/testData/anonMG.ibw")
        inputVec[12345] == -3.404614f-6
    end
    @test begin
        # Regression test 3
        inputVec = readIBW("test/testData/anonZF.ibw")
        length(inputVec) == 254466
    end
    @test begin
        # Regression test 4
        inputVec = readIBW("test/testData/anonZF.ibw")
        inputVec[end] == -2.6649975f-6
    end
    @test begin
        # Regression test 5
        inputVec = readIBW("test/testData/anonPMMA.ibw")
        length(inputVec) == 254961
    end
    @test begin
        # Regression test 6
        inputVec = readIBW("test/testData/anonPMMA.ibw")
        inputVec[1] == -2.3641041f-6
    end
    @test begin
        # Regression test 7
        inputVec = readIBW("test/testData/anonCell.ibw")
        length(inputVec) == 4299111
    end
    @test begin
        # Regression test 8
        inputVec = readIBW("test/testData/anonCell.ibw")
        inputVec[1337] == 1.5439662f-7
    end
end

@testset "Test IBWtoTXT                  " begin
    @test begin
        filename = "Test1.ibw"
        inputData = Float32.([collect(1:100) collect(101:200) collect(201:300)])
        writeIBW(filename,inputData)
        all(inputData[:,3:-1:2] .== IBWtoTXT(filename))
    end
    rm("Test1.ibw")
    @test begin
        filename = "Test1.ibw"
        inputData = Float32.(rand(100,3))
        writeIBW(filename,inputData)
        all(inputData[:,3:-1:2] .== IBWtoTXT(filename))
    end
    rm("Test1.ibw")
end

@testset "Test offsetAndDriftCompensation" begin
    ############################################
    # First data set:
    #   - No noise
    #   - Linear interpolation between keypoints
    #   - 
    sampleRate = 1000.0 # Hz
    #                 TIME   FORCE    DISPLACEMENT
        matOfKeyPoints = [ 0.0    100.0    -200.0        # Starting position
                           5.0    100.0     100.0        # Position of contact
                           6.0   1500.0     110.0        # End of ramp / Start of hold
                          16.0   1500.0     125.0        # End of hold / Start of unload
                          17.0    175.0     105.0        # End of unload / Start of thermal hold
                          47.0    175.0     101.0        # End of thermal hold / Start of retraction
                          48.0    100.0     100.0        # End of retraction / Start of exit
                          50.0    100.0    -100.0]       # End of exit & End of experiment
    xy = generateSyntheticForceDisplacementSignal(matOfKeyPoints,sampleRate) 
    xy , basefit , zeroLine , rampStartIdx, xZero = offsetAndDriftCompensation(Float32.(xy))
    @test begin
        isapprox(zeroLine[1] , (matOfKeyPoints[2,2]-matOfKeyPoints[1,2])/(matOfKeyPoints[2,3]-matOfKeyPoints[1,3]), atol = 1e-6)
    end
    @test begin
        isapprox(zeroLine[2] , 0.5*(matOfKeyPoints[1,2]+matOfKeyPoints[2,2]), atol = 1e-6)
    end
    @test begin
        isapprox(rampStartIdx , sampleRate.*matOfKeyPoints[2,1], rtol = 0.01)
    end
    @test begin
        isapprox(Float32(xZero), matOfKeyPoints[2,3], rtol = 0.01)
    end
    ############################################
    # Second data set:
    #   - Noise
    #   - Linear interpolation between keypoints
    #   - 
    sampleRate = 1000.0 # Hz
    #                 TIME   FORCE    DISPLACEMENT
        matOfKeyPoints = [ 0.0    100.0    -200.0        # Starting position
                           5.0    100.0     100.0        # Position of contact
                           6.0   1500.0     110.0        # End of ramp / Start of hold
                          16.0   1500.0     125.0        # End of hold / Start of unload
                          17.0    175.0     105.0        # End of unload / Start of thermal hold
                          47.0    175.0     101.0        # End of thermal hold / Start of retraction
                          48.0    100.0     100.0        # End of retraction / Start of exit
                          50.0    100.0    -100.0]       # End of exit & End of experiment
    # xy = generateSyntheticForceDisplacementSignal(matOfKeyPoints,sampleRate, stdOfNoise = 2.0,sizeOfNoise =[0.0 0.4 0.4]) 
    # xy , basefit , zeroLine , rampStartIdx, xZero = offsetAndDriftCompensation(Float32.(xy))
    # @test begin
    #     isapprox(zeroLine[1] , (matOfKeyPoints[2,2]-matOfKeyPoints[1,2])/(matOfKeyPoints[2,3]-matOfKeyPoints[1,3]), atol = 0.02)
    # end
    # @test begin
    #     isapprox(zeroLine[2] , 0.5*(matOfKeyPoints[1,2]+matOfKeyPoints[2,2]), atol = 10.0)
    # end
    # @test begin
    #     #println(rampStartIdx)
    #     #println(sampleRate.*matOfKeyPoints[2,1])
    #     isapprox(rampStartIdx, sampleRate.*matOfKeyPoints[2,1], atol = 25)
    # end
    # @test begin
    #     isapprox(Float32(xZero), matOfKeyPoints[2,3], rtol = 0.05)
    # end

end



@testset "Test findStartOfHold           " begin
    sampleRate = 1000.0 # Hz
    #                 TIME   FORCE    DISPLACEMENT
    matOfKeyPoints  = [ 0.0    100.0    -200.0        # Starting position
                        5.0    100.0     100.0        # Position of contact
                        6.0   1500.0     110.0        # End of ramp / Start of hold
                        16.0   1500.0     125.0        # End of hold / Start of unload
                        17.0    175.0     105.0        # End of unload / Start of thermal hold
                        47.0    175.0     101.0        # End of thermal hold / Start of retraction
                        48.0    100.0     100.0        # End of retraction / Start of exit
                        50.0    100.0    -100.0]       # End of exit & End of experiment
    xy = generateSyntheticForceDisplacementSignal(matOfKeyPoints,sampleRate) 
    xy , ~ , ~ , ~, ~ = offsetAndDriftCompensation(Float32.(xy))
    @test begin
        isapprox(findStartOfHold(xy, "first") , sampleRate*matOfKeyPoints[3,1], atol = 10)
    end
    @test begin
        isapprox(findStartOfHold(xy, "last") , sampleRate*matOfKeyPoints[4,1], atol = 10)
    end
    xy = generateSyntheticForceDisplacementSignal(matOfKeyPoints,sampleRate, stdOfNoise = 2.0, sizeOfNoise =[0.0 0.4 0.4]) 
    xy , ~ , ~ , ~, ~ = offsetAndDriftCompensation(Float32.(xy))
    @test begin
        isapprox(findStartOfHold(xy, "first") , sampleRate*matOfKeyPoints[3,1], atol = 30)
    end
    @test begin
        isapprox(findStartOfHold(xy, "last") , sampleRate*matOfKeyPoints[4,1], atol = 30)
    end
end

@testset "determineThermalCreep          " begin
    thermalHoldTime = 10000
    sampleRate = 1000
    holdTimeVals = collect(1:10000) ./ sampleRate
    thermalCreepGenerator(x,y,z) = x .+ y.*holdTimeVals.^z
    @test begin
        xyStart = [collect( 1000.0:-1:500 ) ; 500.0.+10.0*rand(8998) ; collect(500:-1:0.0)]
        xy_hold = thermalCreepGenerator(44.0,4.0,1.1)
        isapprox(-4.0 * 1.1*holdTimeVals[end]^(1.1 - 1.0) , determineThermalCreep(Float32.([xy_hold[end:-1:1] xyStart]),sampleRate, thermalHoldTime, control(false, false), 5.0),rtol = 0.10)
    end
    @test begin
        xyStart = [collect( 1000.0:-1:500 ) ; 500.0.+10.0*rand(8998) ; collect(500:-1:0.0)]
        xy_hold = thermalCreepGenerator(123.0,4.0,1.1)
        isapprox(-4.0 * 1.1*holdTimeVals[end]^(1.1 - 1.0) , determineThermalCreep(Float32.([xy_hold[end:-1:1] xyStart]),sampleRate, thermalHoldTime, control(false, false), 5.0),rtol = 0.10)
    end
    @test begin
        xyStart = [collect( 1000.0:-1:500 ) ; 500.0.+10.0*rand(8998) ; collect(500:-1:0.0)]
        xy_hold = thermalCreepGenerator(123.0,0.20,1.1)
        isapprox(-0.2 * 1.1*holdTimeVals[end]^(1.1 - 1.0) , determineThermalCreep(Float32.([xy_hold[end:-1:1] xyStart]),sampleRate, thermalHoldTime, control(false, false), 5.0),rtol = 0.10)
    end
    @test begin
        xyStart = [collect( 1000.0:-1:500 ) ; 500.0.+10.0*rand(8998) ; collect(500:-1:0.0)]
        xy_hold = thermalCreepGenerator(123.0,4.00,0.8)
        isapprox(-4.0 * 0.8*holdTimeVals[end]^(0.8 - 1.0) , determineThermalCreep(Float32.([xy_hold[end:-1:1] xyStart]),sampleRate, thermalHoldTime, control(false, false), 5.0),rtol = 0.15)
    end
    @test begin
        xyStart = [collect( 1000.0:-1:500 ) ; 500.0.+10.0*rand(8998) ; collect(500:-1:0.0)]
        xy_hold = thermalCreepGenerator(0.0,4.00,0.8)
        isapprox(-4.0 * 0.8*holdTimeVals[end]^(0.8 - 1.0) , determineThermalCreep(Float32.([xy_hold[end:-1:1] xyStart]),sampleRate, thermalHoldTime, control(false, false), 5.0),rtol = 0.15)
    end
    @test begin
        xyStart = [collect( 1000.0:-1:500 ) ; 500.0.+10.0*rand(8998) ; collect(500:-1:0.0)]
        xy_hold = thermalCreepGenerator(123.0,0.00,0.8)
        isapprox(-0.0 * 0.8*holdTimeVals[end]^(0.8 - 1.0) , determineThermalCreep(Float32.([xy_hold[end:-1:1] xyStart]),sampleRate, thermalHoldTime, control(false, false), 5.0),atol = 0.05)
    end
    @test begin
        xyStart = [collect( 1000.0:-1:500 ) ; 500.0.+10.0*rand(8998) ; collect(500:-1:0.0)]
        xy_hold = thermalCreepGenerator(123.0,4.00,0.0)
        isapprox(-4.0 * 0.0*holdTimeVals[end]^(0.0 - 1.0) , determineThermalCreep(Float32.([xy_hold[end:-1:1] xyStart]),sampleRate, thermalHoldTime, control(false, false), 5.0),atol = 0.25)
    end
end

;@testset "determineCreepDuringHold       " begin
    sampleRate = 1000
    holdTimeVals = collect(1:10000) ./ sampleRate
    holdCreepGenerator(x,y,z) = x .+ y.*holdTimeVals.^z
    @test begin
        xy_hold = holdCreepGenerator(44.0,4.0,1.1)
        isapprox(4.0 * 1.1*holdTimeVals[end]^(1.1 - 1.0) , determineCreepDuringHold(xy_hold,sampleRate),rtol = 0.01)
    end
    @test begin
        xy_hold = holdCreepGenerator(123.0,4.0,1.1)
        isapprox(4.0 * 1.1*holdTimeVals[end]^(1.1 - 1.0) , determineCreepDuringHold(xy_hold,sampleRate),rtol = 0.01)
    end
    @test begin
        xy_hold = holdCreepGenerator(123.0,0.20,1.1)
        isapprox(0.2 * 1.1*holdTimeVals[end]^(1.1 - 1.0) , determineCreepDuringHold(xy_hold,sampleRate),rtol = 0.01)
    end
    @test begin
        xy_hold = holdCreepGenerator(123.0,4.00,0.8)
        isapprox(4.0 * 0.8*holdTimeVals[end]^(0.8 - 1.0) , determineCreepDuringHold(xy_hold,sampleRate),rtol = 0.01)
    end
    @test begin
        xy_hold = holdCreepGenerator(0.0,4.00,0.8)
        isapprox(4.0 * 0.8*holdTimeVals[end]^(0.8 - 1.0) , determineCreepDuringHold(xy_hold,sampleRate),rtol = 0.02)
    end
    @test begin
        xy_hold = holdCreepGenerator(123.0,0.00,0.8)
        isapprox(0.0 * 0.8*holdTimeVals[end]^(0.8 - 1.0) , determineCreepDuringHold(xy_hold,sampleRate),atol = 0.01)
    end
    @test begin
        xy_hold = holdCreepGenerator(123.0,4.00,0.0)
        isapprox(4.0 * 0.0*holdTimeVals[end]^(0.0 - 1.0) , determineCreepDuringHold(xy_hold,sampleRate),atol = 0.01)
    end
end








@testset "modulusfitter                  " begin
    ctrl = control(false, false)
    indentSet = metaInfoExperimentalSeries("Test2", 
                                            50.0,
                                        "hemisphere",
                                        "X",
                                        217.51, 
                                        "test/testData/sampleAreaFunction.txt",
                                        "test/testData/",
                                        50000, 
                                        "afm")

    # Regression tests
    @test begin
        hp =  hyperParameters( 2000, 1400, "Oliver-Pharr", false, 0, 0, 0.0 , 5.0)
        isapprox(1.8696, modulusfitter( indentSet, hp, ctrl, "anonMG.ibw"),rtol = 0.05)
    end
    @test begin
        hp =  hyperParameters( 2000, 1400, "Oliver-Pharr", true, 0, 0, 0.0, 5.0)
        isapprox(1.9, modulusfitter( indentSet, hp, ctrl, "anonMG.ibw"),rtol = 0.05)
    end
    @test begin
        hp =  hyperParameters( 2000, 1400, "Feng", false, 0, 0, 0.0, 5.0)
        isapprox( 1.9289, modulusfitter( indentSet, hp, ctrl, "anonMG.ibw"),rtol = 0.05)
    end
    @test begin
        hp =  hyperParameters( 2000, 1400, "Feng", true, 0, 0, 0.0, 5.0)
        isapprox( 1.92197, modulusfitter( indentSet, hp, ctrl, "anonMG.ibw"), rtol = 0.05)
    end
    
    indentSet = metaInfoExperimentalSeries("Test2", 
                                            50.0,
                                        "hemisphere",
                                        "X",
                                        217.51, 
                                        "vickers",
                                        "test/testData/",
                                        50000, 
                                        "ni")
     @test begin
        hp =  hyperParameters( 10, 10, "Oliver-Pharr", false, 0, 0, 0.0, 5.0)
        
        isapprox( 1.92197, modulusfitter( indentSet, hp, ctrl, "anonNI.TXT"), rtol = 5.05)
    end
end
