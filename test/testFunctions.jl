### Packages used
using Distributions

function writeIBW(filename,inputData)

    # Convert matrix input data to a vector with the two matrix columns after one another
    inputData = reshape(inputData, :, 1)

    A = open( filename, "w") 
    write(A, ltoh(Int8(5)))   
    close(A)

    A = open( filename, "w")
    ### Note that the following is just to copy the binary format of the IBW file.
    ### We are filling the bytes with garbage, the only point is that the garbage 
    ### takes the exact amount of space as the header in the IBW file. This allows
    ### us to pull out the data at the end of the function.
    write(A, Int16(0))
    write(A,UInt16(0)) # Does not work           
    
    for xLoop = 1:16
        write(A, Int32(0))
    end

    write(A, UInt32(0))
    write(A, UInt32(0))
    write(A, Int32(length(inputData)))  # npts, works but not yet calculated based on the matrix to write!

    write(A, Int16(2))   # type = single

    write(A, Int16(0))
    for aLoop = 1:6
        write(A, Char('x'))
    end
    write(A, Int16(0))

    for aLoop = 1:32
        write(A, Char('x'))
    end

    for aLoop = 1:6
        write(A, Int32(0))
    end

    write(A,Float64(0.0005)) # sFA[1]
    write(A,Float64(1.0))
    write(A,Float64(1.0))
    write(A,Float64(1.0))

    write(A,Float64(0.0))    # sFB[1]
    write(A,Float64(0.0))
    write(A,Float64(0.0))
    write(A,Float64(0.0))

    for aLoop = 1:4
        write(A, Char('x'))
    end
    
    for aLoop = 1:16
        write(A, Char('x'))
    end
    
    write(A, Int16(0))
    write(A, Int16(0))
    write(A, Float64(0.0))
    write(A, Float64(0.0))
    for aLoop = 1:26
        write(A, Int32(0))
    end
    for aLoop = 1:3
        write(A, Int16(0))
    end
    write(A, Char('x'))
    write(A, Char('x'))
    write(A, Int32(0))
    write(A, Int16(0))
    write(A, Int16(0))
    write(A, Int32(0))
    write(A, Int32(0))
    write(A, Int32(0))

    write(A, Float32.(inputData))
    # Here is where we write the real data.

    close(A)
    # Close file and end operation.
end

function generateSyntheticForceDisplacementSignal(matOfKeyPoints,sampleRate ; stdOfNoise = 0.0,sizeOfNoise =[0.0 0.0 0.0])
    # Should generate a well controlled force-displacement curve for testing.

    enrichedSignal = Array{Float64,2}(undef,0,3)
    for aLoop = 1:(size(matOfKeyPoints,1)-1)       
        numelResults = length(matOfKeyPoints[aLoop,1]:1/sampleRate:matOfKeyPoints[aLoop+1,1])
        tmpMat = collect(matOfKeyPoints[aLoop,1]:1/sampleRate:matOfKeyPoints[aLoop+1,1])
        tmpMat2 = collect(range(matOfKeyPoints[aLoop,2], matOfKeyPoints[aLoop+1,2], length = numelResults))
        tmpMat3 = collect(range(matOfKeyPoints[aLoop,3], matOfKeyPoints[aLoop+1,3], length = numelResults))
        enrichedSignal = [enrichedSignal ; [tmpMat tmpMat2 tmpMat3]]

    end
    noise = sizeOfNoise.*rand(Normal(0,stdOfNoise),size(enrichedSignal))
    enrichedSignal = noise.+enrichedSignal
    return enrichedSignal[:, 3:-1:2]
end

