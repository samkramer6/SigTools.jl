#   edge_detection.jl

using .SigTools

function sobel(N::Int64, M::Int64)

    size_buffer = zeros(N, M)
    sobel_x = [1 0 -1; 2 0 -2; 1 0 -1]
    sobel_y = [1 2 1; 0 0 0; -1 -2 -1]
    sobel = sqrt.(sobel_x .^ 2 .+ sobel_y .^ 2)

    sobel = same_length(size_buffer, sobel_x)

    return sobel[2]::Matrix{Float64}
end

function prewitt(N::Int64, M::Int64)

    size_buffer = zeros(N, M)
    sobel_x = [1 0 -1; 1 0 -1; 1 0 -1]
    sobel_y = [1 1 1; 0 0 0; -1 -1 -1]
    sobel = sqrt.(sobel_x .^ 2 .+ sobel_y .^ 2)

    sobel = same_length(size_buffer, sobel_x)

    return sobel[2]::Matrix{Float64}
end

function roberts(N::Int64, M::Int64)

    size_buffer = zeros(N, M)
    sobel_x = [1 0; 0 -1]
    sobel_y = [0 1; -1 0]
    sobel = sqrt.(sobel_x .^ 2 .+ sobel_y .^ 2)

    sobel = same_length(size_buffer, sobel_x)

    return sobel[2]::Matrix{Float64}
end

"""
SigTools.jl   edge_detection()
This is an edge detection algorithm that can take a variety of different detection kernels.

Inputs:

Outputs:

Sam Kramer
"""
function edge_detection(kernel::AbstractString, image::Matrix{Float64})

    if uppercase(kernel) == "SOBEL"
        det_kernel = sobel(size(image, 1), size(image, 2))
    elseif uppercase(kernel) == "PREWITT"
        det_kernel = prewitt(size(image, 1), size(image, 2))
    elseif uppercase(kernel) == "ROBERTS"
        det_kernel = roberts(size(image, 1), size(image, 2))
    else    # default Sobel filter
        det_kernel = sobel(size(image, 1), size(image, 2))
        println("Defaulting to Sobel Kernel")
    end

    edge_detection_image = convol(det_kernel, image)

    return edge_detection_image::Matrix{Float64}
end
