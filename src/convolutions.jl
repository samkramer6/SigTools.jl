#=  convolutions.jl
    Convolutions Functions
    

    Basic defintions for convolutions for vectors and matrices.
=#

using FFTW
using Test

########################################################################
#                           2D Convolutions                            #
########################################################################

function conv(x::AbstractMatrix, y::AbstractMatrix)

    (x, y) = same_length(x, y);

    f_x = fft(x);
    f_y = fft(y);

    @inbounds conv2_xy = real.(ifft(f_x .* f_y));

    return conv2_xy::Matrix{Float64}
end

function conv(x::Matrix{Float64}, y::Matrix{Float64})

    (x, y) = same_length(x, y);

    f_x = fft(x);
    f_y = fft(y);

    @inbounds conv2_xy = real.(ifft(f_x .* f_y));

    return conv2_xy::Matrix{Float64}
end

########################################################################
#                           1D Convolutions                            #
########################################################################

function conv(x::Vector{Int64}, y::Vector{Int64})

    (x, y) = same_length(x, y);

    f_x = fft(x);
    f_y = fft(y);

    @inbounds conv_xy = real.(ifft(f_x .* f_y));

    return conv_xy::Vector{Float64}
end

function conv(x::Vector{Float64}, y::Vector{Float64})

    (x, y) = same_length(x, y);

    f_x = fft(x);
    f_y = fft(y);

    @inbounds conv_xy = real.(ifft(f_x .* f_y));

    return conv_xy::Vector{Float64}
end

function conv(x::AbstractVector, y::AbstractVector)

    (x, y) = same_length(x, y);

    f_x = fft(x);
    f_y = fft(y);

    @inbounds conv_xy = real.(ifft(f_x .* f_y));

    return conv_xy::AbstractVector
end

########################################################################
#                           same_length()                              #
########################################################################

function same_length(x::AbstractVector, y::AbstractVector)

    if length(x) == length(y)
        return x, y
    end

    leng_diff = abs.(length(x) - length(y));
    pad = vec(zeros(typeof(x[1]), 1, leng_diff));

    if length(x) > length(y)

        @inbounds append!(y, pad);

    elseif length(y) > length(x)

        @inbounds append!(x, pad);

    end

    return x::AbstractVector, y::AbstractVector
end

function same_length(x::Vector{Float64}, y::Vector{Float64})

    if length(x) == length(y)
        return x, y
    end

    leng_diff = abs.(length(x) - length(y));
    pad = vec(zeros(Float64, 1, leng_diff));

    if length(x) > length(y)

        @inbounds append!(y, pad);

    elseif length(y) > length(x)

        @inbounds append!(x, pad);

    end

    return x::Vector{Float64}, y::Vector{Float64}
end

########################################################################
#                            2D same_length()                          #
########################################################################

function same_length(x::Matrix{T} where T <: Real, y::Matrix{T} where T<: Real)
    
	# --Test
		if size(x) == size(y)
        	return x::AbstractMatrix, y::AbstractMatrix
    	end

	# --Plan Height
    	height_plan = size(larger(x, y, 1), 1);
    	width_plan = size(larger(x, y, 2), 2);

	# --Heigth Padding
		height_diff = abs.(height_plan - size(smaller(x, y, 1), 1));
		height_pad = zeros(height_diff, size(smaller(x, y, 1), 2));
	
		if length(height_pad) != 0

            @inbounds a = vcat(smaller(x, y, 1), height_pad);
            b = larger(x, y, 1);

    	elseif length(height_pad) == 0

            a = x;
            b = y;

    	end

	# --Width Padding
		width_diff = abs.(width_plan - size(smaller(x, y, 2), 2));
		width_pad = zeros(width_diff, size(smaller(x, y, 2), 1));
	
		if length(width_pad) != 0
			
			pad = zeros(height_plan, width_diff);
			@inbounds c = hcat(smaller(a, b, 2), pad);
			d = larger(a, b, 2);
			
			
		elseif length(width_pad) == 0

			c = a;
			d = b;

		end

	if size(c) == size(d)
		return c::AbstractMatrix, d::AbstractMatrix
	else
		(c, d) == same_length(c, d);
	end
end

function same_length(x::Matrix{Float64}, y::Matrix{Float64})

	# --Test
		if size(x) == size(y)
            return x::Matrix{Float64}, y::Matrix{Float64}
        end

	# --Plan Height
        height_plan = size(larger(x, y, 1), 1);
        width_plan = size(larger(x, y, 2), 2);

	# --Heigth Padding
        height_diff = abs.(height_plan - size(smaller(x, y, 1), 1));
	    height_pad = zeros(height_diff, size(smaller(x, y, 1), 2));
	
		if length(height_pad) != 0

            @inbounds a = vcat(smaller(x, y, 1), height_pad);
            b = larger(x, y, 1);
    
        elseif length(height_pad) == 0

            a = x;
            b = y;

        end

	# --Width Padding
		width_diff = abs.(width_plan - size(smaller(x, y, 2), 2));
		width_pad = zeros(width_diff, size(smaller(x, y, 2), 1));
	
		if length(width_pad) != 0
			
			pad = zeros(height_plan, width_diff);
		    @inbounds c = hcat(smaller(a, b, 2), pad);
			d = larger(a, b, 2);
			
			
		elseif length(width_pad) == 0

			c = a;
			d = b;

		end

	if size(c) == size(d)
		return c::Matrix{Float64}, d::Matrix{Float64}
	else
		(c, d) == same_length(c, d);
	end
end

function same_length(x::Matrix{Int64}, y::Matrix{Int64})

    x = convert(Matrix{Float64}, x);
    y = convert(Matrix{Float64}, y);

    (x, y) = same_length(x, y);

    return x::Matrix{Float64}, y::Matrix{Float64}
end

########################################################################
#                                larger()                              #
########################################################################

function larger(x::AbstractVector, y::AbstractVector)

    if length(x) > length(y)
        return x::AbstractVector
    elseif length(y) > length(x)
        return y::AbstractVector
    else 
        return x::AbstractVector
    end

end

function larger(x::AbstractMatrix, y::AbstractMatrix, n::Int64)

    if size(x, n) > size(y, n)
        return x::AbstractMatrix
    elseif size(y, n) > size(x, n)
        return y::AbstractMatrix
    else
        return x::AbstractMatrix
    end
    
end

function larger(x::AbstractArray, y::AbstractArray, n::Int64)

    if size(x, n) > size(y, n)
        return x::AbstractArray
    elseif size(y, n) > size(x, n)
        return y::AbstractArray
    else
        return x::AbstractArray
    end
    
end


"""
    Sigtools.jl smaller()
    This function is used to find the smaller of the two matrices in the n dimension.
    This helps the same_length() function that makes matrices the same length.

    Inputs:
        x::AbstractMatrix --> Data matrix 1
        y::AbstractMatrix --> Data matrix 2
        n::Int64 --> Dimension direction (range 1 or 2)

    Outputs:
        x::AbstractMatrix

    Sam Kramer

    Methods:
"""
function smaller(x::AbstractMatrix, y::AbstractMatrix, n::Int64)

    if size(x, n) < size(y, n)
        return x::AbstractMatrix
    elseif size(y, n) < size(x, n)
        return y::AbstractMatrix
    else
        return x::AbstractMatrix
    end

end

function smaller(x::AbstractVecOrMat, y::AbstractVecOrMat, n::Int64)

    if size(x, n) < size(y, n)
        return x::AbstractVecOrMat
    elseif size(y, n) < size(x, n)
        return y::AbstractVecOrMat
    else
        return x::AbstractVecOrMat
    end

end
