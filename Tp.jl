using LinearAlgebra,DelimitedFiles, SparseArrays
using Plots
coordinates = readdlm("coordinates.dat")
coordinates = coordinates[:,2:end]

elements3 = try
    data = readdlm("elements3.dat")
    convert(Array{Int,2}, data[:, 2:end])
catch
    zeros(Int, 0, 0)
end

elements4 = try
    data = readdlm("elements4.dat")
    convert(Array{Int,2}, data[:, 2:end])
catch
    zeros(Int, 0, 0)
end

neumann = try
    data = readdlm("neumann.dat")
    data[:,2:end]
catch
    zeros(0, 0)
end

dirichlet = readdlm("dirichlet.dat")
dirichlet = dirichlet[:,2:end]

FreeNodes = setdiff(1:size(coordinates, 1), unique(dirichlet))

A = spzeros(size(coordinates, 1), size(coordinates, 1))
b = spzeros(size(coordinates, 1), 1)

function stima3(vertices)
    d = size(vertices, 2)
    G = [ones(1, d+1); vertices'] \ [zeros(1, d); I(d)]
    M = det([ones(1, d+1); vertices']) * G * G' / prod(1:d)
    return M
end

function stima4(vertices)
    D_Phi = transpose([vertices[2, :] - vertices[1, :]; vertices[4, :] - vertices[1, :]])
    B = (D_Phi'* D_Phi)
    print(B)
    C1 = [2 -2; -2  2] * B[1  1] + [3 0; 0 -3] * B[1 2] + [2 1; 1 2] * B[2 2]
    C2 = [-1 1; 1 -1] .* B[1 1] .+ [-3 0; 0 3] .* B[1 2] .+ [-1 -2; -2 -1] .* B[2 2]
    M = det(D_Phi) * [C1 C2; C2 C1] / 6
    return M
end

function u_d(x)
    return zeros(size(x, 1))
end

function f(x)
    return ones(size(x, 1))
end

function g(x)
    return zeros(size(x, 1))
end

function show(elements3, elements4, coordinates, u)
    
    triplot(elements3, coordinates[:, 1], coordinates[:, 2], zcolor=u, fill=true)
    plot!(triplot(elements4, coordinates[:, 1], coordinates[:, 2], zcolor=u, fill=true))
    title!("Solution of the Problem")
    view(10, 40)
end

for j in 1:length(elements3[:,1])
    indices = elements3[j, :]
    A[indices, indices] += stima3(coordinates[indices, :])
end

a = [elements4[2, :],:][2,:] 

stima4(coordinates[elements4[2, :],:])
#stima4(coordinates[[2,1], :])

for j in 1:length(elements4[:,1])
    indices = elements4[j, :]
    A[indices, indices] += stima4(coordinates[indices, :])
end

# Volume Forces
for j in 1:size(elements3, 1)
    b[elements3[j, :]] += det([1 1 1; coordinates[elements3[j, :], :]']) * f(sum(coordinates[elements3[j, :], :], dims=1) / 3) / 6
end

for j in 1:size(elements4, 1)
    b[elements4[j, :]] += det([1 1 1; coordinates[elements4[j, 1:3], :]']) * f(sum(coordinates[elements4[j, :], :], dims=1) / 4) / 4
end

# Neumann conditions
for j in 1:size(neumann, 1)
    b[neumann[j, :]] += norm(coordinates[neumann[j, 1], :] - coordinates[neumann[j, 2], :]) * g(sum(coordinates[neumann[j, :], :], dims=1) / 2) / 2
end

# Dirichlet conditions
u = spzeros(size(coordinates, 1))
u[unique(dirichlet)] = u_d(coordinates[unique(dirichlet), :])
b = b - A * u

# <<< method 1
A = Matrix(A)
b = Matrix(b)
b[unique(dirichlet)] = u_d(coordinates[unique(dirichlet), :])
A[unique(dirichlet), :] .= 0.0
A[:, unique(dirichlet)] .= 0.0
temp = unique(dirichlet)

for i in 1:length(temp)
    A[temp[i], temp[i]] = 1.0
end

u = A \ b

show(elements3, elements4, coordinates, u)
