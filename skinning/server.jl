using HttpServer
using WebSockets
using Base.Filesystem
# using .bbw
include("bbw.jl")

"""
    splitpath( path::String ) -> Array{String}

Splits a path into an array of its path components. The inverse of `joinpath()`.
Calling `joinpath( splitpath( path )... )` should produce `path`
(possible with the trailing slash removed).

```jldoctest
julia> splitpath("a/b/c")
("a", "b", "c")
julia> splitpath("/a/b/c")
("/", "a", "b", "c")
```
"""
function splitpath( path::String )
    result = String[]
    
    while path != ""
        path, last = splitdir( path )
        
        ## If path consists of only the path separator, which could happen
        ## when referring to the filesystem root, then last will be empty
        ## and path will be unchanged.
        ## If this is the case, swap them so that we push the root path
        ## marker and then the while loop will terminate.
        if last == ""
            path, last = last, path
        end
        
        push!( result, last )
    end
    
    reverse!( result )
    return result
end

http = HttpHandler() do req::Request, res::Response
    ## If it's not a path inside the jail (current working directory), return a 403.
    jail = realpath( pwd() )
    # println( "Jail path: ", jail )
    
    ## Convert the request path to a real path (remove symlinks, remove . and ..)
    ## UPDATE: realpath() assumes that the input points to something that exists.
    ##         We can't use it, since we don't know if the requested path exists.
    ##         We'll use normpath() instead.
    ## Requests start with "/" referring to the server as the root.
    @assert length(req.resource) > 0 && req.resource[1] == '/'
    ## Drop the leading "/".
    request_path = normpath( req.resource[2:end] )
    ## Convert it to a relative path.
    relative_path = relpath( request_path, jail )
    println( "Requested path: ", request_path )
    ## Return forbidden if the request is above the current working directory.
    if length( relative_path ) > 0 && splitpath( relative_path )[1] == ".."
        @show Response( 403 )
    ## Return not implemented if the request is for a directory.
    ## TODO: If it's a directory, return a file listing.
    elseif isdir( relative_path )
        ## Is there an index.html?
        index = joinpath( relative_path, "index.html" )
        if isfile( index )
            @show HttpServer.FileResponse( index )
        else
            @show Response( 501 )
        end
    ## Return the contents for a file.
    elseif isfile( relative_path )
        @show HttpServer.FileResponse( relative_path )
    else
        @show Response( 404 )
    end
    ## If it's not a file, return a 404.
end

websock = WebSocketHandler() do req::Request, client::WebSocket
    while true
        msg = String( read(client) )
        if startswith( msg, "bbw " )
            parts = split( msg )
            
            nfaces, nvertices, laplacian_mode, solver_mode = parts[2:end]
            nfaces = parse( Int, nfaces )
            nvertices = parse( Int, nvertices )
            
            faces = reinterpret( Int32, read(client) )
            ## These were sent row-major but we want column-major.
            faces = reshape( faces, 3, nfaces ).'
            ## Faces is 0-indexed. Add 1.
            faces .+= 1
            
            vertices = reinterpret( Float32, read(client) )
            ## These were sent row-major but we want column-major.
            vertices = reshape( vertices, :, nvertices ).'
            
            handles = reinterpret( Int32, read(client) )
            ## handles is 0-indexed. Add 1.
            handles .+= 1
            
            # weights = ones( Float32, nvertices, nhandles )./nhandles
            weights = bbw( faces, vertices, handles, laplacian_mode, solver_mode )
            
            weights = convert(Array{Float32}, weights)
            write( client, reinterpret( UInt8, (weights')[:] ) )
        elseif startswith( msg, "linear_blend_skin_2D " )
            parts = split( msg )
            nvertices, = parts[2:end]
            nvertices = parse( Int, nvertices )
            
            vertices = reinterpret( Float32, read(client) )
            weights = reinterpret( Float32, read(client) )
            transforms = JSON.parse( String( read(client) ) )
            
            ## These were sent row-major but we want column-major.
            vertices = reshape( vertices, :, nvertices ).'
            weights = reshape( weights, :, nvertices ).'
            
            ## I can't figure out the clever way to convert transforms into
            ## a 3D array, so I'll do this.
            transformsf = zeros( size(transforms,1), 3, 3 )
            for i in 1:size(transformsf,1)
                for ii in 1:3
                    for ij in 1:3
                        transformsf[i,ii,ij] = transforms[i][ii][ij]
                    end
                end
            end
            transforms = transformsf
            
            deformed = linear_blend_skin_2D( vertices, weights, transforms )
            
            deformed = convert(Array{Float32}, deformed)
            write( client, reinterpret( UInt8, (deformed')[:] ) )
        end
    end
end

server = Server( http, websock )
## Listen on localhost port 8000
run(server, host=IPv4(127,0,0,1), port=9876)
