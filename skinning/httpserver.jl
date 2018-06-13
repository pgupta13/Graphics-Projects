#=
Author: Yotam Gingold <yotam (strudel) yotamgingold.com>
License: Public Domain [CC0](http://creativecommons.org/publicdomain/zero/1.0/)
Description: A simple static file server written in Julia. Listens on localhost port 8000. Doesn't serve files outside current working directory (unless symlinked). Doesn't print directory listings.
URL: https://gist.github.com/yig/f65e86b7730019d4060449f24342fcb4
=#

using HttpServer
using Base.Filesystem

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

server = Server( http )
## Listen on localhost port 8000
run(server, host=IPv4(127,0,0,1), port=8000)
