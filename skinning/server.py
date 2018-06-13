import asyncio
import websockets
import json
import numpy as np
import bbw

async def bbw_server( websocket, path ):
    async for msg in websocket:
        if msg.startswith( "bbw " ):
            parts = msg.split()
            nfaces, nvertices, laplacian_mode, solver_mode = parts[1:]
            nfaces = int( nfaces )
            nvertices = int( nvertices )
            
            faces    = await websocket.recv()
            vertices = await websocket.recv()
            handles  = await websocket.recv()
            
            faces    = np.frombuffer( faces,    dtype = np.int32   ).reshape( nfaces,    -1 ).copy()
            vertices = np.frombuffer( vertices, dtype = np.float32 ).reshape( nvertices, -1 ).copy()
            handles  = np.frombuffer( handles,  dtype = np.int32 ).copy()
            
            weights = bbw.bbw( faces, vertices, handles, laplacian_mode, solver_mode )
            
            await websocket.send( np.ascontiguousarray( weights, np.float32 ).tobytes() )
        
        elif msg.startswith( "linear_blend_skin_2D " ):
            parts = msg.split()
            nvertices, = parts[1:]
            nvertices = int( nvertices )
            
            vertices = await websocket.recv()
            weights  = await websocket.recv()
            transforms  = await websocket.recv()
            
            vertices = np.frombuffer( vertices, dtype = np.float32 ).reshape( nvertices, -1 ).copy()
            transforms  = np.array( json.loads( transforms ) )
            weights  = np.frombuffer( weights,  dtype = np.float32 ).reshape( nvertices, len(transforms) ).copy()
            
            deformed = bbw.linear_blend_skin_2D( vertices, weights, transforms )
            
            await websocket.send( np.ascontiguousarray( deformed, np.float32 ).tobytes() )

port_websocket = 9876
port_http = 8000

## Also start an http server on port 8000
def serve_http( port ):
    import os
    pid = os.fork()
    ## If we are the child, serve files
    if pid == 0:
        import http.server
        import socketserver
        Handler = http.server.SimpleHTTPRequestHandler
        with socketserver.TCPServer(("localhost", port), Handler) as httpd:
            print("Serving HTTP on port", port)
            httpd.serve_forever()

serve_http( port_http )

print("WebSocket server on port", port_websocket )
start_server = websockets.serve( bbw_server, 'localhost', port_websocket )
asyncio.get_event_loop().run_until_complete(start_server)
asyncio.get_event_loop().run_forever()
