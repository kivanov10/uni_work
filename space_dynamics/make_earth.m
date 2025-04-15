% This is a simple script that is used throughout the coursework main
% scripts to generate a realistic looking Earth mesh. It tries to generate
% the custom mesh from a file obtained from the mathworks website. If the
% file is not present it will make a simple gray mesh, in case the user
% does not want to download any additional materials. This script was done
% mainly for convinience and save on space in the main files.

% Earth Mesh taken from  Will Campbell (2024). Earth-sized Sphere with Topography
% (https://www.mathworks.com/matlabcentral/fileexchange/27123-earth-
% sized-sphere-with-topography)

%Try/catch statement in place in case user does not have the custom Earth
%mesh. Instead a gray Earth-sized sphere is generated.
try
earth_sphere('km')
catch
    disp(['Error when trying to create textured Earth mesh. Using generic' ...
        ' 3D mesh to represent the Earth instead. earth_sphere potentially ' ...
        'missing. Please download the file from the link in the make_earth.m ' ...
        'subscript.'])
    earth_radius_km = 6371; % [km]
    [X_mesh, Y_mesh, Z_mesh] = sphere(50); % making mesh
    X_mesh = X_mesh * earth_radius_km; % scaling all axis by Earth radius
    Y_mesh = Y_mesh * earth_radius_km;
    Z_mesh = Z_mesh * earth_radius_km;
    earth_mesh = mesh(X_mesh,Y_mesh,Z_mesh);
    daspect([1 1 1]);
    earth_mesh.FaceColor = [0.5, 0.5, 0.5]; % Grey color for the faces
    earth_mesh.EdgeColor = 'none'; % remove the mesh lines
end

% ECI unit vectors arrows
quiver3(0,0,0,1,0,0,1e4, 'r');
text(1e4,0,0,'i');
quiver3(0,0,0,0,1,0,1e4, 'g');
text(0,1e4,0,'j');
quiver3(0,0,0,0,0,1,1e4, 'b');
text(0,0,1e4,'k');