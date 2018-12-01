# meteor_pbs18

Meteorite crashing into the sea - Physically-Based Simulation HS18 project.

Authors: Silvia Nauer, Mikael Stellio, Sean Bone.

Slack channel: [https://flippingteamworkspace.slack.com](https://flippingteamworkspace.slack.com)

# Cloning with submodules
This Git repository has a submodule for libigl. To clone it correctly use either of the following:

    git clone --recurse-submodules git@gitlab.ethz.ch:bones/meteor_pbs18.git
    git clone --recurse-submodules https://gitlab.ethz.ch/bones/meteor_pbs18.git

# TODO
 - [X] Set up 2D visualization for FLIP particles & cells
 - [X] 2D FLIP working correctly
   - [X] MAC2D data structure
   - [X] 2D FLIP updates
 
 - [ ] Extend to 3D
   - [ ] MAC3D data structure
   - [ ] 3D viz
   - [ ] 3D FLIP updates
   - [ ] Signed distance function

 - [ ] Marching cubes & export mesh at each frame
 - [ ] Import meshes into Maya/Blender for rendering

# Making a video from PNGs
This command will take files numbered `0000.png, 0001.png, 0002.png, ...` in folder `PNG_out` (usually generated in the build directory by WaterSim when "Export PNGs" is enabled) and make an MP4 video of them:
    
    ffmpeg -r 40 -f image2 -s 1280x800 -i PNG_out/%04d.png -vcodec libx264 -crf 20  -pix_fmt yuv420p test.mp4

 - `-r`: frame rate ( = 1/dt, dt is the timestep used in simulation)
 - `-s`: dimensions
 - `-crf`: quality, lower is better
 - `-i`: input images. Change folder name/location if necessary.
