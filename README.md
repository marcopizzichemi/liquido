# **LiquidO**

Forked from https://github.com/marcgrassi/liquido

At commit 59655310ec2700d87c63438685c14405f948ba58


## **Run Modalities**

Simulation can be run in 3 modes

1. Viz mode:
/path/to/build/ofos

This modality will visualize the default detector geometry


2. Batch mode:
/path/to/build/ofos macro.mac

This modality will execute the simulation on the basis of the geometry and primary particles specified in macro.mac


3. Viz macro mode:
/path/to/build/ofos macro.mac 1

This modality will visualize the geometry described in macro.mac. Please comment out the ** /run/beamOn ** command in macro.mac, otherwise it will be executed before visualization.
