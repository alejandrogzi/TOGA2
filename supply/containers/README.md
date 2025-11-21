# Description

## Building container for local execution
If you want to run TOGA2 on a local machine or a multiple CPU cluster node, consider building your container with the `apptainer.def` definition file. This defition contains recipes for pulling the latest public TOGA2 version and installing it alongside with all the necessary dependencies. Building a container from this definition file is straightforward and does not require any additional operations:
```
apptainer build <your_container_name> apptainer.def
```
Note that the resulting container supports only local execution. If you want to use TOGA2 with a cluster job manager as Nextflow executor, proceed to the next section.

## Building container for cluster execution
Configuring containers for cluster job manager compatibility requires additional setup, including that  