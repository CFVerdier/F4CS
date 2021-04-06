# 4CS

Correct-by-Construction Closed-form Controller Synthesis

# Installation

## Install dReal 
Verification of the LF-like conditions is done using an SMT solver. Our tool relies on the SMT solver dReal. The installation guide of dReal 4 can be found [here](https://github.com/dreal/dreal4}{https://github.com/dreal/dreal4).

### Ubuntu / Mac
Currently, we do not support the Python API for dReal. Instead, it is required that dReal is installed on the desktop. Please follow the installation guide as mentioned before. 

In the options passed to the tool, it is required to declare the path of dReal. On Ubuntu, this is e.g.:
```
DREAL_VERSION=4.20.12.1 # your dReal version
/opt/dreal/${DREAL_VERSION}/bin/dreal
```

### Windows
Unfortunately, dReal is build for Linux/Mac. Luckily, dReal is provided as docker image.
* Install Docker https://hub.docker.com/editions/community/docker-ce-desktop-windows/
* Open the Windows power shell and run:
    ```
    docker pull dreal/dreal4
    docker run --rm dreal/dreal4 dreal -h  # Run "dreal -h" 
    ```
    If this displays the help of dReal: great, it is working.
* Create a folder that will be shared with docker: here the SMT files will be written/stored. E.g. the folder that will be shared has the path `D:/docker$\_$connect'.

