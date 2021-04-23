# F4CS

(Formal) Correct-by-Construction Closed-form Controller Synthesis.

# Installation
This tool is still in beta. The project is not yet available on PyPi. Instead, download the project and run
```
pip install -e <path>
```
This will automatically install the required packages. When considering non-polynomial systems or solutions, verification is done via dReal, which requires special installation. See the next section for more details.

## Install dReal 
Verification of the LF-like conditions is done using an SMT solver. In the case of non-polynomial conditions, our tool relies on the SMT solver dReal. The installation guide of dReal 4 can be found [here](https://github.com/dreal/dreal4}{https://github.com/dreal/dreal4).

**To use dReal in our tool, some system specific paths need to be set**. See the [wiki](https://github.com/CFVerdier/F4CS/wiki/SMT-solvers) for more information.

### Ubuntu / Mac
Currently, we do not support the Python API for dReal. Instead, it is required that dReal is installed on the desktop. Please follow the installation guide as mentioned before. 

In the options passed to the tool, it is required to declare the path of dReal. On Ubuntu, this is e.g.:
```
DREAL_VERSION=4.20.12.1 # your dReal version
/opt/dreal/${DREAL_VERSION}/bin/dreal
```

### Windows
Unfortunately, dReal is built for Linux/Mac. Luckily, dReal is provided as docker image.
* Install Docker https://hub.docker.com/editions/community/docker-ce-desktop-windows/
* Open the Windows power shell and run:
    ```
    docker pull dreal/dreal4
    docker run --rm dreal/dreal4 dreal -h  # Run "dreal -h" 
    ```
    If this displays the help of dReal: great, it is working.
* Create a folder that will be shared with docker: here the SMT files will be written/stored. E.g. the folder that will be shared has the path `D:/docker$\_$connect'.

# Examples
Example files are included in the `example` folder. For more details, visit the [example wiki](https://github.com/CFVerdier/F4CS/wiki/Examples)

# Further reading
For more details, please visit the wiki:
https://github.com/CFVerdier/F4CS/wiki

# Non-frequent unimportant questions
- _Why F4CS, and not 4CS?_ Strictly speaking, the 'Formal' is superfluous. However, Python packages should not start with numbers.
- _How do you pronounce F4CS?_ I pronounce it as 'forsees'.

