# F4CS

(Formal) Correct-by-Construction Closed-form Controller Synthesis.

A tool for correct-by-design synthesis of closed-form analytic controllers. It utilizes a counterexample-guided inductive synthesis framework to co-design controllers and certificate functions. Internally, candidate controllers are proposed and verified using an SMT solver. If the candidate controller does not satisfy the specification, a counterexample is extracted, which is then used to refine the proposed controller. This procedure is repeated until a satisfying controller is found, or a maximum number of iterations is met.

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
* Create a folder that will be shared with docker: here the SMT files will be written/stored. E.g. the folder that will be shared has the path `D:/docker_connect'.
* Note: When relying on dReal and running on Windows, make sure Docker is running when performing the synthesis. 

# Examples
Example files are included in the `example` folder. For more details, visit the [example wiki](https://github.com/CFVerdier/F4CS/wiki/Examples)

# Further reading
For more details, please visit the wiki:
https://github.com/CFVerdier/F4CS/wiki

# Related literature
This tool is based on the following work:
  <ul>
      <li> C.F. Verdier and M. Mazo, Jr.
        <i> Formal Controller Synthesis for Hybrid Systems Using Genetic Programming. </i>
        <a href ="https://arxiv.org/abs/2003.14322">arXiv preprint arXiv:2003.14322. </a> (2020)
      <li> C.F. Verdier, R. Babuska, B. Shyrokau and M. Mazo Jr.
          <i> Near Optimal Control with Reachability and Safety Guarantees. </i>
        <a href = "https://doi.org/10.1016/j.ifacol.2019.09.146">5th IFAC Conference on Intelligent Control and Automation Sciences (ICONS). </a> (2019)
      <li> C.F. Verdier and M. Mazo. 
        <i> Formal Controller Synthesis of Analytic Controllers for Sampled-data Systems via Genetic Programming. </i>
        <a href = "https://doi.org/10.1109/CDC.2018.8619121">57th IEEE Conference on Decision and Control (CDC). </a> (2018)
      <li> C.F. Verdier and M. Mazo.
        <i> Formal Controller Synthesis via Genetic Programming. </i> 
        <a href = "https://doi.org/10.1016/j.ifacol.2017.08.1362">IFAC-PapersOnLine. </a> (2017)
    </ul>

# Non-frequent unimportant questions
- _Why F4CS, and not 4CS?_ Strictly speaking, the 'Formal' is superfluous. However, Python packages should not start with numbers.
- _How do you pronounce F4CS?_ I pronounce it as 'forsees'.

