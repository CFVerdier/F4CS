#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 12:40:38 2019

@author: Cees Verdier
"""

import numpy as np
import sympy as sp
import dReal_communication_Windows as dReal
import cma

#import matplotlib.pyplot as plt
#from sympy.plotting import plot as splt



#TODO use the class below for new individuals
class Solution:
    """ Stores and manipulates a candidate solution for RWS
    
    Constructor arguments:
        S (spec): specification from the RWS class
        
    Optimal arguments:
        none
        
    Attributes:
        TODO
    
    Methods:
        none
    """
    #TODO: variable in the entities that it can has: LF, controller, auxillary functions
    def __init__(self,S):
        
        if isinstance(S,RWS):   #Check if it is a RWS spec. 
        #TODO replace with e.g. a dictionary to switch modes 
            #HARDCODED CONTROLLER AND CLBF
            self.S = S
            self.par = [-1,-2,-78]
            self.par_len = len(self.par)
            self.p = sp.symbols('p0:{}'.format(self.par_len))
            p = self.p
            x = S.var
            self.k_sym_p = sp.Array([p[0]*x[0]+p[1]*x[1]])
            self.V_sym_p = x[0]**2+x[1]**2+x[0]*x[1]+p[2]
            
            self.k_sym = None
            self.V_sym = None
            self.k_fun = None
            self.V_fun = None
            self.dV_sym = None
            self.dVfun = None
            self.dtV_sym = None
            
            #TODO: store these values automatically
            self.sample_fitness = 0
            self.fitness = 0
            
            #substitute the parameter vector and create the functions
            self.substitute_parameters(self.par)

            
            #TODO: Remove spec specific data from the solution.

    def substitute_parameters(self,z):
        self.par = z
        self.k_sym = self.k_sym_p.subs([(self.p[i],self.par[i]) for i in range(self.par_len)])
        self.V_sym = self.V_sym_p.subs([(self.p[i],self.par[i]) for i in range(self.par_len)])
        self.make_functions()
        
    def make_functions(self):    
        #Make a Controller and LBF function
        self.k_fun = sp.lambdify([self.S.var],self.k_sym,"numpy")
        self.V_fun =  sp.lambdify([self.S.var],self.V_sym,"numpy")

        #Compute the derivatives and make into functions
        self.dV_sym = sp.diff(self.V_sym,[self.S.var])
        self.dV_fun =  sp.lambdify([self.S.var],self.dV_sym,"numpy")
        self.dtV_sym = sp.Matrix(self.dV_sym).dot(self.S.f_sym).subs(zip(self.S.input,self.k_sym))
               

class Spec:
    """ Base class of a specification.
    
    Constructor arguments:
        variables: tuple of symbolic state variables.
            Example: sympy.symbols('x1, x2')
        inputs: tuple of symbolic input variables.
            Example: sympy.symbols('u1')
        options (dictionary): dictionary with all the settings relating to the specification. 
        
    Attributes:
        var (list): User defined state variables
        input (list): User defined system inputs
        n (int): number system variables
        m (int): number of input variables

        
    Methods
        sample_set(interval, numper of samples): sample a number of random samples from a hyperrectangle, spanned by the interval. Returns an array of float64
        sample_set_complement(OuterList,InnerList,numsamp): Samples a set Outerset\Innerset, where numsamp denotes the number of samples, Outerlist and Innerlist the interval of the Outerset and Innerset respectively
    
    Child classes:
        RWS: Reach while stay specification for continuous time systems
    
    """
    #constructor
    def __init__(self,variables,inputs,f_sym,options):
        """Constructor"""
        
        self._number_conditions = 1 #number of conditions     
        self.var = variables
        self.input = inputs
        self.options = options
        self.dprecision = self.options.get('dprecision',0.01) #Defaulit 0.01
        self.numsamp = self.options.get('numsamp',100)  #Default: 100
        self.max_samples = self.options.get('max_samp',1000) #Defaulit 0.01
        self.epsilon = self.options.get('epsilon',0.1) #Robustness buffer for the sample-based filtness. Default 0.1
        
        #Make functions of the dynamics
        self.f_sym = f_sym
        self.f_fun =  sp.lambdify([self.var,self.input],self.f_sym,"numpy")
        
        self.n = len(self.var)
        self.m = len(self.input)
        
        self.conditions =  None
        self.condition_set = (True,)
        self.verification_result = [None]*self._number_conditions
        
        self.rng = np.random.default_rng() #Initialize random generator

    def sample_set(self,setList,numsamp =None):
        #If there is no argument of numpsamp, use atribute
        if numsamp is None:
            numsamp = self.numsamp
        Array = np.array(setList)
        samples = (Array[:,1]-Array[:,0])*self.rng.random((numsamp,self.n))+Array[:,0]
        return samples
 
    def sample_set_complement(self,OuterList,InnerList,numsamp = None):
        """ Sample a set of the form Outer\Inner"""
        #If there is no argument of numpsamp, use atribute
        if numsamp is None:
            numsamp = self.numsamp
        OutArray = np.array(OuterList)
        InArray = np.array(InnerList)
        samples = (OutArray[:,1]-OutArray[:,0])*self.rng.random((numsamp,self.n))+OutArray[:,0]
        
        # Cut out a set
        #select a dimension
        sdims = self.rng.integers(self.n,size=numsamp)
        #percentage
        perc = 2*self.rng.random(numsamp)-1
        for i in range(numsamp):
            sdim = sdims[i]
            if perc[i]>=0:
                samples[i,sdim] = InArray[sdim,1]+ perc[i]*( OutArray[sdim,1]- InArray[sdim,1])
            else:
                samples[i,sdim] = InArray[sdim,0]-perc[i]*( OutArray[sdim,0]- InArray[sdim,0])
        return samples
    
    
    def parameter_fitness(self,x,solution):
        """Given a parameter vector x, compute the fitness """
        solution.substitute_parameters(x) #Substitute parameters
        return -1*self.sample_fitness(solution)
    
    def sample_fitness(self,solution):        
        return 1.0
    
    def SMT_fitness(self,solution):
        #Call SMT solver to verify the conditions.
        self.verify(solution)
        #Total SMT fitness is equal to the number of true statements, normalized to number of conditions
        return (np.sum([int(n['sat']) for n in self.verification_result]))/self._number_conditions
    
    #TODO: if we store the sample fitness of a candidate after parameter optimization
    def fitness(self,solution):
        return (self.sample_fitness(solution)+self.SMT_fitness(solution))/2
          
    def create_conditions(self,solution):     
        self.conditions = (True,)

    def verify(self,solution):
        #Create the conditions to verify
        self.create_conditions(solution)   
        #For each condition, verify the condition with dReal
        for i in range(0,self._number_conditions):
            self.verification_result[i]=dReal.dReal_verify(self.conditions[i],self.condition_set[i],self.var,self.dprecision,self.options['path'],"con{}".format(i+1))
            #If there is a counter example, sample the data and append.
            if isinstance(self.verification_result[i]['violation'],tuple):
                self.data_sets[i] =np.append(self.data_sets[i],[self.verification_result[i]['violation']],axis =0)
                #Saturate w.r.t. the maximum number of samples. (FIFO)
                if len(self.data_sets[i]) > self.max_samples:
                    self.data_sets[i] = self.data_sets[i][-self.max_samples:]
    
 
class RWS(Spec):
    """Represents the RWS specification and implements the corresponding fitness function and verification.
    Subclass of Spec.
    
    Constructor arguments:
    variables: tuple of symbolic state variables.
        Example: sympy.symbols('x1, x2')
    inputs: tuple of symbolic input variables.
        Example: sympy.symbols('u1,')
    f_sym (sympy expression): Symbolic expression of the system dynamics
    options (dictionary): dictionary with all the settings relating to the specification. 
        Required options:
            Slist, Ilist, Olist: a list of the lower and upperbounds of the safe set, initial set and goal set
            path: path where the SMT files are stored.
        Optional
            numsamp:        Number of samples. Default 100
            rdelta:         Inflation of the boundary. Default: 0.01
            gamma:          (Arbitrary) decrease of the LF. Default :0.01, 
            c:              (Arbitrary) nonnegative parameter (see manual). Default: 0.01
            dprecision:     Precision of dReal. Default: 0.01
            

    """
    def __init__(self,variables,inputs,f_sym,options):
        
        #Call the __init__ function of the Spec parent class first.
        Spec.__init__(self,variables,inputs,f_sym,options)
        
        self._number_conditions = 3 #number of RWS conditions
        
        Slist = self.options['Slist']
        Ilist = self.options['Ilist']
        Olist = self.options['Olist']
        self.c = self.options.get('c',0.01)  #arbitrary small constant. Default 0.01
        self.gamma = self.options.get('gamma',0.01)  #decrease of the LBF. Default 0.01
        
        #Create an inflated safe set to create a conservative boundary set
        rdelta = self.options.get('rdelta',0.01) #Default =0.01
        Rlist = [[Slist[i][0]-rdelta,Slist[i][1]+rdelta] for i in range(0,self.n)]
        
        #Create sample sets
        Idata = self.sample_set(Ilist)
        dSdata = self.sample_set_complement(Rlist, Slist)
        SnOdata = self.sample_set_complement(Slist, Olist)
        self.data_sets = [Idata,dSdata,SnOdata]
            
        #Create symbolic domains for SMT solver
        Sset =sp.And()
        Rset = sp.And()
        Iset = sp.And()
        Oset =sp.And()
        for i in range(0,self.n): 
            Sset = sp.And(Sset,sp.And(self.var[i]>= Slist[i][0],self.var[i]<=Slist[i][1]))
            Iset = sp.And(Iset,sp.And(self.var[i]>= Ilist[i][0],self.var[i]<=Ilist[i][1]))
            Oset = sp.And(Oset,sp.And(self.var[i]>= Olist[i][0], self.var[i]<=Olist[i][1]))
            SnOset =sp.And(Sset,sp.Not(Oset))
            #Create the closure of R\S to create a conservative boundary of S. 
            Rset = sp.And(Rset,sp.And(self.var[i]>= Rlist[i][0],self.var[i]<=Rlist[i][1]))
            Sopenset = sp.And(Sset,sp.And(self.var[i]> Slist[i][0],self.var[i]<Slist[i][1]))
            clRnSset = sp.And(Rset,sp.Not(Sopenset))

        self.condition_set = (Iset,clRnSset,SnOset) 
        self.conditions = None
        self.verification_result = [None]*self._number_conditions
    
    #Define the fitness function for RWS
    def sample_fitness(self,solution):
        #Define pointwise functions for each LBF condition
        def fit1(x):
            return np.minimum(-solution.V_fun(x)-self.epsilon,0.0)
        def fit2(x):
            return np.minimum(solution.V_fun(x)-self.c-self.epsilon,0.0)
        def fit3(x):
            dtV =np.dot(solution.dV_fun(x),self.f_fun(x,solution.k_fun(x)))[0]
            return np.minimum(np.maximum(solution.V_fun(x)-self.c+self.epsilon,-dtV-self.gamma-self.epsilon),0.0)
        
        #TODO: optimize
        fit1_data = np.array([fit1(point) for point in self.data_sets[0]])
        fit2_data = np.array([fit2(point) for point in self.data_sets[1]])
        fit3_data = np.array([fit3(point) for point in self.data_sets[2]])
        
        fit1_val = 1/(1+np.linalg.norm(fit1_data))
        fit2_val = 1/(1+np.linalg.norm(fit2_data))
        fit3_val = 1/(1+np.linalg.norm(fit3_data))
        
        w2 = np.floor(fit1_val)
        w3 = np.floor(w2*fit2_val)
        
        return (fit1_val + w2*fit2_val + w3*fit3_val)/3
        
    def create_conditions(self,solution):     

        #create the conditions to verify
        con1 = solution.V_sym <= 0
        con2 = solution.V_sym > 0
        con3 = sp.Or(solution.V_sym >0, solution.dtV_sym <= -self.gamma)
        
        self.conditions = (con1, con2, con3)

class Synthesis:
    """ Synthesizes a controller and certificate function for a given specification
    
    Constructor arguments:
        specification (spec): specification
        options: dictionary with the desired specifications.
        
    Optimal arguments:
        TODO
        
    Attributes:
        TODO
    
    Methods:
        TODO
    """
    #TODO: variable in the entities that it can has: LF, controller, auxillary functions
    def __init__(self,specification,options):
        self.method = self.options.get('Method','Template')
        self.spec = specification
        
        
class TBS(Synthesis):
    """Template based controller and certificate function synthesis. Subclass of synthesis
    
        Constructor arguments:
        specification (spec): specification
        options: dictionary with the desired specifications.
        
    Optimal arguments:
        TODO
        
    Attributes:
        TODO
    
    Methods:
        TODO
    
    """

### DEMO ################
        
def demo():
    var_list = x1,x2= sp.symbols('x1,x2')
    input_list = u1, =sp.symbols('u1,') 
      
    #Dynamics
    f_sym =  sp.Matrix([x2,u1])  #Column vector
    
    Slist = [[-15,15],[-15,15]]
    Ilist = [[-5,5],[-5,5]]
    Olist = [[-1,1],[-1,1]]
    
    #path where the SMT files will be stored
    path = 'e:/docker_connect/data'
    
    options = {'Slist':Slist,   #Interval list of the safe set
               'Ilist':Ilist,   #Interval list of the initial set
               'Olist':Olist,   #Interval list of the goal set
               'numsamp':100,  #Number of (initial) samples
               'max_samp': 300, #Maximum number of samples (when adding violations)
               'rdelta':0.01,   #Inflation of the boundary
               'gamma':0.01,    #(arbitrary) decrease of the LF
               'c':0.01,        #(arbitrary) nonnegative parameter (see manual)
               'path': path,
               'epsilon':0.1,  #Robustness buffer for the sample-based fitness
               'dprecision': 0.01}    #Path where the SMT files will be stored
    
    #Initialize specification
    spec = RWS(var_list,input_list,f_sym,options)
    #Create an (hardcoded) individual
    ind = Solution(spec)
    
    #Give the individual arbitrary new parameters for testing
    ind.par = [1,1,1]
    ind.substitute_parameters(ind.par)
    
    sigma0 = 0.5
    cma.fmin(spec.parameter_fitness,ind.par,sigma0,args={ind,},options = {'verbose':-9})
    spec.verify(ind)
        
  
#demo()


## debug

var_list = x1,x2= sp.symbols('x1,x2')
input_list = u1, =sp.symbols('u1,') 
  
#Dynamics
f_sym =  sp.Matrix([x2,u1])  #Column vector

Slist = [[-15,15],[-15,15]]
Ilist = [[-5,5],[-5,5]]
Olist = [[-1,1],[-1,1]]

#path where the SMT files will be stored
path = 'e:/docker_connect/data'

options = {'Slist':Slist,   #Interval list of the safe set
           'Ilist':Ilist,   #Interval list of the initial set
           'Olist':Olist,   #Interval list of the goal set
           'numsamp':100,  #Number of (initial) samples
           'max_samp': 300, #Maximum number of samples (when adding violations)
           'rdelta':0.01,   #Inflation of the boundary
           'gamma':0.01,    #(arbitrary) decrease of the LF
           'c':0.01,        #(arbitrary) nonnegative parameter (see manual)
           'path': path,
           'epsilon':0.1,  #Robustness buffer for the sample-based fitness
           'dprecision': 0.01}    #Path where the SMT files will be stored

#Initialize specification
spec = RWS(var_list,input_list,f_sym,options)
#Create an (hardcoded) individual
ind = Solution(spec)

#Give the individual arbitrary new parameters for testing
ind.par = [1,1,1]
ind.substitute_parameters(ind.par)

sigma0 = 0.5
cma.fmin(spec.parameter_fitness,ind.par,sigma0,args={ind,},options = {'verbose':1,'seed': 1})
#spec.verify(ind)

## Timing experiment 
import timeit
ind.par = [1,1,1]
ind.substitute_parameters(ind.par)

time = timeit.timeit('cma.fmin(spec.parameter_fitness,ind.par,sigma0,args={ind,},options = {"verbose":-9,"seed": 1})','from __main__ import spec, ind, sigma0, cma',number =10)
print(time)