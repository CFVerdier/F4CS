#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 15 12:40:38 2019

@author: Cees Verdier
"""

import numpy as np
import sympy as sp
import dReal_communication_Windows as dReal

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
    def __init__(self,S):
        
        if isinstance(S,RWS):   #Check if it is a RWS spec. TODO replace with e.g. a dictionary to switch modes 
            #HARDCODED CONTROLLER AND CLBF
            x1, x2 = S.var
            self.k_sym = sp.Array([-1*x1-2*x2])
            self.V_sym = x1**2+x2**2+x1*x2-78
            
            ##TODO:make dependent on active mode
            
            #Make a Controller and LBF function
            self.k_fun = sp.lambdify([S.var],self.k_sym,"numpy")
            self.V_fun =  sp.lambdify([S.var],self.V_sym,"numpy")
    
            #Compute the derivatives and make into functions
            self.dV_sym = sp.diff(self.V_sym,[S.var])
            self.dV_fun =  sp.lambdify([S.var],self.dV_sym,"numpy")
            self.dtV_sym = sp.Matrix(self.dV_sym).dot(S.f_sym).subs(zip(S.input,self.k_sym))
        

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
        self.var = variables
        self.input = inputs
        self.options = options
        
        #Make functions of the dynamics
        self.f_sym = f_sym
        self.f_fun =  sp.lambdify([self.var,self.input],self.f_sym,"numpy")
        
        self.n = len(self.var)
        self.m = len(self.input)
        
        self.rng = np.random.default_rng() #Initialize random generator

    def sample_set(self,setList,numsamp):
        Array = np.array(setList)
        samples = (Array[:,1]-Array[:,0])*self.rng.random((numsamp,self.n))+Array[:,0]
        return samples
 
    def sample_set_complement(self,OuterList,InnerList,numsamp):
        """ Sample a set of the form Outer\Inner"""
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
        
        Slist = self.options['Slist']
        Ilist = self.options['Ilist']
        Olist = self.options['Olist']
        numsamp = self.options.get('numsamp',100)  #Default: 100
        self.c = self.options.get('c',0.01)  #arbitrary small constant. Default 0.01
        self.gamma = self.options.get('gamma',0.01)  #decrease of the LBF. Default 0.01
        
        #Create an inflated safe set to create a conservative boundary set
        rdelta = self.options.get('rdelta',0.01) #Default =0.01
        Rlist = [[Slist[i][0]-rdelta,Slist[i][1]+rdelta] for i in range(0,self.n)]
        
        #Create sample sets
        self.Idata = self.sample_set(Ilist,numsamp)
        self.dSdata = self.sample_set_complement(Rlist, Slist, numsamp)
        self.SnOdata = self.sample_set_complement(Slist, Olist,numsamp)
     
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
        self.con1set =  Iset
        self.con2set = clRnSset
        self.con3set = SnOset
    
    #Define the fitness function for RWS
    def sample_fitness(self,solution):
        #Define pointwise functions for each LBF condition
        def fit1(x):
            return np.minimum(-solution.V_fun(x),0.0)
        def fit2(x):
            return np.minimum(solution.V_fun(x)-self.c,0.0)
        def fit3(x):
            dtV =np.dot(solution.dV_fun(x),self.f_fun(x,solution.k_fun(x)))[0]
            return np.minimum(np.maximum(solution.V_fun(x)-self.c,-dtV-self.gamma),0.0)
        
        fit1_data = np.array([fit1(point) for point in self.Idata])
        fit2_data = np.array([fit2(point) for point in self.dSdata])
        fit3_data = np.array([fit3(point) for point in self.SnOdata])
        
        fit1_val = 1/(1+np.linalg.norm(fit1_data))
        fit2_val = 1/(1+np.linalg.norm(fit2_data))
        fit3_val = 1/(1+np.linalg.norm(fit3_data))
        
        w2 = np.floor(fit1_val)
        w3 = np.floor(w2*fit2_val)
        
        return (fit1_val + w2*fit2_val + w3*fit3_val)/3
        
    def verify(self,solution):     
        # call dReal
        path = self.options['path']
        dprecision = self.options.get('dprecision',0.01) #Defaulit 0.01
        
        #create the conditions to verify
        con1 = solution.V_sym <= 0
        con2 = solution.V_sym > 0
        con3 = sp.Or(solution.V_sym >0, solution.dtV_sym <= -self.gamma)

        #TODO: clean up the result
        result1 =dReal.dReal_verify(con1,self.con1set,self.var,dprecision,path,'con1')
        result2 =dReal.dReal_verify(con2,self.con2set,self.var,dprecision,path,'con2')
        result3 =dReal.dReal_verify(con3,self.con3set,self.var,dprecision,path,'con3')
        
        #ToDo: extract counter examples and append to set
       
  
### DEMO ################
        
def demo():
    
    #Variable Declaration
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
               'numsamp':100,  #Number of samples
               'rdelta':0.01,   #Inflation of the boundary
               'gamma':0.01,    #(arbitrary) decrease of the LF
               'c':0.01,        #(arbitrary) nonnegative parameter (see manual)
               'path': path,
               'dprecision': 0.01}    #Path where the SMT files will be stored
    
    #Initialize specification
    S = RWS(var_list,input_list,f_sym,options)
    #Create an (hardcoded) individual
    ind = Solution(S)

    fitness = S.sample_fitness(ind)
    print(fitness)
    #verify
    S.verify(ind)
    
demo()