# -*- coding: utf-8 -*-


import sympy
import numpy as np
import collections
import os, os.path
import subprocess

sym_dict = collections.OrderedDict({"Add":"+","Mul":"*","Pow":"^",
           "StrictGreaterThan":">","GreaterThan":">=",
           "StrictLessThan":"<","LessThan":"<=",
           "And":"and", "Or":"or","Not":"not","Max":"max","Abs":"abs","Half":"0.5"})

sym_dict_exceptions = collections.OrderedDict({"Max":"max","Min":"min"})

num_dict =["Float","Integer","Zero","One","Symbol","NegativeOne"]


#TODO Rational numbers not supported
#TODO: error in dReal not captured: everything in violation key. This is because the powershell will run well anyway, hence no timeout errors like usual.
#TODO: booleans (for example, if true, this is not translated) 

def symname2lisp(fun): 
    """Translator function between symbolic python and SMT2 function names.
    Part of the sym2lisp translator""" 
    expr = fun
    for word in sym_dict:
        expr = expr.replace(word,sym_dict[word])
    return expr

def sym2lisp(expr):
    """Translator function between symbolic python expression and SMT2 expression"""
    name =expr.func.__name__
    if name in num_dict:
        sform = str(expr)
    elif name in sym_dict_exceptions:
        sform = "("+ symname2lisp(name)
        sform = sform + " " + sym2lisp(expr.args[0])
        for arg in expr.args[1:-1]:
            sform = sform + " (" + symname2lisp(name) +" "+ sym2lisp(arg)
        sform = sform + " " + sym2lisp(expr.args[-1])
        for arg in expr.args[1:]:
            sform = sform + ")"
    else:
        sform = "("+ symname2lisp(name)
        for arg in expr.args:
            sform = sform + " " +sym2lisp(arg)
        sform = sform + ")"
    return sform


def writeSMT2file(symbolic_expr,domain,symbol_list,dReal_precision,file_path,file_name="file.smt2"):
    """ Create and SMT file to check whether inequality symbolic_expr is satisfied using dReal."""
    
    #Check file suffix
    if not file_name.endswith('.smt2'):
        file_name = file_name + '.smt2'
    #Write settings in a string
    text = "(set-logic QF_NRA)\n"
    #text = text + "(set-info :precision "+ str(dReal_precision)+" )\n"
    #add declaration of constants
    for var in symbol_list:
        text = text + "(declare-fun " + str(var) +  " () Real)\n"
    #define domain
    text = text + "(assert "+ sym2lisp(domain) + ")\n"
    #Add inequality
    text =  text + "(assert (not  "+ sym2lisp(symbolic_expr) + "))\n"
    #check satisfiability
    text =  text + "(check-sat)\n"
    text = text + "(exit)"
    
    with open(os.path.join(file_path,file_name),'w+') as f:
        f.write(text)
        
    print("SMT2 File exported at " + os.path.join(file_path,file_name))
    return 
        
def call_dReal(file_path,file_name ="file.smt2",dreal_precision = 0.001,t_max= None):
    """ Call dReal and return a dictionary containing whether the inequality
    is satisfied, if not, the violations stored as a string, and a time-out flag
    """
    # Check file suffix
    if not file_name.endswith('.smt2'):
        file_name = file_name + '.smt2'
    
    #Initialize results
    result ={}
    print("Calling dReal")
    try:
        #dReal
        
        #Old method: broken 
        #client = docker.from_env()
        #vol = {file_path:{'bind':"/data"}}
        #outputdReal = client.containers.run("dreal/dreal4", "dreal data/"+file_name +" --model",volumes = vol,auto_remove = True).decode("utf-8")
            
        if t_max == None:
            outputdReal=subprocess.check_output(['powershell.exe',"docker run -v " + file_path + ":/data --rm dreal/dreal4 dreal data/" +file_name +" --model"],shell =True,stderr=subprocess.STDOUT).decode("utf-8")
        else:
            outputdReal=subprocess.check_output(['powershell.exe',"docker run -v " +file_path + ":/data --rm dreal/dreal4 dreal data/" +file_name +" --model"],shell =True,stderr=subprocess.STDOUT,timeout=t_max).decode("utf-8")
    except Exception:
        outputdReal = 'time-out'
        result['time-out']= True
        print("dReal time-out ({} seconds or other unexpected result.)".format(t_max))
    #Process data: unsat = 1, delta-sat = 0 (unsat proofs the inequality)
    if outputdReal == 'time-out':
        result['sat'] = False
        result['violation'] =  None
    elif outputdReal == 'unsat\n':
        result['sat']=True
        result['time-out'] = False
        result['violation'] = None
        print('Inequality Satisfied')
    else:
        result['sat'] = False
        result['time-out']= False
        result['violation'] =  outputdReal
        print('Inequality not verified')
    return result

def readViolations(dreal_result,symbol_list):
    """ Read violations. Note that accuracy might be lost due to the use of eval()"""
    string = dreal_result
    for var in symbol_list:
        string = string.replace(str(var),'"'+str(var)+'"')
    array = string.splitlines()
    array = array[1:]
    #Store point as a dictionary
    viol_dict ={}
    for entry in array:
        d =eval('{'+entry +'}')
        viol_dict.update(d)
    #Take the first vertex as counter example
    new_sample = ()
    for var in symbol_list:
#        new_sample = new_sample +(viol_dict[str(var)][0],)
        new_sample = new_sample +(np.mean(viol_dict[str(var)]),)
    return new_sample

def dReal_verify(symbolic_expr,domain,symbol_list,dReal_precision,file_path,file_name ="file.smt2",time_out = None):
    """ Main wrapper function"""
    # Check file suffix
    if not file_name.endswith('.smt2'):
        file_name = file_name + '.smt2'
    #Write SMT file
    writeSMT2file(symbolic_expr,domain,symbol_list,dReal_precision,file_path,file_name)
    #Call dReal
    result = call_dReal(file_path,file_name,dreal_precision = dReal_precision,t_max= time_out)
    #Read violations
    if result['violation']!= None and result['time-out'] == False:
        result['violation'] = readViolations(result['violation'],symbol_list)
        print("Violation at"+ str(result['violation']))
    return result

def demo():
    # Define path variables
    path = 'e:/docker_connect'
    
    #Declare symbolic variables
    varlist = x, y, z = sympy.symbols("x y z")
    #Specify domain
    dom = sympy.And(x>= -1, x<= 1, y>= -1, y<= 1,z>= -1, z<= 1)
    #Declare symbolic inequality
    F = x*y+1 >= sympy.cos(x/3.)+x/2.
    
    #print(sym2lisp(F))
    
    #varlist =  x1, x2, u1, t1, a1, a2 = sympy.symbols("x1 x2 u1 t1 a1 a2")
    #textstring ='delta-sat with delta = 0.001\nx1 : [-0.5, -0.4999977863608696116]\nx2 : [-0.4342024765383196705, -0.4341866642512871022]\nu1 : [-3.371923588897039803, -3.371893035553346962]\nt1 : [0, 8.382643932090520633e-07]\na1 : [-0.01000000000000000021, -0.009997748906802021371]\na2 : [-0.1600000000000000033, -0.1599857142871917715]'
    #res = readViolations(textstring,varlist)
    
    result =dReal_verify(F,dom,varlist,0.01,path)
    return result

#result = demo()