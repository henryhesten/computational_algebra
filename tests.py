# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 16:58:16 2018

@author: henry
"""


import numpy as np
from numpy import *
from matplotlib import *
from matplotlib.pyplot import *
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, mark_inset,InsetPosition
import numpy as np
import scipy.optimize as spo
import scipy.integrate as integrate
from scipy import constants, pi, exp , sqrt
#from scipy.special import factorial
import copy
import math
import scipy.signal
import re
import itertools

import sys,os

import abc

if(os.environ["PATH"][0]=="C" or os.environ["PATH"][0]=="D"):
    pathh = "D:/hesten-bec/"
    sys.path.insert(0,pathh)
    import code_.cubehelix
else:
    import cubehelix
    pathh = "/home/henry/storage/hesten-bec/"
    sys.path.insert(0,pathh)

    
sys.path.insert(0,pathh)
from code_.models.myCommon import *
from classes import *
#%%
if not op("q",["w"]) == op("q",["w"]):
    raise Exception("n1")

#%% addition test
term = addi([numb(10),numb(20)])
term2 = addi([numb(20),numb(10)])
if term == term2:
    raise Exception("ad1")
if not term == addi([numb(10),numb(20)]):
    raise Exception("ad2")

term3 = addi([addi([sym("q"),sym("w")]),addi([sym("e"),sym("r")]),sym("t")])
if term3.simplify() != addi([sym("e"),sym("q"),sym("r"),sym("t"),sym("w")]):
    raise Exception("ad3")
    
term4 = addi([sym("q"),sym("w"),mult([sym("q"),numb(10)]),sym("w")])
if term4.simplify() != addi([mult([numb(11),sym('q',[])]),mult([numb(2),sym('w',[])])]):
    raise Exception("ad4")
    
arr = [numb(10),sym("a"),sym("a",["l"]),sym("b"),op("q"),op("q",["l"]),op("w"),summ(sym("q",["i"]),["i",numb(0),sym("N")]),summ(sym("q",["j"]),["j",numb(0),sym("N")]),summ(sym("w",["i"]),["i",numb(0),sym("N")]),summ(op("z",["i"]),["i",numb(0),sym("N")]),trac(sym("h")),trac(sym("g")),trac(op("z")),trac(op("y")),mult([sym("a"),sym("b")]),mult([sym("c"),sym("d")]),mult([sym("r"),op("t")]),mult([op("t"),op("y")]),mult([op("u"),op("v")]),addi([sym("h"),sym("j")]),addi([sym("k"),sym("l")]),addi([op("y"),op("u")]),addi([op("u"),op("v")]), frac(sym("q"),sym("w")),frac(sym("y"),sym("z"))]
np.random.shuffle(arr)
term5 = addi(arr)
term5.doSort()
if term5 != addi([numb(10),sym('a',['l']),sym('a',[]),sym('b',[]),summ(sym('q',['i']),['i',numb(0),sym('N',[])]),summ(sym('w',['i']),['i',numb(0),sym('N',[])]),summ(sym('q',['j']),['j',numb(0),sym('N',[])]),frac(sym('q',[]),sym('w',[])),frac(sym('y',[]),sym('z',[])),addi([sym('h',[]),sym('j',[])]),addi([sym('k',[]),sym('l',[])]),mult([sym('a',[]),sym('b',[])]),mult([sym('c',[]),sym('d',[])]),op('q',['l']),op('q',[]),op('w',[]),trac(sym('g',[])),trac(sym('h',[])),trac(op('y',[])),trac(op('z',[])),summ(op('z',['i']),['i',numb(0),sym('N',[])]),addi([op('u',[]),op('v',[])]),addi([op('y',[]),op('u',[])]),mult([sym('r',[]),op('t',[])]),mult([op('t',[]),op('y',[])]),mult([op('u',[]),op('v',[])])]):
    raise Exception("ad5")
    
term6 = addi([mult([sym('a',[]),sym('x',[])]),mult([sym('a',[]),sym('y',[])]),mult([sym('b',[]),sym('x',[])]),mult([sym('b',[]),sym('y',[])])])
fac = addi([sym("a"),sym("b")])
if term6.factorise(fac) != addi([sym('x',[]),sym('y',[])]):
    raise Exception("ad6")
    
term7 = mult([fac,addi([frac(sym("x"),sym("a")),frac(sym("y"),sym("b"))])])
term7b = term7.simplify()
if term7b != addi([sym('x',[]),sym('y',[]),frac(mult([sym('a',[]),sym('y',[])]),sym('b',[])),frac(mult([sym('b',[]),sym('x',[])]),sym('a',[]))]):
    raise Exception("ad7")
term7F = term7b.factorise(fac)
if term7F != addi([frac(sym("x"),sym("a")),frac(sym("y"),sym("b"))]):
    raise Exception("ad7B")
    
term8 = addi([frac(sym("a"),sym("b")),frac(sym("q"),sym("w"))])
if term8.mergeFrac() != frac(addi([mult([sym('a',[]),sym('w',[])]),mult([sym('b',[]),sym('q',[])])]),mult([sym('b',[]),sym('w',[])])):
    raise Exception("ad8")
#%%
termM1 = mult([op("a",["w"]),op("a^T",["w"])])
if termM1.simplify() != addi([numb(1),mult([op('a^T',['w']),op('a',['w'])])]):
    raise Exception("m1")

termM2 = mult([op("b"),addi([op("c1"),op("c2")]),addi([op("d1"),op("d2")])])
if termM2.simplify() != addi([mult([op('b',[]),op('c1',[]),op('d1',[])]),mult([op('b',[]),op('c1',[]),op('d2',[])]),mult([op('b',[]),op('c2',[]),op('d1',[])]),mult([op('b',[]),op('c2',[]),op('d2',[])])]):
    raise Exception("m2")

termM3 = mult([  op("a",["w"]), op("sp^T",["m"]), op("a^T",["w"]), op("sp",["n"]), op("a^T",["h"])  ])
if termM3.simplify() != addi([mult([numb(-2),delta(['m', 'n']),op('a^T',['h']),op('a^T',['w']),op('a',['w']),op('f',['m'])]),mult([numb(-2),delta(['m', 'n']),delta(['w', 'h']),op('a^T',['w']),op('f',['m'])]),mult([numb(-2),delta(['m', 'n']),op('a^T',['h']),op('f',['m'])]),mult([delta(['m', 'n']),op('a^T',['h']),op('a^T',['w']),op('a',['w'])]),mult([delta(['m', 'n']),delta(['w', 'h']),op('a^T',['w'])]),mult([delta(['m', 'n']),op('a^T',['h'])]),mult([delta(['w', 'h']),op('a^T',['w']),op('sp',['n']),op('sp^T',['m'])]),mult([op('a^T',['h']),op('a^T',['w']),op('a',['w']),op('sp',['n']),op('sp^T',['m'])]),mult([op('a^T',['h']),op('sp',['n']),op('sp^T',['m'])])]):
    raise Exception("m3")
    
termM4 = mult([op("a"),sym("b"),summ(sym("c",["m"]),["m",numb(0),sym("M")]),summ(op("d",["n"]),["n",numb(0),sym("N")])])
if termM4.simplify() != mult([sym('b',[]),summ(mult([sym('c',['m']),summ(mult([op('a',[]),op('d',['n'])]),['n',numb(0),sym('N',[])])]),['m',numb(0),sym('M',[])])]):
    raise Exception("m4")
    
termM5 = mult([op("b"),op("sp",["l"]),op("sp^T",["l"]),op("c"),op("sp^T",["l"]),op("sp",["l"]),op("d"),op("sp",["l"]),op("sp^T",["k"])])

termM6 = mult([addi([numb(1),mult([numb(-2),op('f',['l'])])]),op('d',[]),op('sp',['l'])])
if termM6.simplify() != addi([mult([numb(-2),op('f',['l']),op('d',[]),op('sp',['l'])]),mult([op('d',[]),op('sp',['l'])])]):
    raise Excpation("m5")

#%%
termD1 = delta(["w","w"])
if termD1.simplify() != numb(1):
    raise Exception("d1")
    
#%%
termT1 = trac(mult([addi([op("b"),op("c")]),op("d")]))
if termT1.simplify() != addi([trac(mult([op('b',[]),op('d',[])])),trac(mult([op('c',[]),op('d',[])]))]):
    raise Exception("t1")
    
termT2 = trac(termM4)
if termT2.simplify() != mult([sym('b',[]),summ(sym('c',['m']),['m',numb(0),sym('M',[])]),summ(trac(mult([op('a',[]),op('d',['n'])])),['n',numb(0),sym('N',[])])]):
    raise Exception("t2")
    
termT3 = trac(mult([op("a"),rho(),op("b")]))
if termT3.simplify() != trac(mult([op('b',[]),op('a',[]),rho()])):
    raise Exception("t3")


#%%
termS1 = summ(addi([sym("q",["m"]),sym("w",["m"])]),["m",numb(0),sym("N")])
if termS1.simplify() != addi([summ(sym('q',['m']),['m',numb(0),sym('N',[])]),summ(sym('w',['m']),['m',numb(0),sym('N',[])])]):
    raise Exception("S1")
    
termS2 = summ(mult([sym("q",["m"]),sym("w")]),["m",numb(0),sym("N")])
if termS2.simplify() != mult([sym('w',[]),summ(sym('q',['m']),['m',numb(0),sym('N',[])])]):
    raise Exception("S2")
    
termS3 = summ(sym("e",["r"]),["m",numb(1),sym("N")])
if termS3.simplify() != addi([mult([numb(-1),sym('e',['r'])]),mult([sym('N',[]),sym('e',['r'])])]):
    raise Exception("S3")

termS4 = mult([sym("q"),op("a"),sym("w"),summ(mult([op("b"),op("c")]),["m",numb(0),sym("M")]),sym("e"),op("d"),sym("r")])
if termS4.simplify() != mult([sym('e',[]),sym('q',[]),sym('r',[]),sym('w',[]),summ(mult([op('a',[]),op('b',[]),op('c',[]),op('d',[])]),['m',numb(0),sym('M',[])])]):
    raise Exception("S4")
#%%
termF1 = frac(frac(sym("q"),sym("w")),frac(sym("e"),sym("r")))
if termF1.simplify() != frac(mult([sym('q',[]),sym('r',[])]),mult([sym('e',[]),sym('w',[])])):
    raise Exception("f1")

termF2 = frac(mult([sym("q"),frac(sym("a"),sym("b"))]),mult([sym("w"),frac(sym("c"),sym("d"))]))
if termF2.simplify() != frac(mult([sym('a',[]),sym('d',[]),sym('q',[])]),mult([sym('b',[]),sym('c',[]),sym('w',[])])):
    raise Exception("f2")

termF3 = frac(mult([sym("q"),sym("w")]),mult([sym("e"),sym("w")]))
if termF3.simplify() != frac(sym('q',[]),sym('e',[])):
    raise Exception("f3")

termF4 = frac(mult([sym("q"),sym("w"),sym("e")]),mult([sym("e"),sym("w")]))
if termF4.simplify() != sym("q"):
    raise Exception("f4")
       
termF6 = mult([frac(op("q"),sym("w")),op("z"),frac(op("a"),sym("b"))])
if termF6.simplify() != frac(mult([op('q',[]),op('z',[]),op('a',[])]),mult([sym('b',[]),sym('w',[])])):
    raise Exception("f6")
    
termF7 = frac( addi([ frac( mult([sym("a"),sym("b")]),sym("c")), frac(sym("e"),sym("f"))]), addi([ frac( mult([sym("x"),sym("y")]),sym("z")), frac(sym("v"),sym("u"))]))
if termF7.mergeFrac() != frac(mult([addi([mult([mult([sym('a',[]),sym('b',[])]),sym('f',[])]),mult([sym('c',[]),sym('e',[])])]),mult([sym('z',[]),sym('u',[])])]),mult([addi([mult([mult([sym('x',[]),sym('y',[])]),sym('u',[])]),mult([sym('z',[]),sym('v',[])])]),mult([sym('c',[]),sym('f',[])])])):
    raise Exception("f7")
#%%
termE1 = equality(sym("q"),sym("w"))
if termE1.simplify({"preserveEq":False}) != equality(addi([sym('q',[]),mult([numb(-1),sym('w',[])])]),numb(0)):
    raise Exception("E1")

termE2 = equality(frac(sym("a"),sym("b")),frac(sym("c"),sym("d")))
if termE2.simplify({"preserveEq":False}) != equality(addi([mult([numb(-1),sym('b',[]),sym('c',[])]),mult([sym('a',[]),sym('d',[])])]),numb(0)):
    raise Exception("e2")

termE3 = equality(addi([frac(sym("a"),sym("b")),numb(1)]),addi([frac(sym("c"),sym("d")),numb(1)]))
if termE3.removeFrac() != equality(addi([mult([numb(-1),sym('b',[]),sym('c',[])]),mult([sym('a',[]),sym('d',[])])]),numb(0)):
    raise Exception("e3")
    
termE4 = equality(addi([mult([sym("a"),sym("x")]),mult([sym("a"),sym("y")])]),numb(0))
if termE4.factoriseNest(sym("a")) != equality(addi([sym('x',[]),sym('y',[])]),numb(0)):
    raise Exception("e4")

termE5 = summ(trac(mult([op("a",["m"]),op("b",["n"]),addi([sym("c",["k"]),sym("d",["m"])])])),["k",numb(0),sym("M")])
termE5b = summ(trac(mult([op("a",["z"]),op("b",["y"]),addi([sym("c",["x"]),sym("d",["w"])])])),["u",numb(0),sym("M")])
if termE5.equalMinusInds(termE5b,["u","v","w","x","y","z"]) != {'u': 'k', 'w': 'm', 'x': 'k', 'y': 'n', 'z': 'm'}:
    raise Exception("E5b")
termE5c = summ(trac(mult([op("a",["z"]),op("b",["u"]),addi([sym("c",["x"]),sym("d",["w"])])])),["u",numb(0),sym("M")])
if termE5.equalMinusInds(termE5c,["u","v","w","x","y","z"]) != False:
    raise Exception("E5c")

if term5.equalMinusInds(trac(op("a",["m"])),["m"]) != False:
    raise Exception("E5d")

#%%
termH1 = equality(addi([numb(1),sym("q"),sym("w^*"),op("a"),op("b^T"),mult([op("q^T"),sym("w"),op("e")])]), summ(trac(mult([op("q^T"),sym("w"),op("e")])),["m",numb(0),sym("M")]))
if termH1.HermConj() != equality(addi([numb(1),sym('q^*',[]),sym('w',[]),op('a^T',[]),op('b',[]),mult([op('e^T',[]),sym('w^*',[]),op('q',[])])]),summ(trac(mult([op('e^T',[]),sym('w^*',[]),op('q',[])])),['m',numb(0),sym('M',[])])):
    raise Exception("h1")

termH2 = frac(mult([op("q^T"),sym("w"),op("e")]),mult([sym("w"),sym("e")]))
if termH2.HermConj() != frac(mult([op('e^T',[]),sym('w^*',[]),op('q',[])]),mult([sym('e^*',[]),sym('w^*',[])])):
    raise Exception("h2")
#%%
termD1 = deriv( addi([sym("a"),sym("b")]),sym("t"))
if termD1.simplify() != addi([deriv(sym('a',[]),sym('t',[])),deriv(sym('b',[]),sym('t',[]))]):
    raise Exception("D1")

termD2 = deriv( mult([sym("a"),sym("b")]),sym("t"))
if termD2.simplify() != addi([mult([sym('a',[]),deriv(sym('b',[]),sym('t',[]))]),mult([sym('b',[]),deriv(sym('a',[]),sym('t',[]))])]):
    raise Exception("D2")

termD3 = deriv( frac(addi([mult([sym("a"),sym("b")]),sym("c")]),sym("e")),sym("t"))
if termD3.simplify() != frac(addi([mult([numb(-1),sym('a',[]),sym('b',[]),deriv(sym('e',[]),sym('t',[]))]),mult([numb(-1),sym('c',[]),deriv(sym('e',[]),sym('t',[]))]),mult([sym('a',[]),sym('e',[]),deriv(sym('b',[]),sym('t',[]))]),mult([sym('b',[]),sym('e',[]),deriv(sym('a',[]),sym('t',[]))]),mult([sym('e',[]),deriv(sym('c',[]),sym('t',[]))])]),mult([sym('e',[]),sym('e',[])])):
    raise Exception("D3")

#%%
termB1 = trac(op("a",["m"]))
termB1B = termB1.substituteMinusInds(trac(op("a",["z"])),trac(op("b",["z"])))
if termB1B != trac(op('b',['m'])):
    raise Exception("B1")

termB2 = summ( trac(mult([op("a",["m"]),op("b",["n"])])), ["m",numb(0),sym("M")])
termB2B = termB2.substituteMinusInds( trac(mult([op("a",["z"]),op("b",["x"])])), addi([sym("Q"),sym("W")]), ["z","x"] )
if termB2B != summ(addi([sym('Q',[]),sym('W',[])]),['m',numb(0),sym('M',[])]):
    raise Exception("B2")

termB3 = summ(trac(mult([op('a^T',['p']),op('a',['p']),op('f',['j']),rho()])), ["m",numb(0),sym("M")])
termB3B = termB3.substituteMinusInds(mult([op('a^T',['z']),op('a',['z']),op('f',['y']),rho()]),sym("Q"),["z","y"])
if termB3B != summ(trac(sym('Q',[])),['m',numb(0),sym('M',[])]):
    raise Exception("B3")
#%%
print("DONE TESTS")