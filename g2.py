# -*- coding: utf-8 -*-
"""
Created on Mon Jan 21 15:57:36 2019

@author: Henry Hesten
"""

from numpy import *
    
from classes import *
exec(open("tests.py").read())

latexLen = 60
#%%
rr = rho()
am = op("a",["m"])
aTm = op("a^T",["m"])
ap = op("a",["p"])
aTp = op("a^T",["p"])
aq = op("a",["q"])
aTq = op("a^T",["q"])
spj = op("sp",["j"])
smj = op("sp^T",["j"])
spl = op("sp",["l"])
sml = op("sp^T",["l"])
spk = op("sp",["k"])
smk = op("sp^T",["k"])
fl = mult([spl,sml])
fk = mult([spk,smk])
def lind(opp):
    return addi([
            mult([opp.HermConj(),opp,rr]),
            mult([rr,opp.HermConj(),opp]),
            mult([numb(-2),opp,rr,opp.HermConj()]),
            ])
masterA = [
        summ(mult([numb(-0.5),sym("kap"), lind(am)]),["m",numb(0),sym("M")]),
        summ(mult([numb(-0.5),sym("gu",["j"]), lind(spj)]),["j",numb(0),sym("J")]),
        summ(mult([numb(-0.5),sym("gd",["j"]), lind(smj)]),["j",numb(0),sym("J")]),
        summ(summ(mult([numb(-0.5),sym("ld",["m","j"]), lind(mult([aTm,smj]))]),["j",numb(0),sym("J")]),["m",numb(0),sym("M")]),
        summ(summ(mult([numb(-0.5),sym("lu",["m","j"]), lind(mult([am,spj]))]),["j",numb(0),sym("J")]),["m",numb(0),sym("M")]),
        ]

master = addi(masterA)

az = op("a",["z"])
aTz = op("a^T",["z"])
ay = op("a",["y"])
aTy = op("a^T",["y"])
fu = op("f",["u"])
fv = op("f",["v"])
corrA = [
        [trac(mult([aTz, az, fu, rho()])), addi([ mult([ trac(mult([aTz, az, rho()])) , trac(mult([fu, rho()])) ]), sym("Cnf",["z","u"])]) ],
        [trac(mult([aTz, aTy, az, ay, rho()])), addi([ mult([ trac(mult([aTz, az, rho()])) , trac(mult([aTy,ay, rho()])) ]), sym("Cnn",["z","y"]), mult([numb(-1),delta(["z","y"]),trac(mult([aTz,az,rho()]))])]) ],
        [trac(mult([fv,fu, rho()])), addi([ mult([ trac(mult([fv, rho()])) , trac(mult([fu, rho()])) ]), sym("Cff",["v","u"])]) ],
        [trac(mult([aTz, aTy, az, ay, fu, rho()])), 
             addi([ sym("Cnnf",["z","y","u"]),
                    mult([trac(mult([aTz,az,rho()])), sym("Cnf",["y","u"])]),
                    mult([trac(mult([aTy,ay,rho()])), sym("Cnf",["z","u"])]),
                    mult([trac(mult([fu,rho()])), sym("Cnn",["z","y"])]),
                    mult([trac(mult([aTz, az, rho()])), trac(mult([aTy, ay, rho()])) ,trac(mult([fu, rho()])) ]),
                    mult([numb(-1),delta(["z","y"]),trac(mult([aTy, ay, rho()])) ,trac(mult([fu, rho()]))]),
                    mult([numb(-1),delta(["z","y"]),sym("Cnf",["y","u"])])
                ])],
        [trac(mult([aTz, az, fv, fu, rho()])), 
             addi([ sym("Cnff",["z","v","u"]),
                    mult([trac(mult([aTz,az,rho()])), sym("Cff",["v","u"])]),
                    mult([trac(mult([fv,rho()])), sym("Cnf",["z","u"])]),
                    mult([trac(mult([fu,rho()])), sym("Cnf",["z","v"])]),
                    mult([trac(mult([aTz, az, rho()])), trac(mult([fv, rho()])) ,trac(mult([fu, rho()])) ]),
                ])],
        ]
def subCorr(term):
    for fromm,to in corrA:
        term = term.substituteMinusInds(fromm,to)
    return term.simplify().collectByTrace()
#%%
npEq = equality(deriv(trac(mult([aTp,ap,rr])),sym("t")),trac(mult([aTp,ap,master])) )
npSim = npEq.simplify()
print(mySplitLatex( npSim.collectByTrace().toLatex(), latexLen))

#%%
flEq = equality(deriv(trac(mult([fl,rr])),sym("t")),trac(mult([fl,master])) )
flSim = flEq.simplify()
print(mySplitLatex(flSim.collectByTrace().toLatex(), latexLen))
#%%
AtpAtqApAq = equality(deriv(trac(mult([aTp,aTq,ap,aq,rr])),sym("t")),trac(mult([aTp,aTq,ap,aq,master])) )
AtpAtqApAqSim = AtpAtqApAq.simplify({"doRep":True})
print(mySplitLatex(AtpAtqApAqSim.collectByTrace().toLatex(), latexLen))
#%%
AtpApFl = equality(deriv(trac(mult([aTp,ap,fl,rr])),sym("t")),trac(mult([aTp,ap,fl,master])) )
AtpApFlSim = AtpApFl.simplify()
print(mySplitLatex(AtpApFlSim.collectByTrace().toLatex(), latexLen))
#%%
FlFk = equality(deriv(trac(mult([fl,fk,rr])),sym("t")),trac(mult([fl,fk,master])) )
FlFkSim = FlFk.simplify()
print(mySplitLatex(FlFkSim.collectByTrace().toLatex(), latexLen))

#%%
npSim2 = subCorr(npSim)
print(mySplitLatex(npSim2.toLatex(), latexLen))
#%%
flSim2 = subCorr(flSim)
print(mySplitLatex(flSim2.toLatex(), latexLen))
#%%
AtpAtqApAqSim2 = subCorr(AtpAtqApAqSim).isolateTerm(deriv(sym('Cnn',['z', 'x']),sym('t',[])))
print(mySplitLatex(AtpAtqApAqSim2.toLatex(), latexLen))
#%%
AtpApFlSim2 = subCorr(AtpApFlSim).isolateTerm(deriv(sym('Cnf',['z', 'x']),sym('t',[])))
print(mySplitLatex(AtpApFlSim2.toLatex(), latexLen))
#%%
FlFkSim2 = subCorr(FlFkSim).isolateTerm(deriv(sym('Cff',['k', 'l']),sym('t',[])))
print(mySplitLatex(FlFkSim2.toLatex(), latexLen))
#%%
AtpAtqApAqSim3 = AtpAtqApAqSim2.substituteMinusInds(npSim2.lhs,npSim2.rhs).simplify().collectByTrace(True)
print(mySplitLatex(AtpAtqApAqSim3.toLatex(), 55))
#%%
tmp = AtpApFlSim2.substituteMinusInds(npSim2.lhs,npSim2.rhs).simplify()
AtpApFlSim3 = tmp.substituteMinusInds(flSim2.lhs,flSim2.rhs).simplify().collectByTrace(True)
print(mySplitLatex(AtpApFlSim3.toLatex(), 55))
#%%
FlFkSim3 = FlFkSim2.substituteMinusInds(flSim2.lhs,flSim2.rhs).simplify().collectByTrace(True)
print(mySplitLatex(FlFkSim3.toLatex(), 53))

#%%
cnfA = []
otherA = []
for tt in AtpApFlSim3.rhs.termA:
    tA = tt.getTrTermA(True)
    match=False
    for t2 in tA:
        if t2.equalMinusInds(sym("Cnf",["z","x"]),["z","x"]):
            match = True
    if match == True:
        cnfA.append(tt)
    else:
        otherA.append(tt)
#cnfS = equality(minus(None,addi(cnfA)),addi(otherA))
#print(mySplitLatex(cnfS.simplify().collectByTrace(True).toLatex(),55))

nBar = equality(sym("nb"), sum(addi(otherA)))
dBar = equality(mult([sym("db"),sym('Cnf',['p', 'l'])]), sum(addi(cnfA)))
dBar2 = equality(sym("db"),dBar.rhs.factoriseNest(sym('Cnf',['p', 'l'])))
print(mySplitLatex(nBar.simplify().collectByTrace(True).toLatex(),52))
#print(mySplitLatex(dBar.simplify().collectByTrace(True).toLatex(),55))
print(mySplitLatex(dBar2.simplify().collectByTrace(True).toLatex(),65))
