# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 15:35:35 2018

@author: Henry Hesten
"""

import numpy as np
import copy
import re

import abc

#%%
def myFlatten(lst):
    ret = []
    for el in lst:
        if isListOrArray(el):
            ret += myFlatten(el)
        else:
            ret.append(el)
    return ret
def mySet(lst, fnc=None):
    if fnc==None:
        fnc = lambda x:x
    dd = {}
    for el in lst:
        dd[fnc(el)] = el
    return list(dd.values())
def isListOrArray(aa):
    if isinstance(aa,np.ndarray):
        if len(aa.shape) == 0:
            return False
        else:
            return True
    return isinstance(aa,list) or isinstance(aa,tuple) or isinstance(aa,np.matrix)
def getDigitsVarBase(num,baseArr,offsets=None):
    if isListOrArray(num):
        return map(lambda x: getDigitsVarBase(x,baseArr,offsets=offsets), num)
    #print(num,baseArr)
    if isinstance(baseArr, int):
        num2 = num
        if(num2==0):
            num2 = 1
        baseArr = [baseArr]*int(np.log(num2)/np.log(baseArr)+1)
    if(offsets==None):
        offsets=[0]*len(baseArr)
    ret = [0]*len(baseArr)
    for i,b in list(enumerate(baseArr))[::-1]:
        tmp = int(num/b)
        ret[i] = int(num-tmp*b+offsets[i])
        num=tmp
    if num>0:
        return None;
    return ret

def inds2Str(inds):
    if inds == []:
        return ""
    elif len(inds) == 1:
        return "_{}".format(str(inds[0]))
    else:
        return "_{{{}}}".format(",".join(map(str,inds)))

class term(abc.ABC):
    simplifiedBool = False
    sortd = False
    
    @abc.abstractproperty
    def containsTerms(self):
        pass
    @abc.abstractmethod
    def toStr(self):
        pass
    def __str__(self):
        return self.toStr()
    @abc.abstractmethod
    def toLatex(self, top=False):
        pass
    @abc.abstractmethod
    def toCodeP(self, dicc=None):
        pass
    def toCode(self, dicc=None):
        code = self.toCodeP(dicc)
        if dicc == None:
            return code
        else:
            ret = code
            if "simp" in dicc:
                if dicc["simp"] == True:
                    if self.isSimplified():
                        ret = "[# {} ]".format(ret)
                    else:
                        ret = "[ {} ]".format(ret)
            return ret
    @abc.abstractmethod
    def containsOp(self):
        pass
    @abc.abstractmethod
    def HermConj(self):
        pass
    def containsOpForSort(self):
        if isinstance(self,trac):
            return True
        else:
            return self.containsOp()
    @abc.abstractmethod
    def getInds(self):
        pass
    @abc.abstractmethod
    def simplifyTypeP(self,dicc):
        pass
    
    ## returns a dictionary of {otherI: selfI} required to make self match other
    ## Or returns False
    def equalMinusInds(self,other,freeInds):
        freeInds = list(set(freeInds))
        if type(other) != type(self):
            return False
        else:
            return self.equalMinusIndsP(other,freeInds)
    @abc.abstractmethod
    def equalMinusIndsP(self,other,freeInds):
        pass
    @abc.abstractmethod
    def eqP(self,other):
        pass
    def __eq__(self,other):
        if isinstance(other, self.__class__):
            return self.eqP(other)
        else:
            return False
    def __neq__(self,other):
        return not self.__eq__(other)
    def __lt__(self,other):
        sf = sortFunc(self,other)
        if sf==-1:
            return True
        else:
            return False
    def __gt__(self,other):
        sf = sortFunc(self,other)
        if sf==1:
            return True
        else:
            return False
    def __le__(self,other):
        sf = sortFunc(self,other)
        if sf<=0:
            return True
        else:
            return False
    def __ge__(self,other):
        sf = sortFunc(self,other)
        if sf>=0:
            return True
        else:
            return False
    @abc.abstractmethod
    def sortStr(self):
        pass
    @abc.abstractmethod
    def doSort(self):
        pass
    @abc.abstractmethod
    def isSimplified(self):
        pass
    @abc.abstractmethod
    def factoriseNest(self,term1):
        pass
    @abc.abstractmethod
    def substitute(self,fromm,to):
        pass
    ## substitutes from as to, but the indices don't  have to match
    def substituteMinusInds(self,fromm,to,freeInds=None):
        if freeInds == None:
            freeInds = to.getInds()
        return self.substituteMinusIndP(fromm,to,freeInds)
    @abc.abstractmethod
    def substituteMinusIndP(self,fromm,to,freeInds):
        pass
    @abc.abstractmethod
    def changeInd(self,fromI,toI): ## returns True is an index was changed
        pass
    def sameTermSymOp(self,other): # [bool, mult(syms_ops), num1, num2]
        arr = [self,other]
        arr2 = []
        for i,term in enumerate(arr):
            if isinstance(term,mult):
                if isinstance(term.termA[0],numb):
                    numm = term.termA[0]
                    rest = term.termA[1:]
                else:
                    numm = numb(1)
                    rest = term.termA
                if len(rest)==1:
                    arr2.append([numm,rest[0]])
                else:
                    arr2.append([numm,mult(rest)])
            else:
                arr2.append([numb(1),term])
        if arr2[0][1] == arr2[1][1]:
            return [True,arr2[0][1],arr2[0][0],arr2[1][0]]
        else:
            return [False,None,None,None]
    
    def simplify(self, dicc=None):
        dd = {
                "mergeFrac":False,
                "preserveEq":True,
                "doRep":True,
                }
        if dicc != None:
            dd.update(dicc)
        return self.simplifyType(dd)
    def simplifyType(self,dicc):
        if self.isSimplified() and dicc["mergeFrac"]:
            return self
        else:
            ret = copy.deepcopy(self)
            while(True):
                old = copy.deepcopy(ret)
                ret = ret.simplifyTypeP(dicc)
                if dicc["mergeFrac"]:
                    ret = ret.mergeFrac()
                if old==ret:
                    break
            ret.simplifiedBool = True
            return ret
    def mergeFrac(self):
        return self
    def collectByTrace(self, incCorr):
        raise Exception("Not implemented {}".format(type(self)))
    
    ## returns [i,z] where self = z*term**i
    def powerOf(self,term):
        return [0,self]

class numb(term):
    containsTerms = False;
    sortd = True
    def __init__(self,val):
        self.val = val
    def toStr(self):
        v = self.val
        if int(v) == v:
            return str(int(v))
        else:
            return str(v)
    def toLatex(self, top=False):
        v = self.val
        if int(v) == v:
            return str(int(v))
        else:
            return str(v)
    def toCodeP(self,dicc):
        return "numb({})".format(self.val)
    def containsOp(self):
        return False
    def HermConj(self):
        return self
    def getInds(self):
        return []
    def simplifyTypeP(self,dicc):
        return self
    def equalMinusIndsP(self,other,freeInds):
        if self == other:
            return {}
        else:
            return False
    def eqP(self,other):
        if self.val==other.val:
            return True
        else:
            return False
    def sortStr(self):
        #return "a{}".format(self.val)
        return str(self.val)
    def doSort(self):
        pass
    def isSimplified(self):
        return True
    def changeInd(self,fromI,toI):
        return False
    def factoriseNest(self,term1):
        raise Exception("Cannot Factorise a Number")
    def substitute(self,fromm,to):
        if self == fromm:
            return to
        else:
            return self
    def substituteMinusIndP(self,fromm,to,freeInds):
        return self.substitute(fromm,to)
    def getTrTermA(self, incCorr=False):
        return []
    
class sym(term):
    containsTerms = False;
    sortd = True
    def __init__(self,name,inds=[]):
        if not isListOrArray(inds):
            raise Exception("Indices must be list")
        self.name = name
        self.inds = inds # [], ["m"], ["m","n"]
    def toStr(self):
        return self.name+inds2Str(self.inds)
    def toLatex(self, top=False):
        if self.name == "al":
            return "\\alpha"+inds2Str(self.inds)
        if self.name == "ld":
            return "\\Lambda^{\\downarrow}"+inds2Str(self.inds)
        if self.name == "lu":
            return "\\Lambda^{\\uparrow}"+inds2Str(self.inds)
        if self.name == "kap":
            return "\\kappa"+inds2Str(self.inds)
        if self.name == "gu":
            return "\\Gamma^{\\uparrow}"+inds2Str(self.inds)
        if self.name == "gd":
            return "\\Gamma^{\\downarrow}"+inds2Str(self.inds)
        if self.name == "Cnf":
            return "C(n_{},f_{})".format(self.inds[0],self.inds[1])
        if self.name == "Cnn":
            return "C(n_{},n_{})".format(self.inds[0],self.inds[1])
        if self.name == "Cff":
            return "C(f_{},f_{})".format(self.inds[0],self.inds[1])
        if self.name == "Cnnf":
            return "C(n_{},n_{},f_{})".format(self.inds[0],self.inds[1],self.inds[2])
        if self.name == "Cnff":
            return "C(n_{},f_{},f_{})".format(self.inds[0],self.inds[1],self.inds[2])
        return self.name+inds2Str(self.inds)
    def toCodeP(self,dicc):
        return "sym('{}',{})".format(self.name,self.inds)
    def containsOp(self):
        return False
    def HermConj(self):
        m = re.search('(.*)\^\*$',self.name)
        if m == None:
            return sym(self.name+"^*",self.inds)
        else:
            return sym(m.group(1),self.inds)
    def getInds(self):
        return self.inds
    def simplifyTypeP(self,dicc):
        return self
    def equalMinusIndsP(self,other,freeInds):
        if self.name != other.name:
            return False
        if len(self.inds) != len(other.inds):
            return False
        dd = {}
        for j in range(len(self.inds)):
            si = self.inds[j]
            oi = other.inds[j]
            if oi in freeInds:
                dd[oi] = si
            elif oi==si:
                continue
            else:
                return False
        return dd
    def eqP(self,other):
        if self.name==other.name and self.inds==other.inds:
            return True
        else:
            return False
    def sortStr(self):
        return "{}{}".format(self.name,self.inds)
    def doSort(self):
        pass
    def isSimplified(self):
        return True
    def changeInd(self,fromI,toI):
        chg = False
        for j,ind in enumerate(self.inds):
            if ind == fromI:
                self.inds[j] = toI
                chg = True
        return chg
    def factoriseNest(self,term1):
        raise Exception("Cannot Factorise a Symbol")
    def substitute(self,fromm,to):
        if self == fromm:
            return to
        else:
            return self
    def substituteMinusIndP(self,fromm,to,freeInds):
        dd = self.equalMinusInds(fromm,freeInds)
        if dd==False:
            return self
        else:
            to = copy.deepcopy(to)
            for k in dd.keys():
                to.changeInd(k,dd[k])
            return to
    def getTrTermA(self, incCorr=False):
        if incCorr:
            namA = ["Cnn","Cnf","Cff","Cnnf","Cnff"]
            if self.name in namA:
                return [self]
            else:
                return []
        else:
            return []
    def collectByTrace(self, incCorr):
        return self

class op(term):
    containsTerms = False;
    sortd = True
    def __init__(self,name,inds=[]):
        if not isListOrArray(inds):
            raise Exception("Indices must be list")
        self.name = name
        self.inds = inds # [], ["m"], ["m","n"]
    def toStr(self):
        return self.name+inds2Str(self.inds)
    def toLatex(self, top=False):
        return self.name+inds2Str(self.inds)
    def toCodeP(self,dicc):
        return "op('{}',{})".format(self.name,self.inds)
    def containsOp(self):
        return True
    def HermConj(self):
        m = re.search('(.*)\^T$',self.name)
        if m == None:
            return op(self.name+"^T",self.inds)
        else:
            return op(m.group(1),self.inds)
    def getInds(self):
        return self.inds
    def simplifyTypeP(self,dicc):
        return self
    def equalMinusIndsP(self,other,freeInds):
        if self.name != other.name:
            return False
        if len(self.inds) != len(other.inds):
            return False
        dd = {}
        for j in range(len(self.inds)):
            si = self.inds[j]
            oi = other.inds[j]
            if oi in freeInds:
                dd[oi] = si
            elif oi==si:
                continue
            else:
                return False
        return dd
    def eqP(self,other):
        if self.name==other.name and self.inds==other.inds:
            return True
        else:
            return False
    def sortStr(self):
        return "{}{}".format(self.name,self.inds)
    def doSort(self):
        self.sortd = True
        pass
    def isSimplified(self):
        return True
    def changeInd(self,fromI,toI):
        chg = False
        for j,ind in enumerate(self.inds):
            if ind == fromI:
                self.inds[j] = toI
                chg = True
        return chg
    def factoriseNest(self,term1):
        raise Exception("Cannot Factorise an op")
    def substitute(self,fromm,to):
        if self == fromm:
            return to
        else:
            return self
    def substituteMinusIndP(self,fromm,to,freeInds):
        dd = self.equalMinusInds(fromm,freeInds)
        if dd==False:
            return self
        else:
            to = copy.deepcopy(to)
            for k in dd.keys():
                to.changeInd(k,dd[k])
            return to
    def getTrTermA(self, incCorr=False):
        return []
class rho(op):
    name = "rho"
    inds = []
    def __init__(self):
        pass
    def sortStr(self):
        return "zzz"
    def toCodeP(self,dicc):
        return "rho()"
    def HermConj(self):
        return self
    def toLatex(self, top=False):
        return "\\rho"
    def getTrTermA(self, incCorr=False):
        return []
    
class delta(term):
    containsTerms = False;
    def __init__(self,inds):
        if len(inds) != 2:
            raise Exception("Wrong number of inds: Delta")
        for i in [0,1]:
            if not isinstance(inds[i],str):
                raise Exception("Delta contains list as ind")
        self.inds = inds
    def toStr(self):
        return "d"+inds2Str(self.inds)
    def toLatex(self, top=False):
        return "\\delta"+inds2Str(self.inds)
    def toCodeP(self,dicc):
        return "delta({})".format(self.inds)
    def containsOp(self):
        return False
    def HermConj(self):
        return self
    def getInds(self):
        return self.inds
    def simplifyTypeP(self,dicc):
        if self.inds[0] == self.inds[1]:
            return numb(1)
        return self
    def equalMinusIndsP(self,other,freeInds):
        dd = {}
        for j in range(2):
            si = self.inds[j]
            oi = other.inds[j]
            if oi in freeInds:
                dd[oi] = si
            elif oi==si:
                continue
            else:
                return False
        return dd
    def eqP(self,other):
        if self.inds != other.inds:
            return False
        return True
    def sortStr(self):
        #return "c{}".format(self.inds)
        return str(self.inds)
    def doSort(self):
        pass
    def isSimplified(self):
        return True
    def changeInd(self,fromI,toI):
        chg = False
        for j,ind in enumerate(self.inds):
            if ind == fromI:
                self.inds[j] = toI
                chg = True
        return chg
    def factoriseNest(self,term1):
        raise Exception("Cannot Factorise a delta")
    def substitute(self,fromm,to):
        if self == fromm:
            return to
        else:
            return self
    def substituteMinusIndP(self,fromm,to,freeInds):
        dd = self.equalMinusInds(fromm,freeInds)
        if dd==False:
            return self
        else:
            to = copy.deepcopy(to)
            for k in dd.keys():
                to.changeInd(k,dd[k])
            return to
    def getTrTermA(self, incCorr=False):
        return []

class trac(term):
    containsTerms = True;
    sortd = True
    def __init__(self,iner):
        self.iner = iner
    def toStr(self):
        return "Tr_{{{}}}".format(self.iner.toStr())
    def toLatex(self, top=False):
        if isinstance(self.iner,mult):
            tA = self.iner.termA
            if tA[-1] == rho():
                if len(tA) == 2:
                    if isinstance(tA[0],op):
                        if tA[0].name == "f":
                            if len(tA[0].inds) == 1:
                                return "\\braket{{ f_{{ {} }} }}".format(tA[0].inds[0])
                elif len(tA) == 3:
                    good = True
                    namA = ["a^T","a"]
                    ind = None
                    for i in [0,1]:
                        if isinstance(tA[i],op):
                            if tA[i].name == namA[i]:
                                if ind == None:
                                    ind = tA[i].inds
                                elif tA[i].inds != ind:
                                    good = False
                            else:
                                good = False
                        else:
                            good = False
                    if good:
                        if len(ind) == 1:
                            return "\\braket{{ n_{{ {} }} }}".format(ind[0])
        return "\\tr\\{{{}\\}}".format(self.iner.toLatex())
    def toCodeP(self,dicc):
        return "trac({})".format(self.iner.toCode(dicc))
    def containsOp(self):
        return False
    def HermConj(self):
        return trac(self.iner.HermConj())
    def getInds(self):
        return self.iner.getInds()
    def simplifyTypeP(self,dicc):
        self.iner = self.iner.simplifyType(dicc)
        
        ## expand tr( + )
        if isinstance(self.iner, addi):
            ntA = []
            for tt in self.iner.termA:
                ntA.append(trac(tt))
            return addi(ntA)
        
        ## take non-op out
        if isinstance(self.iner, mult):
            out = []
            i = 0
            '''
            while i< len(self.iner.termA):
                if self.iner.termA[i].containsOp():
                    out.append(self.iner.termA.pop(i))
                else:
                    i+=1
            '''
            while not self.iner.termA[0].containsOp():
                out.append(self.iner.termA.pop(0))
                if len(self.iner.termA) == 0:
                    break
            if len(out)>0:
                if len(self.iner.termA ) ==0 :
                    raise Exception("trace of nothing")
                nt = mult(out+[trac(self.iner)])
                nt = nt.simplifyType(dicc)
                return nt
            
        ## cycle
        if isinstance(self.iner, mult):
            strA = list(map(lambda x: x.sortStr(), self.iner.termA))
            maxI = max(list(enumerate(strA)), key=lambda x:x[1])[0]
            if maxI<len(self.iner.termA)-1:
                ntA = self.iner.termA[maxI+1:] + self.iner.termA[:maxI+1]
                nt = trac(mult(ntA))
                nt = nt.simplifyType(dicc)
                return nt
            
        ## move sums outside
        if isinstance(self.iner, summ):
            inn = self.iner
            nt = summ(trac(inn.iner),[inn.sumI,inn.sumMin,inn.sumMax])
            nt = nt.simplifyType(dicc)
            return nt
        
        if isinstance(self.iner, rho):
            return numb(1)
        
        return self
    def equalMinusIndsP(self,other,freeInds):
        return self.iner.equalMinusInds(other.iner,freeInds)
    def eqP(self,other):
        if self.iner != other.iner:
            return False
        return True
    def sortStr(self):
        #if self.iner.containsOp():
        #    ret = "u"
        #else:
        #    ret = "d"
        ret = self.iner.sortStr()
        return ret
    def doSort(self):
        self.iner.doSort()
    def isSimplified(self):
        return self.iner.isSimplified()
    def changeInd(self,fromI,toI):
        return self.iner.changeInd(fromI,toI)
    def factoriseNest(self,term1):
        raise Exception("Cannot Factorise a Trace")
    def substitute(self,fromm,to):
        if self == fromm:
            return to
        else:
            return trac(self.iner.substitute(fromm,to))
    def substituteMinusIndP(self,fromm,to,freeInds):
        dd = self.equalMinusInds(fromm,freeInds)
        if dd==False:
            return trac(self.iner.substituteMinusIndP(fromm,to,freeInds))
        else:
            to = copy.deepcopy(to)
            for k in dd.keys():
                to.changeInd(k,dd[k])
            return to
    def collectByTrace(self, incCorr=False):
        return self
    def getTrTermA(self, incCorr=False):
        return [self]
    
class summ(term):
    containsTerms = True;
    def __init__(self,iner,sumData):
        self.iner = iner
        self.sumI = sumData[0]
        self.sumMin = sumData[1]
        self.sumMax = sumData[2]
    def getSumData(self):
        return [self.sumI,self.sumMin,self.sumMax]
    def toStr(self):
        return "SUM_{{{}={}}}^{{{}}}{{{}}}".format(self.sumI,self.sumMin.toStr(),self.sumMax.toStr(),self.iner.toStr())
    def toLatex(self, top=False):
        return "\sum_{{{}={}}}^{{{}}}{{{}}}".format(self.sumI,self.sumMin.toLatex(),self.sumMax.toLatex(),self.iner.toLatex())
    def toCodeP(self,dicc):
        return "summ({},['{}',{},{}])".format(self.iner.toCode(dicc),self.sumI,self.sumMin.toCode(dicc),self.sumMax.toCode(dicc))
    def containsOp(self):
        return self.iner.containsOp()
    def HermConj(self):
        return summ(self.iner.HermConj(),self.getSumData())
    def getInds(self):
        return self.iner.getInds()
    def simplifyTypeP(self,dicc): # Do not move ops out of sums it breaks trace
        self.iner = self.iner.simplifyType(dicc)
        
        ## expand sum +
        if isinstance(self.iner, addi):
            ntA = []
            for tt in self.iner.termA:
                ntA.append(summ(tt,[self.sumI,self.sumMin,self.sumMax]))
            nt = addi(ntA)
            nt.simplifyType(dicc)
            return nt
        
        ## check deltas
        if isinstance(self.iner, mult):
            for i,tt in enumerate(self.iner.termA):
                if isinstance(tt,delta):
                    for j,ind in enumerate(tt.inds):
                        if ind == self.sumI:
                            toI = tt.inds[1-j]
                            nt = self.iner
                            nt.changeInd(self.sumI,toI)
                            nt = nt.simplifyType(dicc)
                            return nt
                        
        ## move not summed over, not op outside
        if isinstance(self.iner, mult):
            out = []
            i = 0
            while i<len(self.iner.termA):
                if self.iner.termA[i].containsOp() or (self.sumI in self.iner.termA[i].getInds()):
                    i += 1
                else:
                    out.append(self.iner.termA.pop(i))
                    self.iner.simplifiedBool = False
            if len(self.iner.termA) == 0:
                nt = mult(out+[addi([self.sumMax, mult([numb(-1),self.sumMin])])])
                nt.simplifyType(dicc)
                return nt
            if len(out)>0:
                nt = mult(out+[self])
                nt.simplifyType(dicc)
                return nt
        else:
            if not (self.sumI in self.iner.getInds()):
                nt = mult( [self.iner, addi([self.sumMax, mult([numb(-1),self.sumMin])])] )
                nt = nt.simplifyType(dicc)
                return nt
        
        return self
    def equalMinusIndsP(self,other,freeInds):
        dd = {}
        oi = other.sumI
        si = self.sumI
        if oi in freeInds:
            dd[oi] = si
        elif oi!=si:
            return False
        freeInds = copy.deepcopy(freeInds)
        if oi in freeInds:
            freeInds.remove(oi)
        other = copy.deepcopy(other)
        other.changeInd(oi,si)
        tmp = self.iner.equalMinusInds(other.iner,freeInds)
        if tmp == False:
            return tmp
        else:
            dd.update(tmp)
            return dd
    def eqP(self,other):
        if self.iner != other.iner:
            return False
        if self.sumI != other.sumI:
            return False
        if self.sumMin != other.sumMin:
            return False
        if self.sumMax != other.sumMax:
            return False
        return True
    def sortStr(self):
        #if self.containsOp():
        #    ret = "w"
        #else:
        #    ret = "h"
        ret = str(self.sumI)+self.sumMin.sortStr()+self.sumMax.sortStr()
        ret += self.iner.sortStr()
        return ret
    def doSort(self):
        self.sortd = True
        self.iner.doSort()
    def isSimplified(self):
        return self.iner.isSimplified()
    def changeInd(self,fromI,toI):
        return self.iner.changeInd(fromI,toI)
    def factoriseNest(self,term1):
        raise Exception("Cannot Factorise a Sum")
    def substitute(self,fromm,to):
        if self == fromm:
            return to
        else:
            return summ(self.iner.substitute(fromm,to),self.getSumData())
    def substituteMinusIndP(self,fromm,to,freeInds):
        dd = self.equalMinusInds(fromm,freeInds)
        if dd==False:
            return summ(self.iner.substituteMinusIndP(fromm,to,freeInds), self.getSumData())
        else:
            to = copy.deepcopy(to)
            for k in dd.keys():
                to.changeInd(k,dd[k])
            return to
    def getTrTermA(self, incCorr=False):
        return self.iner.getTrTermA( incCorr)
    
class addi(term):
    containsTerms = True;
    def __init__(self,termA):
        self.termA = termA
    def toStr(self):
        return "({})".format(" + ".join(map(lambda x: x.toStr(), self.termA)) )
    def toLatex(self, top=False):
        ret = ""
        if not top:
            ret += "\\left("
        for i,tt in enumerate(self.termA):
            ll = tt.toLatex()
            if i>0 and ll[0] != "-":
                ret += " + "
            else:
                ret += " "
            ret += ll
        if not top:
            ret += "\\right)"
        return ret
        #return "({})".format(" + ".join(map(lambda x: x.toLatex(), self.termA)) )
    def toCodeP(self,dicc):
        return "addi([{}])".format(",".join(map(lambda x: x.toCode(dicc), self.termA)))
    def containsOp(self):
        return any(list(map(lambda x: x.containsOp(), self.termA)))
    def HermConj(self):
        return addi(list(map(lambda x: x.HermConj(), self.termA)))
    def getInds(self):
        return myFlatten(list(map(lambda x: x.getInds(), self.termA)))
    def writeTermA(self,newTermA):
        self.termA = newTermA
        self.sortd = False
        self.simplifiedBool = False
    def writeTermAElem(self,newElem,ind):
        self.termA[ind] = newElem
        self.sortd = False
        self.simplifiedBool = False
    def doSort(self):
        if not self.sortd:
            #self.termA = sorted(self.termA, key=lambda x: x.sortStr())
            self.termA = sorted(self.termA)
            self.sortd = True
    def simplifyTypeP(self,dicc):
        self.writeTermA( list(map(lambda x: x.simplifyType(dicc), self.termA)) )
        self.doSort()
        if len(self.termA) == 1:
            return self.termA[0]
        if len(self.termA) == 0:
            raise Exception("Addition contains no terms")
        if self.termA[0] == numb(0):
            self.writeTermA(self.termA[1:])
        
        i=0
        while(i<len(self.termA)): ## collapse additions of additions
            if isinstance(self.termA[i], addi):
                self.writeTermA( self.termA + self.termA.pop(i).termA )
            else:
                i += 1
        
        i=0
        while(i<len(self.termA)): ## combine terms
            j = i+1
            while(j<len(self.termA)):
                st = self.termA[i].sameTermSymOp(self.termA[j])
                if st[0]:
                    numm = numb( st[2].val+st[3].val )
                    self.termA.pop(j) ################################
                    self.writeTermAElem(mult([numm,st[1]]), i)
                else:
                    j += 1
            i += 1
            
        i=0
        while(i<=len(self.termA)): ## combine numbers
            cond = False
            if i<len(self.termA):
                if isinstance(self.termA[i],numb):
                    i += 1
                else:
                    cond = True
            else:
                cond = True
            if cond:
                if i>1:
                    newN = numb(sum(list(map(lambda x: x.val, self.termA[:i]))))
                    self.writeTermA([newN]+self.termA[i:])
                break
                    
        if not self.sortd:
            self.doSort()
        
        return self
    def equalMinusIndsP(self,other,freeInds):
        dd = {}
        oA = other.termA
        sA = self.termA
        dd = {}
        if len(oA) != len(sA):
            return False
        for i in range(len(sA)):
            tmp = sA[i].equalMinusInds(oA[i],freeInds)
            if tmp == False:
                return False
            dd.update(tmp)
            for k in tmp.keys():
                freeInds = copy.deepcopy(freeInds)
                freeInds.remove(k)
                other = copy.deepcopy(other)
                other.changeInd(k,tmp[k])
                oA = other.termA
        return dd
    def eqP(self,other):
        if len(self.termA) != len(other.termA):
            return False
        if list(self.termA) != list(other.termA):
            return False
        return True
    def sortStr(self):
        #if self.containsOp():
        #    ret = "v"
        #else:
        #    ret = "f"
        ret = self.termA[0].sortStr()
        return ret
    def isSimplified(self):
        return self.simplifiedBool and np.prod(list(map(lambda x: x.isSimplified(), self.termA)))
    def changeInd(self,fromI,toI): ## THIS OPERATES STRANGELY!!
        boolA = list(map(lambda x: x.changeInd(fromI,toI), self.termA)) ## list of bools
        if any(boolA):
            self.simplifiedBool = False
            return True
        else:
            return False
        
    def factorise(self,term1):
        if not isinstance(term1, addi):
            term1 = addi([term1])
        divisA = []
        for t in self.termA:
            tmp = []
            for t2 in term1.termA:
                tmp.append(frac(t,t2).simplify())
            divisA.append(tmp)
        
        retA = []
        #i = 0
        while len(divisA)>0:
            factorisedTerm = False
            for j,pontT in enumerate(divisA[0]):
                indA = [0]
                posA = [j]
                for k,div in list(enumerate(divisA))[1:]:
                    for l,t in enumerate(div):
                        if t == pontT:
                            if not (l in posA):
                                posA.append(l)
                                indA.append(k)
                valid = True
                for l in range(len(term1.termA)):
                    if not (l in posA):
                        valid = False
                if valid:
                    retA.append(pontT)
                    for k in sorted(indA)[::-1]:
                        divisA.pop(k)
                    factorisedTerm = True
            if not factorisedTerm:
                raise Exception("Connot factorise {}\nterm {} missing".format(self,self.termA[0]))
        return addi(retA)
    def factoriseNest(self,term1):
        try:
            return self.factorise(term1)
        except:
            newTA = list(map(lambda x: x.factoriseNest(term1), self.termA))
            return addi(newTA).simplifyType(dicc)
    def substitute(self,fromm,to):
        if self == fromm:
            return to
        else:
            tA = list(map(lambda x: x.substitute(fromm,to), self.termA))
            return addi(tA)
    def substituteMinusIndP(self,fromm,to,freeInds):
        dd = self.equalMinusInds(fromm,freeInds)
        if dd==False:
            tA = list(map(lambda x: x.substituteMinusIndP(fromm,to,freeInds), self.termA))
            return addi(tA)
        else:
            to = copy.deepcopy(to)
            for k in dd.keys():
                to.changeInd(k,dd[k])
            return to
    
    ## returns an array where self = arr[0] + arr[1]*term + arr[2]*term**2
    def groupPower(self,term):
        arr = []
        for t in self.termA:
            i,rest = t.powerOf(term)
            while(len(arr)<=i):
                arr.append([])
            arr[i].append(rest)
        for i in range(len(arr)):
            arr[i] = addi(arr[i])
        return arr
    def mergeFrac(self):
        tA = list(map(lambda x: x.mergeFrac(),self.termA))
        if frac in map(type,tA):
            numA = []
            denA = []
            for t in tA:
                if isinstance(t, frac):
                    if t.den in denA:
                        for i,t2 in enumerate(denA):
                            if t.den == t2:
                                numA.append(denA[:i]+denA[i+1:]+[t.num])
                    else:
                        for arr in numA:
                            arr.append(t.den)
                        numA.append(denA+[t.num])
                        denA.append(t.den)
                else:
                    numA.append(denA+[t])
            num = addi(list(map(lambda x: mult(x),numA)))
            den = mult(denA)
            return frac(num,den)
        else:
            return addi(tA)
    def collectByTrace(self, incCorr=False):
        dd = {}
        for t in self.termA:
            trTerms = t.getTrTermA( incCorr)
            key = "_".join(sorted(list(map(lambda x: x.toCode(), trTerms))))
            if not key in dd:
                dd[key] = [trTerms,[]]
            dd[key][1].append(t)
        lenA = []
        keyA = dd.keys()
        for k in keyA:
            ll = 0
            for t in dd[k][0]:
                if isinstance(t,trac):
                    inn = t.iner
                    if isinstance(inn,mult):
                        ll += len(inn.termA)
                    else:
                        ll += 1
                elif isinstance(t,sym):
                    ll += len(t.name)
            lenA.append([k,ll])
        lenA = sorted(lenA, key=lambda x: x[1] )
        ret = []
        for k,l in lenA:
            ret += dd[k][1]
        return addi(ret)
            
    def getTrTermA(self, incCorr=False):
        return mySet(myFlatten(list(map(lambda t: t.getTrTermA( incCorr), self.termA))), lambda x: x.toCode())
        
        
class mult(term):
    containsTerms = True;
    
    ## using f = sp sm = (sz+1)/2
    # [ op1_i, op2_i, lambda i: [op1,op2]_i ]
    sComA = [
            ["a","a^T", lambda x: numb(1)],
            ["sp^T","sp", lambda x: addi([numb(1), mult([numb(-2),op("f",x)])]) ],
            ["sp^T","f", lambda x: op("sp^T",x)],
            ["sp","f", lambda x: mult([numb(-1),op("sp",x)]) ]
        ]
    for op1 in ["a","a^T"]:
        for op2 in ["sp","sp^T","f"]:
            sComA.append([op2,op1,lambda x: numb(0)])
    
    multRepA = [
            [ ["sp","sp^T"], op("f",["j"]), "j" ],
        ]
    
    def __init__(self,termA):
        self.termA = termA
    def toStr(self):
        if self.termA[0] == numb(-1):
            return "-{}".format(" ".join(map(lambda x: x.toStr(), self.termA[1:])) )
        return "{}".format(" ".join(map(lambda x: x.toStr(), self.termA)) )
    def toLatex(self, top=False):
        if self.termA[0] == numb(-1):
            return "-{}".format(" ".join(map(lambda x: x.toLatex(), self.termA[1:])) )
        return "{}".format(" ".join(map(lambda x: x.toLatex(), self.termA)) )
    def toCodeP(self,dicc):
        return "mult([{}])".format(",".join(map(lambda x: x.toCode(dicc), self.termA)))
    def containsOp(self):
        return any(list(map(lambda x: x.containsOp(), self.termA)))
    def HermConj(self):
        return mult(list(map(lambda x: x.HermConj(), self.termA[::-1])))
    def getInds(self):
        return myFlatten(list(map(lambda x: x.getInds(), self.termA)))
    def writeTermA(self,newTermA):
        self.termA = newTermA
        self.sortd = False
        self.simplifiedBool = False
    def writeTermAElem(self,newElem,ind):
        self.termA[ind] = newElem
        self.sortd = False
        self.simplifiedBool = False
    def doSort(self): ## Does Not Commute!
        if not self.sortd:
            notOpA = []
            opA = []
            for tt in self.termA:
                if tt.containsOp():
                    opA.append(tt)
                else:
                    notOpA.append(tt)
            notOpA = sorted(notOpA)
            
            for i,tt in enumerate(opA[:-1]): ## same operator diff ind
                tt2 = opA[i+1]
                if isinstance(tt,op) and isinstance(tt2,op):
                    if tt.name == tt2.name:
                        if len(tt.inds) == 1 and len(tt2.inds) == 1:
                            if tt.inds[0] > tt2.inds[0]:
                                opA[i] = tt2
                                opA[i+1] = tt
                                
            self.writeTermA(notOpA+opA)
            self.sortd = True
    def simplifyTypeP(self,dicc):
        forDebugTermA = self.termA
        self.writeTermA( list(map(lambda x: x.simplifyType(dicc), self.termA)) )
        self.doSort()
        if len(self.termA) == 1:
            return self.termA[0]
        if len(self.termA) == 0:
            raise Exception("Multiplication contains no terms")
        if self.termA[0] == numb(0):
            return numb(0)
        if self.termA[0] == numb(1):
            self.writeTermA(self.termA[1:])
        
        i=0
        while(i<len(self.termA)): ## collapse mult of mult
            if isinstance(self.termA[i], mult):
                self.writeTermA( self.termA[:i] + self.termA[i].termA + self.termA[i+1:] )
            else:
                i += 1
        
        i=0
        while(i<=len(self.termA)): ## combine numbers
            cond = False
            if i<len(self.termA):
                if isinstance(self.termA[i],numb):
                    i += 1
                else:
                    cond = True
            else:
                cond = True
            if cond:
                if i>1:
                    newN = numb(np.prod(list(map(lambda x: x.val, self.termA[:i]))))
                    self.writeTermA([newN]+self.termA[i:])
                break
        
        if not self.sortd:
            self.doSort()
        
        ## commutations
        for i in range(0,len(self.termA)-1):
            t0 = self.termA[i]
            if not t0.containsOp():
                continue
            t1 = self.termA[i+1]
            if (not isinstance(t0,op)) or (not isinstance(t1,op)):
                continue
            for com in self.sComA:
                if t0.name == com[0] and t1.name == com[1]:
                    ind0 = t0.getInds()
                    ind1 = t1.getInds()
                    if len(ind0) != 1 or len(ind1) != 1:
                        raise Exception("Wrong number of indices in commutator {}, {}".format(t0.toStr(),t1.toStr()))
                    commutator = com[2](ind0)
                    nt1 = mult(self.termA[:i]+[t1,t0]+self.termA[i+2:])
                    if ind0 == ind1:
                        nt2 = mult(self.termA[:i]+[commutator]+self.termA[i+2:])
                    else:
                        nt2 = mult([delta(ind0+ind1)]+ self.termA[:i]+[commutator]+self.termA[i+2:])
                    nt = addi([nt1,nt2])
                    nt = nt.simplifyType(dicc)
                    return nt
                
        ## replacements
        if dicc["doRep"]:
            for i in range(0,len(self.termA)-1):
                t0 = self.termA[i]
                t1 = self.termA[i+1]
                if (not isinstance(t0,op)) or (not isinstance(t1,op)) or t0.inds!=t1.inds:
                    continue
                for rep in self.multRepA:
                    if t0.name == rep[0][0] and t1.name == rep[0][1]:
                        replacement = copy.deepcopy(rep[1])
                        replacement.changeInd(rep[2],t0.inds[0])
                        ntA = self.termA[:i] + [replacement] + self.termA[i+2:]
                        nt = mult(ntA)
                        nt = nt.simplifyType(dicc)
                        return nt
        
        ## expansions
        if any(list(map(lambda x: isinstance(x,addi),self.termA))):
            tAA = []
            for tt in self.termA:
                if isinstance(tt, addi):
                    tAA.append(tt.termA)
                else:
                    tAA.append([tt])
            baseA = list(map(len, tAA))
            ntA = []
            for ind in range( np.prod(baseA) ):
                indA = getDigitsVarBase(ind,baseA)
                tmp = []
                for i in range(len(tAA)):
                    tmp.append( tAA[i][indA[i]] )
                ntA.append( mult(tmp) )
            nt = addi(ntA)
            nt = nt.simplifyType(dicc)
            return nt
        
        ## change inds of terms multiplied by delta
        for i,tt in enumerate(self.termA):
            if isinstance(tt, delta):
                changeI = tt.inds[1]
                changeTo = tt.inds[0]
                for j,tt2 in enumerate(self.termA):
                    if i==j:
                        continue
                    tt2.changeInd(changeI,changeTo)
        
        ## grow fractions
        for i,tt in enumerate(self.termA):
            if isinstance(tt, frac):
                num = mult(self.termA[:i] + [tt.num] + self.termA[i+1:])
                nt = frac(num,tt.den).simplifyType(dicc)
                return nt
        
        ## Put ops inside summs, for trace
        if summ in map(type, self.termA):
            befA = []
            aftA = []
            notOpA = []
            foundSum = False
            sumT = None
            for tt in self.termA:
                if type(tt) == summ and foundSum == False:
                    foundSum = True
                    sumT = tt
                else:
                    if tt.containsOp():
                        if foundSum:
                            aftA.append(tt)
                        else:
                            befA.append(tt)
                    else:
                        notOpA.append(tt)
            if len(befA)+len(aftA)>0:
                newIner = mult(befA+[sumT.iner]+aftA).simplifyType(dicc)
                newS = summ(newIner, sumT.getSumData()).simplifyType(dicc)
                ret = mult(notOpA+[newS]).simplifyType(dicc)
                #print(ret.toCode())
                return ret
                    
        if not self.sortd:
            self.doSort()
        return self
    def equalMinusIndsP(self,other,freeInds):
        dd = {}
        oA = other.termA
        sA = self.termA
        if len(oA) != len(sA):
            return False
        for i in range(len(sA)):
            tmp = sA[i].equalMinusInds(oA[i],freeInds)
            if tmp == False:
                return False
            dd.update(tmp)
            for k in tmp.keys():
                freeInds = copy.deepcopy(freeInds)
                freeInds.remove(k)
                other = copy.deepcopy(other)
                other.changeInd(k,tmp[k])
                oA = other.termA
        return dd
    def eqP(self,other):
        if len(self.termA) != len(other.termA):
            return False
        for i in range(len(self.termA)):
            if list(self.termA) != list(other.termA):
                return False
        return True
    def sortStr(self):
        #if self.containsOp():
        #    ret = "u"
        #else:
        #    ret = "e"
        ret = self.termA[0].sortStr()
        return ret
    def isSimplified(self):
        return self.simplifiedBool and np.prod(list(map(lambda x: x.isSimplified(), self.termA)))
    def changeInd(self,fromI,toI): ## THIS OPERATES STRANGELY!!
        boolA = list(map(lambda x: x.changeInd(fromI,toI), self.termA)) ## list of bools
        if any(boolA):
            self.simplifiedBool = False
            return True
        else:
            return False
    def factoriseNest(self,term1):
        for i,t in enumerate(self.termA):
            try:
                ret = t.factoriseNest(term1)
                tA = copy.copy(self.termA)
                tA[i] = ret
                return mult(tA).simplifyType(dicc)
            except:
                continue
    def substitute(self,fromm,to):
        if self == fromm:
            return to
        else:
            tA = list(map(lambda x: x.substitute(fromm,to), self.termA))
            return mult(tA)
    def substituteMinusIndP(self,fromm,to,freeInds):
        dd = self.equalMinusInds(fromm,freeInds)
        if dd==False:
            tA = list(map(lambda x: x.substituteMinusIndP(fromm,to,freeInds), self.termA))
            return mult(tA)
        else:
            to = copy.deepcopy(to)
            for k in dd.keys():
                to.changeInd(k,dd[k])
            return to
    
    ## returns [i,z] where self = z*term**i
    def powerOf(self,term):
        if term.containsOp():
            raise Exception("Power of cannot contain Op")
        tA = copy.copy(self.termA)
        poww = 0
        i=0
        while(i<len(tA)):
            if tA[i] == term:
                poww += 1
                tA.pop(i)
            else:
                i += 1
        return [poww,mult(tA)]
    def mergeFrac(self):
        tA = list(map(lambda x: x.mergeFrac(),self.termA))
        return mult(tA)
    def getTrTermA(self, incCorr=False):
        return mySet(myFlatten(list(map(lambda t: t.getTrTermA( incCorr), self.termA))), lambda x: x.toCode())
    def collectByTrace(self, incCorr=False):
        return self

class frac(term):
    containsTerms = True;
    sortd = False
    def __init__(self,num,den):
        if den.containsOp():
            raise Exception("operator in fraction")
        self.num = num
        self.den = den
    def toStr(self):
        return "({})/({})".format( self.num.toStr(), self.den.toStr() )
    def toLatex(self, top=False):
        return "\\frac{{{}}}{{{}}}".format( self.num.toLatex(), self.den.toLatex() )
    def toCodeP(self,dicc):
        return "frac({},{})".format( self.num.toCode(dicc), self.den.toCode(dicc) )
    def containsOp(self):
        return self.num.containsOp()
    def HermConj(self):
        return frac(self.num.HermConj(),self.den.HermConj())
    def getInds(self):
        return myFlatten( self.num.getInds() + self.den.getInds() ) ########## use set? !! #####
    def doSort(self):
        if not self.sortd:
            self.num.doSort()
            self.sortd = True
    def simplifyTypeP(self,dicc):
        self.num = self.num.simplifyType(dicc)
        self.den = self.den.simplifyType(dicc)
        
        self = self.flattenDoubleDiv()
        
        if not isinstance(self.num,mult):
            nmA = [self.num]
        else:
            nmA = self.num.termA
        if not isinstance(self.den,mult):
            dnA = [self.den]
        else:
            dnA = self.den.termA
        
        ## flatten double division
        #for i,tt in enumerate(nmA):
        #    if isinstance(tt, frac):
        #        nmA[i] = tt.num
        #        dnA.append(tt.den)
        #for i,tt in enumerate(dnA):
        #    if isinstance(tt, frac):
        #        dnA[i] = tt.num
        #        nmA.append(tt.den)
                
        ## cancel
        i = 0
        while(i<len(nmA)):
            j = 0
            while(j<len(dnA)):
                if nmA[i] == dnA[j]:
                    nmA.pop(i)
                    dnA.pop(j)
                    i -= 1
                    break
                else:
                    j += 1
            i += 1
        
        if len(nmA) == 0:
            self.num = numb(1)
        elif len(nmA) == 1:
            self.num = nmA[0].simplifyType(dicc)
        else:
            self.num = mult(nmA).simplifyType(dicc)
            
        if len(dnA) == 0:
            return self.num
        elif len(dnA) == 1:
            self.den = dnA[0].simplifyType(dicc)
            if self.den == numb(1):
                return self.num
        else:
            self.den = mult(dnA).simplifyType(dicc)
                        
        return self
    def equalMinusIndsP(self,other,freeInds):
        oA = [other.num,other.den]
        sA = [self.num,self.den]
        dd = {}
        for i in range(len(sA)):
            tmp = sA[i].equalMinusInds(oA[i],freeInds)
            if tmp == False:
                return False
            dd.update(tmp)
            for k in tmp.keys():
                freeInds = copy.deepcopy(freeInds)
                freeInds.remove(k)
                other = copy.deepcopy(other)
                other.changeInd(k,tmp[k])
                oA = [other.num,other.den]
        return dd
    def eqP(self,other):
        if self.num == other.num and self.den == other.den:
            return True
        else:
            return False
    def sortStr(self):
        ret = self.num.sortStr()+self.den.sortStr()
        return ret
    def isSimplified(self):
        return self.simplifiedBool and self.num.isSimplified() and self.den.isSimplified()
    def changeInd(self,fromI,toI): ## THIS OPERATES STRANGELY!!
        boolA = [self.num.changeInd(fromI,toI), self.den.changeInd(fromI,toI)]
        if any(boolA):
            self.simplifiedBool = False
            return True
        else:
            return False
    def factoriseNest(self,term1):
        num = self.num.factoriseNest(term1)
        return frac(num,self.den).simplifyType(dicc)
    def substitute(self,fromm,to):
        if self == fromm:
            return to
        else:
            return frac(self.num.substitute(fromm,to),self.den.substitute(fromm,to))
    def substituteMinusIndP(self,fromm,to,freeInds):
        dd = self.equalMinusInds(fromm,freeInds)
        if dd==False:
            return frac(self.num.substituteMinusIndP(fromm,to,freeInds),self.den.substituteMinusIndP(fromm,to,freeInds))
        else:
            to = copy.deepcopy(to)
            for k in dd.keys():
                to.changeInd(k,dd[k])
            return to
    def flattenDoubleDiv(self):
        slf = copy.deepcopy(self)
        if not isinstance(slf.num,mult):
            nmA = [slf.num]
        else:
            nmA = slf.num.termA
        if not isinstance(slf.den,mult):
            dnA = [slf.den]
        else:
            dnA = slf.den.termA
        
        ## flatten double division
        for i,tt in enumerate(nmA):
            if isinstance(tt, frac):
                nmA[i] = tt.num
                dnA.append(tt.den)
        for i,tt in enumerate(dnA):
            if isinstance(tt, frac):
                dnA[i] = tt.num
                nmA.append(tt.den)
        
        if len(nmA) == 1:
            nm = nmA[0]
        else:
            nm = mult(nmA)
            
        if len(dnA) == 1:
            dn = dnA[0]
        else:
            dn = mult(dnA)
        return frac(nm,dn)
    def mergeFrac(self):
        #nm = self.num.mergeFrac()
        #dn = self.den.mergeFrac()
        return frac(self.num.mergeFrac(),self.den.mergeFrac()).flattenDoubleDiv()
    def collectByTrace(self, incCorr=False):
        return frac(self.num.collectByTrace(incCorr),self.den)
    

class equality(term):
    simplifiedBool = False
    sortd = False
    containsTerms = True
    
    def __init__(self,lhs,rhs):
        self.lhs = lhs
        self.rhs = rhs
    def toStr(self):
        return "{} = {}".format(self.lhs.toStr(),self.rhs.toStr())
    def toLatex(self, top=False):
        return "{} =&\\; {}\\\\".format(self.lhs.toLatex(True),self.rhs.toLatex(True))
    def toCodeP(self,dicc):
        return "equality({},{})".format(self.lhs.toCode(dicc),self.rhs.toCode(dicc))
    def containsOp(self):
        return self.lhs.containsOp() or self.rhs.containsOp()
    def HermConj(self):
        return equality(self.lhs.HermConj(),self.rhs.HermConj())
    def getInds(self):
        return self.lhs.getInds() + self.rhs.getInds()
    def simplifySides(self):
        nlhs = self.lhs.simplifyType(dicc)
        nrhs = self.rhs.simplifyType(dicc)
        return equality(nlhs,nrhs)
    def simplifyTypeP(self,dicc):
        self.lhs = self.lhs.simplifyType(dicc)
        self.rhs = self.rhs.simplifyType(dicc)
        
        if dicc["preserveEq"]:
            return self
        else:
            if not isinstance(self.lhs,mult):
                lhA = [self.lhs]
            else:
                lhA = self.lhs.termA
            if not isinstance(self.rhs,mult):
                rhA = [self.rhs]
            else:
                rhA = self.rhs.termA
                
            for i in [0,1]:
                ar1 = [lhA,rhA][i]
                ar2 = [rhA,lhA][i]
                for j,tt in enumerate(ar1):
                    if isinstance(tt,frac):
                        ar1[j] = tt.num
                        ar2.append(tt.den)
            
            if len(lhA) == 1:
                nlh = lhA[0]
            else:
                nlh = mult(lhA)
            if len(rhA) == 1:
                nrh = rhA[0]
            else:
                nrh = mult(rhA)
            
            nlhs = addi([nlh,mult([numb(-1),nrh])])
            
            return equality(nlhs.simplifyType(dicc),numb(0))
    
    def removeFrac(self):
        ret = self.simplify({"preserveEq":False})
        if ret.rhs != numb(0):
            raise Exception("Remove Frac RHS != 0")
        if not isinstance(ret.lhs, addi):
            print("RemoveFrac is not an addi")
            return ret
        
        multA = []
        for i,tt in enumerate(ret.lhs.termA):
            nt = mult(multA+[tt]).simplify()
            if not isinstance(nt,mult):
                nt = mult([nt])
            for tt2 in nt.termA:
                if not isinstance(tt2,frac):
                    continue
                multA.append(tt2.den)
        nlhs = mult(multA+[ret.lhs]).simplify()
        return equality(nlhs,numb(0))
    
    def equalMinusIndsP(self,other,freeInds):
        oA = [other.lhs,other.rhs]
        sA = [self.lhs,self.rhs]
        dd = {}
        for i in range(len(sA)):
            tmp = sA[i].equalMinusInds(oA[i],freeInds)
            if tmp == False:
                return False
            dd.update(tmp)
            for k in tmp.keys():
                freeInds = copy.deepcopy(freeInds)
                freeInds.remove(k)
                other = copy.deepcopy(other)
                other.changeInd(k,tmp[k])
                oA = [other.lhs,other.rhs]
        return dd
    def eqP(self,other):
        return self.lhs==other.lhs and self.rhs==other.rhs
    
    def sortStr(self):
        pass
    
    def doSort(self):
        self.lhs.doSort()
        self.rhs.doSort()
    
    def isSimplified(self):
        return self.simplifiedBool and self.rhs.isSimplified() and self.lhs.isSimplified()
    
    def changeInd(self,fromI,toI): ## returns True is an index was changed
        return self.lhs.changeInd(fromI,toI) or self.rhs.changeInd(fromI,toI)
    def factoriseNest(self,term1,doSimp=True):
        if self.lhs != numb(0):
            nlhs = self.lhs.factoriseNest(term1)
        else:
            nlhs = numb(0)
        if self.rhs != numb(0):
            nrhs = self.rhs.factoriseNest(term1)
        else:
            nrhs = numb(0)
        ret = equality(nlhs,nrhs)
        if doSimp:
            ret = ret.simplify()
        return ret
    def substitute(self,fromm,to):
        if self == fromm:
            return to
        else:
            return equality(self.lhs.substitute(fromm,to),self.rhs.substitute(fromm,to))
    def substituteMinusIndP(self,fromm,to,freeInds):
        dd = self.equalMinusInds(fromm,freeInds)
        if dd==False:
            return equality(self.lhs.substituteMinusIndP(fromm,to,freeInds),self.rhs.substituteMinusIndP(fromm,to,freeInds))
        else:
            to = copy.deepcopy(to)
            for k in dd.keys():
                to.changeInd(k,dd[k])
            return to
    def solveFor(self,term):
        print("####################\nOnly works for fully expanded equations, no fractions")
        ret = self.simplify({"preserveEq":False})
        if not isinstance(ret.lhs,addi):
            raise Exception("LHS must be addition")
        powA = ret.lhs.groupPower(term)
        if len(powA)>2:
            raise Exception("Only implemented linear")
        if len(powA)<2:
            raise Exception("Term not in equation")
        ret = equality(term,  frac(minus(None,powA[0]),powA[1]).simplify()  )
        return ret
    
    ## combines everything into a single fraction, required for simplification
    def mergeFrac(self):
        return equality(self.lhs.mergeFrac(),self.rhs.mergeFrac())
    
    def collectByTrace(self, incCorr=False):
        return equality( self.lhs.collectByTrace(incCorr), self.rhs.collectByTrace(incCorr) )
    
    def isolateTerm(self,trm):
        if isinstance(self.lhs,addi):
            lhsA = self.lhs.termA
        else:
            lhsA = [self.lhs]
        if isinstance(self.rhs,addi):
            rhsA = self.rhs.termA
        else:
            rhsA = [self.rhs]
        
        nlhsA = []
        nrhsA = []
        
        freeI = trm.getInds()
        for tt in lhsA:
            if tt.equalMinusInds(trm,freeI):
                nlhsA.append(tt)
            else:
                nrhsA.append(mult([numb(-1),tt]))
        for tt in rhsA:
            if tt.equalMinusInds(trm,freeI):
                nlhsA.append(mult([numb(-1),tt]))
            else:
                nrhsA.append(tt)
        
        if len(nlhsA) == 0:
            nlhs = numb(0)
        else:
            nlhs = addi(nlhsA).simplify()
        if len(nrhsA) == 0:
            nrhs = numb(0)
        else:
            nrhs = addi(nrhsA).simplify()
        
        return equality( nlhs, nrhs ).simplify()
            

class deriv(term):
    containsTerms = True;
    sortd = False
    def __init__(self,num,den):
        if den.containsOp():
            raise Exception("operator in derivative")
        self.num = num
        self.den = den
    def toStr(self):
        return "(d {})/(d {})".format( self.num.toStr(), self.den.toStr() )
    def toLatex(self, top=False):
        return "\\frac{{\dd {}}}{{\dd {}}}".format( self.num.toLatex(), self.den.toLatex() )
    def toCodeP(self,dicc):
        return "deriv({},{})".format( self.num.toCode(dicc), self.den.toCode(dicc) )
    def containsOp(self):
        return self.num.containsOp()
    def HermConj(self):
        return frac(self.num.HermConj(),self.den.HermConj())
    def getInds(self):
        return myFlatten( self.num.getInds() + self.den.getInds() ) ########## use set? !! #####
    def doSort(self):
        if not self.sortd:
            self.num.doSort()
            self.sortd = True
    def simplifyTypeP(self,dicc):
        self.num = self.num.simplifyType(dicc)
        self.den = self.den.simplifyType(dicc)
        
        if isinstance(self.num,numb) or isinstance(self.num,delta):
            return numb(0)
        
        if isinstance(self.num, addi):
            return addi(list(map(lambda t: deriv(t,self.den),self.num.termA)))
        if isinstance(self.num, mult):
            pA = self.num.termA
            tA = []
            for i in range(len(pA)):
                nt = mult(pA[:i]+[deriv(pA[i],self.den)]+pA[i+1:])
                tA.append(nt)
            return addi(tA)
        if isinstance(self.num, frac):
            old = self.num
            nn = minus( mult([deriv(old.num,self.den),old.den]),  mult([old.num,deriv(old.den,self.den)]) )
            nd = mult([ old.den, old.den ])
            return frac(nn,nd)
                        
        return self
    def equalMinusIndsP(self,other,freeInds):
        oA = [other.num,other.den]
        sA = [self.num,self.den]
        dd = {}
        for i in range(len(sA)):
            tmp = sA[i].equalMinusInds(oA[i],freeInds)
            if tmp == False:
                return False
            dd.update(tmp)
            for k in tmp.keys():
                freeInds = copy.deepcopy(freeInds)
                freeInds.remove(k)
                other = copy.deepcopy(other)
                other.changeInd(k,tmp[k])
                oA = [other.num,other.den]
        return dd
    def eqP(self,other):
        if self.num == other.num and self.den == other.den:
            return True
        else:
            return False
    def sortStr(self):
        ret = self.num.sortStr()+self.den.sortStr()
        return ret
    def isSimplified(self):
        return self.simplifiedBool and self.num.isSimplified() and self.den.isSimplified()
    def changeInd(self,fromI,toI): ## THIS OPERATES STRANGELY!!
        boolA = [self.num.changeInd(fromI,toI), self.den.changeInd(fromI,toI)]
        if any(boolA):
            self.simplifiedBool = False
            return True
        else:
            return False
    def factoriseNest(self,term1):
        num = self.num.factoriseNest(term1)
        return deriv(num,self.den).simplifyType(dicc)
    def substitute(self,fromm,to):
        if self == fromm:
            return to
        else:
            return deriv(self.num.substitute(fromm,to),self.den.substitute(fromm,to))
    def substituteMinusIndP(self,fromm,to,freeInds):
        dd = self.equalMinusInds(fromm,freeInds)
        if dd==False:
            return deriv(self.num.substituteMinusIndP(fromm,to,freeInds),self.den.substituteMinusIndP(fromm,to,freeInds))
        else:
            to = copy.deepcopy(to)
            for k in dd.keys():
                to.changeInd(k,dd[k])
            return to
    def collectByTrace(self, incCorr=False):
        return deriv( self.num.collectByTrace(incCorr), self.den )
    def getTrTermA(self, incCorr=False):
        return self.num.getTrTermA( incCorr)
    
    

def sortFunc(term1,term2):
    sortArr = [numb,sym,op,rho,delta,deriv,trac,summ,frac,addi,mult]

    if term1.containsOpForSort() and (not term2.containsOpForSort()):
        return 1
    if (not term1.containsOpForSort()) and term2.containsOpForSort():
        return -1
    try:
        ind1 = sortArr.index(type(term1))
    except Exception:
        print("term1: {}, type: {} not in {}".format(term1.toStr(),type(term1),sortArr))
        raise Exception()
    try:
        ind2 = sortArr.index(type(term2))
    except Exception:
        print("term2: {}, type: {} not in {}".format(term2,type(term2),sortArr))
        raise Exception()
    
    if ind1<ind2:
        return -1
    elif ind1>ind2:
        return 1
    else:
        if term1.sortStr() < term2.sortStr():
            return -1
        elif term1.sortStr() > term2.sortStr():
            return 1
        else:
            return 0

def minus(a,b):
    if a == None:
        return mult([numb(-1),b])
    else:
        return addi([a,mult([numb(-1),b])])
