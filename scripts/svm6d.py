#! /usr/bin/env python
#$ -S /gs/software/python-2.4.2/bin/python
import assess
import datafunc
import svm
import libsvmWrap
import os
from ext.libsvm import LINEAR, POLY, RBF, SIGMOID, LINEARRIDGE
from ext import libsvm

class DifferentCSVM (svm.SVM) :
    """ An extension of the SVM so i can manipulate Cpos and Cneg """
    Cpos=1.0
    Cneg=1.0

    def __init__(self, arg = None, **args):
        svm.SVM.__init__(self, arg, **args)


    def trainLibsvm(self, data, **args) :
        Cpos=self.Cpos
        Cneg=self.Cneg

        print 'Cpos, Cneg: ', Cpos,Cneg

        # prepare data for the libsvm wrapper :
        # set kernel:, plotStrings=plotStrings
        if hasattr(self, 'kernel') and self.kernel is not None :
            kernel = self.kernel
        else :
            kernel = data.kernel
        kernelType = kernel.__class__.__name__
        #param = libsvm.svm_parameter(kernel_type = LINEAR)
        param = libsvmWrap.svmParameter(kernel_type = LINEAR,
                                        svm_type = self.svm_type)
        param.C = self.C
        param.nu = self.nu
        param.shrinking = 1
        cridge = None
        if kernelType == "Polynomial" :
            # (gamma x' y + coef0)^degree
            param.kernel_type = POLY
            param.degree = kernel.degree
            param.coef0 = kernel.additiveConst
            param.gamma = 1
        elif kernelType == "Gaussian":
            # exp(-gamma * |x - y|^2), plotStrings=plotStrings
            param.kernel_type = RBF
            param.gamma = kernel.gamma
        elif kernelType == "Cosine" :
            # i'm using the sigmoid kernel as the cosine kernel
            param.kernel_type = SIGMOID
            param.kernel_type = SIGMOID
        elif kernelType == "LinearRidge" :
            print 'using linear ridge kernel'
            param.kernel_type = LINEARRIDGE
            if not hasattr(self, 'ridge') :
                RidgeList = [1.0/Cneg, 1.0/Cpos]
                #Cpos = Cneg = 1e10
                self.ridge = [RidgeList[data.labels.Y[i]]
                              for i in range(len(data))]
            cridge = arrayWrap.doubleArray(self.ridge)
            param.ridge = cridge.cast()

        s=libsvm.DecisionFunction()
        if data.__class__.__name__ == 'SparseCDataSet' :
            data.convert2libsvm()
            prob = data.prob
            libsvm.svm_train_one4(prob, param.this, Cpos, Cneg, s.this)
            data.libsvmDestroy()
        else :
            prob=libsvmWrap.svmProblem (data)
            libsvm.svm_train_one4(prob.this, param.this, Cpos, Cneg, s.this)
            del prob

        if cridge is not None :
            del cridge

        b = - s.rho

        numSV = s.numSV
        alphaArray = libsvm.doubleArray_frompointer(s.alpha)
        svArray = libsvm.intArray_frompointer(s.svID)
        svID = [svArray[i] for i in range(numSV)]
        alpha = [alphaArray[i] for i in range(numSV)]

        return alpha, b, svID

    def getC(self, data) :
        c = [self.Cneg, self.Cpos]
        C = [c[data.labels.Y[i]] for i in range(len(data))]
        return C

def scoreLinear(X,lab,w,fn):
  scores=[]
  for ix in range(len(X)):
    score = 0.0
    for i in range(len(w)):
      score += X[ix][i]*w[i]
    scores += [(score,lab[ix])]
  scores.sort(reverse=True)
  labels = [sc[1] for sc in scores]
  f=open(fn + '.res',"w")
  for l in labels:
    f.write(str(l)+'\n')
  f.close

d = datafunc.DataSet("../6d.pyml", labelsColumn = -1)
s = DifferentCSVM()
s.Cpos = 10**float(os.environ.get("CPOS"))
s.Cneg = 10**float(os.environ.get("CNEG"))
fn = os.environ.get("CPOS") + '-' + os.environ.get("CNEG")
s.train(d)
r=s.test(d)
print("Model parameters:")
print(s.model.b)
print(s.model.w)
scoreLinear(d.X,d.labels.Y,s.model.w,fn)
