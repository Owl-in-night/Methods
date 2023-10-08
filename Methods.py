from cProfile import label
import numpy as np
from sympy import*
import matplotlib.pyplot as plt
class Method:
    def __init__(self,x,y,a):
        self.x=x
        self.y=y
        self.a=a
    def lagrange(self):
        result=0
        for i in range(len(self.y)):
            value=self.y[i]
            
            for j in range(len(self.y)):
                if i!=j:
                    value*=(self.a-self.x[j])/(self.x[i]-self.x[j])
            result+=value
        return result
    def PolLagrange(self):
        x=Symbol("x")
        polynom=0
        for i in range(len(self.x)):
            term=1
            for j in range(len(self.x)):
                if j!=i:
                    term*=(x-self.x[j])/(self.x[i]-self.x[j])
            polynom+=term*self.y[i]
        return polynom
    def table(self):
        matrix=[[0]*len(self.x) for i in range(len(self.x))]
        for i in range(len(self.x)):
            matrix[i][0]=self.y[i]
        for i in range(1,len(self.x)):
            for j in range(i,len(self.x)):
                matrix[j][i]=(matrix[j][i-1]-matrix[j-1][i-1])/(self.x[j]-self.x[j-i])
        return matrix
    def diagonal(self):
        matrix=self.table()
        list=[matrix[i][i] for i in range(len(matrix))]
        return list
    def newton(self):
        result=self.diagonal()[0]
        for i in range(1,len(self.diagonal())):
            value=self.diagonal()[i]
            for j in range(i):
                value*=(self.a-self.x[j])
            result+=value
        return result
    def PolNewton(self):
        x=Symbol("x")
        polynom=self.y[0]
        for i in range(1,len(self.x)):
            factor=self.diagonal()[i]
            term=1
            for j in range(i):
                term*=(x-self.x[j])
            polynom+=term*factor
        return polynom
    
    def evaluation(self,f,x):
        return eval(f)
    
    def graph(self,method,f):
        x=Symbol("x")
        methods={"1":["Interpolación de Lagrange", self.PolLagrange()],"2":["Interpolación de Newton",self.PolNewton()]}
        plt.subplots(figsize=(15,10))
        plt.grid(1)
        if method=="1":
            plt.title(methods[method][0])
            px=lambdify(x,methods[method][1])
        elif method=="2":
            plt.title(methods[method][0])
            px=lambdify(x,methods[method][1])
        xi=np.linspace(min(self.x),max(self.x),100)
        fi=px(xi)
        plt.ylabel("y")
        plt.xlabel("x")
        plt.plot(xi,fi,color="b",label="Interpolation polynom")
        plt.plot(xi,self.evaluation(f,xi),color="r", label=f"${f}$")
        for i in range(len(self.x)):
            plt.scatter(self.x[i],self.y[i],color="purple",linewidths=4)
        plt.scatter(self.a,px(self.a),color="g",label="Point evaluated")
        plt.legend(loc="best")
        plt.show()